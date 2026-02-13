# =============================================================================
# SGCB differential expression analysis core function
# Unified pipeline: information-geometric natural gradient + hierarchical Bayes
# + manifold distance test
# Optional: calibrated bootstrap confidence intervals
# Rcpp-accelerated
# =============================================================================

#' @useDynLib SGCB, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats p.adjust rnorm median
NULL


#' SGCB Differential Expression Analysis
#'
#' Information-geometric differential expression analysis based on the
#' generalized gamma (GG) distribution.
#' Fits GG(\eqn{\alpha, \beta, \gamma}) parameters via natural gradient +
#' hierarchical Bayesian MAP estimation, limma-style fitFDist variance
#' shrinkage, Firth penalty (n <= 10), and manifold distance Wald
#' \eqn{\chi^2} test for distributional shift detection (active when
#' n >= 6 per group; contributes to SGCB_Score and pvalue_dd).
#' Output includes independent DD (Differential Distribution) p-values and
#' SGCB_Score omnibus score (Cauchy combination of DE + DD + DV).
#'
#' @param counts Count matrix (genes x samples)
#' @param group Grouping vector (two-level factor or character)
#' @param alpha FDR threshold, default 0.1
#' @param use_manifold_test Whether to enable manifold distance test, default TRUE
#' @param bootstrap Whether to compute calibrated bootstrap CIs, default FALSE
#' @param n_boot Number of bootstrap replicates, default 200 (only used when bootstrap=TRUE)
#' @param min_count Minimum count threshold, default 10
#' @param min_samples Minimum number of samples, default 2
#' @param config SGCBConfig configuration object (optional)
#' @param ... Backward-compatible parameters (method, n_dropout, dropout_rate etc. are deprecated; silently ignored)
#' @return SGCBResults object
#' @export
#' @examples
#' \donttest{
#' set.seed(42)
#' counts <- matrix(rnbinom(500 * 6, mu = 200, size = 10), nrow = 500, ncol = 6)
#' rownames(counts) <- paste0("gene_", seq_len(500))
#' counts[1:50, 4:6] <- counts[1:50, 4:6] * 3
#' group <- rep(c("ctrl", "treat"), each = 3)
#' res <- sgcbDE(counts, group)
#' print(res)
#' sig <- significantGenes(res, padj_cutoff = 0.05)
#' }
sgcbDE <- function(counts, group, alpha = 0.1,
                   use_manifold_test = TRUE,
                   bootstrap = FALSE, n_boot = 200,
                   min_count = 10, min_samples = 2,
                   config = NULL, ...) {

    cfg <- if (is.null(config)) defaultSGCBConfig() else config
    if (is.null(config)) cfg@bootB <- as.integer(n_boot)

    # Input validation
    stopifnot(is.matrix(counts) || is.data.frame(counts))
    counts <- as.matrix(counts)
    stopifnot(ncol(counts) == length(group))

    # Gene filtering
    keep <- rowSums(counts >= min_count) >= min_samples
    counts_filt <- counts[keep, , drop = FALSE]
    gene_names <- rownames(counts_filt)
    n_genes <- nrow(counts_filt)
    if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(n_genes))

    # Group assignment
    group <- factor(group)
    stopifnot(nlevels(group) == 2)
    group_int <- as.integer(group) - 1L
    ctrl_idx <- which(group_int == 0)
    treat_idx <- which(group_int == 1)
    n_ctrl <- length(ctrl_idx)
    n_treat <- length(treat_idx)

    # ==========================================================================
    # Normalization (DESeq2-style median-of-ratios, Anders & Huber 2010)
    # Robust to DE genes: takes the median of per-gene ratios to geometric mean
    # ==========================================================================
    log_counts <- log(counts_filt + 1)
    geo_mean_per_gene <- exp(rowMeans(log_counts))
    usable <- geo_mean_per_gene > 0
    ratios <- counts_filt[usable, , drop = FALSE] / geo_mean_per_gene[usable]
    sf <- apply(ratios, 2, median)
    sf <- sf / exp(mean(log(sf)))
    sf[sf < 0.01] <- 1
    norm_counts <- t(t(counts_filt) / sf)

    eps <- 1e-8

    # ==========================================================================
    # Core inference: information geometry (C++)
    # - Natural gradient + Fisher matrix GG parameter fitting (hierarchical Bayes MAP)
    # - Firth/Jeffreys penalty (n <= 10) + gamma=1 Gamma reduction (n <= 5)
    # - limma-style fitFDist variance shrinkage (Smyth 2004)
    # - Moderated t-test + manifold distance Wald chi-squared test
    # - min-Bonferroni combination
    # ==========================================================================
    result <- sgcb_info_geom_inference_cpp(
        norm_counts, group_int,
        n_dropout = 0L,
        dropout_rate = 0.0,
        use_natural_grad = TRUE,
        use_hierarchical = TRUE,
        eps = eps
    )

    # ==========================================================================
    # p-value combination and multiple testing correction
    # ==========================================================================
    n_min_grp <- min(n_ctrl, n_treat)
    manifold_active <- use_manifold_test && n_min_grp >= 6

    padj_t <- p.adjust(result$pvalue_t, method = "BH")
    padj_mu_wald <- p.adjust(result$pvalue_mu_wald, method = "BH")
    padj_manifold <- p.adjust(result$pvalue_manifold, method = "BH")
    padj_manifold[!manifold_active] <- NA_real_

    # ==========================================================================
    # Cauchy combination function (Liu & Xie 2020 JASA)
    # T = Σ wᵢ·tan(π(0.5 - pᵢ)), p = 0.5 - arctan(T/Σwᵢ)/π
    # ==========================================================================
    cauchy_combine <- function(p_mat, weights = NULL) {
        pm <- base::as.matrix(p_mat)
        pm[pm < eps] <- eps
        pm[pm > 1 - eps] <- 1 - eps
        n_p <- base::ncol(pm)
        w <- if (is.null(weights)) rep(1, n_p) else weights
        tan_vals <- tan(pi * (0.5 - pm))
        T_vec <- base::as.vector(tan_vals %*% w)
        p_out <- 0.5 - atan(T_vec / sum(w)) / pi
        p_out[p_out < eps] <- eps
        p_out[p_out > 1] <- 1
        p_out
    }

    # ==========================================================================
    # Primary DE test: moderated t-test (robust to arbitrary data-generating process)
    # pvalue_mu_wald retained as independent channel, only participates in SGCB_Score omnibus
    # ==========================================================================
    pvalue <- result$pvalue_t
    padj <- p.adjust(pvalue, method = "BH")

    # ==========================================================================
    # LFC shrinkage (Cauchy prior, apeglm-style)
    # ==========================================================================
    lfc_se <- sqrt(result$var_shrunk * (1.0 / n_ctrl + 1.0 / n_treat) + eps)
    lfc_shrunk <- shrink_lfc_cauchy(result$log2FC, lfc_se, prior_scale = 1.0)

    # ==========================================================================
    # DV/DG tests (differential variance / differential dynamics)
    # ==========================================================================
    dv_result <- .sgcbDVTest(
        result$ctrl_alpha, result$ctrl_beta, result$ctrl_gamma,
        result$treat_alpha, result$treat_beta, result$treat_gamma,
        n_ctrl, n_treat, eps
    )
    dg_result <- .sgcbDGTest(
        result$ctrl_alpha, result$ctrl_gamma,
        result$treat_alpha, result$treat_gamma,
        n_ctrl, n_treat, eps
    )
    dv_padj <- p.adjust(dv_result$pvalue, method = "BH")
    dg_alpha_padj <- p.adjust(dg_result$alpha_pvalue, method = "BH")
    dg_gamma_padj <- p.adjust(dg_result$gamma_pvalue, method = "BH")

    # ==========================================================================
    # Optional: GG-calibrated bootstrap (S2)
    # Posterior predictive parameter perturbation + variance-calibrated statistic
    # Perturb GG params ~ N(theta_hat, I(theta_hat)^{-1}/n) -> wider null -> better FDR control
    # ==========================================================================
    boot_cols <- NULL
    if (bootstrap) {
        ctrl_mat <- norm_counts[, ctrl_idx, drop = FALSE]
        treat_mat <- norm_counts[, treat_idx, drop = FALSE]
        treat_vec <- as.vector(treat_mat) + eps
        alpha_rep <- rep(result$ctrl_alpha, n_treat)
        beta_rep <- rep(result$ctrl_beta, n_treat)
        gamma_rep <- rep(result$ctrl_gamma, n_treat)
        ll_treat_vec <- gg_loglik_vec(treat_vec, alpha_rep, beta_rep, gamma_rep, eps)
        ll_treat_mat <- matrix(ll_treat_vec, nrow = n_genes, ncol = n_treat)
        T_obs <- -rowMeans(ll_treat_mat)

        # GG-calibrated: posterior parameter perturbation bootstrap
        # Approximate Fisher information -> parameter SE: sigma(log alpha) ~ 1/sqrt(n * trigamma(alpha) * alpha^2)
        B <- cfg@bootB
        a_ctrl <- pmax(result$ctrl_alpha, 0.05)
        b_ctrl <- pmax(result$ctrl_beta, eps)
        g_ctrl <- pmax(result$ctrl_gamma, 0.05)
        se_log_a <- 1 / sqrt(n_ctrl * trigamma(a_ctrl) * a_ctrl^2 + eps)
        se_log_b <- 1 / sqrt(n_ctrl * a_ctrl * g_ctrl + eps)
        se_log_g <- 1 / sqrt(n_ctrl * (a_ctrl / g_ctrl^2 + eps) + eps)

        T_null_all <- matrix(NA_real_, n_genes, B)
        sapply(seq_len(B), function(b) {
            a_pert <- a_ctrl * exp(rnorm(n_genes, 0, se_log_a))
            b_pert <- b_ctrl * exp(rnorm(n_genes, 0, se_log_b))
            g_pert <- g_ctrl * exp(rnorm(n_genes, 0, se_log_g))
            a_pert <- pmax(a_pert, 0.05)
            b_pert <- pmax(b_pert, eps)
            g_pert <- pmax(g_pert, 0.05)
            null_samples <- gg_sample_mat(n_genes, as.integer(n_treat),
                                          a_pert, b_pert, g_pert, eps)
            null_vec <- as.vector(null_samples) + eps
            a_rep_b <- rep(a_ctrl, n_treat)
            b_rep_b <- rep(b_ctrl, n_treat)
            g_rep_b <- rep(g_ctrl, n_treat)
            ll_null_vec <- gg_loglik_vec(null_vec, a_rep_b, b_rep_b, g_rep_b, eps)
            ll_null_mat <- matrix(ll_null_vec, nrow = n_genes, ncol = n_treat)
            T_null_all[, b] <<- -rowMeans(ll_null_mat)
            NULL
        })

        pvalue_bootstrap <- compute_empirical_pvalues(T_null_all, T_obs)
        quantiles_boot <- row_quantiles_cpp(T_null_all, c(0.025, 0.5, 0.975), 0L)

        boot_cols <- data.frame(
            T_obs = T_obs,
            T_null_median = quantiles_boot[, 2],
            T_CI_lo = quantiles_boot[, 1],
            T_CI_hi = quantiles_boot[, 3],
            pvalue_bootstrap = pvalue_bootstrap
        )
    }

    # ==========================================================================
    # DD independent output (Differential Distribution)
    # Cauchy combination of manifold + DV -> pvalue_dd
    # ==========================================================================
    pvalue_dd <- rep(NA_real_, n_genes)
    if (manifold_active) {
        p_dd_mat <- cbind(result$pvalue_manifold, dv_result$pvalue)
        pvalue_dd <- cauchy_combine(p_dd_mat)
    }
    padj_dd <- rep(NA_real_, n_genes)
    padj_dd[!is.na(pvalue_dd)] <- p.adjust(pvalue_dd[!is.na(pvalue_dd)], method = "BH")

    # ==========================================================================
    # SGCB_Score: Cauchy combination of DE(t) + GG(mu_wald) + DD(manifold) + DV [+ bootstrap]
    # Four-channel omnibus: each channel independently calibrated before combination
    # Score = -log10(p_omnibus), higher = stronger evidence
    # ==========================================================================
    p_channels <- list(result$pvalue_t, result$pvalue_mu_wald, dv_result$pvalue)
    if (manifold_active) p_channels <- c(list(result$pvalue_t, result$pvalue_mu_wald, result$pvalue_manifold), list(dv_result$pvalue))
    if (!is.null(boot_cols)) p_channels <- c(p_channels, list(boot_cols$pvalue_bootstrap))
    p_omni_mat <- do.call(cbind, p_channels)
    sgcb_score_p <- cauchy_combine(p_omni_mat)
    sgcb_score_padj <- p.adjust(sgcb_score_p, method = "BH")
    sgcb_score <- -log10(pmax(sgcb_score_p, eps))

    # ==========================================================================
    # Build results
    # ==========================================================================
    p_manifold_out <- result$pvalue_manifold
    p_manifold_out[!manifold_active] <- NA_real_

    results <- data.frame(
        gene_id = gene_names,
        baseMean = (result$mean_ctrl + result$mean_treat) / 2,
        log2FoldChange = result$log2FC,
        log2FC_gg = result$log2FC_gg,
        log2FC_shrunk = lfc_shrunk,
        lfcSE = lfc_se,
        mean_ctrl = result$mean_ctrl,
        mean_treat = result$mean_treat,
        geodesic_dist = result$geodesic_dist,
        ctrl_alpha = result$ctrl_alpha,
        ctrl_beta = result$ctrl_beta,
        ctrl_gamma = result$ctrl_gamma,
        treat_alpha = result$treat_alpha,
        treat_beta = result$treat_beta,
        treat_gamma = result$treat_gamma,
        dv_log2_var_ratio = dv_result$log2_var_ratio,
        dv_stat = dv_result$stat,
        dv_pvalue = dv_result$pvalue,
        dv_padj = dv_padj,
        dv_var_ctrl = dv_result$var_ctrl,
        dv_var_treat = dv_result$var_treat,
        dg_alpha_log2_ratio = dg_result$alpha_log2_ratio,
        dg_alpha_stat = dg_result$alpha_stat,
        dg_alpha_pvalue = dg_result$alpha_pvalue,
        dg_alpha_padj = dg_alpha_padj,
        dg_gamma_log2_ratio = dg_result$gamma_log2_ratio,
        dg_gamma_stat = dg_result$gamma_stat,
        dg_gamma_pvalue = dg_result$gamma_pvalue,
        dg_gamma_padj = dg_gamma_padj,
        stat = result$t_stat,
        pvalue = pvalue,
        pvalue_t = result$pvalue_t,
        pvalue_mu_wald = result$pvalue_mu_wald,
        pvalue_manifold = p_manifold_out,
        padj = padj,
        padj_t = padj_t,
        padj_mu_wald = padj_mu_wald,
        padj_manifold = padj_manifold,
        pvalue_dd = pvalue_dd,
        padj_dd = padj_dd,
        SGCB_Score = sgcb_score,
        SGCB_Score_p = sgcb_score_p,
        SGCB_Score_padj = sgcb_score_padj,
        row.names = gene_names,
        stringsAsFactors = FALSE
    )

    if (!is.null(boot_cols)) {
        results <- cbind(results, boot_cols)
    }

    # Metadata
    attr(results, "alpha") <- alpha
    attr(results, "method") <- "infogeom"
    attr(results, "df_prior") <- result$df_total
    attr(results, "n_ctrl") <- n_ctrl
    attr(results, "n_treat") <- n_treat
    attr(results, "n_filtered") <- sum(!keep)
    attr(results, "bootstrap") <- bootstrap
    attr(results, "var_prior_trend") <- result$var_prior_trend
    attr(results, "var_shrunk") <- result$var_shrunk
    attr(results, "gg_weight") <- result$gg_weight
    attr(results, "d_gg_eff") <- result$d_gg_eff

    class(results) <- c("SGCBResults", "data.frame")
    results
}

# =============================================================================
# Backward compatibility (deprecated, delegates to unified entry point)
# =============================================================================

#' @rdname sgcbDE-deprecated
#' @title Deprecated SGCB Functions
#' @description These functions are deprecated. Use \code{\link{sgcbDE}} instead.
#' @param counts Count matrix (genes x samples)
#' @param group Group vector
#' @param ... Passed to \code{sgcbDE}
#' @export
sgcbDE_fast <- function(counts, group, ...) {
    .Deprecated("sgcbDE", msg = "sgcbDE_fast is deprecated. Use sgcbDE() instead.")
    sgcbDE(counts, group, ...)
}

#' @rdname sgcbDE-deprecated
#' @export
sgcbDE_infogeom <- function(counts, group, ...) {
    .Deprecated("sgcbDE", msg = "sgcbDE_infogeom is deprecated. Use sgcbDE() instead.")
    sgcbDE(counts, group, ...)
}

#' Print SGCBResults
#' @param x SGCBResults object
#' @param ... Additional arguments
#' @export
print.SGCBResults <- function(x, ...) {
    alpha <- attr(x, "alpha")
    if (is.null(alpha)) { print.data.frame(x); return(invisible(x)) }
    n_sig <- sum(x$padj < alpha, na.rm = TRUE)
    cat("SGCB Differential Expression Results\n")
    cat("-------------------------------------\n")
    cat("Genes tested:", nrow(x), "\n")
    cat("Significant (padj <", alpha, "):", n_sig, "\n")
    df_prior <- attr(x, "df_prior")
    if (!is.null(df_prior)) cat("Prior df:", round(df_prior, 1), "\n")
    if (isTRUE(attr(x, "bootstrap"))) cat("Bootstrap CI: available\n")
    cat("\nTop genes:\n")
    top <- head(x[order(x$pvalue), c("log2FoldChange", "pvalue", "padj")], 6)
    print.data.frame(top)
    invisible(x)
}

#' Extract significant genes
#' @param x SGCBResults object
#' @param padj_cutoff Adjusted p-value threshold
#' @param lfc_cutoff Absolute log2FC threshold
#' @return Subset of significant genes (data.frame)
#' @export
#' @examples
#' \donttest{
#' set.seed(42)
#' counts <- matrix(rnbinom(500 * 6, mu = 200, size = 10), nrow = 500, ncol = 6)
#' rownames(counts) <- paste0("gene_", seq_len(500))
#' counts[1:50, 4:6] <- counts[1:50, 4:6] * 3
#' group <- rep(c("ctrl", "treat"), each = 3)
#' res <- sgcbDE(counts, group)
#' sig <- significantGenes(res, padj_cutoff = 0.05, lfc_cutoff = 1)
#' nrow(sig)
#' }
significantGenes <- function(x, padj_cutoff = 0.1, lfc_cutoff = 0) {
    stopifnot(inherits(x, "SGCBResults"))
    sig <- x[!is.na(x$padj) & x$padj < padj_cutoff & abs(x$log2FoldChange) > lfc_cutoff, ]
    sig <- sig[order(sig$pvalue), ]
    class(sig) <- "data.frame"
    sig
}
