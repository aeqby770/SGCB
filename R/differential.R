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


#' SGCB differential analysis for a two-group comparison
#'
#' Fits gene-wise generalized gamma (GG) models for a two-group RNA-seq
#' contrast and returns mean-shift and distribution-shift summaries in one
#' result table. The pipeline uses natural-gradient hierarchical MAP
#' estimation, limma-style variance moderation, and optional diagnostic
#' channels derived from the fitted GG parameters.
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
#' @param prior_de Prior inclusion probability for DE channel in joint posterior;
#'   numeric in (0,1) or \code{"auto"}.
#' @param prior_dv Prior inclusion probability for DV channel in joint posterior;
#'   numeric in (0,1) or \code{"auto"}.
#' @param prior_dg Prior inclusion probability for DG channel in joint posterior;
#'   numeric in (0,1) or \code{"auto"}.
#' @param prior_sd_de Prior SD of DE log2 effect for ABF evidence
#' @param prior_sd_dv Prior SD of DV log2 effect for ABF evidence
#' @param prior_sd_dg_alpha Prior SD of DG alpha log2 effect for ABF evidence
#' @param prior_sd_dg_gamma Prior SD of DG gamma log2 effect for ABF evidence
#' @param ... Backward-compatible parameters that are currently accepted but not used by the current implementation
#' @return An object of class \code{SGCBResults}
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
                   config = NULL,
                   prior_de = "auto",
                   prior_dv = "auto",
                   prior_dg = "auto",
                   prior_sd_de = 1.0,
                   prior_sd_dv = 0.5,
                   prior_sd_dg_alpha = 0.5,
                   prior_sd_dg_gamma = 0.5,
                   ...) {

    cfg <- if (is.null(config)) defaultSGCBConfig() else config
    if (is.null(config)) cfg@bootB <- as.integer(n_boot)
    dots <- list(...)
    use_natural_grad_opt <- if ("use_natural_grad" %in% names(dots)) isTRUE(dots$use_natural_grad) else TRUE
    use_hierarchical_opt <- if ("use_hierarchical" %in% names(dots)) isTRUE(dots$use_hierarchical) else TRUE
    use_firth_opt <- if ("use_firth" %in% names(dots)) isTRUE(dots$use_firth) else TRUE
    use_gamma_submodel_opt <- if ("use_gamma_submodel" %in% names(dots)) isTRUE(dots$use_gamma_submodel) else TRUE
    use_gg_variance_opt <- if ("use_gg_variance" %in% names(dots)) isTRUE(dots$use_gg_variance) else TRUE

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
    eps <- 1e-8

    # ==========================================================================
    # Normalization (edgeR-style TMM, Robinson & Oshlack 2010)
    # C++ path: tmm_factors_cpp(counts, ref_idx)
    # ==========================================================================
    lib_sizes <- colSums(counts_filt)
    ref_idx <- which.min(abs(lib_sizes - median(lib_sizes))) - 1L
    tmm <- tmm_factors_cpp(counts_filt, as.integer(ref_idx),
                           trim_m = 0.3, trim_a = 0.05,
                           eps = eps, n_threads = 0L)
    eff_lib <- lib_sizes * tmm
    sf <- eff_lib / exp(mean(log(eff_lib)))
    sf[sf < 0.01] <- 1
    norm_counts <- t(t(counts_filt) / sf)

    # Continuous relaxation for zero counts under continuous GG likelihood
    norm_counts <- norm_counts + 0.5

    # ==========================================================================
    # Core inference: information geometry (C++)
    # - Natural gradient + Fisher matrix GG parameter fitting (hierarchical Bayes MAP)
    # - Firth/Jeffreys penalty (n <= 10) + gamma=1 Gamma reduction (n <= 5)
    # - limma-style fitFDist variance shrinkage (Smyth 2004)
    # - Moderated t-test + manifold distance Wald chi-squared test
    # - Channel aggregation is handled in R (DD and SGCB_Score via Cauchy combination)
    # ==========================================================================
    result <- sgcb_info_geom_inference_cpp(
        norm_counts, group_int,
        n_dropout = 0L,
        dropout_rate = 0.0,
        use_natural_grad = use_natural_grad_opt,
        use_hierarchical = use_hierarchical_opt,
        use_firth = use_firth_opt,
        use_gamma_submodel = use_gamma_submodel_opt,
        use_gg_variance = use_gg_variance_opt,
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
    # pvalue_mu_wald retained as diagnostic column (not used in primary/omnibus testing)
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
        n_ctrl, n_treat, eps,
        norm_counts = norm_counts, ctrl_idx = ctrl_idx, treat_idx = treat_idx
    )
    dg_result <- .sgcbDGTest(
        result$ctrl_alpha, result$ctrl_gamma,
        result$treat_alpha, result$treat_gamma,
        n_ctrl, n_treat, eps
    )
    dv_padj <- p.adjust(dv_result$pvalue, method = "BH")
    dv_padj_model <- p.adjust(dv_result$p_model, method = "BH")
    dv_padj_model_raw <- p.adjust(dv_result$p_model_raw, method = "BH")
    dv_padj_bf <- p.adjust(dv_result$p_bf, method = "BH")
    dv_padj_mix <- p.adjust(dv_result$p_mix, method = "BH")
    dg_alpha_padj <- p.adjust(dg_result$alpha_pvalue, method = "BH")
    dg_gamma_padj <- p.adjust(dg_result$gamma_pvalue, method = "BH")
    dg_entropy_padj <- p.adjust(dg_result$entropy_pvalue, method = "BH")
    dg_shape_padj <- p.adjust(dg_result$shape_pvalue, method = "BH")

    # ==========================================================================
    # Empirical-null SE calibration for ABF channels
    # ==========================================================================
    de_cal <- .sgcbCalibratedSE(result$log2FC, lfc_se, eps)
    dv_cal <- .sgcbCalibratedSE(dv_result$log2_cv2_ratio, dv_result$se_log2_cv2, eps)
    dg_alpha_cal <- .sgcbCalibratedSE(dg_result$alpha_log2_ratio, dg_result$alpha_se_log2, eps)
    dg_gamma_cal <- .sgcbCalibratedSE(dg_result$gamma_log2_ratio, dg_result$gamma_se_log2, eps)

    # ==========================================================================
    # ABF evidence from effect-size + uncertainty channels
    # Wakefield-style approximation (normal likelihood + normal prior)
    # ==========================================================================
    bf_de_abf <- .sgcbABF1(
        effect = result$log2FC,
        se = de_cal$se_cal,
        prior_sd = prior_sd_de
    )
    bf_dv_abf <- .sgcbABF1(
        effect = dv_result$log2_cv2_ratio,
        se = dv_cal$se_cal,
        prior_sd = prior_sd_dv
    )
    bf_dg_abf <- .sgcbABFDiag(
        effect_mat = cbind(dg_result$alpha_log2_ratio, dg_result$gamma_log2_ratio),
        se_mat = cbind(dg_alpha_cal$se_cal, dg_gamma_cal$se_cal),
        prior_sd = c(prior_sd_dg_alpha, prior_sd_dg_gamma)
    )

    # LR-BF proxy from model-based DV channel (less conservative than bootstrap-only p)
    lr_proxy_dv <- stats::qchisq(.sgcbClipProb(dv_result$pvalue, eps), df = 1, lower.tail = FALSE)
    lr_proxy_dg <- stats::qchisq(.sgcbClipProb(dg_result$shape_pvalue, eps), df = 2, lower.tail = FALSE)
    bf_dv_lr <- .sgcbLRBF(lr_proxy_dv, n_ctrl + n_treat, 1)
    bf_dg_lr <- .sgcbLRBF(lr_proxy_dg, n_ctrl + n_treat, 2)
    bf_de <- bf_de_abf
    bf_dv <- sqrt(bf_dv_abf) * bf_dv_lr
    bf_dg <- sqrt(bf_dg_abf) * bf_dg_lr

    prior_de_used <- .sgcbResolvePrior(
        prior_de,
        p = pvalue,
        bf = bf_de,
        prior_alpha = 1.0,
        prior_beta = 9.0
    )
    prior_dv_beta <- ifelse(min(n_ctrl, n_treat) <= 6, 999.0,
                            ifelse(min(n_ctrl, n_treat) <= 10, 199.0, 19.0))
    prior_dv_upper <- ifelse(min(n_ctrl, n_treat) <= 10, 0.10, 0.25)
    prior_dv_used <- .sgcbResolvePrior(
        prior_dv,
        p = dv_result$pvalue,
        bf = bf_dv,
        upper = prior_dv_upper,
        prior_alpha = 1.0,
        prior_beta = prior_dv_beta
    )
    prior_dg_used <- .sgcbResolvePrior(
        prior_dg,
        p = dg_result$shape_pvalue,
        bf = bf_dg,
        prior_alpha = 1.0,
        prior_beta = 9.0
    )

    # ==========================================================================
    # Joint DE/DV/DG posterior and model probabilities (8 models)
    # ==========================================================================
    joint_post <- .sgcbModelPosterior3(
        p_de = pvalue,
        p_dv = dv_result$pvalue,
        p_dg = dg_result$shape_pvalue,
        prior_de = prior_de_used,
        prior_dv = prior_dv_used,
        prior_dg = prior_dg_used,
        bf_de = bf_de,
        bf_dv = bf_dv,
        bf_dg = bf_dg
    )

    de_call <- .sgcbPosteriorCalls(joint_post$p_de_post, target_fdr = alpha)
    dv_call <- .sgcbPosteriorCalls(joint_post$p_dv_post, target_fdr = alpha)
    dg_call <- .sgcbPosteriorCalls(joint_post$p_dg_post, target_fdr = alpha)
    de_fdr_call <- padj <= alpha
    dv_fdr_call <- dv_padj <= alpha
    dv_model_fdr_call <- dv_padj_model <= alpha
    dv_model_raw_fdr_call <- dv_padj_model_raw <= alpha
    dv_bf_fdr_call <- dv_padj_bf <= alpha
    dv_mix_fdr_call <- dv_padj_mix <= alpha
    dg_fdr_call <- dg_shape_padj <= alpha

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
        se_log_b <- 1 / sqrt(n_ctrl * a_ctrl * g_ctrl^2 + eps)
        se_log_g <- 1 / sqrt(n_ctrl * (1 + a_ctrl * digamma(a_ctrl)^2 + 2 * digamma(a_ctrl) + a_ctrl * trigamma(a_ctrl)) + eps)

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

        quantiles_boot <- row_quantiles_cpp(T_null_all, c(0.025, 0.5, 0.975), 0L)
        T_null_median <- quantiles_boot[, 2]
        dev_obs <- abs(T_obs - T_null_median)
        dev_null <- abs(T_null_all - T_null_median)
        pvalue_bootstrap <- (1 + rowSums(dev_null >= dev_obs)) / (1 + B)

        boot_cols <- data.frame(
            T_obs = T_obs,
            T_null_median = T_null_median,
            T_CI_lo = quantiles_boot[, 1],
            T_CI_hi = quantiles_boot[, 3],
            pvalue_bootstrap = pvalue_bootstrap
        )
    }

    # ==========================================================================
    # DD independent output (Differential Distribution)
    # Cauchy combination of DV(CV²) + DG(gamma), avoiding manifold double-counting
    # ==========================================================================
    p_dd_mat <- cbind(dv_result$pvalue, dg_result$shape_pvalue)
    pvalue_dd <- cauchy_combine(p_dd_mat)
    padj_dd <- p.adjust(pvalue_dd, method = "BH")

    # ==========================================================================
    # SGCB_Score: Cauchy combination of DE(t) + DD [+ bootstrap]
    # Keep channels as orthogonal as possible; manifold/mu_wald are diagnostic columns
    # Score = -log10(p_omnibus), higher = stronger evidence
    # ==========================================================================
    p_channels <- list(result$pvalue_t, pvalue_dd)
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
        dv_log2_cv2_ratio = dv_result$log2_cv2_ratio,
        dv_se_log2_cv2 = dv_result$se_log2_cv2,
        dv_stat = dv_result$stat,
        dv_pvalue = dv_result$pvalue,
        dv_padj = dv_padj,
        dv_pvalue_model = dv_result$p_model,
        dv_padj_model = dv_padj_model,
        dv_pvalue_model_raw = dv_result$p_model_raw,
        dv_padj_model_raw = dv_padj_model_raw,
        dv_pvalue_bf = dv_result$p_bf,
        dv_padj_bf = dv_padj_bf,
        dv_pvalue_mix = dv_result$p_mix,
        dv_padj_mix = dv_padj_mix,
        dv_var_ctrl = dv_result$var_ctrl,
        dv_var_treat = dv_result$var_treat,
        dv_cv2_ctrl = dv_result$cv2_ctrl,
        dv_cv2_treat = dv_result$cv2_treat,
        dg_alpha_log2_ratio = dg_result$alpha_log2_ratio,
        dg_alpha_se_log2 = dg_result$alpha_se_log2,
        dg_alpha_stat = dg_result$alpha_stat,
        dg_alpha_pvalue = dg_result$alpha_pvalue,
        dg_alpha_padj = dg_alpha_padj,
        dg_gamma_log2_ratio = dg_result$gamma_log2_ratio,
        dg_gamma_se_log2 = dg_result$gamma_se_log2,
        dg_gamma_stat = dg_result$gamma_stat,
        dg_gamma_pvalue = dg_result$gamma_pvalue,
        dg_gamma_padj = dg_gamma_padj,
        dg_entropy_diff = dg_result$entropy_diff,
        dg_entropy_stat = dg_result$entropy_stat,
        dg_entropy_pvalue = dg_result$entropy_pvalue,
        dg_entropy_padj = dg_entropy_padj,
        dg_shape_stat = dg_result$shape_stat,
        dg_shape_pvalue = dg_result$shape_pvalue,
        dg_shape_padj = dg_shape_padj,
        stat = result$t_stat,
        pvalue = pvalue,
        pvalue_t = result$pvalue_t,
        pvalue_mu_wald = result$pvalue_mu_wald,
        pvalue_manifold = p_manifold_out,
        padj = padj,
        padj_t = padj_t,
        padj_mu_wald = padj_mu_wald,
        padj_manifold = padj_manifold,
        de_fdr_call = de_fdr_call,
        dv_fdr_call = dv_fdr_call,
        dv_model_fdr_call = dv_model_fdr_call,
        dv_model_raw_fdr_call = dv_model_raw_fdr_call,
        dv_bf_fdr_call = dv_bf_fdr_call,
        dv_mix_fdr_call = dv_mix_fdr_call,
        dg_fdr_call = dg_fdr_call,
        p_de_post = joint_post$p_de_post,
        de_post_call = de_call$call,
        p_dv_post = joint_post$p_dv_post,
        dv_post_call = dv_call$call,
        p_dg_post = joint_post$p_dg_post,
        dg_post_call = dg_call$call,
        p_de_only_post = joint_post$p_de_only_post,
        p_dv_only_post = joint_post$p_dv_only_post,
        p_dg_only_post = joint_post$p_dg_only_post,
        p_de_dv_post = joint_post$p_de_dv_post,
        p_de_dg_post = joint_post$p_de_dg_post,
        p_dv_dg_post = joint_post$p_dv_dg_post,
        p_de_dv_dg_post = joint_post$p_de_dv_dg_post,
        bf_de = joint_post$bf_de,
        bf_dv = joint_post$bf_dv,
        bf_dg = joint_post$bf_dg,
        model_prob_null = joint_post$model_prob[, "model_prob_null"],
        model_prob_de = joint_post$model_prob[, "model_prob_de"],
        model_prob_dv = joint_post$model_prob[, "model_prob_dv"],
        model_prob_dg = joint_post$model_prob[, "model_prob_dg"],
        model_prob_de_dv = joint_post$model_prob[, "model_prob_de_dv"],
        model_prob_de_dg = joint_post$model_prob[, "model_prob_de_dg"],
        model_prob_dv_dg = joint_post$model_prob[, "model_prob_dv_dg"],
        model_prob_de_dv_dg = joint_post$model_prob[, "model_prob_de_dv_dg"],
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
    attr(results, "df_prior") <- result$df_total - (n_ctrl + n_treat - 2)
    attr(results, "n_ctrl") <- n_ctrl
    attr(results, "n_treat") <- n_treat
    attr(results, "n_filtered") <- sum(!keep)
    attr(results, "bootstrap") <- bootstrap
    attr(results, "var_prior_trend") <- result$var_prior_trend
    attr(results, "var_shrunk") <- result$var_shrunk
    attr(results, "gg_weight") <- result$gg_weight
    attr(results, "d_gg_eff") <- result$d_gg_eff
    attr(results, "use_natural_grad") <- use_natural_grad_opt
    attr(results, "use_hierarchical") <- use_hierarchical_opt
    attr(results, "use_firth") <- use_firth_opt
    attr(results, "use_gamma_submodel") <- use_gamma_submodel_opt
    attr(results, "use_gg_variance") <- use_gg_variance_opt
    attr(results, "prior_de") <- prior_de
    attr(results, "prior_dv") <- prior_dv
    attr(results, "prior_dg") <- prior_dg
    attr(results, "prior_de_used") <- prior_de_used
    attr(results, "prior_dv_used") <- prior_dv_used
    attr(results, "prior_dv_upper") <- prior_dv_upper
    attr(results, "prior_dg_used") <- prior_dg_used
    attr(results, "prior_sd_de") <- prior_sd_de
    attr(results, "prior_sd_dv") <- prior_sd_dv
    attr(results, "prior_sd_dg_alpha") <- prior_sd_dg_alpha
    attr(results, "prior_sd_dg_gamma") <- prior_sd_dg_gamma
    attr(results, "abf_inflation_de") <- de_cal$inflation
    attr(results, "abf_inflation_dv") <- dv_cal$inflation
    attr(results, "abf_inflation_dg_alpha") <- dg_alpha_cal$inflation
    attr(results, "abf_inflation_dg_gamma") <- dg_gamma_cal$inflation
    attr(results, "decision_mode") <- "frequentist_primary"
    attr(results, "primary_de_channel") <- "padj"
    attr(results, "primary_dv_channel") <- "dv_padj"
    attr(results, "primary_dg_channel") <- "dg_shape_padj"
    attr(results, "post_call_fdr") <- alpha
    attr(results, "post_de_threshold") <- de_call$threshold
    attr(results, "post_dv_threshold") <- dv_call$threshold
    attr(results, "post_dg_threshold") <- dg_call$threshold
    attr(results, "post_de_realized_fdr") <- de_call$realized_fdr
    attr(results, "post_dv_realized_fdr") <- dv_call$realized_fdr
    attr(results, "post_dg_realized_fdr") <- dg_call$realized_fdr

    class(results) <- c("SGCBResults", "data.frame")
    results
}

#' @noRd
.align_sgcb_sample_data <- function(counts, sample_data) {
    counts <- as.matrix(counts)
    sample_data <- as.data.frame(sample_data, stringsAsFactors = FALSE)
    stopifnot(ncol(counts) == nrow(sample_data))
    if (!is.null(colnames(counts)) && !is.null(rownames(sample_data)) && !identical(colnames(counts), rownames(sample_data))) {
        stopifnot(all(colnames(counts) %in% rownames(sample_data)))
        sample_data <- sample_data[colnames(counts), , drop = FALSE]
    }
    sample_data
}

#' @noRd
.resolve_sgcb_subset <- function(sample_data, counts, sample_subset = NULL) {
    idx <- rep(TRUE, nrow(sample_data))
    if (is.null(sample_subset)) {
        return(idx)
    }
    if (is.logical(sample_subset)) {
        stopifnot(length(sample_subset) == nrow(sample_data))
        idx <- sample_subset
    } else if (is.numeric(sample_subset)) {
        idx <- rep(FALSE, nrow(sample_data))
        idx[as.integer(sample_subset)] <- TRUE
    } else {
        sample_ids <- rownames(sample_data)
        if (is.null(sample_ids)) {
            sample_ids <- colnames(counts)
        }
        idx <- sample_ids %in% as.character(sample_subset)
    }
    idx
}

#' Contrast-defined deployment wrapper for SGCB
#'
#' Subsets a user-specified two-level comparison from \code{sample_data} and
#' forwards the resulting count matrix and group vector to \code{sgcbDE()}.
#' This helper does not fit covariates and does not alter the two-group
#' inferential core.
#'
#' @param counts Count matrix (genes x samples)
#' @param sample_data Sample metadata with one row per sample
#' @param group_col Column in \code{sample_data} defining the grouping factor
#' @param contrast_levels Character vector of length 2 giving the levels to compare, in control-then-treatment order
#' @param sample_subset Optional logical, numeric, or character subset applied before the contrast is formed
#' @inheritParams sgcbDE
#' @return An object of class \code{SGCBResults} with additional attributes describing the wrapper call
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' counts <- matrix(rnbinom(300 * 9, mu = 150, size = 8), nrow = 300, ncol = 9)
#' colnames(counts) <- paste0("sample_", seq_len(9))
#' sample_data <- data.frame(
#'   condition = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#'   row.names = colnames(counts)
#' )
#' res_ab <- sgcbContrast(counts, sample_data, group_col = "condition", contrast_levels = c("A", "B"))
#' }
sgcbContrast <- function(counts, sample_data, group_col, contrast_levels,
                         sample_subset = NULL, alpha = 0.1,
                         use_manifold_test = TRUE,
                         bootstrap = FALSE, n_boot = 200,
                         min_count = 10, min_samples = 2,
                        config = NULL, ...) {
    counts <- as.matrix(counts)
    sample_data <- .align_sgcb_sample_data(counts, sample_data)
    contrast_levels <- as.character(contrast_levels)
    stopifnot(length(contrast_levels) == 2)
    idx_subset <- .resolve_sgcb_subset(sample_data, counts, sample_subset)
    idx_contrast <- idx_subset & sample_data[[group_col]] %in% contrast_levels
    counts_sub <- counts[, idx_contrast, drop = FALSE]
    group_sub <- factor(sample_data[[group_col]][idx_contrast], levels = contrast_levels)
    stopifnot(nlevels(group_sub) == 2)
    res <- sgcbDE(
        counts = counts_sub,
        group = group_sub,
        alpha = alpha,
        use_manifold_test = use_manifold_test,
        bootstrap = bootstrap,
        n_boot = n_boot,
        min_count = min_count,
        min_samples = min_samples,
        config = config,
        ...
    )
    attr(res, "deployment_mode") <- "contrast_defined_wrapper"
    attr(res, "group_col") <- group_col
    attr(res, "contrast_levels") <- contrast_levels
    attr(res, "n_wrapper_samples") <- ncol(counts_sub)
    attr(res, "sample_subset_applied") <- !is.null(sample_subset)
    res
}

#' Pairwise deployment wrapper for SGCB
#'
#' Enumerates all requested two-level comparisons from \code{group_col} and
#' applies \code{sgcbContrast()} to each pair. Each element of the returned
#' list is still a separate two-group SGCB analysis.
#'
#' @param counts Count matrix (genes x samples)
#' @param sample_data Sample metadata with one row per sample
#' @param group_col Column in \code{sample_data} defining the grouping factor
#' @param sample_subset Optional logical, numeric, or character subset applied before pairwise contrasts are formed
#' @param contrast_levels Optional character vector restricting the levels to be paired
#' @inheritParams sgcbDE
#' @return A named list of \code{SGCBResults} objects with class \code{SGCBContrastList}
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' counts <- matrix(rnbinom(300 * 9, mu = 150, size = 8), nrow = 300, ncol = 9)
#' colnames(counts) <- paste0("sample_", seq_len(9))
#' sample_data <- data.frame(
#'   condition = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#'   row.names = colnames(counts)
#' )
#' res_list <- sgcbPairwise(counts, sample_data, group_col = "condition")
#' names(res_list)
#' }
sgcbPairwise <- function(counts, sample_data, group_col,
                         sample_subset = NULL, contrast_levels = NULL,
                         alpha = 0.1, use_manifold_test = TRUE,
                         bootstrap = FALSE, n_boot = 200,
                         min_count = 10, min_samples = 2,
                        config = NULL, ...) {
    counts <- as.matrix(counts)
    sample_data <- .align_sgcb_sample_data(counts, sample_data)
    idx_subset <- .resolve_sgcb_subset(sample_data, counts, sample_subset)
    available_levels <- unique(as.character(sample_data[[group_col]][idx_subset]))
    target_levels <- if (is.null(contrast_levels)) available_levels else as.character(contrast_levels)
    contrast_list <- utils::combn(target_levels, 2, simplify = FALSE)
    out <- lapply(contrast_list, function(level_pair) {
        sgcbContrast(
            counts = counts,
            sample_data = sample_data,
            group_col = group_col,
            contrast_levels = level_pair,
            sample_subset = idx_subset,
            alpha = alpha,
            use_manifold_test = use_manifold_test,
            bootstrap = bootstrap,
            n_boot = n_boot,
            min_count = min_count,
            min_samples = min_samples,
            config = config,
            ...
        )
    })
    names(out) <- vapply(contrast_list, function(level_pair) paste(level_pair, collapse = "_vs_"), character(1))
    attr(out, "deployment_mode") <- "pairwise_wrapper"
    attr(out, "group_col") <- group_col
    attr(out, "contrast_levels") <- target_levels
    attr(out, "sample_subset_applied") <- !is.null(sample_subset)
    class(out) <- c("SGCBContrastList", "list")
    out
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
    top <- utils::head(x[order(x$pvalue), c("log2FoldChange", "pvalue", "padj")], 6)
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
