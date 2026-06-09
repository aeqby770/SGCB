#' GG Distributional Regression
#'
#' Fits a Generalized Gamma GLM with design matrix support for differential
#' expression (DE) and optionally differential variability (DV) analysis.
#'
#' V1 (mean-only): \code{log(beta_gi) = x_i' b_g}, alpha/gamma gene-specific.
#' V2 (mean + dispersion): additionally fits \code{log(alpha_gi) = w_i' d_g}
#' and tests DV via likelihood ratio.
#'
#' @param counts Integer count matrix (genes x samples), with rownames as gene IDs.
#' @param design Design matrix (samples x p). If NULL, constructed from \code{group}.
#' @param group Factor or character vector of group labels.
#' @param contrast Numeric contrast vector of length p for DE test.
#' @param design_disp Dispersion design matrix (samples x q) for DV test.
#'   If NULL, only V1 (DE) is run. If \code{"auto"}, uses the same design as
#'   \code{design}.
#' @param contrast_disp Dispersion contrast vector of length q for DV output.
#' @param alpha Significance level for \code{padj}/\code{padj_dv}.
#' @param min_count Minimum count threshold for gene filtering.
#' @param min_samples Minimum number of samples exceeding \code{min_count}.
#' @param prior_de Prior inclusion probability for DE channel in joint posterior;
#'   numeric in (0,1) or \code{"auto"}.
#' @param prior_dv Prior inclusion probability for DV channel in joint posterior;
#'   numeric in (0,1) or \code{"auto"}.
#' @param prior_sd_de Prior SD for DE effect (log2 scale) in ABF.
#' @param prior_sd_dv Prior SD for DV effect (log2 scale) in ABF.
#' @param ... Additional arguments forwarded to the C++ engines.
#'
#' @return A \code{data.frame} (class \code{SGCBDRegResults}) with DE columns
#'   and optionally DV columns (\code{pvalue_dv}, \code{padj_dv},
#'   \code{dv_log2ratio}, \code{lr_stat_dv}).
#'
#' @export
sgcbDReg <- function(counts, design = NULL, group = NULL, contrast = NULL,
                     design_disp = NULL, contrast_disp = NULL,
                     alpha = 0.05, min_count = 10, min_samples = 2,
                     prior_de = "auto", prior_dv = "auto",
                     prior_sd_de = 1.0, prior_sd_dv = 0.5, ...) {

    if (is.null(design) && is.null(group))
        stop("Either 'design' or 'group' must be provided.")

    if (is.null(design)) {
        grp <- factor(group)
        design <- stats::model.matrix(~ grp)
    }

    stopifnot(ncol(counts) == nrow(design))

    P <- ncol(design)
    if (is.null(contrast)) {
        contrast <- rep(0, P)
        contrast[P] <- 1
    }
    stopifnot(length(contrast) == P)

    # ---- Gene filtering ----
    keep <- rowSums(counts >= min_count) >= min_samples
    y_filt <- counts[keep, , drop = FALSE]
    gene_ids <- rownames(y_filt)
    if (is.null(gene_ids)) gene_ids <- paste0("gene_", seq_len(nrow(y_filt)))

    # ---- TMM normalization (edgeR-style, Robinson & Oshlack 2010) ----
    lib_raw <- colSums(y_filt)
    ref_idx <- which.min(abs(lib_raw - stats::median(lib_raw))) - 1L
    tmm <- tmm_factors_cpp(y_filt, as.integer(ref_idx),
                           trim_m = 0.3, trim_a = 0.05,
                           eps = 1e-8, n_threads = 0L)
    eff <- lib_raw * tmm
    eff <- eff / exp(mean(log(eff)))
    y_norm <- sweep(y_filt, 2, eff, "/")

    # Continuous relaxation for zero counts under continuous GG likelihood
    # (consistent with sgcbDE; GG density requires x > 0)
    y_norm <- y_norm + 0.5

    # ---- V1 Fit (DE) ----
    dots <- list(...)
    v1_args <- dots[intersect(names(dots),
        c("max_outer", "max_inner_ag", "ag_lr", "tol", "eps"))]
    res <- do.call(sgcb_dreg_v1_cpp,
                   c(list(y = y_norm, design = design,
                          contrast = as.numeric(contrast)), v1_args))

    # ---- Output (DE) ----
    out <- data.frame(
        gene_id            = gene_ids,
        baseMean           = res$baseMean,
        log2FoldChange     = res$log2FoldChange,
        lfcSE              = res$lfcSE,
        stat               = res$stat,
        pvalue             = res$pvalue,
        padj               = stats::p.adjust(res$pvalue, method = "BH"),
        alpha              = res$alpha,
        gamma              = res$gamma,
        dispersion         = res$dispersion,
        dispersion_shrunk  = res$dispersion_shrunk,
        stringsAsFactors   = FALSE
    )
    out$de_fdr_call <- out$padj <= alpha

    de_cal <- .sgcbCalibratedSE(out$log2FoldChange, out$lfcSE)
    bf_de_abf <- .sgcbABF1(
        effect = out$log2FoldChange,
        se = de_cal$se_cal,
        prior_sd = prior_sd_de
    )
    prior_de_used <- .sgcbResolvePrior(
        prior_de,
        p = out$pvalue,
        bf = bf_de_abf,
        prior_alpha = 1.0,
        prior_beta = 9.0
    )
    de_post <- .sgcbModelPosterior1(
        p_de = out$pvalue,
        prior_de = prior_de_used,
        bf_de = bf_de_abf
    )
    de_call <- .sgcbPosteriorCalls(de_post$p_de_post, target_fdr = alpha)
    out$p_de_post <- de_post$p_de_post
    out$de_post_call <- de_call$call
    out$p_de_only_post <- de_post$p_de_only_post
    out$bf_de <- de_post$bf_de
    out$model_prob_null <- de_post$model_prob[, "model_prob_null"]
    out$model_prob_de <- de_post$model_prob[, "model_prob_de"]

    # ---- V2 Fit (DV) ----
    run_dv <- !is.null(design_disp)
    if (identical(design_disp, "auto")) {
        design_disp <- design
        run_dv <- TRUE
    }
    if (run_dv) {
        Q <- ncol(design_disp)
        if (is.null(contrast_disp)) {
            contrast_disp <- rep(0, Q)
            contrast_disp[Q] <- 1
        }
        dv_args <- dots[intersect(names(dots),
            c("max_outer_dv", "max_inner_d", "d_lr", "eps"))]
        dv <- do.call(sgcb_dreg_dv_cpp,
                      c(list(y = y_norm, design_mean = design,
                             design_disp = design_disp,
                             contrast_disp = as.numeric(contrast_disp),
                             coef_v1 = res$coefficients,
                             alpha_v1 = res$alpha,
                             gamma_v1 = res$gamma), dv_args))
        out$pvalue_dv    <- dv$pvalue_dv
        out$padj_dv      <- stats::p.adjust(dv$pvalue_dv, method = "BH")
        out$dv_fdr_call  <- out$padj_dv <= alpha
        out$dv_log2ratio <- dv$dv_log2ratio
        out$dv_se_log2   <- dv$dv_se_log2
        out$lr_stat_dv   <- dv$lr_stat_dv
        out$ll_null_dv   <- dv$ll_null
        out$ll_full_dv   <- dv$ll_full

        dv_cal <- .sgcbCalibratedSE(out$dv_log2ratio, out$dv_se_log2)
        bf_dv_abf <- .sgcbABF1(
            effect = out$dv_log2ratio,
            se = dv_cal$se_cal,
            prior_sd = prior_sd_dv
        )
        bf_dv_lr <- .sgcbLRBF(out$lr_stat_dv, ncol(y_norm), ncol(design_disp) - 1)
        coef_test_dv <- dv$d_coefficients[, seq.int(2, ncol(design_disp)), drop = FALSE] / log(2)
        bf_dv_laplace <- .sgcbLaplaceBF(
            ll_full = dv$ll_full,
            ll_null = dv$ll_null,
            logdet_full = dv$logdet_full_d,
            logdet_null = dv$logdet_null_d,
            df_test = pmax(ncol(design_disp) - 1, 1)
        )
        bf_dv_laplace_prior <- .sgcbLaplaceBF(
            ll_full = dv$ll_full,
            ll_null = dv$ll_null,
            logdet_full = dv$logdet_full_d,
            logdet_null = dv$logdet_null_d,
            df_test = pmax(ncol(design_disp) - 1, 1),
            coef_test = coef_test_dv,
            prior_sd = rep(prior_sd_dv, ncol(coef_test_dv))
        )
        n_eff_dv <- ncol(y_norm) / pmax(ncol(design_disp), 1)
        prior_dv_beta <- ifelse(n_eff_dv <= 6, 999.0,
                                ifelse(n_eff_dv <= 10, 199.0, 19.0))
        prior_dv_upper <- ifelse(n_eff_dv <= 10, 0.10, 0.25)
        prior_dv_used <- .sgcbResolvePrior(
            prior_dv,
            p = out$pvalue_dv,
            bf = bf_dv_laplace,
            lower = 0.001,
            upper = prior_dv_upper,
            prior_alpha = 1.0,
            prior_beta = prior_dv_beta
        )

        joint_post <- .sgcbModelPosterior2(
            p_de = out$pvalue,
            p_dv = out$pvalue_dv,
            prior_de = prior_de_used,
            prior_dv = prior_dv_used,
            bf_de = bf_de_abf,
            bf_dv = pmax(bf_dv_laplace, 1e-12)
        )
        de_call <- .sgcbPosteriorCalls(joint_post$p_de_post, target_fdr = alpha)
        dv_call <- .sgcbPosteriorCalls(joint_post$p_dv_post, target_fdr = alpha)
        out$p_de_post <- joint_post$p_de_post
        out$de_post_call <- de_call$call
        out$p_dv_post <- joint_post$p_dv_post
        out$dv_post_call <- dv_call$call
        out$p_de_only_post <- joint_post$p_de_only_post
        out$p_dv_only_post <- joint_post$p_dv_only_post
        out$p_de_dv_post <- joint_post$p_de_dv_post
        out$bf_de <- joint_post$bf_de
        out$bf_dv <- joint_post$bf_dv
        out$bf_dv_lr <- bf_dv_lr
        out$bf_dv_laplace <- bf_dv_laplace
        out$bf_dv_laplace_prior <- bf_dv_laplace_prior
        out$bf_dv_abf <- bf_dv_abf
        out$model_prob_null <- joint_post$model_prob[, "model_prob_null"]
        out$model_prob_de <- joint_post$model_prob[, "model_prob_de"]
        out$model_prob_dv <- joint_post$model_prob[, "model_prob_dv"]
        out$model_prob_de_dv <- joint_post$model_prob[, "model_prob_de_dv"]
    }

    rownames(out) <- gene_ids
    attr(out, "design")      <- design
    attr(out, "contrast")    <- contrast
    attr(out, "df_prior")    <- res$df_prior
    attr(out, "var_prior")   <- res$var_prior
    attr(out, "df_total")    <- res$df_total
    attr(out, "df_residual") <- res$df_residual
    attr(out, "n_outer")     <- res$n_outer
    attr(out, "loglik")      <- res$loglik
    attr(out, "prior_de")    <- prior_de
    attr(out, "prior_dv")    <- prior_dv
    attr(out, "prior_de_used") <- prior_de_used
    if (run_dv) attr(out, "prior_dv_used") <- prior_dv_used
    if (run_dv) attr(out, "prior_dv_upper") <- prior_dv_upper
    attr(out, "prior_sd_de") <- prior_sd_de
    attr(out, "prior_sd_dv") <- prior_sd_dv
    attr(out, "abf_inflation_de") <- de_cal$inflation
    if (run_dv) attr(out, "abf_inflation_dv") <- dv_cal$inflation
    attr(out, "post_call_fdr") <- alpha
    attr(out, "post_de_threshold") <- de_call$threshold
    attr(out, "post_de_realized_fdr") <- de_call$realized_fdr
    if (run_dv) attr(out, "post_dv_threshold") <- dv_call$threshold
    if (run_dv) attr(out, "post_dv_realized_fdr") <- dv_call$realized_fdr
    attr(out, "decision_mode") <- "frequentist_primary"
    attr(out, "primary_de_channel") <- "padj"
    if (run_dv) attr(out, "primary_dv_channel") <- "padj_dv"
    attr(out, "method")      <- if (run_dv) "GG-DReg-V2" else "GG-DReg-V1"
    class(out) <- c("SGCBDRegResults", "data.frame")
    out
}

#' @export
print.SGCBDRegResults <- function(x, ...) {
    n <- nrow(x)
    n_de <- sum(x$padj < 0.05, na.rm = TRUE)
    has_dv <- "padj_dv" %in% names(x)
    n_dv <- if (has_dv) sum(x$padj_dv < 0.05, na.rm = TRUE) else 0L
    cat(sprintf("SGCBDRegResults: %d genes, %d DE (padj<0.05)", n, n_de))
    if (has_dv) cat(sprintf(", %d DV (padj_dv<0.05)", n_dv))
    cat(sprintf("\nMethod: %s | df_total: %.1f\n",
                attr(x, "method"), attr(x, "df_total")))
    invisible(x)
}

#' GG-DReg for proteomics data
#'
#' Thin wrapper around \code{\link{sgcbDReg}} with defaults tuned for
#' quantitative proteomics (e.g. TMT, LFQ intensities). Uses lower
#' \code{min_count} since protein abundances are typically smaller than
#' RNA-seq counts.
#'
#' @inheritParams sgcbDReg
#' @export
sgcbProtein <- function(counts, design = NULL, group = NULL, contrast = NULL,
                        design_disp = NULL, contrast_disp = NULL,
                        min_count = 1, min_samples = 2, ...) {
    sgcbDReg(counts, design = design, group = group, contrast = contrast,
             design_disp = design_disp, contrast_disp = contrast_disp,
             min_count = min_count, min_samples = min_samples, ...)
}

#' GG-DReg for pseudo-bulk single-cell data
#'
#' Thin wrapper around \code{\link{sgcbDReg}} with defaults tuned for
#' pseudo-bulk aggregated counts from single-cell RNA-seq.
#'
#' @inheritParams sgcbDReg
#' @export
sgcbPseudobulk <- function(counts, design = NULL, group = NULL,
                           contrast = NULL, design_disp = NULL,
                           contrast_disp = NULL,
                           min_count = 10, min_samples = 3, ...) {
    sgcbDReg(counts, design = design, group = group, contrast = contrast,
             design_disp = design_disp, contrast_disp = contrast_disp,
             min_count = min_count, min_samples = min_samples, ...)
}
