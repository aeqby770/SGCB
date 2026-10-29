# SGCB: conditional arithmetic-mean contrasts and zero-aware GG relative variability.
#' @useDynLib SGCB, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats p.adjust median plogis
NULL

#' Two-group mean and relative-variability analysis
#'
#' Mean contrasts test unconditional arithmetic means, including zero values.
#' Exact generalized-gamma likelihood calculations are reused by both endpoints.
#' The mean model uses pure generalized-gamma positive components and
#' nonzero-frequency uncertainty. Its log-location law uses a two-normal
#' approximation matching the first four cumulants; log-frequency uncertainty
#' enters through its conditional first two moments. Mean inference requires
#' finite first moments. A shared prior describes the control-group arithmetic
#' mean and paired GG shapes under both hypotheses. A global point-null fraction
#' and Gaussian nonzero effects are independent of that reference baseline.
#' The first group level defines the control; changing that reference changes
#' this modeling assumption. Posterior summaries condition on fitted population
#' parameters and the learned finite nuisance catalogue.
#' Strictly positive inputs retain the paired positive-law posterior. With observed zeros,
#' relative variability fits complete positive-law and nonzero-frequency
#' states within the finite-second-moment domain. The paired population prior
#' is symmetrized under group exchange. Sample-aligned influence values carry
#' positive-parameter and frequency covariance into each log-CV-squared contrast.
#' A Studentized reference matches degrees of freedom to variance-contribution
#' heterogeneity across sample columns within each group.
#' Gene-wise regularized GG fits remain descriptive distribution diagnostics.

#' @param counts Nonnegative gene-by-sample counts or normalized abundances.
#' @param group A factor or vector defining two independent groups.
#' @param alpha Benjamini-Hochberg threshold, default 0.05.
#' @param min_count Minimum value used by the expression filter.
#' @param min_samples Required number of samples meeting min_count.
#' @param normalization Either "tmm" for counts or "none" for normalized abundance.
#' @param offset Offset used only for descriptive log-variance and log-spread summaries.
#'   Zero values remain zero in the GG distribution model.
#' @param shape_prior_sd Prior standard deviation of Q for the descriptive gene-wise GG fit; default 1.
#' @param threads Number of OpenMP threads. Defaults to SGCB_THREADS when set,
#'   otherwise to all available cores.
#' @return An SGCBResults data frame. pvalue is cumulative posterior null mass
#'   ordered by posterior null probability and normalized by total null mass;
#'   padj is its complete-family BH adjustment. Under the population model this
#'   controls posterior expected false discovery proportion. log2FoldChange and
#'   lfcSE are the mean and standard deviation of the joint log2 mean-ratio
#'   posterior, including the zero-effect class; log2FC_shrunk is the same posterior
#'   mean. stat is their standardized ratio, not a t-test statistic.
#'   de_inference_available records whether both groups have positive data.
#'   Unavailable mean tests remain in the family with pvalue=1.
#'   Variance-moderation attributes remain descriptive log-variance diagnostics.
#'   dv_log2_cv2_ratio is the posterior mean log2 CV-squared contrast;
#'   dv_cv2_ctrl/treat are exponentiated posterior mean log CV-squared values.
#'   dv_pvalue contains a two-sided reference-tail bound with observed zeros
#'   and a paired-posterior tail bound for strictly positive inputs. dv_padj
#'   is its complete-family BH screening adjustment. dv_ci_lower/upper give
#'   the corresponding 95 percent interval. The observed-zero reference
#'   conditions on fitted state weights and uses a local influence approximation.
#'   dv_se_log2_cv2 combines fitted-state variation and the local uncertainty scale.
#'   DE_joint_q/DV_joint_q are exploratory mixed-score
#'   adjustments across both channels. All input genes passing the expression
#'   filter remain in each adjustment family, including unavailable results.
#'   ctrl_mu/sigma/Q, treat_mu/sigma/Q, fit_converged, fit_covariance_available,
#'   and shape-Q tests describe the separate gene-wise GG diagnostic.
#'   logspread_pvalue is a separate median-centered Brown-Forsythe summary.
#'   Population controls, convergence, and posterior summaries are attributes.
#'   When no gene has at least four positive observations in both groups,
#'   DV scores are unavailable; mean inference can still borrow nuisance information.
#' @export
#' @examples
#' \donttest{
#' source(system.file("examples", "minimal_workflow.R", package="SGCB"))
#' }
sgcbDE <- function(counts, group, alpha=0.05, min_count=10, min_samples=2,
                   normalization="tmm", offset=0.5, shape_prior_sd=1,
                   threads=as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))) {
    eps <- 1e-8
    counts <- as.matrix(counts)
    group <- factor(group)
    stopifnot(ncol(counts)==length(group),nlevels(group)==2,all(!is.na(group)),
              all(is.finite(counts)),all(counts>=0),
              is.finite(offset),offset>=0,shape_prior_sd>0)
    keep <- rowSums(counts>=min_count)>=min_samples
    x <- counts[keep,,drop=FALSE]
    gene_id <- rownames(x)
    if (is.null(gene_id)) gene_id <- paste0("gene_",seq_len(nrow(x)))
    n_genes <- nrow(x)
    stopifnot(n_genes>0)
    group_int <- as.integer(group)-1L
    ctrl <- which(group_int==0)
    treat <- which(group_int==1)
    n0 <- length(ctrl)
    n1 <- length(treat)
    stopifnot(n0>=2,n1>=2)
    time_normalization <- proc.time()[[3]]
    normalization <- match.arg(normalization,c("tmm","none"))
    sf <- rep(1,ncol(x))
    if (normalization=="tmm") {
        lib <- colSums(x)
        ref <- which.min(abs(lib-median(lib)))-1L
        tmm <- tmm_factors_cpp(x,as.integer(ref),0.3,0.05,eps,as.integer(threads))
        effective <- lib*tmm
        sf <- effective/exp(mean(log(effective)))
        sf[sf<0.01] <- 1
    }
    x <- sweep(x,2,sf,"/")
    time_mean <- proc.time()[[3]]

    # Moderated log-variances are descriptive; the DE endpoint is an arithmetic mean.
    logx <- log2(x+offset)
    stopifnot(all(is.finite(logx)))
    mean0 <- rowMeans(x[,ctrl,drop=FALSE])
    mean1 <- rowMeans(x[,treat,drop=FALSE])
    base_mean <- (mean0+mean1)/2
    var0 <- matrixStats::rowVars(logx[,ctrl,drop=FALSE])+eps
    var1 <- matrixStats::rowVars(logx[,treat,drop=FALSE])+eps
    pooled <- ((n0-1)*var0+(n1-1)*var1)/(n0+n1-2)
    moderated <- adaptive_squeeze_var_hierarchical_cpp(pooled,base_mean+offset+eps,n0+n1,FALSE,eps)
    kernel_cache <- new.env(parent=emptyenv())
    mean_population <- .populationMean(x,group_int,as.integer(threads),kernel_cache)
    mean_posterior <- mean_population$genes
    logfc <- mean_posterior$effect
    lfc_se <- mean_posterior$se
    statistic <- logfc/(lfc_se+eps)
    p_de <- mean_posterior$p
    time_spread <- proc.time()[[3]]
    p_spread <- .brownForsytheDV(x+offset,ctrl,treat,n0,n1,eps)
    time_distribution <- proc.time()[[3]]

    native_fit <- hurdle_gg_cpp(x,group_int,1/shape_prior_sd^2,as.integer(threads))
    fit <- as.data.frame(native_fit)
    converged <- as.logical(fit$converged0 & fit$converged1)
    covariance <- as.logical(fit$covariance0 & fit$covariance1)
    variance_q <- fit$shape_var0_model+fit$shape_var1_model
    df_q <- variance_q^2/(fit$shape_var0_model^2/(fit$n_positive0-3)+
                         fit$shape_var1_model^2/(fit$n_positive1-3))
    p_q <- 2*stats::pt(-abs(fit$Q1-fit$Q0)/sqrt(variance_q),df_q)
    shape_available <- converged & covariance & is.finite(p_q)
    p_q[!shape_available] <- NA_real_
    time_population <- proc.time()[[3]]
    population <- .populationGG(x,group_int,as.integer(threads),alpha,kernel_cache=kernel_cache)
    posterior <- population$genes
    effect <- posterior$effect
    effect_se <- posterior$se
    p_dv <- posterior$p
    available <- posterior$available==1
    p_dv[!available] <- NA_real_
    moments <- is.finite(posterior$log_cv2_0_mean) & is.finite(posterior$log_cv2_1_mean)
    df_dv <- rep(NA_real_,n_genes)
    time_post <- proc.time()[[3]]
    q_de <- p.adjust(p_de,"BH",n=n_genes)
    q_dv <- p.adjust(p_dv,"BH",n=n_genes)
    q_spread <- p.adjust(p_spread,"BH",n=n_genes)
    joint <- c(p_de,p_dv)
    joint[!is.finite(joint)] <- 1
    joint <- matrix(p.adjust(joint,"BH"),nrow=n_genes)
    de_call <- mean_posterior$available==1 & q_de<alpha
    dv_call <- !is.na(q_dv) & q_dv<alpha
    result <- data.frame(gene_id=gene_id,baseMean=base_mean,
        log2FoldChange=logfc,lfcSE=lfc_se,stat=statistic,pvalue=p_de,padj=q_de,
        log2FC_shrunk=logfc,
        mean_ctrl=mean0,mean_treat=mean1,
        ctrl_mu=fit$mu0,ctrl_sigma=fit$sigma0,ctrl_Q=fit$Q0,
        treat_mu=fit$mu1,treat_sigma=fit$sigma1,treat_Q=fit$Q1,
        ctrl_nonzero_fraction=fit$positive0,treat_nonzero_fraction=fit$positive1,
        ctrl_positive_n=fit$n_positive0,treat_positive_n=fit$n_positive1,
        fitted_mean_ctrl=fit$fitted_mean0,fitted_mean_treat=fit$fitted_mean1,
        dv_cv2_ctrl=exp(posterior$log_cv2_0_mean),dv_cv2_treat=exp(posterior$log_cv2_1_mean),
        dv_log2_cv2_ratio=effect,dv_se_log2_cv2=effect_se,dv_df=df_dv,
        dv_ci_lower=posterior$lower,dv_ci_upper=posterior$upper,
        dv_pvalue=p_dv,dv_padj=q_dv,
        logspread_pvalue=p_spread,logspread_padj=q_spread,
        positive_log_variance_ctrl=exp(fit$logvar_positive0),
        positive_log_variance_treat=exp(fit$logvar_positive1),
        shape_Q_difference=fit$Q1-fit$Q0,shape_Q_se=sqrt(variance_q),
        shape_Q_pvalue=p_q,shape_Q_padj=p.adjust(p_q,"BH",n=n_genes),
        fit_converged=converged,fit_covariance_available=covariance,
        de_inference_available=mean_posterior$available==1,
        dv_moments_finite=moments,dv_inference_available=available,
        shape_inference_available=shape_available,
        DE_joint_q=joint[,1],DV_joint_q=joint[,2],DE=de_call,DV=dv_call,
        channel=ifelse(de_call & dv_call,"DE_and_DV",ifelse(de_call,"DE_only",
                       ifelse(dv_call,"DV_only","neither"))),
        row.names=gene_id,stringsAsFactors=FALSE)
    attr(result,"alpha") <- alpha
    attr(result,"normalization") <- normalization
    attr(result,"size_factors") <- sf
    attr(result,"offset") <- offset
    attr(result,"shape_prior_sd") <- shape_prior_sd
    attr(result,"de_endpoint") <- mean_population$endpoint
    attr(result,"de_score_type") <- mean_population$score_type
    attr(result,"mean_model") <- mean_population[c("atoms","pairs","control","convergence",
        "A_grid","location_prior","effect_profile","posterior_convergence","rng_state")]
    attr(result,"mean_posterior") <- mean_posterior
    attr(result,"mean_timing") <- mean_population$timing
    attr(result,"mean_integration_checks") <- mean_population$integration_checks
    attr(result,"dv_covariance") <- population$covariance
    attr(result,"n_ctrl") <- n0
    attr(result,"n_treat") <- n1
    attr(result,"n_filtered") <- sum(!keep)
    attr(result,"df_prior") <- moderated$df_prior
    attr(result,"var_shrunk") <- moderated$var_shrunk
    attr(result,"var_prior_trend") <- moderated$var_prior_trend
    attr(result,"fit_diagnostics") <- list(
        ctrl=data.frame(converged=as.logical(fit$converged0),covariance=as.logical(fit$covariance0),
            iterations=fit$iterations0,gradient=fit$score0),
        treat=data.frame(converged=as.logical(fit$converged1),covariance=as.logical(fit$covariance1),
            iterations=fit$iterations1,gradient=fit$score1))
    attr(result,"population_model") <- population[c("atoms","positive_atoms","pairs","reference","control","convergence")]
    attr(result,"population_posterior") <- posterior
    attr(result,"population_timing") <- population$timing
    attr(result,"dv_score_type") <- population$score_type
    attr(result,"distribution_fits") <- fit
    attr(result,"native_source_id") <- attr(native_fit,"native_source_id")
    attr(result,"fit_started") <- as.character(Sys.time()- (proc.time()[[3]]-time_normalization))
    attr(result,"run_id") <- Sys.getenv("SGCB_RUN_ID")
    attr(result,"elapsed_components") <- c(normalization=time_mean-time_normalization,
        arithmetic_mean=time_spread-time_mean,logspread=time_distribution-time_spread,
        gg_distribution=time_population-time_distribution,
        population_inference=time_post-time_population,
        postprocessing=proc.time()[[3]]-time_post)
    attr(result,"method") <- "population_zero_aware_generalized_gamma"
    attr(result,"primary_de_channel") <- "padj"
    attr(result,"primary_dv_channel") <- "dv_padj"
    class(result) <- c("SGCBResults","data.frame")
    result
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
#' @param ... Additional arguments forwarded to \code{sgcbDE()}.
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
                         sample_subset = NULL, alpha = 0.05,
                         min_count = 10, min_samples = 2, ...) {
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
        min_count = min_count,
        min_samples = min_samples,
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
#' @param ... Additional arguments forwarded to \code{sgcbDE()}.
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
                         alpha = 0.05, min_count = 10, min_samples = 2, ...) {
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
            min_count = min_count,
            min_samples = min_samples,
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
    cat("SGCB mean and relative-variability results\n")
    cat("-------------------------------------\n")
    cat("Genes tested:", nrow(x), "\n")
    cat("Mean discoveries (BH <", alpha, "):", n_sig, "\n")
    cat("CV-squared discoveries (BH <", alpha, "):",sum(x$DV,na.rm=TRUE),"\n")
    cat("CV-squared inference available:",sum(x$dv_inference_available),"\n")
    df_prior <- attr(x, "df_prior")
    if (!is.null(df_prior)) cat("Prior df:", round(df_prior, 1), "\n")
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
significantGenes <- function(x, padj_cutoff = 0.05, lfc_cutoff = 0) {
    stopifnot(inherits(x, "SGCBResults"))
    sig <- x[!is.na(x$padj) & x$padj < padj_cutoff & abs(x$log2FoldChange) > lfc_cutoff, ]
    sig <- sig[order(sig$pvalue), ]
    class(sig) <- "data.frame"
    sig
}
