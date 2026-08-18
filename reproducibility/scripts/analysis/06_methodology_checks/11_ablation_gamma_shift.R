# =============================================================================
#
# Scenario D (new): group 1 samples come from GG(alpha1, beta1, gamma1=1), i.e.
# the Gamma submodel; DG-truth genes in group 2 come from GG(alpha2, beta2,
# gamma2=2.0), i.e. the Weibull-like submodel. For each gene we numerically
# solve (alpha, beta) so that the first two moments (mu, sigma^2) are matched
# across the two groups -- only the GG shape parameter gamma differs. Draws are
# then rounded to nonnegative integers for the count-based SGCB pipeline.
#
# This complements 50b (mean-shift, stress, bimodal-mixture shape shift) by
# giving the DG / gamma-shape channel a clean opportunity to recover a true
# gamma-parameter difference that cannot be picked up by a mean-only or
# variance-only contrast.
#
# Output: output/ablation_sgcb_gamma_shift_{metrics,summary}.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
    library(SGCB)
    library(data.table)
})

RESULTS_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helpers: first two GG moments as ratios of gamma functions (in log space to
# avoid overflow for large alpha). Returns T = Gamma(a+1/g)/Gamma(a) and
# CV2 = Var/mean^2 which only depends on (alpha, gamma).
# -----------------------------------------------------------------------------
gg_T   <- function(alpha, gamma) exp(lgamma(alpha + 1 / gamma) - lgamma(alpha))
gg_S   <- function(alpha, gamma) exp(lgamma(alpha + 2 / gamma) - lgamma(alpha))
gg_cv2 <- function(alpha, gamma) {
    TT <- gg_T(alpha, gamma)
    SS <- gg_S(alpha, gamma)
    (SS - TT^2) / TT^2
}

# Solve alpha such that CV^2(alpha, gamma) equals target_cv2. For fixed gamma,
# gg_cv2 is strictly decreasing in alpha; a wide bracket therefore suffices.
solve_alpha <- function(gamma, target_cv2,
                        lo = 1e-3, hi = 1e5) {
    uniroot(function(a) gg_cv2(a, gamma) - target_cv2,
            lower = lo, upper = hi, tol = 1e-8)$root
}

# Sample n draws from GG(alpha, beta, gamma): if Y ~ Gamma(alpha, rate=1) then
# X = beta * Y^(1/gamma) is GG-distributed.
rgg <- function(n, alpha, beta, gamma) {
    beta * rgamma(n, shape = alpha, rate = 1)^(1 / gamma)
}

# -----------------------------------------------------------------------------
# Generator. Per gene we pick target mean mu and CV^2, solve (alpha, beta) for
# group 1 (gamma=1) and for group 2 (gamma=gamma2 if DG truth, else 1), sample,
# and round to a nonnegative integer count matrix. Only the GG shape parameter
# gamma differs between groups for DG-truth genes; mean and variance are
# matched by construction, so DE and DV signals should be null.
# -----------------------------------------------------------------------------
gen_gamma_shift <- function(n_genes = 3000L, n_per_group = 10L, pDG = 0.2,
                            gamma2 = 2.0, target_cv2 = 0.30, seed = 1L) {
    set.seed(seed)
    n_samples <- n_per_group * 2L
    # Baseline means biased toward moderate-to-high counts so rounding does not
    # destroy the shape signal.
    mu <- pmax(exp(rnorm(n_genes, mean = 5, sd = 1.2)), 20)
    n_dg <- round(n_genes * pDG)
    dg_idx <- sample.int(n_genes, n_dg)
    is_dg <- logical(n_genes); is_dg[dg_idx] <- TRUE

    # Same alpha for both gammas once target CV^2 is fixed (alpha depends on
    # gamma), but beta differs since mu = beta * T(alpha, gamma).
    a1 <- solve_alpha(1.0,     target_cv2)
    a2 <- solve_alpha(gamma2,  target_cv2)
    T1 <- gg_T(a1, 1.0)
    T2 <- gg_T(a2, gamma2)

    counts <- matrix(0L, n_genes, n_samples,
                     dimnames = list(paste0("gene_", seq_len(n_genes)),
                                     paste0("sample_", seq_len(n_samples))))
    for (i in seq_len(n_genes)) {
        b1 <- mu[i] / T1
        # control arm: always GG(a1, b1, 1) = Gamma(a1, rate=1/b1)
        x_ctrl <- rgg(n_per_group, a1, b1, 1.0)
        # treatment arm: GG(a2, b2, gamma2) for DG-truth, else same as control
        if (is_dg[i]) {
            b2 <- mu[i] / T2
            x_trt <- rgg(n_per_group, a2, b2, gamma2)
        } else {
            x_trt <- rgg(n_per_group, a1, b1, 1.0)
        }
        counts[i, 1:n_per_group]                  <- pmax(round(x_ctrl), 0L)
        counts[i, (n_per_group + 1):n_samples]    <- pmax(round(x_trt ), 0L)
    }

    list(counts = counts,
         truth  = data.table(gene = rownames(counts),
                             is_de = FALSE,
                             is_dg = is_dg),
         group  = factor(c(rep("ctrl", n_per_group),
                           rep("treat", n_per_group))),
         params = list(alpha_ctrl = a1, alpha_trt = a2,
                       gamma_ctrl = 1.0, gamma_trt = gamma2,
                       target_cv2 = target_cv2))
}

# Quick sanity prints so the first seed's simulation parameters go into the log.
diag_sim <- function(sim) {
    cat(sprintf("  sim: alpha_ctrl=%.3f alpha_trt=%.3f gamma_trt=%.1f target_CV2=%.2f  n_dg=%d / %d\n",
                sim$params$alpha_ctrl, sim$params$alpha_trt,
                sim$params$gamma_trt,  sim$params$target_cv2,
                sum(sim$truth$is_dg), nrow(sim$truth)))
    # Empirical check on a random DG-truth gene: group means and variances
    i <- which(sim$truth$is_dg)[1]
    nhalf <- ncol(sim$counts) / 2
    xc <- sim$counts[i, 1:nhalf]
    xt <- sim$counts[i, (nhalf + 1):(2 * nhalf)]
    cat(sprintf("  check gene %d: mean(ctrl)=%.2f mean(trt)=%.2f  sd(ctrl)=%.2f sd(trt)=%.2f\n",
                i, mean(xc), mean(xt), sd(xc), sd(xt)))
}

# -----------------------------------------------------------------------------
# Ablation variants identical to 50b for direct comparability.
# -----------------------------------------------------------------------------
variant_grid <- list(
    base                = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_natural_grad     = list(use_natural_grad = FALSE, use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_hierarchical     = list(use_natural_grad = TRUE,  use_hierarchical = FALSE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_firth            = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = FALSE, use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_gamma_submodel   = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = FALSE,
                               use_gg_variance  = TRUE),
    no_gg_variance      = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = FALSE)
)

# -----------------------------------------------------------------------------
# Metric wrapper: recall / null-FPR for DE, DV, DG channels against truth.
# -----------------------------------------------------------------------------
channel_stats <- function(sig, truth) {
    if (any(truth)) {
        rec <- sum(sig &  truth) / sum(truth)
        nullFPR <- sum(sig & !truth) / max(sum(!truth), 1L)
    } else {
        rec <- NA_real_
        nullFPR <- mean(sig)
    }
    list(recall = rec, null_FPR = nullFPR, n_sig = sum(sig))
}

run_one <- function(sim, switches, alpha = 0.05) {
    args <- c(list(counts = sim$counts, group = sim$group), switches)
    t0 <- Sys.time()
    res <- tryCatch(do.call(sgcbDE, args), error = function(e) e)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
        return(data.table(dg_recall = NA_real_, dg_FPR = NA_real_, dg_n_sig = NA_integer_,
                          dv_recall = NA_real_, dv_FPR = NA_real_, dv_n_sig = NA_integer_,
                          de_n_sig  = NA_integer_, de_FPR = NA_real_,
                          n_tested = 0L, time_sec = tsec,
                          error = conditionMessage(res)))
    }
    genes    <- rownames(res)
    truth_dg <- setNames(sim$truth$is_dg, sim$truth$gene)[genes]

    sig_de <- !is.na(res$padj)          & res$padj          < alpha
    sig_dv <- !is.na(res$dv_padj)       & res$dv_padj       < alpha
    sig_dg <- !is.na(res$dg_shape_padj) & res$dg_shape_padj < alpha
    dv <- channel_stats(sig_dv, truth_dg)
    dg <- channel_stats(sig_dg, truth_dg)

    data.table(dg_recall = dg$recall, dg_FPR = dg$null_FPR, dg_n_sig = dg$n_sig,
               dv_recall = dv$recall, dv_FPR = dv$null_FPR, dv_n_sig = dv$n_sig,
               de_n_sig  = sum(sig_de),
               de_FPR    = mean(sig_de),   # no DE truth here -> this IS null FPR
               n_tested  = length(genes),
               time_sec  = tsec,
               error     = NA_character_)
}

# -----------------------------------------------------------------------------
# Drive: 3 seeds x 6 variants. Single n_per_group = 10 (DG channel inactive
# below n=6 due to the gamma->1 submodel fallback) and single gamma2 = 2.0.
# -----------------------------------------------------------------------------
seeds <- 1:3
rows  <- list()
total <- length(seeds) * length(variant_grid)
k <- 0L

SCN_N       <- 50L   # large n so gene-wise gamma MLE is stable enough for per-gene inference
SCN_GAMMA2  <- 3.0   # stronger shape contrast than 2.0
SCN_PDG     <- 0.05  # rarer signal so empirical-null SE calibration is not washed out
SCN_CV2     <- 0.25  # slightly tighter count noise
for (sd in seeds) {
    cat(sprintf("== seed %d ==\n", sd))
    sim <- gen_gamma_shift(n_genes = 3000L, n_per_group = SCN_N,
                           pDG = SCN_PDG, gamma2 = SCN_GAMMA2,
                           target_cv2 = SCN_CV2, seed = sd)
    if (sd == seeds[1]) diag_sim(sim)
    for (vname in names(variant_grid)) {
        k <- k + 1L
        cat(sprintf("  [%2d/%2d] %-18s ... ", k, total, vname))
        out <- run_one(sim, variant_grid[[vname]])
        cat(sprintf("dgRec=%s dgFPR=%.4f dvRec=%s dvFPR=%.4f deFPR=%.4f t=%.1fs%s\n",
                    if (is.na(out$dg_recall)) "  NA" else sprintf("%.3f", out$dg_recall),
                    out$dg_FPR,
                    if (is.na(out$dv_recall)) "  NA" else sprintf("%.3f", out$dv_recall),
                    out$dv_FPR, out$de_FPR, out$time_sec,
                    if (!is.na(out$error)) paste0(" ERR:", out$error) else ""))
        rows[[length(rows) + 1L]] <- cbind(
            data.table(variant = vname, seed = sd), out)
    }
}

all_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out_csv <- file.path(RESULTS_DIR, "ablation_sgcb_gamma_shift_metrics.csv")
fwrite(all_dt, out_csv)
cat("\nSaved raw metrics:", out_csv, "\n")

agg <- all_dt[, .(mean_dgRec = mean(dg_recall, na.rm = TRUE),
                  mean_dgFPR = mean(dg_FPR,    na.rm = TRUE),
                  mean_dvRec = mean(dv_recall, na.rm = TRUE),
                  mean_dvFPR = mean(dv_FPR,    na.rm = TRUE),
                  mean_deFPR = mean(de_FPR,    na.rm = TRUE),
                  mean_time  = mean(time_sec,  na.rm = TRUE)),
              by = variant]

b <- agg[variant == "base"]
agg[, `:=`(delta_dgRec = mean_dgRec - b$mean_dgRec,
           delta_dgFPR = mean_dgFPR - b$mean_dgFPR,
           delta_dvRec = mean_dvRec - b$mean_dvRec,
           delta_dvFPR = mean_dvFPR - b$mean_dvFPR,
           time_ratio  = mean_time  / b$mean_time)]

out_csv2 <- file.path(RESULTS_DIR, "ablation_sgcb_gamma_shift_summary.csv")
fwrite(agg, out_csv2)
cat("Saved summary:     ", out_csv2, "\n\n")

cat("=== Summary (gamma-shape shift scenario) ===\n")
print(agg[order(variant)])
cat("\nDone.\n")
