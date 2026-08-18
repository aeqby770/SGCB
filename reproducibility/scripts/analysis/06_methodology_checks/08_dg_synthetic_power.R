#!/usr/bin/env Rscript
# =============================================================================
# Q2c: DG synthetic ground-truth validation
#
# then test whether DG channel detects them.
#
# Design:
#   - Sample sizes per group: 5, 10, 20, 50
#   - Metrics: sensitivity, specificity, FPR (null genes), AUC
# =============================================================================
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})
devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

set.seed(2026)

N_GENES   <- 10000
N_DE      <- 500
ALPHA_BASE <- 2.0
BETA_BASE  <- 100.0
GAMMA_BASE <- 1.0

gamma_fcs <- c(1.5, 2.0, 3.0)
sample_sizes <- c(5L, 10L, 20L, 50L)

results <- list()

for (n_per_group in sample_sizes) {
  for (gfc in gamma_fcs) {
    cat(sprintf("\n=== n=%d, gamma_fc=%.1f ===\n", n_per_group, gfc))

    alpha_vec <- rep(ALPHA_BASE, N_GENES)
    beta_vec  <- rep(BETA_BASE, N_GENES)
    gamma_ctrl_vec <- rep(GAMMA_BASE, N_GENES)
    gamma_treat_vec <- rep(GAMMA_BASE, N_GENES)

    de_idx <- sample(N_GENES, N_DE)
    gamma_treat_vec[de_idx] <- GAMMA_BASE * gfc

    # Generate counts via gg_sample_mat (returns continuous GG values)
    ctrl_mat <- gg_sample_mat(N_GENES, n_per_group, alpha_vec, beta_vec, gamma_ctrl_vec)
    treat_mat <- gg_sample_mat(N_GENES, n_per_group, alpha_vec, beta_vec, gamma_treat_vec)

    # Round to integer counts (SGCB expects count-like data)
    counts <- cbind(round(pmax(ctrl_mat, 0)), round(pmax(treat_mat, 0)))
    storage.mode(counts) <- "integer"
    group <- factor(rep(c("ctrl", "treat"), each = n_per_group))

    # Run SGCB
    res <- tryCatch(
      sgcbDE(counts, group),
      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(res)) next

    df <- as.data.frame(res)
    truth <- rep(FALSE, N_GENES)
    truth[de_idx] <- TRUE

    # Extract DG p-values
    for (ch in c("dg_shape_pvalue", "dg_gamma_pvalue", "dg_alpha_pvalue")) {
      pv <- df[[ch]]
      if (is.null(pv) || all(is.na(pv))) next

      valid <- !is.na(pv)
      pv <- pv[valid]
      tr <- truth[valid]

      # FPR on null genes
      fpr_05 <- mean(pv[!tr] < 0.05, na.rm = TRUE)
      # Sensitivity on DE genes
      sens_05 <- mean(pv[tr] < 0.05, na.rm = TRUE)
      # Specificity
      spec_05 <- 1 - fpr_05

      # AUC (simple Mann-Whitney)
      if (sum(tr) > 0 && sum(!tr) > 0) {
        auc <- tryCatch({
          r <- rank(pv)
          n1 <- sum(tr); n0 <- sum(!tr)
          u <- sum(r[tr]) - n1 * (n1 + 1) / 2
          1 - u / (n1 * n0)
        }, error = function(e) NA_real_)
      } else {
        auc <- NA_real_
      }

      results[[length(results) + 1]] <- data.table(
        n_per_group = n_per_group,
        gamma_fc = gfc,
        channel = ch,
        fpr_05 = round(fpr_05, 4),
        sensitivity_05 = round(sens_05, 4),
        specificity_05 = round(spec_05, 4),
        auc = round(auc, 4),
        n_valid = sum(valid),
        n_de = sum(tr)
      )

      cat(sprintf("  %s: FPR=%.3f, Sens=%.3f, AUC=%.3f\n", ch, fpr_05, sens_05, auc))
    }
  }
}

dt <- rbindlist(results)
cat("\n\n========== FULL RESULTS ==========\n")
print(dt, nrows = 100)

# Summary pivot: AUC by (n, gamma_fc, channel)
cat("\n\n========== AUC PIVOT ==========\n")
auc_pivot <- dcast(dt, n_per_group + gamma_fc ~ channel, value.var = "auc")
print(auc_pivot)

cat("\n\n========== SENSITIVITY@5% PIVOT ==========\n")
sens_pivot <- dcast(dt, n_per_group + gamma_fc ~ channel, value.var = "sensitivity_05")
print(sens_pivot)

cat("\n\n========== FPR@5% PIVOT (should be <= 0.05) ==========\n")
fpr_pivot <- dcast(dt, n_per_group + gamma_fc ~ channel, value.var = "fpr_05")
print(fpr_pivot)

outdir <- paste0(PROJECT_ROOT, "/benchmark/output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
fwrite(dt, file.path(outdir, "dg_synthetic_power.csv"))
cat("\nSaved to", file.path(outdir, "dg_synthetic_power.csv"), "\n")
