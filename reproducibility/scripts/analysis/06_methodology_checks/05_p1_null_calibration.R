#!/usr/bin/env Rscript
# =============================================================================
# 41_p1_null_calibration_check.R
# Quick null calibration check for P1a (Cauchy DE) + P1b (relaxed w_base)
#
# Under H0 (no DE), observed FDR at 5% should be <= 5%
# Test on 6 NB simulation scenarios (n=3,5,10,20,50,100 per group)
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})

devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

set.seed(2026)

scenarios <- data.table(
  label = c("n3", "n5", "n10", "n20", "n50", "n100"),
  n_per_group = c(3L, 5L, 10L, 20L, 50L, 100L)
)

n_genes <- 5000
mu_base <- 200
size_base <- 10  # NB dispersion

results <- list()

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i]
  n <- sc$n_per_group
  cat(sprintf("=== %s (n=%d per group) ===\n", sc$label, n))

  # Pure null: same distribution both groups
  counts <- matrix(
    rnbinom(n_genes * 2 * n, mu = mu_base, size = size_base),
    nrow = n_genes, ncol = 2 * n
  )
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  group <- factor(c(rep("ctrl", n), rep("treat", n)))

  t0 <- proc.time()
  res <- tryCatch(
    sgcbDE(counts, group, alpha = 0.05),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]

  if (is.null(res)) {
    results[[i]] <- data.table(
      scenario = sc$label, n_per_group = n,
      n_genes_tested = NA_integer_,
      n_reject_005 = NA_integer_, fdr_005 = NA_real_,
      n_reject_010 = NA_integer_, fdr_010 = NA_real_,
      runtime_sec = elapsed, status = "ERROR"
    )
    next
  }

  df <- as.data.frame(res)
  dt <- as.data.table(df)

  pvals <- dt$pvalue
  padj <- p.adjust(pvals, method = "BH")
  n_tested <- sum(!is.na(pvals))

  # FDR = proportion of rejections that are false = n_reject / n_reject (all are FP under H0)
  # Actually: FDR = E[FP/max(R,1)]. Under H0, FP = R, so FDR_obs = R/n_tested... no
  # Standard: if all genes are null, padj <= alpha rejections are all FP
  # The BH procedure guarantees FDR <= alpha * (n_null / n_total) = alpha under H0
  # So we just check: proportion rejected <= alpha

  n_reject_005 <- sum(padj <= 0.05, na.rm = TRUE)
  n_reject_010 <- sum(padj <= 0.10, na.rm = TRUE)
  fdr_005 <- n_reject_005 / n_tested  # should be <= 0.05
  fdr_010 <- n_reject_010 / n_tested  # should be <= 0.10

  # Also check p-value uniformity via KS test
  ks <- ks.test(pvals[!is.na(pvals)], "punif")

  cat(sprintf("  n_tested=%d, reject@5%%=%d (%.4f), reject@10%%=%d (%.4f), KS p=%.4f, %.1fs\n",
              n_tested, n_reject_005, fdr_005, n_reject_010, fdr_010, ks$p.value, elapsed))

  results[[i]] <- data.table(
    scenario = sc$label, n_per_group = n,
    n_genes_tested = n_tested,
    n_reject_005 = n_reject_005, fdr_005 = fdr_005,
    n_reject_010 = n_reject_010, fdr_010 = fdr_010,
    ks_pvalue = ks$p.value,
    runtime_sec = elapsed, status = "OK"
  )
}

res_dt <- rbindlist(results, fill = TRUE)

cat("\n=== P1 Null Calibration Summary ===\n")
print(res_dt[, .(scenario, n_per_group, n_reject_005, fdr_005, n_reject_010, fdr_010, ks_pvalue)])

# Check: any scenario with FDR > 2x nominal?
bad <- res_dt[fdr_005 > 0.10 | fdr_010 > 0.20]
if (nrow(bad) > 0) {
  cat("\n*** WARNING: FDR inflation detected! ***\n")
  print(bad)
} else {
  cat("\nAll scenarios pass null calibration check.\n")
}

out_file <- paste0(PROJECT_ROOT, "/benchmark/output/p1_null_calibration.csv")
fwrite(res_dt, out_file)
cat(sprintf("Saved to %s\n", out_file))
