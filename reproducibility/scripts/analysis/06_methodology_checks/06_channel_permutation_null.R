#!/usr/bin/env Rscript
# =============================================================================
# 43_channel_permutation_null.R
# Q2a: Permutation null calibration for ALL SGCB channels (DE, DV, DG, DD)
#
# Strategy:
#   Use real TCGA-KIRC tumor expression data (homogeneous tissue)
#   are all approximately Uniform(0,1) under permutation null
#
# Output:
#   benchmark/output/channel_permutation_null.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})

devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

set.seed(2026)

OUT_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")

# =============================================================================
# 1. Load TCGA-KIRC tumor expression (all tumor, no real group difference)
# =============================================================================
tcga <- readRDS(paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc.rds"))
expr_full <- tcga$expression
cat(sprintf("TCGA-KIRC: %d genes x %d tumor samples\n", nrow(expr_full), ncol(expr_full)))

# Subsample to manageable size: 50 per group (100 total)
n_per_group <- 50
n_total <- 2 * n_per_group
idx <- sample(ncol(expr_full), n_total)
counts <- expr_full[, idx]

# Filter low-expression genes
keep <- rowSums(counts >= 10) >= 5
counts <- counts[keep, ]
cat(sprintf("After filtering: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# =============================================================================
# 2. Permutation null: random group assignment (no real difference)
# =============================================================================
n_perm <- 5
results <- list()

for (perm_i in seq_len(n_perm)) {
  cat(sprintf("\n=== Permutation %d/%d ===\n", perm_i, n_perm))
  
  group <- factor(sample(rep(c("ctrl", "treat"), each = n_per_group)))
  
  t0 <- proc.time()
  res <- tryCatch(
    sgcbDE(counts, group, alpha = 0.05),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  Runtime: %.1f sec\n", elapsed))
  
  if (is.null(res)) next
  
  df <- as.data.frame(res)
  dt <- as.data.table(df)
  
  # Extract all channel p-values (use actual column names from sgcbDE output)
  channels <- c("pvalue", "dv_pvalue", "dg_shape_pvalue", "dg_gamma_pvalue", "pvalue_dd", "pvalue_manifold", "SGCB_Score_p")
  avail <- intersect(channels, names(dt))
  cat(sprintf("  Available channels: %s\n", paste(avail, collapse = ", ")))
  
  for (ch in avail) {
    pvals <- dt[[ch]]
    pvals <- pvals[!is.na(pvals)]
    n_valid <- length(pvals)
    
    if (n_valid < 100) {
      cat(sprintf("    %s: only %d valid p-values, skip\n", ch, n_valid))
      next
    }
    
    # KS test for uniformity
    ks <- ks.test(pvals, "punif")
    
    # BH rejection counts at nominal levels
    padj <- p.adjust(pvals, "BH")
    n_reject_005 <- sum(padj <= 0.05)
    n_reject_010 <- sum(padj <= 0.10)
    
    # Observed FPR (proportion rejected under null)
    fpr_005 <- n_reject_005 / n_valid
    fpr_010 <- n_reject_010 / n_valid
    
    # Median p-value (should be ~0.5 under null)
    med_p <- median(pvals)
    
    cat(sprintf("    %s: n=%d, med_p=%.3f, reject@5%%=%d (%.4f), KS p=%.4f\n",
                ch, n_valid, med_p, n_reject_005, fpr_005, ks$p.value))
    
    results[[length(results) + 1]] <- data.table(
      permutation = perm_i,
      channel = ch,
      n_valid = n_valid,
      median_p = med_p,
      n_reject_005 = n_reject_005,
      fpr_005 = fpr_005,
      n_reject_010 = n_reject_010,
      fpr_010 = fpr_010,
      ks_stat = ks$statistic,
      ks_pvalue = ks$p.value,
      runtime_sec = elapsed
    )
  }
}

# =============================================================================
# 3. Summary
# =============================================================================
res_dt <- rbindlist(results)

cat("\n\n=== Channel Permutation Null Summary ===\n")
summary_dt <- res_dt[, .(
  mean_median_p = mean(median_p),
  mean_fpr_005 = mean(fpr_005),
  max_fpr_005 = max(fpr_005),
  mean_ks_pvalue = mean(ks_pvalue),
  min_ks_pvalue = min(ks_pvalue),
  n_perms = .N
), by = channel]

print(summary_dt)

# Flag problematic channels
bad <- summary_dt[mean_fpr_005 > 0.10 | min_ks_pvalue < 0.01]
if (nrow(bad) > 0) {
  cat("\n*** WARNING: These channels show potential miscalibration ***\n")
  print(bad)
} else {
  cat("\nAll channels pass permutation null calibration.\n")
}

fwrite(res_dt, file.path(OUT_DIR, "channel_permutation_null.csv"))
cat(sprintf("Saved to %s\n", file.path(OUT_DIR, "channel_permutation_null.csv")))
