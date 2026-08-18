PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages(library(data.table))
devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

set.seed(2026)
tcga <- readRDS(paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc.rds"))
expr <- tcga$expression
idx <- sample(ncol(expr), 100)
counts <- expr[, idx]
keep <- rowSums(counts >= 10) >= 5
counts <- counts[keep, ]
cat(sprintf("%d genes x %d samples\n", nrow(counts), ncol(counts)))

group <- factor(sample(rep(c("ctrl", "treat"), each = 50)))
res <- sgcbDE(counts, group)
df <- as.data.frame(res)

# Check gamma estimates
cat("\n=== Gamma estimates ===\n")
cat("ctrl_gamma: min=", round(min(df$ctrl_gamma),4), " median=", round(median(df$ctrl_gamma),4),
    " max=", round(max(df$ctrl_gamma),4), "\n")
cat("treat_gamma: min=", round(min(df$treat_gamma),4), " median=", round(median(df$treat_gamma),4),
    " max=", round(max(df$treat_gamma),4), "\n")
cat("ctrl_gamma==1:", sum(df$ctrl_gamma == 1), "/", nrow(df), "\n")
cat("treat_gamma==1:", sum(df$treat_gamma == 1), "/", nrow(df), "\n")

# Manually compute DG z-stats to diagnose
log_gamma_ratio <- log(df$treat_gamma / df$ctrl_gamma)
log_alpha_ratio <- log(df$treat_alpha / df$ctrl_alpha)

cat("\n=== log(gamma_ratio) ===\n")
cat("fraction exactly 0:", mean(log_gamma_ratio == 0), "\n")
cat("sd:", round(sd(log_gamma_ratio), 5), "\n")
cat("MAD/0.6745:", round(median(abs(log_gamma_ratio)) / 0.6745, 5), "\n")
cat("quantiles:\n")
print(round(quantile(log_gamma_ratio, c(0,.01,.05,.25,.5,.75,.95,.99,1)), 5))

cat("\n=== log(alpha_ratio) ===\n")
cat("fraction exactly 0:", mean(log_alpha_ratio == 0), "\n")
cat("sd:", round(sd(log_alpha_ratio), 5), "\n")
cat("MAD/0.6745:", round(median(abs(log_alpha_ratio)) / 0.6745, 5), "\n")
cat("quantiles:\n")
print(round(quantile(log_alpha_ratio, c(0,.01,.05,.25,.5,.75,.95,.99,1)), 5))

# Check individual channel p-values
cat("\n=== dg_gamma_pvalue ===\n")
pgam <- df$dg_gamma_pvalue[!is.na(df$dg_gamma_pvalue)]
cat("fraction < 0.05:", round(mean(pgam < 0.05), 4), "\n")
cat("fraction < 0.01:", round(mean(pgam < 0.01), 4), "\n")
cat("fraction > 0.95:", round(mean(pgam > 0.95), 4), "\n")

cat("\n=== dg_alpha_pvalue ===\n")
palp <- df$dg_alpha_pvalue[!is.na(df$dg_alpha_pvalue)]
cat("fraction < 0.05:", round(mean(palp < 0.05), 4), "\n")
cat("fraction < 0.01:", round(mean(palp < 0.01), 4), "\n")
cat("fraction > 0.95:", round(mean(palp > 0.95), 4), "\n")

cat("\n=== dg_entropy_pvalue ===\n")
pent <- df$dg_entropy_pvalue[!is.na(df$dg_entropy_pvalue)]
cat("fraction < 0.05:", round(mean(pent < 0.05), 4), "\n")

cat("\n=== z-stat SDs ===\n")
cat("z_alpha sd:", round(sd(df$dg_alpha_stat, na.rm=TRUE), 3), "\n")
cat("z_gamma sd:", round(sd(df$dg_gamma_stat, na.rm=TRUE), 3), "\n")
cat("z_entropy sd:", round(sd(df$dg_entropy_stat, na.rm=TRUE), 3), "\n")

# Check: for genes with dg_shape_pvalue < 0.05, what are their gamma values?
cat("\n=== Rejected genes (shape p < 0.05) ===\n")
rej <- df$dg_shape_pvalue < 0.05
cat("n_rejected:", sum(rej), "/", nrow(df), "\n")
if (sum(rej) > 0) {
  cat("median ctrl_gamma (rejected):", round(median(df$ctrl_gamma[rej]), 4), "\n")
  cat("median ctrl_gamma (not rejected):", round(median(df$ctrl_gamma[!rej]), 4), "\n")
  cat("median ctrl_alpha (rejected):", round(median(df$ctrl_alpha[rej]), 4), "\n")
  cat("median ctrl_alpha (not rejected):", round(median(df$ctrl_alpha[!rej]), 4), "\n")
  cat("median baseMean (rejected):", round(median(df$baseMean[rej]), 1), "\n")
  cat("median baseMean (not rejected):", round(median(df$baseMean[!rej]), 1), "\n")
}
