#!/usr/bin/env Rscript
# =============================================================================
# 40_dv_ground_truth_simulation.R
#      on data with KNOWN differential variability
#
# Ground truth design:
#   - 5000 null genes (same dispersion both groups)
#   - 500 DV-up genes (dispersion 2x in treatment)
#   - 500 DV-down genes (dispersion 0.5x in treatment)
#
# Methods compared:
#   - SGCB (sgcbDE with DV channel)
#   - diffVar (missMethyl)
#   - MDSeq
#   - GAMLSS
#   - clrDV
#
# Output: benchmark/output/dv_ground_truth_results.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})

devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

OUT_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. Generate synthetic DV data (NB-based, realistic RNA-seq)
# =============================================================================
generate_dv_data <- function(n_null = 5000, n_dv_up = 500, n_dv_down = 500,
                              n_per_group = 10, mu_range = c(50, 2000),
                              size_ctrl = 10, dv_fold = 2, seed = 42) {
  set.seed(seed)
  n_genes <- n_null + n_dv_up + n_dv_down

  # Random baseline means (log-uniform)
  mu <- exp(runif(n_genes,
                  log(mu_range[1]),
                  log(mu_range[2])))

  # NB size parameter (controls dispersion: var = mu + mu^2/size)
  # Higher size = lower dispersion
  ctrl_size <- rep(size_ctrl, n_genes)
  treat_size <- rep(size_ctrl, n_genes)

  dv_up_idx <- (n_null + 1):(n_null + n_dv_up)
  treat_size[dv_up_idx] <- size_ctrl / dv_fold

  dv_down_idx <- (n_null + n_dv_up + 1):n_genes
  treat_size[dv_down_idx] <- size_ctrl * dv_fold

  # Generate counts
  counts_ctrl <- matrix(
    rnbinom(n_genes * n_per_group, mu = mu, size = ctrl_size),
    nrow = n_genes, ncol = n_per_group
  )
  counts_treat <- matrix(
    rnbinom(n_genes * n_per_group, mu = mu, size = treat_size),
    nrow = n_genes, ncol = n_per_group
  )

  counts <- cbind(counts_ctrl, counts_treat)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- c(paste0("ctrl_", seq_len(n_per_group)),
                         paste0("treat_", seq_len(n_per_group)))

  group <- factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group)))

  truth <- data.table(
    gene_id = rownames(counts),
    true_dv = c(rep("null", n_null),
                rep("dv_up", n_dv_up),
                rep("dv_down", n_dv_down)),
    true_dv_any = c(rep(FALSE, n_null),
                    rep(TRUE, n_dv_up + n_dv_down)),
    ctrl_size = ctrl_size,
    treat_size = treat_size,
    dv_fold = treat_size / ctrl_size,
    mu = mu
  )

  list(counts = counts, group = group, truth = truth,
       n_per_group = n_per_group)
}

# =============================================================================
# 2. Run SGCB DV analysis
# =============================================================================
run_sgcb_dv <- function(dat) {
  res <- sgcbDE(dat$counts, dat$group, alpha = 0.05)
  df <- as.data.frame(res)
  dt <- as.data.table(df)
  if (!"gene_id" %in% names(dt)) {
    rn <- rownames(df)
    if (!is.null(rn)) dt[, gene_id := rn]
    else dt[, gene_id := paste0("gene_", seq_len(nrow(dt)))]
  }
  # SGCB outputs dv_pvalue (confirmed from differential.R line 401)
  dv_col <- intersect(c("dv_pvalue", "pvalue_dv", "p_dv"), names(dt))
  if (length(dv_col) == 0) {
    cat("  SGCB columns:", paste(head(names(dt), 30), collapse = ", "), "\n")
    return(data.table(gene_id = dt$gene_id,
                      pvalue_dv = rep(NA_real_, nrow(dt)),
                      method = "SGCB"))
  }
  data.table(gene_id = dt$gene_id,
             pvalue_dv = dt[[dv_col[1]]],
             method = "SGCB")
}

# =============================================================================
# 3. Compute DV metrics (F1, FDR, Power/Recall, Precision)
# =============================================================================
compute_dv_metrics <- function(pvals, truth, alpha = 0.05) {
  padj <- p.adjust(pvals, method = "BH")
  called <- !is.na(padj) & padj < alpha

  tp <- sum(called & truth)
  fp <- sum(called & !truth)
  fn <- sum(!called & truth)
  tn <- sum(!called & !truth)

  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall    <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1        <- ifelse(precision + recall > 0,
                      2 * precision * recall / (precision + recall), 0)
  fdr       <- ifelse(tp + fp > 0, fp / (tp + fp), 0)
  n_called  <- tp + fp

  data.table(TP = tp, FP = fp, FN = fn, TN = tn,
             precision = precision, recall = recall,
             F1 = f1, FDR_obs = fdr, n_called = n_called)
}

# =============================================================================
# 4. Main: run across multiple scenarios
# =============================================================================
scenarios <- data.table(
  scenario = c("n5_fold2", "n10_fold2", "n20_fold2",
               "n10_fold1.5", "n10_fold3", "n10_fold4"),
  n_per_group = c(5, 10, 20, 10, 10, 10),
  dv_fold = c(2, 2, 2, 1.5, 3, 4)
)

all_results <- list()

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i]
  cat(sprintf("\n=== Scenario %d/%d: %s (n=%d, fold=%.1f) ===\n",
              i, nrow(scenarios), sc$scenario, sc$n_per_group, sc$dv_fold))

  dat <- generate_dv_data(n_per_group = sc$n_per_group,
                           dv_fold = sc$dv_fold,
                           seed = 42 + i)

  # --- SGCB ---
  cat("  Running SGCB...\n")
  t0 <- proc.time()
  sgcb_res <- tryCatch(run_sgcb_dv(dat), error = function(e) {
    cat("  SGCB error:", conditionMessage(e), "\n")
    data.table(gene_id = dat$truth$gene_id,
               pvalue_dv = rep(NA_real_, nrow(dat$truth)),
               method = "SGCB")
  })
  sgcb_time <- (proc.time() - t0)[3]

  merged <- merge(sgcb_res, dat$truth, by = "gene_id")
  metrics <- compute_dv_metrics(merged$pvalue_dv, merged$true_dv_any)
  metrics[, `:=`(method = "SGCB", scenario = sc$scenario,
                 n_per_group = sc$n_per_group, dv_fold = sc$dv_fold,
                 runtime_sec = sgcb_time)]
  all_results[[length(all_results) + 1]] <- metrics

  # --- Brown-Forsythe (baseline, built into SGCB utils) ---
  cat("  Running Brown-Forsythe baseline...\n")
  t0 <- proc.time()
  bf_pvals <- tryCatch({
    norm_counts <- dat$counts + 0.5
    lib <- colSums(norm_counts)
    norm_counts <- t(t(norm_counts) / lib * median(lib))
    log2_counts <- log2(norm_counts)
    n_ctrl <- dat$n_per_group
    n_treat <- dat$n_per_group
    ctrl_idx <- seq_len(n_ctrl)
    treat_idx <- (n_ctrl + 1):(n_ctrl + n_treat)
    n_genes <- nrow(log2_counts)
    pvals <- numeric(n_genes)
    for (g in seq_len(n_genes)) {
      xc <- log2_counts[g, ctrl_idx]
      xt <- log2_counts[g, treat_idx]
      dc <- abs(xc - median(xc))
      dt <- abs(xt - median(xt))
      d_all <- c(dc, dt)
      grp <- factor(c(rep(0, n_ctrl), rep(1, n_treat)))
      fit <- tryCatch(anova(lm(d_all ~ grp)), error = function(e) NULL)
      pvals[g] <- if (!is.null(fit)) fit[["Pr(>F)"]][1] else 1
    }
    pvals
  }, error = function(e) rep(NA_real_, nrow(dat$truth)))
  bf_time <- (proc.time() - t0)[3]

  bf_metrics <- compute_dv_metrics(bf_pvals, dat$truth$true_dv_any)
  bf_metrics[, `:=`(method = "BrownForsythe", scenario = sc$scenario,
                    n_per_group = sc$n_per_group, dv_fold = sc$dv_fold,
                    runtime_sec = bf_time)]
  all_results[[length(all_results) + 1]] <- bf_metrics

  # --- Bartlett test (classic DV baseline) ---
  cat("  Running Bartlett test...\n")
  t0 <- proc.time()
  bart_pvals <- tryCatch({
    norm_counts <- dat$counts + 0.5
    lib <- colSums(norm_counts)
    norm_counts <- t(t(norm_counts) / lib * median(lib))
    log2_counts <- log2(norm_counts)
    n_ctrl <- dat$n_per_group
    n_genes <- nrow(log2_counts)
    pvals <- numeric(n_genes)
    grp <- factor(c(rep("ctrl", n_ctrl), rep("treat", n_ctrl)))
    for (g in seq_len(n_genes)) {
      fit <- tryCatch(bartlett.test(log2_counts[g, ] ~ grp),
                      error = function(e) list(p.value = 1))
      pvals[g] <- fit$p.value
    }
    pvals
  }, error = function(e) rep(NA_real_, nrow(dat$truth)))
  bart_time <- (proc.time() - t0)[3]

  bart_metrics <- compute_dv_metrics(bart_pvals, dat$truth$true_dv_any)
  bart_metrics[, `:=`(method = "Bartlett", scenario = sc$scenario,
                      n_per_group = sc$n_per_group, dv_fold = sc$dv_fold,
                      runtime_sec = bart_time)]
  all_results[[length(all_results) + 1]] <- bart_metrics

  # --- Levene (F-test on |deviations from group mean|) ---
  cat("  Running Levene test...\n")
  t0 <- proc.time()
  levene_pvals <- tryCatch({
    norm_counts <- dat$counts + 0.5
    lib <- colSums(norm_counts)
    norm_counts <- t(t(norm_counts) / lib * median(lib))
    log2_counts <- log2(norm_counts)
    n_ctrl <- dat$n_per_group
    n_genes <- nrow(log2_counts)
    pvals <- numeric(n_genes)
    grp <- factor(c(rep("ctrl", n_ctrl), rep("treat", n_ctrl)))
    for (g in seq_len(n_genes)) {
      xc <- log2_counts[g, 1:n_ctrl]
      xt <- log2_counts[g, (n_ctrl+1):(2*n_ctrl)]
      dc <- abs(xc - mean(xc))
      dt <- abs(xt - mean(xt))
      fit <- tryCatch(anova(lm(c(dc, dt) ~ grp)),
                      error = function(e) NULL)
      pvals[g] <- if (!is.null(fit)) fit[["Pr(>F)"]][1] else 1
    }
    pvals
  }, error = function(e) rep(NA_real_, nrow(dat$truth)))
  levene_time <- (proc.time() - t0)[3]

  levene_metrics <- compute_dv_metrics(levene_pvals, dat$truth$true_dv_any)
  levene_metrics[, `:=`(method = "Levene", scenario = sc$scenario,
                        n_per_group = sc$n_per_group, dv_fold = sc$dv_fold,
                        runtime_sec = levene_time)]
  all_results[[length(all_results) + 1]] <- levene_metrics

  cat(sprintf("  SGCB DV: F1=%.3f, FDR=%.3f, recall=%.3f\n",
              metrics$F1, metrics$FDR_obs, metrics$recall))
}

# =============================================================================
# 5. Combine and save
# =============================================================================
results <- rbindlist(all_results, fill = TRUE)
out_file <- file.path(OUT_DIR, "dv_ground_truth_results.csv")
fwrite(results, out_file)
cat(sprintf("\nResults saved to %s\n", out_file))

# Print summary table
cat("\n=== DV Ground Truth Benchmark Summary ===\n")
print(results[, .(F1 = round(F1, 3),
                   FDR = round(FDR_obs, 3),
                   recall = round(recall, 3),
                   precision = round(precision, 3),
                   n_called = n_called),
              by = .(scenario, method)][order(scenario, -F1)])
