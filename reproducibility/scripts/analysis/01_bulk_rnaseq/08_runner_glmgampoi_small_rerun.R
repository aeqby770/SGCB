PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(glmGamPoi)

CONFIG <- list(
  SEED       = 12345L,
  N_REPS_SIM = 10L,
  OUT_DIR    = paste0(PROJECT_ROOT, "/benchmark/output/paper_figures")
)

generate_nb <- function(n_genes, n_per_group, pDE, effect_size, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_samples <- 2L * n_per_group
  base_mean <- pmax(exp(rnorm(n_genes, 4, 2)), 1)
  disp <- 0.1 + 1 / sqrt(base_mean)
  counts <- matrix(0L, n_genes, n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("s", seq_len(n_samples))
  n_de <- round(n_genes * pDE)
  de_idx <- sample(n_genes, n_de)
  de_dir <- sample(c(-1, 1), n_de, replace = TRUE)
  for (i in seq_len(n_genes)) {
    mu <- base_mean[i]; sz <- 1 / disp[i]
    counts[i, 1:n_per_group] <- rnbinom(n_per_group, size = sz, mu = mu)
    mu_t <- mu
    idx_in <- match(i, de_idx)
    if (!is.na(idx_in)) mu_t <- mu * ifelse(de_dir[idx_in] == 1, effect_size, 1/effect_size)
    counts[i, (n_per_group+1):n_samples] <- rnbinom(n_per_group, size = sz, mu = mu_t)
  }
  group <- factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group)))
  truth <- rownames(counts)[de_idx]
  list(counts = counts, group = group, truth = truth)
}

eval_de <- function(dt, truth, fdr = 0.05) {
  dt <- dt[!is.na(padj)]
  sig <- dt$gene[dt$padj < fdr]
  TP <- sum(sig %in% truth); FP <- sum(!(sig %in% truth))
  precision <- TP / max(1, TP + FP)
  recall <- TP / max(1, length(truth))
  f1 <- 2 * precision * recall / max(0.001, precision + recall)
  actual_fdr <- FP / max(1, TP + FP)
  data.table(TP = TP, FP = FP, precision = round(precision, 4),
             recall = round(recall, 4), F1 = round(f1, 4),
             actual_FDR = round(actual_fdr, 4), n_sig = length(sig))
}

run_glmgampoi <- function(counts, group) {
  col_data <- data.frame(condition = group, row.names = colnames(counts))
  fit <- glm_gp(counts, design = ~condition, col_data = col_data, size_factors = "normed_sum")
  contrast_str <- paste0("condition", levels(group)[2])
  res <- test_de(fit, contrast = contrast_str)
  data.table(gene = res$name, pvalue = res$pval, padj = res$adj_pval,
             log2FC = res$lfc, method = "glmGamPoi", time = 0)
}

cat("=== Rerun glmGamPoi small-sample benchmark ===\n")
small_results <- list()
for (n_pg in c(3L, 4L, 5L, 6L, 8L, 10L, 15L, 20L)) {
  cat("n=", n_pg, "... ")
  for (rep_i in seq_len(CONFIG$N_REPS_SIM)) {
    sim_data <- generate_nb(3000, n_pg, 0.1, 2.0, seed = CONFIG$SEED + 500 + 50*n_pg + rep_i)
    keep <- rowSums(sim_data$counts > 5) >= 2
    counts_sm <- sim_data$counts[keep, ]
    truth_sm <- intersect(sim_data$truth, rownames(counts_sm))
    res <- run_glmgampoi(counts_sm, sim_data$group)
    ev <- eval_de(res, truth_sm)
    ev[, `:=`(method = "glmGamPoi", n_per_group = n_pg, rep = rep_i)]
    small_results[[length(small_results) + 1]] <- ev
  }
  cat("OK\n")
}
small_dt <- rbindlist(small_results)
small_new <- small_dt[, .(mean_F1 = mean(F1), sd_F1 = sd(F1),
                           mean_recall = mean(recall), mean_FDR = mean(actual_FDR)),
                       by = .(method, n_per_group)]
cat("\n=== New glmGamPoi results ===\n")
print(small_new)

old_8 <- fread(file.path(CONFIG$OUT_DIR, "test6_small_sample_8methods.csv"))
updated <- rbind(old_8[method != "glmGamPoi"], small_new, fill = TRUE)
fwrite(updated, file.path(CONFIG$OUT_DIR, "test6_small_sample_8methods.csv"))
cat("\nUpdated test6_small_sample_8methods.csv\n")
