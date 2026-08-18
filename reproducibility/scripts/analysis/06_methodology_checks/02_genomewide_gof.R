PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(edgeR)
library(parallel)
devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

CONFIG <- list(
  SEED = 12345L,
  BULK_DIR = paste0(PROJECT_ROOT, "/benchmark/data/bulk"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/methodology_checks"),
  MIN_COUNT = 10L,
  MIN_SAMPLES = 5L,
  GTEX_SUBSET = 80L,
  CHUNK_SIZE = 250L,
  N_WORKERS = max(1L, parallel::detectCores() - 1L),
  EPS = 1e-8
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

brca_obj <- readRDS(file.path(CONFIG$BULK_DIR, "tcga_brca_tumor_vs_normal.rds"))
gtex_obj <- readRDS(file.path(CONFIG$BULK_DIR, "gtex_liver.rds"))

dataset_meta <- data.table(
  dataset = c("TCGA_BRCA", "GTEx_Liver"),
  source = c("BRCA tumor vs normal", "GTEx liver random split")
)

raw_list <- lapply(dataset_meta$dataset, function(ds_name) {
  if (identical(ds_name, "TCGA_BRCA")) {
    counts <- brca_obj$counts
    group <- factor(brca_obj$group)
    sample_names <- colnames(counts)
  } else {
    counts_full <- round(gtex_obj$counts)
    keep_cols <- sample(seq_len(ncol(counts_full)), CONFIG$GTEX_SUBSET)
    counts <- counts_full[, keep_cols, drop = FALSE]
    sample_names <- colnames(counts)
    group <- factor(rep(c("split_A", "split_B"), each = length(keep_cols) / 2L))
  }

  keep_rows <- rowSums(counts >= CONFIG$MIN_COUNT) >= CONFIG$MIN_SAMPLES
  counts <- counts[keep_rows, , drop = FALSE]
  gene_ids <- rownames(counts)
  counts_fit <- counts[gene_ids, , drop = FALSE]

  dge <- DGEList(counts = counts_fit, group = group)
  dge <- calcNormFactors(dge)
  sf <- dge$samples$lib.size * dge$samples$norm.factors
  sf <- sf / exp(mean(log(sf)))
  norm_counts <- t(t(counts_fit) / sf) + 0.5

  idx_ctrl <- which(group == levels(group)[1])
  idx_treat <- which(group == levels(group)[2])
  mle_ctrl <- fit_gg_profile_init_cpp(norm_counts[, idx_ctrl, drop = FALSE])
  mle_treat <- fit_gg_profile_init_cpp(norm_counts[, idx_treat, drop = FALSE])
  base_mean <- rowMeans(norm_counts)

  split_idx <- split(seq_along(gene_ids), ceiling(seq_along(gene_ids) / CONFIG$CHUNK_SIZE))

  out <- rbindlist(lapply(split_idx, function(idx_block) {
    rbindlist(lapply(idx_block, function(i) {
      x_ctrl <- as.numeric(norm_counts[i, idx_ctrl])
      x_treat <- as.numeric(norm_counts[i, idx_treat])

      x_ctrl <- sort(pmax(x_ctrl, CONFIG$EPS))
      x_treat <- sort(pmax(x_treat, CONFIG$EPS))

      n_ctrl <- length(x_ctrl)
      n_treat <- length(x_treat)

      emp_ctrl <- seq_len(n_ctrl) / n_ctrl
      emp_treat <- seq_len(n_treat) / n_treat

      a_ctrl <- pmax(mle_ctrl$alpha[i], 0.05)
      b_ctrl <- pmax(mle_ctrl$beta[i], CONFIG$EPS)
      g_ctrl <- pmax(mle_ctrl$gamma[i], 0.05)
      a_treat <- pmax(mle_treat$alpha[i], 0.05)
      b_treat <- pmax(mle_treat$beta[i], CONFIG$EPS)
      g_treat <- pmax(mle_treat$gamma[i], 0.05)

      mu_ctrl <- mean(x_ctrl)
      var_ctrl <- max(stats::var(x_ctrl), CONFIG$EPS)
      shape_ctrl <- max(mu_ctrl * mu_ctrl / var_ctrl, 0.05)
      scale_ctrl <- max(var_ctrl / mu_ctrl, CONFIG$EPS)
      meanlog_ctrl <- mean(log(x_ctrl))
      sdlog_ctrl <- max(stats::sd(log(x_ctrl)), 1e-4)

      mu_treat <- mean(x_treat)
      var_treat <- max(stats::var(x_treat), CONFIG$EPS)
      shape_treat <- max(mu_treat * mu_treat / var_treat, 0.05)
      scale_treat <- max(var_treat / mu_treat, CONFIG$EPS)
      meanlog_treat <- mean(log(x_treat))
      sdlog_treat <- max(stats::sd(log(x_treat)), 1e-4)

      cdf_gg_ctrl <- pgamma((x_ctrl / b_ctrl)^g_ctrl, shape = a_ctrl, lower.tail = TRUE)
      cdf_gamma_ctrl <- pgamma(x_ctrl, shape = shape_ctrl, scale = scale_ctrl, lower.tail = TRUE)
      cdf_lnorm_ctrl <- plnorm(x_ctrl, meanlog = meanlog_ctrl, sdlog = sdlog_ctrl)

      cdf_gg_treat <- pgamma((x_treat / b_treat)^g_treat, shape = a_treat, lower.tail = TRUE)
      cdf_gamma_treat <- pgamma(x_treat, shape = shape_treat, scale = scale_treat, lower.tail = TRUE)
      cdf_lnorm_treat <- plnorm(x_treat, meanlog = meanlog_treat, sdlog = sdlog_treat)

      ll_gg_ctrl <- sum(log(g_ctrl) - g_ctrl * a_ctrl * log(b_ctrl) - lgamma(a_ctrl) + (g_ctrl * a_ctrl - 1) * log(x_ctrl) - (x_ctrl / b_ctrl)^g_ctrl)
      ll_gamma_ctrl <- sum(dgamma(x_ctrl, shape = shape_ctrl, scale = scale_ctrl, log = TRUE))
      ll_lnorm_ctrl <- sum(dlnorm(x_ctrl, meanlog = meanlog_ctrl, sdlog = sdlog_ctrl, log = TRUE))

      ll_gg_treat <- sum(log(g_treat) - g_treat * a_treat * log(b_treat) - lgamma(a_treat) + (g_treat * a_treat - 1) * log(x_treat) - (x_treat / b_treat)^g_treat)
      ll_gamma_treat <- sum(dgamma(x_treat, shape = shape_treat, scale = scale_treat, log = TRUE))
      ll_lnorm_treat <- sum(dlnorm(x_treat, meanlog = meanlog_treat, sdlog = sdlog_treat, log = TRUE))

      aic_ctrl <- c(GG = 6 - 2 * ll_gg_ctrl, Gamma = 4 - 2 * ll_gamma_ctrl, LogNormal = 4 - 2 * ll_lnorm_ctrl)
      bic_ctrl <- c(GG = log(n_ctrl) * 3 - 2 * ll_gg_ctrl, Gamma = log(n_ctrl) * 2 - 2 * ll_gamma_ctrl, LogNormal = log(n_ctrl) * 2 - 2 * ll_lnorm_ctrl)
      aic_treat <- c(GG = 6 - 2 * ll_gg_treat, Gamma = 4 - 2 * ll_gamma_treat, LogNormal = 4 - 2 * ll_lnorm_treat)
      bic_treat <- c(GG = log(n_treat) * 3 - 2 * ll_gg_treat, Gamma = log(n_treat) * 2 - 2 * ll_gamma_treat, LogNormal = log(n_treat) * 2 - 2 * ll_lnorm_treat)

      rbind(
        data.table(
          dataset = ds_name,
          gene_id = gene_ids[i],
          group_label = levels(factor(c("ctrl", "treat")))[1],
          n_obs = n_ctrl,
          baseMean = base_mean[i],
          ks_gg = max(abs(emp_ctrl - cdf_gg_ctrl)),
          ks_gamma = max(abs(emp_ctrl - cdf_gamma_ctrl)),
          ks_lognorm = max(abs(emp_ctrl - cdf_lnorm_ctrl)),
          madcdf_gg = mean(abs(emp_ctrl - cdf_gg_ctrl)),
          madcdf_gamma = mean(abs(emp_ctrl - cdf_gamma_ctrl)),
          madcdf_lognorm = mean(abs(emp_ctrl - cdf_lnorm_ctrl)),
          aic_gg = aic_ctrl[["GG"]],
          aic_gamma = aic_ctrl[["Gamma"]],
          aic_lognorm = aic_ctrl[["LogNormal"]],
          bic_gg = bic_ctrl[["GG"]],
          bic_gamma = bic_ctrl[["Gamma"]],
          bic_lognorm = bic_ctrl[["LogNormal"]],
          best_aic = names(which.min(aic_ctrl)),
          best_bic = names(which.min(bic_ctrl)),
          delta_aic_gamma_vs_gg = aic_ctrl[["Gamma"]] - aic_ctrl[["GG"]],
          delta_aic_lognorm_vs_gg = aic_ctrl[["LogNormal"]] - aic_ctrl[["GG"]],
          delta_bic_gamma_vs_gg = bic_ctrl[["Gamma"]] - bic_ctrl[["GG"]],
          delta_bic_lognorm_vs_gg = bic_ctrl[["LogNormal"]] - bic_ctrl[["GG"]]
        ),
        data.table(
          dataset = ds_name,
          gene_id = gene_ids[i],
          group_label = levels(factor(c("ctrl", "treat")))[2],
          n_obs = n_treat,
          baseMean = base_mean[i],
          ks_gg = max(abs(emp_treat - cdf_gg_treat)),
          ks_gamma = max(abs(emp_treat - cdf_gamma_treat)),
          ks_lognorm = max(abs(emp_treat - cdf_lnorm_treat)),
          madcdf_gg = mean(abs(emp_treat - cdf_gg_treat)),
          madcdf_gamma = mean(abs(emp_treat - cdf_gamma_treat)),
          madcdf_lognorm = mean(abs(emp_treat - cdf_lnorm_treat)),
          aic_gg = aic_treat[["GG"]],
          aic_gamma = aic_treat[["Gamma"]],
          aic_lognorm = aic_treat[["LogNormal"]],
          bic_gg = bic_treat[["GG"]],
          bic_gamma = bic_treat[["Gamma"]],
          bic_lognorm = bic_treat[["LogNormal"]],
          best_aic = names(which.min(aic_treat)),
          best_bic = names(which.min(bic_treat)),
          delta_aic_gamma_vs_gg = aic_treat[["Gamma"]] - aic_treat[["GG"]],
          delta_aic_lognorm_vs_gg = aic_treat[["LogNormal"]] - aic_treat[["GG"]],
          delta_bic_gamma_vs_gg = bic_treat[["Gamma"]] - bic_treat[["GG"]],
          delta_bic_lognorm_vs_gg = bic_treat[["LogNormal"]] - bic_treat[["GG"]]
        )
      )
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)

  out[, group_label := fifelse(group_label == "ctrl", as.character(levels(group)[1]), as.character(levels(group)[2]))]
  out[, base_bin := cut(
    log10(baseMean + 1),
    breaks = unique(quantile(log10(baseMean + 1), probs = seq(0, 1, 0.1), na.rm = TRUE)),
    include.lowest = TRUE
  )]
  out
})

gof_raw <- rbindlist(raw_list, use.names = TRUE, fill = TRUE)

gof_summary <- gof_raw[, .(
  n_genes = .N,
  median_ks_gg = median(ks_gg, na.rm = TRUE),
  median_ks_gamma = median(ks_gamma, na.rm = TRUE),
  median_ks_lognorm = median(ks_lognorm, na.rm = TRUE),
  median_madcdf_gg = median(madcdf_gg, na.rm = TRUE),
  median_madcdf_gamma = median(madcdf_gamma, na.rm = TRUE),
  median_madcdf_lognorm = median(madcdf_lognorm, na.rm = TRUE),
  frac_best_aic_gg = mean(best_aic == "GG", na.rm = TRUE),
  frac_best_aic_gamma = mean(best_aic == "Gamma", na.rm = TRUE),
  frac_best_aic_lognorm = mean(best_aic == "LogNormal", na.rm = TRUE),
  frac_best_bic_gg = mean(best_bic == "GG", na.rm = TRUE),
  median_delta_aic_gamma_vs_gg = median(delta_aic_gamma_vs_gg, na.rm = TRUE),
  median_delta_aic_lognorm_vs_gg = median(delta_aic_lognorm_vs_gg, na.rm = TRUE),
  median_delta_bic_gamma_vs_gg = median(delta_bic_gamma_vs_gg, na.rm = TRUE),
  median_delta_bic_lognorm_vs_gg = median(delta_bic_lognorm_vs_gg, na.rm = TRUE)
), by = .(dataset, group_label, base_bin)]

gof_overall <- gof_raw[, .(
  n_genes = .N,
  median_ks_gg = median(ks_gg, na.rm = TRUE),
  median_ks_gamma = median(ks_gamma, na.rm = TRUE),
  median_ks_lognorm = median(ks_lognorm, na.rm = TRUE),
  frac_best_aic_gg = mean(best_aic == "GG", na.rm = TRUE),
  frac_best_bic_gg = mean(best_bic == "GG", na.rm = TRUE),
  median_delta_aic_gamma_vs_gg = median(delta_aic_gamma_vs_gg, na.rm = TRUE),
  median_delta_aic_lognorm_vs_gg = median(delta_aic_lognorm_vs_gg, na.rm = TRUE)
), by = .(dataset, group_label)]

fwrite(gof_raw, file.path(CONFIG$OUT_DIR, "genomewide_gof_raw.csv.gz"))
fwrite(gof_summary, file.path(CONFIG$OUT_DIR, "genomewide_gof_binned_summary.csv"))
fwrite(gof_overall, file.path(CONFIG$OUT_DIR, "genomewide_gof_overall_summary.csv"))
