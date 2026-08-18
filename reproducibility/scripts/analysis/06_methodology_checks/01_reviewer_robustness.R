# =============================================================================
# 35_reviewer_concern_tests.R
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  SEED = 12345L,
  BULK_DIR = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/reviewer_checks"),
  FDR = 0.05,
  N_GENES_NULL = 3000L,
  N_SAMPLES_NULL = 20L,
  N_REP_NULL = 12L,
  N_GENES_3L = 3500L,
  N_SAMPLES_3L = 20L,
  N_REP_3L = 10L,
  N_REP_DOSE = 8L,
  DV_DOSE_LEVELS = c(1.2, 1.6, 2.0, 2.5, 3.0, 4.0),
  DD_DOSE_LEVELS = c(0.5, 1.0, 1.5, 2.0, 2.5),
  N_REP_TCGA_SPLIT = 10L,
  N_REP_TCGA_PERM = 20L,
  PAN_CANCERS = c("brca", "coad", "hnsc", "kirc", "lihc", "luad", "stad"),
  N_GENES_CONFOUND = 3500L,
  N_SAMPLES_CONFOUND = 80L,
  N_REP_CONFOUND = 12L,
  CONFOUND_P_DV = 0.15,
  CONFOUND_P_DD = 0.15,
  CONFOUND_DV_RATIO = 2.5,
  CONFOUND_DD_ALPHA_SHIFT = 0.8,
  CONFOUND_DD_GAMMA_SHIFT = -0.35,
  CONFOUND_SCENARIO = c("balanced", "imbalanced"),
  CONFOUND_PROP_CTRL_S1 = c(0.5, 0.8),
  CONFOUND_PROP_TREAT_S1 = c(0.5, 0.2)
)

dir.create(CONFIG$OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)
devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))
library(SGCB)

metric_dt <- function(sig, truth) {
  tp <- sum(sig %in% truth)
  fp <- sum(!(sig %in% truth))
  fn <- sum(!(truth %in% sig))
  prec <- tp / max(1, tp + fp)
  rec <- tp / max(1, length(truth))
  f1 <- 2 * prec * rec / max(1e-8, prec + rec)
  data.table(
    TP = tp,
    FP = fp,
    FN = fn,
    precision = round(prec, 4),
    recall = round(rec, 4),
    F1 = round(f1, 4),
    actual_FDR = round(fp / max(1, tp + fp), 4)
  )
}

extract_sets <- function(res, fdr = 0.05) {
  sig_de <- rownames(res)[which(!is.na(res$padj) & res$padj < fdr)]
  sig_dv <- rownames(res)[which(!is.na(res$dv_padj) & res$dv_padj < fdr)]
  sig_dd <- rownames(res)[which(!is.na(res$padj_dd) & res$padj_dd < fdr)]
  sig_omni <- rownames(res)[which(!is.na(res$SGCB_Score_padj) & res$SGCB_Score_padj < fdr)]
  dv_only <- setdiff(sig_dv, sig_de)
  dd_only <- setdiff(sig_dd, sig_de)
  shape_proxy <- setdiff(dd_only, sig_dv)
  list(sig_de = sig_de, sig_dv = sig_dv, sig_dd = sig_dd, sig_omni = sig_omni,
       dv_only = dv_only, dd_only = dd_only, shape_proxy = shape_proxy)
}

sample_gg <- function(n, alpha, beta, gamma_p) {
  y <- rgamma(n, shape = alpha, rate = 1)
  pmax(round(beta * y^(1 / gamma_p)), 0L)
}

gg_mean <- function(alpha, beta, gamma_p) {
  beta * gamma(alpha + 1 / gamma_p) / gamma(alpha)
}

jaccard <- function(a, b) length(intersect(a, b)) / max(1, length(union(a, b)))

generate_three_layer_data <- function(
  n_g,
  n_s,
  p_de,
  p_dv,
  p_dd,
  de_effect,
  dv_ratio,
  dd_alpha_shift,
  dd_gamma_shift
) {
  n_pg <- n_s / 2L

  alpha_base <- runif(n_g, 2, 6)
  gamma_base <- runif(n_g, 0.8, 2.2)
  base_mean <- pmax(exp(rnorm(n_g, 4, 1.2)), 1)
  beta_base <- base_mean / pmax(gg_mean(alpha_base, 1, gamma_base), 1e-6)

  de_idx <- sample(n_g, round(n_g * p_de))
  rest1 <- setdiff(seq_len(n_g), de_idx)
  dv_idx <- sample(rest1, round(n_g * p_dv))
  rest2 <- setdiff(rest1, dv_idx)
  dd_idx <- sample(rest2, round(n_g * p_dd))

  counts <- matrix(0L, n_g, n_s)
  rownames(counts) <- paste0("Gene", seq_len(n_g))

  lapply(seq_len(n_g), function(i) {
    a <- alpha_base[i]
    b <- beta_base[i]
    g <- gamma_base[i]

    ctrl <- sample_gg(n_pg, a, b, g)

    a2 <- a
    b2 <- b * ifelse(i %in% de_idx, de_effect, 1)
    g2 <- g

    a2 <- ifelse(i %in% dv_idx, a / dv_ratio, a2)
    m0 <- gg_mean(a, b, g)
    m_dv <- gg_mean(a2, 1, g2)
    b2 <- ifelse(i %in% dv_idx, m0 / pmax(m_dv, 1e-6), b2)

    a2 <- ifelse(i %in% dd_idx, a + dd_alpha_shift, a2)
    g2 <- ifelse(i %in% dd_idx, pmax(g + dd_gamma_shift, 0.2), g2)
    m_dd <- gg_mean(a2, 1, g2)
    b2 <- ifelse(i %in% dd_idx, m0 / pmax(m_dd, 1e-6), b2)

    counts[i, ] <<- c(ctrl, sample_gg(n_pg, a2, b2, g2))
    NULL
  })

  keep <- rowSums(counts > 5) >= 2
  counts <- counts[keep, ]
  group <- factor(c(rep("ctrl", n_pg), rep("treat", n_pg)))

  truth_de <- intersect(paste0("Gene", de_idx), rownames(counts))
  truth_dv <- intersect(paste0("Gene", dv_idx), rownames(counts))
  truth_dd <- intersect(paste0("Gene", dd_idx), rownames(counts))
  truth_dist <- union(truth_dv, truth_dd)
  truth_omni <- union(truth_de, truth_dist)

  list(
    counts = counts,
    group = group,
    truth_de = truth_de,
    truth_dv = truth_dv,
    truth_dd = truth_dd,
    truth_dist = truth_dist,
    truth_omni = truth_omni
  )
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
null_res <- rbindlist(lapply(seq_len(CONFIG$N_REP_NULL), function(i) {
  set.seed(CONFIG$SEED + i)
  n_g <- CONFIG$N_GENES_NULL
  n_s <- CONFIG$N_SAMPLES_NULL
  n_pg <- n_s / 2L
  mu <- pmax(exp(rnorm(n_g, 4, 1.2)), 1)
  size <- pmax(mu / runif(n_g, 0.2, 0.6), 1)
  counts <- sapply(seq_len(n_s), function(j) rnbinom(n_g, mu = mu, size = size))
  counts <- matrix(as.integer(counts), nrow = n_g, ncol = n_s)
  rownames(counts) <- paste0("Gene", seq_len(n_g))
  group <- factor(c(rep("ctrl", n_pg), rep("treat", n_pg)))

  fit <- sgcbDE(counts, group)
  data.table(
    replicate = i,
    channel = c("DE", "DV", "DG_gamma", "DD", "SGCB_Score"),
    observed_fpr = c(
      mean(fit$padj < CONFIG$FDR, na.rm = TRUE),
      mean(fit$dv_padj < CONFIG$FDR, na.rm = TRUE),
      mean(fit$dg_gamma_padj < CONFIG$FDR, na.rm = TRUE),
      mean(fit$padj_dd < CONFIG$FDR, na.rm = TRUE),
      mean(fit$SGCB_Score_padj < CONFIG$FDR, na.rm = TRUE)
    )
  )
}), fill = TRUE)

null_summary <- null_res[, .(
  mean_fpr = round(mean(observed_fpr, na.rm = TRUE), 4),
  sd_fpr = round(sd(observed_fpr, na.rm = TRUE), 4),
  q95_fpr = round(as.numeric(quantile(observed_fpr, 0.95, na.rm = TRUE)), 4)
), by = channel]

fwrite(null_res, file.path(CONFIG$OUT_DIR, "A_null_channel_fpr_replicates.csv"))
fwrite(null_summary, file.path(CONFIG$OUT_DIR, "A_null_channel_fpr_summary.csv"))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
three_layer_res <- rbindlist(lapply(seq_len(CONFIG$N_REP_3L), function(rep_i) {
  set.seed(CONFIG$SEED + 1000 + rep_i)
  sim <- generate_three_layer_data(
    n_g = CONFIG$N_GENES_3L,
    n_s = CONFIG$N_SAMPLES_3L,
    p_de = 0.05,
    p_dv = 0.05,
    p_dd = 0.05,
    de_effect = 2.0,
    dv_ratio = 3.0,
    dd_alpha_shift = 0.8,
    dd_gamma_shift = -0.35
  )

  fit <- sgcbDE(sim$counts, sim$group)
  sets <- extract_sets(fit, CONFIG$FDR)

  rbindlist(list(
    cbind(data.table(replicate = rep_i, family = "channel", endpoint = "DE", target = "truth_DE"), metric_dt(sets$sig_de, sim$truth_de)),
    cbind(data.table(replicate = rep_i, family = "channel", endpoint = "DV", target = "truth_DV"), metric_dt(sets$sig_dv, sim$truth_dv)),
    cbind(data.table(replicate = rep_i, family = "channel", endpoint = "DD", target = "truth_DV_or_DD"), metric_dt(sets$sig_dd, sim$truth_dist)),
    cbind(data.table(replicate = rep_i, family = "channel", endpoint = "SGCB_Score", target = "truth_DE_or_DV_or_DD"), metric_dt(sets$sig_omni, sim$truth_omni)),
    cbind(data.table(replicate = rep_i, family = "category", endpoint = "DV_only_not_DE", target = "truth_DV"), metric_dt(sets$dv_only, sim$truth_dv)),
    cbind(data.table(replicate = rep_i, family = "category", endpoint = "DD_only_not_DE", target = "truth_DV_or_DD"), metric_dt(sets$dd_only, sim$truth_dist)),
    cbind(data.table(replicate = rep_i, family = "category", endpoint = "shape_proxy_DD_not_DE_not_DV", target = "truth_DD"), metric_dt(sets$shape_proxy, sim$truth_dd))
  ), fill = TRUE)
}), fill = TRUE)

three_layer_summary <- three_layer_res[, .(
  mean_precision = round(mean(precision, na.rm = TRUE), 4),
  mean_recall = round(mean(recall, na.rm = TRUE), 4),
  mean_F1 = round(mean(F1, na.rm = TRUE), 4),
  mean_actual_FDR = round(mean(actual_FDR, na.rm = TRUE), 4)
), by = .(family, endpoint, target)]

fwrite(three_layer_res, file.path(CONFIG$OUT_DIR, "B_threelayer_metrics_replicates.csv"))
fwrite(three_layer_summary, file.path(CONFIG$OUT_DIR, "B_threelayer_metrics_summary.csv"))

dose_dv_res <- rbindlist(lapply(seq_along(CONFIG$DV_DOSE_LEVELS), function(i_level) {
  dose_level <- CONFIG$DV_DOSE_LEVELS[i_level]
  rbindlist(lapply(seq_len(CONFIG$N_REP_DOSE), function(rep_i) {
    set.seed(CONFIG$SEED + 5000 + i_level * 100 + rep_i)
    sim <- generate_three_layer_data(
      n_g = CONFIG$N_GENES_3L,
      n_s = CONFIG$N_SAMPLES_3L,
      p_de = 0,
      p_dv = 0.10,
      p_dd = 0,
      de_effect = 1,
      dv_ratio = dose_level,
      dd_alpha_shift = 0,
      dd_gamma_shift = 0
    )

    fit <- sgcbDE(sim$counts, sim$group)
    sets <- extract_sets(fit, CONFIG$FDR)

    rbindlist(list(
      cbind(data.table(replicate = rep_i, dose_type = "DV", dose_level = dose_level, endpoint = "DV_channel", target = "truth_DV"), metric_dt(sets$sig_dv, sim$truth_dv)),
      cbind(data.table(replicate = rep_i, dose_type = "DV", dose_level = dose_level, endpoint = "DV_only_not_DE", target = "truth_DV"), metric_dt(sets$dv_only, sim$truth_dv)),
      cbind(data.table(replicate = rep_i, dose_type = "DV", dose_level = dose_level, endpoint = "DD_channel", target = "truth_DV"), metric_dt(sets$sig_dd, sim$truth_dv))
    ), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

dose_dd_res <- rbindlist(lapply(seq_along(CONFIG$DD_DOSE_LEVELS), function(i_level) {
  dose_level <- CONFIG$DD_DOSE_LEVELS[i_level]
  rbindlist(lapply(seq_len(CONFIG$N_REP_DOSE), function(rep_i) {
    set.seed(CONFIG$SEED + 6000 + i_level * 100 + rep_i)
    sim <- generate_three_layer_data(
      n_g = CONFIG$N_GENES_3L,
      n_s = CONFIG$N_SAMPLES_3L,
      p_de = 0,
      p_dv = 0.05,
      p_dd = 0.05,
      de_effect = 1,
      dv_ratio = 1 + dose_level,
      dd_alpha_shift = 0.8 * dose_level,
      dd_gamma_shift = -0.35 * dose_level
    )

    fit <- sgcbDE(sim$counts, sim$group)
    sets <- extract_sets(fit, CONFIG$FDR)

    rbindlist(list(
      cbind(data.table(replicate = rep_i, dose_type = "DD", dose_level = dose_level, endpoint = "DD_channel", target = "truth_DV_or_DD"), metric_dt(sets$sig_dd, sim$truth_dist)),
      cbind(data.table(replicate = rep_i, dose_type = "DD", dose_level = dose_level, endpoint = "DD_only_not_DE", target = "truth_DV_or_DD"), metric_dt(sets$dd_only, sim$truth_dist)),
      cbind(data.table(replicate = rep_i, dose_type = "DD", dose_level = dose_level, endpoint = "shape_proxy_DD_not_DE_not_DV", target = "truth_DD"), metric_dt(sets$shape_proxy, sim$truth_dd))
    ), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

dose_res <- rbindlist(list(dose_dv_res, dose_dd_res), fill = TRUE)

dose_summary <- dose_res[, .(
  mean_precision = round(mean(precision, na.rm = TRUE), 4),
  mean_recall = round(mean(recall, na.rm = TRUE), 4),
  mean_F1 = round(mean(F1, na.rm = TRUE), 4),
  mean_actual_FDR = round(mean(actual_FDR, na.rm = TRUE), 4)
), by = .(dose_type, dose_level, endpoint, target)]

dose_monotonic <- dose_summary[, .(
  spearman_rho_recall = round(cor(dose_level, mean_recall, method = "spearman", use = "pairwise.complete.obs"), 4),
  pearson_r_recall = round(cor(dose_level, mean_recall, method = "pearson", use = "pairwise.complete.obs"), 4)
), by = .(dose_type, endpoint, target)]

fwrite(dose_res, file.path(CONFIG$OUT_DIR, "B_dose_response_replicates.csv"))
fwrite(dose_summary, file.path(CONFIG$OUT_DIR, "B_dose_response_summary.csv"))
fwrite(dose_monotonic, file.path(CONFIG$OUT_DIR, "B_dose_response_monotonicity.csv"))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
tcga <- readRDS(file.path(CONFIG$BULK_DIR, "tcga_brca_tumor_vs_normal.rds"))
counts_tcga <- tcga$counts
group_tcga <- factor(tcga$group, levels = c("normal", "tumor"))
keep_tcga <- rowSums(counts_tcga >= 10) >= 5L
counts_tcga <- counts_tcga[keep_tcga, ]

fit_obs <- sgcbDE(counts_tcga, group_tcga)
sets_obs <- extract_sets(fit_obs, CONFIG$FDR)

idx_normal <- which(group_tcga == "normal")
idx_tumor <- which(group_tcga == "tumor")

split_repro <- rbindlist(lapply(seq_len(CONFIG$N_REP_TCGA_SPLIT), function(i) {
  set.seed(CONFIG$SEED + 2000 + i)
  n1 <- length(idx_normal)
  n2 <- length(idx_tumor)
  n1h <- floor(n1 / 2)
  n2h <- floor(n2 / 2)

  n1_a <- sample(idx_normal, n1h)
  n2_a <- sample(idx_tumor, n2h)
  idx_a <- c(n1_a, n2_a)
  idx_b <- setdiff(seq_len(ncol(counts_tcga)), idx_a)

  fit_a <- sgcbDE(counts_tcga[, idx_a], group_tcga[idx_a])
  fit_b <- sgcbDE(counts_tcga[, idx_b], group_tcga[idx_b])
  sets_a <- extract_sets(fit_a, CONFIG$FDR)
  sets_b <- extract_sets(fit_b, CONFIG$FDR)

  data.table(
    replicate = i,
    metric = c("DV_only_jaccard", "DD_only_jaccard", "shape_proxy_jaccard"),
    value = c(
      jaccard(sets_a$dv_only, sets_b$dv_only),
      jaccard(sets_a$dd_only, sets_b$dd_only),
      jaccard(sets_a$shape_proxy, sets_b$shape_proxy)
    )
  )
}), fill = TRUE)

perm_baseline <- rbindlist(lapply(seq_len(CONFIG$N_REP_TCGA_PERM), function(i) {
  set.seed(CONFIG$SEED + 3000 + i)
  fit_p <- sgcbDE(counts_tcga, sample(group_tcga, length(group_tcga), replace = FALSE))
  sets_p <- extract_sets(fit_p, CONFIG$FDR)
  data.table(
    replicate = i,
    n_dv_only = length(sets_p$dv_only),
    n_dd_only = length(sets_p$dd_only),
    n_shape_proxy = length(sets_p$shape_proxy)
  )
}), fill = TRUE)

obs_counts <- data.table(
  n_dv_only = length(sets_obs$dv_only),
  n_dd_only = length(sets_obs$dd_only),
  n_shape_proxy = length(sets_obs$shape_proxy)
)

tcga_perm_test <- data.table(
  endpoint = c("DV_only", "DD_only", "shape_proxy"),
  observed_n = c(obs_counts$n_dv_only, obs_counts$n_dd_only, obs_counts$n_shape_proxy),
  perm_mean = c(mean(perm_baseline$n_dv_only), mean(perm_baseline$n_dd_only), mean(perm_baseline$n_shape_proxy)),
  perm_p = c(
    (1 + sum(perm_baseline$n_dv_only >= obs_counts$n_dv_only)) / (1 + nrow(perm_baseline)),
    (1 + sum(perm_baseline$n_dd_only >= obs_counts$n_dd_only)) / (1 + nrow(perm_baseline)),
    (1 + sum(perm_baseline$n_shape_proxy >= obs_counts$n_shape_proxy)) / (1 + nrow(perm_baseline))
  )
)

split_summary <- split_repro[, .(
  mean_value = round(mean(value, na.rm = TRUE), 4),
  sd_value = round(sd(value, na.rm = TRUE), 4)
), by = metric]

pan_files <- file.path(
  CONFIG$BULK_DIR,
  paste0("tcga_", CONFIG$PAN_CANCERS, "_tumor_vs_normal.rds")
)

pan_results <- lapply(seq_along(CONFIG$PAN_CANCERS), function(i) {
  cancer <- toupper(CONFIG$PAN_CANCERS[i])
  dat_i <- readRDS(pan_files[i])
  counts_i <- dat_i$counts
  group_i <- factor(dat_i$group, levels = c("normal", "tumor"))
  keep_i <- rowSums(counts_i >= 10) >= 5L
  counts_i <- counts_i[keep_i, ]

  fit_i <- sgcbDE(counts_i, group_i)
  sets_i <- extract_sets(fit_i, CONFIG$FDR)

  dt_i <- data.table(
    gene = rownames(fit_i),
    padj_dv = fit_i$dv_padj,
    padj_dd = fit_i$padj_dd
  )
  dt_i[, dv_only := gene %in% sets_i$dv_only]
  dt_i[, dd_only := gene %in% sets_i$dd_only]

  dv_rank <- dt_i[dv_only == TRUE, .(gene, padj_dv)]
  setorder(dv_rank, padj_dv, gene)
  dv_rank[, rank := seq_len(.N)]

  dd_rank <- dt_i[dd_only == TRUE, .(gene, padj_dd)]
  setorder(dd_rank, padj_dd, gene)
  dd_rank[, rank := seq_len(.N)]

  list(
    cancer = cancer,
    dv_rank = dv_rank[, .(gene, rank)],
    dd_rank = dd_rank[, .(gene, rank)]
  )
})

pan_cancers <- unlist(lapply(pan_results, function(x) x$cancer), use.names = FALSE)
pan_dv_rank_map <- setNames(lapply(pan_results, function(x) x$dv_rank), pan_cancers)
pan_dd_rank_map <- setNames(lapply(pan_results, function(x) x$dd_rank), pan_cancers)

pan_dv_long <- rbindlist(lapply(seq_along(pan_cancers), function(i) {
  dt_i <- pan_dv_rank_map[[pan_cancers[i]]]
  data.table(cancer = rep(pan_cancers[i], nrow(dt_i)), gene = dt_i$gene, rank = dt_i$rank)
}), fill = TRUE)
pan_dd_long <- rbindlist(lapply(seq_along(pan_cancers), function(i) {
  dt_i <- pan_dd_rank_map[[pan_cancers[i]]]
  data.table(cancer = rep(pan_cancers[i], nrow(dt_i)), gene = dt_i$gene, rank = dt_i$rank)
}), fill = TRUE)

pan_dv_recur <- pan_dv_long[, .(
  n_cancers = uniqueN(cancer),
  mean_rank = round(mean(rank), 2),
  sd_rank = round(sd(rank), 2)
), by = gene][order(-n_cancers, mean_rank, gene)]
pan_dd_recur <- pan_dd_long[, .(
  n_cancers = uniqueN(cancer),
  mean_rank = round(mean(rank), 2),
  sd_rank = round(sd(rank), 2)
), by = gene][order(-n_cancers, mean_rank, gene)]

pan_pairs <- CJ(cancer1 = pan_cancers, cancer2 = pan_cancers, sorted = FALSE)[cancer1 < cancer2]

dv_rank_consistency <- rbindlist(lapply(seq_len(nrow(pan_pairs)), function(i) {
  c1 <- pan_pairs$cancer1[i]
  c2 <- pan_pairs$cancer2[i]
  overlap <- merge(pan_dv_rank_map[[c1]], pan_dv_rank_map[[c2]], by = "gene", suffixes = c("_1", "_2"))
  top1 <- head(pan_dv_rank_map[[c1]]$gene, 200)
  top2 <- head(pan_dv_rank_map[[c2]]$gene, 200)
  data.table(
    cancer1 = c1,
    cancer2 = c2,
    n_overlap = nrow(overlap),
    spearman_rho = cor(overlap$rank_1, overlap$rank_2, method = "spearman", use = "pairwise.complete.obs"),
    top200_jaccard = jaccard(top1, top2)
  )
}), fill = TRUE)

dd_rank_consistency <- rbindlist(lapply(seq_len(nrow(pan_pairs)), function(i) {
  c1 <- pan_pairs$cancer1[i]
  c2 <- pan_pairs$cancer2[i]
  overlap <- merge(pan_dd_rank_map[[c1]], pan_dd_rank_map[[c2]], by = "gene", suffixes = c("_1", "_2"))
  top1 <- head(pan_dd_rank_map[[c1]]$gene, 200)
  top2 <- head(pan_dd_rank_map[[c2]]$gene, 200)
  data.table(
    cancer1 = c1,
    cancer2 = c2,
    n_overlap = nrow(overlap),
    spearman_rho = cor(overlap$rank_1, overlap$rank_2, method = "spearman", use = "pairwise.complete.obs"),
    top200_jaccard = jaccard(top1, top2)
  )
}), fill = TRUE)

rank_consistency_summary <- rbindlist(list(
  dv_rank_consistency[, .(
    domain = "DV_only",
    mean_overlap = round(mean(n_overlap, na.rm = TRUE), 2),
    mean_spearman_rho = round(mean(spearman_rho, na.rm = TRUE), 4),
    mean_top200_jaccard = round(mean(top200_jaccard, na.rm = TRUE), 4)
  )],
  dd_rank_consistency[, .(
    domain = "DD_only",
    mean_overlap = round(mean(n_overlap, na.rm = TRUE), 2),
    mean_spearman_rho = round(mean(spearman_rho, na.rm = TRUE), 4),
    mean_top200_jaccard = round(mean(top200_jaccard, na.rm = TRUE), 4)
  )]
), fill = TRUE)

fwrite(split_repro, file.path(CONFIG$OUT_DIR, "C_tcga_split_half_repro.csv"))
fwrite(split_summary, file.path(CONFIG$OUT_DIR, "C_tcga_split_half_repro_summary.csv"))
fwrite(perm_baseline, file.path(CONFIG$OUT_DIR, "C_tcga_permutation_baseline.csv"))
fwrite(tcga_perm_test, file.path(CONFIG$OUT_DIR, "C_tcga_permutation_test_summary.csv"))
fwrite(pan_dv_recur, file.path(CONFIG$OUT_DIR, "C_pan_cancer_dv_only_recurrent_rank.csv"))
fwrite(pan_dd_recur, file.path(CONFIG$OUT_DIR, "C_pan_cancer_dd_only_recurrent_rank.csv"))
fwrite(dv_rank_consistency, file.path(CONFIG$OUT_DIR, "C_pan_cancer_dv_rank_consistency.csv"))
fwrite(dd_rank_consistency, file.path(CONFIG$OUT_DIR, "C_pan_cancer_dd_rank_consistency.csv"))
fwrite(rank_consistency_summary, file.path(CONFIG$OUT_DIR, "C_pan_cancer_rank_consistency_summary.csv"))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
confound_grid <- CJ(
  scenario = CONFIG$CONFOUND_SCENARIO,
  replicate = seq_len(CONFIG$N_REP_CONFOUND),
  sorted = FALSE
)

confound_res <- rbindlist(lapply(seq_len(nrow(confound_grid)), function(i) {
  scenario_i <- confound_grid$scenario[i]
  rep_i <- confound_grid$replicate[i]
  scenario_idx <- match(scenario_i, CONFIG$CONFOUND_SCENARIO)

  set.seed(CONFIG$SEED + 7000 + scenario_idx * 100 + rep_i)
  n_g <- CONFIG$N_GENES_CONFOUND
  n_s <- CONFIG$N_SAMPLES_CONFOUND
  n_c <- n_s / 2L

  prop_ctrl_s1 <- CONFIG$CONFOUND_PROP_CTRL_S1[scenario_idx]
  prop_treat_s1 <- CONFIG$CONFOUND_PROP_TREAT_S1[scenario_idx]
  n_ctrl_s1 <- round(n_c * prop_ctrl_s1)
  n_treat_s1 <- round(n_c * prop_treat_s1)

  subtype <- c(
    sample(c(rep("S1", n_ctrl_s1), rep("S2", n_c - n_ctrl_s1))),
    sample(c(rep("S1", n_treat_s1), rep("S2", n_c - n_treat_s1)))
  )
  subtype_s2 <- subtype == "S2"
  group <- factor(c(rep("ctrl", n_c), rep("treat", n_c)))

  alpha_base <- runif(n_g, 2, 6)
  gamma_base <- runif(n_g, 0.8, 2.2)
  base_mean <- pmax(exp(rnorm(n_g, 4, 1.2)), 1)
  beta_base <- base_mean / pmax(gg_mean(alpha_base, 1, gamma_base), 1e-6)

  dv_idx <- sample(n_g, round(n_g * CONFIG$CONFOUND_P_DV))
  rest_idx <- setdiff(seq_len(n_g), dv_idx)
  dd_idx <- sample(rest_idx, round(n_g * CONFIG$CONFOUND_P_DD))

  counts <- matrix(0L, nrow = n_g, ncol = n_s)
  rownames(counts) <- paste0("Gene", seq_len(n_g))

  lapply(seq_len(n_g), function(g_i) {
    a <- alpha_base[g_i]
    b <- beta_base[g_i]
    g <- gamma_base[g_i]

    m0 <- gg_mean(a, b, g)
    a_dv <- a / CONFIG$CONFOUND_DV_RATIO
    b_dv <- m0 / pmax(gg_mean(a_dv, 1, g), 1e-6)

    a_dd <- a + CONFIG$CONFOUND_DD_ALPHA_SHIFT
    g_dd <- pmax(g + CONFIG$CONFOUND_DD_GAMMA_SHIFT, 0.2)
    b_dd <- m0 / pmax(gg_mean(a_dd, 1, g_dd), 1e-6)

    is_dv <- g_i %in% dv_idx
    is_dd <- g_i %in% dd_idx

    a_vec <- rep(a, n_s)
    g_vec <- rep(g, n_s)
    b_vec <- rep(b, n_s)

    a_vec <- ifelse(subtype_s2 & is_dv, a_dv, a_vec)
    b_vec <- ifelse(subtype_s2 & is_dv, b_dv, b_vec)

    a_vec <- ifelse(subtype_s2 & is_dd, a_dd, a_vec)
    g_vec <- ifelse(subtype_s2 & is_dd, g_dd, g_vec)
    b_vec <- ifelse(subtype_s2 & is_dd, b_dd, b_vec)

    y <- rgamma(n_s, shape = a_vec, rate = 1)
    counts[g_i, ] <<- pmax(round(b_vec * y^(1 / g_vec)), 0L)
    NULL
  })

  fit <- sgcbDE(counts, group)
  sets <- extract_sets(fit, CONFIG$FDR)

  truth_dv <- paste0("Gene", dv_idx)
  truth_dd <- paste0("Gene", dd_idx)

  data.table(
    scenario = scenario_i,
    replicate = rep_i,
    de_fpr = mean(fit$padj < CONFIG$FDR, na.rm = TRUE),
    dv_fpr = mean(fit$dv_padj < CONFIG$FDR, na.rm = TRUE),
    dd_fpr = mean(fit$padj_dd < CONFIG$FDR, na.rm = TRUE),
    score_fpr = mean(fit$SGCB_Score_padj < CONFIG$FDR, na.rm = TRUE),
    n_dv_only = length(sets$dv_only),
    n_dd_only = length(sets$dd_only),
    n_shape_proxy = length(sets$shape_proxy),
    dv_confounded_capture = sum(sets$sig_dv %in% truth_dv) / max(1, length(truth_dv)),
    dd_confounded_capture = sum(sets$sig_dd %in% truth_dd) / max(1, length(truth_dd))
  )
}), fill = TRUE)

confound_summary <- confound_res[, .(
  mean_de_fpr = round(mean(de_fpr, na.rm = TRUE), 4),
  mean_dv_fpr = round(mean(dv_fpr, na.rm = TRUE), 4),
  mean_dd_fpr = round(mean(dd_fpr, na.rm = TRUE), 4),
  mean_score_fpr = round(mean(score_fpr, na.rm = TRUE), 4),
  mean_n_dv_only = round(mean(n_dv_only), 1),
  mean_n_dd_only = round(mean(n_dd_only), 1),
  mean_n_shape_proxy = round(mean(n_shape_proxy), 1),
  mean_dv_confounded_capture = round(mean(dv_confounded_capture, na.rm = TRUE), 4),
  mean_dd_confounded_capture = round(mean(dd_confounded_capture, na.rm = TRUE), 4)
), by = scenario]

confound_long <- melt(
  confound_summary,
  id.vars = "scenario",
  variable.name = "metric",
  value.name = "value"
)
confound_contrast <- dcast(confound_long, metric ~ scenario, value.var = "value")
confound_contrast[, delta_imbalanced_minus_balanced := imbalanced - balanced]

fwrite(confound_res, file.path(CONFIG$OUT_DIR, "D_shape_variance_confound_replicates.csv"))
fwrite(confound_summary, file.path(CONFIG$OUT_DIR, "D_shape_variance_confound_summary.csv"))
fwrite(confound_contrast, file.path(CONFIG$OUT_DIR, "D_shape_variance_confound_contrast.csv"))
fwrite(confound_res, file.path(CONFIG$OUT_DIR, "D_subtype_imbalance_replicates.csv"))
fwrite(confound_summary, file.path(CONFIG$OUT_DIR, "D_subtype_imbalance_summary.csv"))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
summary_out <- rbindlist(list(
  cbind(section = "A_null_channel_fpr", null_summary),
  cbind(section = "B_threelayer", three_layer_summary),
  cbind(section = "B_dose_monotonicity", dose_monotonic),
  cbind(section = "C_tcga_permutation", tcga_perm_test),
  cbind(section = "C_pan_cancer_rank_consistency", rank_consistency_summary),
  cbind(section = "D_shape_variance_confound", confound_summary),
  cbind(section = "D_shape_variance_confound_delta", confound_contrast)
), fill = TRUE)

fwrite(summary_out, file.path(CONFIG$OUT_DIR, "00_reviewer_checks_summary.csv"))

cat("\nDone. Outputs in:", CONFIG$OUT_DIR, "\n")
