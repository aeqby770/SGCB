# =============================================================================
#   dataset: cptac | ecoli_lfq | sgsds_ratio2 | sgsds_ratio2.5 | tmt_mir
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
THREADS <- as.integer(Sys.getenv("BENCHMARK_THREADS", "1"))
THREADS <- ifelse(is.na(THREADS) | THREADS < 1L, 1L, THREADS)
set.seed(12345)
Sys.setenv(
  OMP_NUM_THREADS = THREADS,
  OPENBLAS_NUM_THREADS = THREADS,
  MKL_NUM_THREADS = THREADS,
  VECLIB_MAXIMUM_THREADS = THREADS,
  NUMEXPR_NUM_THREADS = THREADS,
  BLIS_NUM_THREADS = THREADS,
  R_DATATABLE_NUM_THREADS = THREADS
)
data.table::setDTthreads(THREADS)

args      <- commandArgs(trailingOnly = TRUE)
DATASET   <- args[1]
NULL_SEED <- as.integer(args[2])

BASE     <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics")
DATA_DIR <- file.path(BASE, "data")
OUT_DIR  <- file.path(BASE, "sgcb/results")
IS_NULL  <- !is.na(NULL_SEED)
SGCB_ARGS <- list(
  prior_de = "auto",
  prior_dv = "auto",
  prior_sd_de = 1.0,
  prior_sd_dv = 0.5
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
library(data.table)
library(SGCB)

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, "SGCB__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, "SGCB__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== Prot SGCB: %s %s=====\n",
            Sys.time(), DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

switch(as.character(file.exists(OUT_CSV)),
  "TRUE" = { cat("  [SKIP] Already exists:", OUT_CSV, "\n"); q("no", status = 0) },
  "FALSE" = NULL)

source(file.path(BASE, "load_proteomics_data.R"))

prot_log2  <- loaded$prot_log2
group      <- loaded$group
pep_count  <- loaded$pep_count
is_spike   <- loaded$is_spike
sim_truth  <- loaded$sim_truth  # NULL for non-simulated datasets

# Null FDR
switch(as.character(IS_NULL),
  "TRUE" = {
    set.seed(NULL_SEED)
    group <- factor(sample(as.character(group)))
    cat("  [NULL FDR] Labels permuted\n")
    NULL
  },
  "FALSE" = NULL
)

cat(sprintf("  Data: %d proteins x %d samples, groups: %s\n",
            nrow(prot_log2), ncol(prot_log2),
            paste(paste0(levels(group), "=", table(group)), collapse=", ")))

# ===== SGCB =====
t0 <- proc.time()
med_log2 <- median(prot_log2, na.rm = TRUE)
shifted  <- prot_log2 - med_log2 + log2(1000)
pseudo_counts <- round(2^shifted)
pseudo_counts[is.na(pseudo_counts) | pseudo_counts < 0] <- 0L
storage.mode(pseudo_counts) <- "integer"
res <- do.call(sgcbProtein, c(list(counts = pseudo_counts, group = group), SGCB_ARGS))
rr  <- as.data.frame(res)
prior_de_used <- c(as.numeric(attr(res, "prior_de_used")), NA_real_)[1]
prior_dv_used <- c(as.numeric(attr(res, "prior_dv_used")), NA_real_)[1]

res_df <- data.table(
  gene           = rr$gene_id,
  pvalue         = rr$pvalue,
  padj           = rr$padj,
  log2FoldChange = rr$log2FoldChange,
  p_de_post      = rr$p_de_post,
  p_de_only_post = rr$p_de_only_post,
  de_post_call   = rr$de_post_call,
  de_fdr_call    = rr$de_fdr_call
)
has_dv_col <- "pvalue_dv" %in% names(rr)
switch(as.character(has_dv_col),
  "TRUE" = {
    res_df[, dv_pvalue := rr$pvalue_dv]
    res_df[, dv_padj := rr$padj_dv]
    res_df[, dv_pvalue_model := rr$pvalue_dv]
    res_df[, dv_padj_model := rr$padj_dv]
    res_df[, p_dv_post := rr$p_dv_post]
    res_df[, p_dv_only_post := rr$p_dv_only_post]
    res_df[, p_de_dv_post := rr$p_de_dv_post]
    res_df[, dv_post_call := rr$dv_post_call]
    res_df[, dv_fdr_call := rr$dv_fdr_call]
    NULL
  },
  "FALSE" = NULL
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  SGCB done: %.2f sec, %d proteins\n", time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]

has_truth <- !is.null(is_spike) && !IS_NULL && any(is_spike)
metrics <- data.table(method = "SGCB", n_features = nrow(res_df), time_sec = time_sec,
                       n_sig_fdr5 = sum(res_df$padj < 0.05),
                       n_post_de_0.8 = sum(res_df$p_de_post > 0.8, na.rm = TRUE),
                       n_post_deonly_0.8 = sum(res_df$p_de_only_post > 0.8, na.rm = TRUE),
                       n_post_de_call = sum(res_df$de_post_call, na.rm = TRUE),
                       n_de_fdr_call = sum(res_df$de_fdr_call, na.rm = TRUE),
                       n_post_dv_call = ifelse(has_dv_col, sum(res_df$dv_post_call, na.rm = TRUE), NA_integer_),
                       n_dv_fdr_call = ifelse(has_dv_col, sum(res_df$dv_fdr_call, na.rm = TRUE), NA_integer_),
                       prior_de_used = prior_de_used,
                       prior_dv_used = prior_dv_used,
                       prior_dv_upper = as.numeric(attr(res, "prior_dv_upper")),
                       post_de_threshold = as.numeric(attr(res, "post_de_threshold")),
                       post_dv_threshold = as.numeric(attr(res, "post_dv_threshold")),
                       post_de_realized_fdr = as.numeric(attr(res, "post_de_realized_fdr")),
                       post_dv_realized_fdr = as.numeric(attr(res, "post_dv_realized_fdr")),
                       dataset = DATASET, is_null = IS_NULL)
switch(as.character(has_truth),
  "TRUE" = {
    matched_spike <- is_spike[match(res_df$gene, rownames(prot_log2))]
    matched_spike[is.na(matched_spike)] <- FALSE
    sig <- res_df$padj < 0.05
    tp <- sum(sig & matched_spike); fp <- sum(sig & !matched_spike)
    n_sp <- sum(matched_spike)
    pr <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rc <- ifelse(n_sp == 0, 0, tp / n_sp)
    f1 <- ifelse(pr + rc == 0, 0, 2 * pr * rc / (pr + rc))
    metrics[, `:=`(n_spike = n_sp, TP_fdr5 = tp, FP_fdr5 = fp,
                   precision = round(pr, 4), recall = round(rc, 4), F1 = round(f1, 4))]
    cat(sprintf("  Spike-in: TP=%d FP=%d F1=%.4f\n", tp, fp, f1))
    NULL
  },
  "FALSE" = NULL
)

# ===== Simulation ground truth (sim_proteomics) =====
has_sim_truth <- !is.null(sim_truth) && !IS_NULL
has_dv <- "dv_pvalue" %in% names(res_df)

sim_truth_cols <- data.table(
  n_true_de = NA_integer_, n_true_dv = NA_integer_,
  sim_TP_de = NA_integer_, sim_FP_de = NA_integer_,
  sim_prec_de = NA_real_, sim_recall_de = NA_real_, sim_F1_de = NA_real_,
  sim_TP_de_post08 = NA_integer_, sim_FP_de_post08 = NA_integer_,
  sim_prec_de_post08 = NA_real_, sim_recall_de_post08 = NA_real_, sim_F1_de_post08 = NA_real_,
  sim_TP_dv = NA_integer_, sim_FP_dv = NA_integer_,
  sim_prec_dv = NA_real_, sim_recall_dv = NA_real_, sim_F1_dv = NA_real_
)

switch(as.character(has_sim_truth),
  "TRUE" = {
    mt <- sim_truth[match(res_df$gene, sim_truth$gene), ]
    td <- ifelse(is.na(mt$is_de), FALSE, mt$is_de)
    tv <- ifelse(is.na(mt$is_dv), FALSE, mt$is_dv)
    sig <- res_df$padj < 0.05
    sig_post <- res_df$p_de_only_post > 0.8
    tp_d <- sum(sig & td); fp_d <- sum(sig & !td); n_d <- sum(td)
    tp_d_post <- sum(sig_post & td, na.rm = TRUE)
    fp_d_post <- sum(sig_post & !td, na.rm = TRUE)
    pr_d <- ifelse(tp_d + fp_d == 0, 0, tp_d / (tp_d + fp_d))
    rc_d <- ifelse(n_d == 0, 0, tp_d / n_d)
    f1_d <- ifelse(pr_d + rc_d == 0, 0, 2 * pr_d * rc_d / (pr_d + rc_d))
    pr_d_post <- ifelse(tp_d_post + fp_d_post == 0, 0, tp_d_post / (tp_d_post + fp_d_post))
    rc_d_post <- ifelse(n_d == 0, 0, tp_d_post / n_d)
    f1_d_post <- ifelse(pr_d_post + rc_d_post == 0, 0, 2 * pr_d_post * rc_d_post / (pr_d_post + rc_d_post))

    n_v <- sum(tv)
    tp_v <- NA_integer_; fp_v <- NA_integer_
    pr_v <- NA_real_; rc_v <- NA_real_; f1_v <- NA_real_
    switch(as.character(has_dv),
      "TRUE" = {
        sig_v <- res_df$dv_padj < 0.05
        tp_v <- sum(sig_v & tv, na.rm = TRUE)
        fp_v <- sum(sig_v & !tv, na.rm = TRUE)
        pr_v <- ifelse(tp_v + fp_v == 0, 0, tp_v / (tp_v + fp_v))
        rc_v <- ifelse(n_v == 0, 0, tp_v / n_v)
        f1_v <- ifelse(pr_v + rc_v == 0, 0, 2 * pr_v * rc_v / (pr_v + rc_v))
        NULL
      },
      "FALSE" = NULL
    )
    sim_truth_cols <- data.table(
      n_true_de = n_d, n_true_dv = n_v,
      sim_TP_de = tp_d, sim_FP_de = fp_d,
      sim_prec_de = round(pr_d, 4), sim_recall_de = round(rc_d, 4), sim_F1_de = round(f1_d, 4),
      sim_TP_de_post08 = tp_d_post, sim_FP_de_post08 = fp_d_post,
      sim_prec_de_post08 = round(pr_d_post, 4), sim_recall_de_post08 = round(rc_d_post, 4), sim_F1_de_post08 = round(f1_d_post, 4),
      sim_TP_dv = tp_v, sim_FP_dv = fp_v,
      sim_prec_dv = round(pr_v, 4), sim_recall_dv = round(rc_v, 4), sim_F1_dv = round(f1_v, 4)
    )
    cat(sprintf("  Sim truth DE: TP=%d FP=%d F1=%.4f\n", tp_d, fp_d, f1_d))
    NULL
  },
  "FALSE" = NULL
)

metrics <- cbind(metrics, sim_truth_cols)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
