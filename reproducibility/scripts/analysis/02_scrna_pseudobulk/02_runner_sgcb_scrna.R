# =============================================================================
#   dataset: kang18 | segerstolpe | xin | zilionis | jakel | muscat_sim
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

BASE     <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk")
PB_DIR   <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/05_pseudobulk")
OUT_DIR  <- file.path(BASE, "sgcb/results")
IS_NULL  <- !is.na(NULL_SEED)
SGCB_ARGS <- list(
  prior_de = "auto",
  prior_dv = "auto",
  prior_dg = "auto",
  prior_sd_de = 1.0,
  prior_sd_dv = 0.5,
  prior_sd_dg_alpha = 0.5,
  prior_sd_dg_gamma = 0.5
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
library(data.table)
library(SGCB)

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, "SGCB__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, "SGCB__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== scRNA SGCB: %s %s=====\n",
            Sys.time(), DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

switch(as.character(file.exists(OUT_CSV)),
  "TRUE" = { cat("  [SKIP] Already exists:", OUT_CSV, "\n"); q("no", status = 0) },
  "FALSE" = NULL)

pb_map <- c(
  kang18       = "kang18_pseudobulk.rds",
  segerstolpe  = "segerstolpe_t2d_pseudobulk.rds",
  xin          = "xin_t2d_pseudobulk.rds",
  zilionis     = "zilionis_lung_pseudobulk.rds",
  jakel        = "jakel_pseudobulk.rds",
  muscat_sim   = "muscat_sim_pseudobulk.rds"
)

ds     <- readRDS(file.path(PB_DIR, pb_map[DATASET]))
counts <- as.matrix(ds$counts)
storage.mode(counts) <- "integer"
rownames(counts) <- make.unique(rownames(counts))
group  <- factor(gsub("[^A-Za-z0-9._]", "_", as.character(ds$group)))
sim_truth <- ds$truth

switch(as.character(IS_NULL),
  "TRUE" = {
    set.seed(NULL_SEED)
    group <- factor(sample(as.character(group)))
    cat("  [NULL FDR] Labels permuted\n")
    NULL
  },
  "FALSE" = NULL
)

cat(sprintf("  Data: %d genes x %d samples, groups: %s\n",
            nrow(counts), ncol(counts),
            paste(paste0(levels(group), "=", table(group)), collapse=", ")))

t0  <- proc.time()
res <- do.call(sgcbDE, c(list(counts = counts, group = group), SGCB_ARGS))
rr  <- as.data.frame(res)

res_df <- data.table(
  gene           = rr$gene_id,
  pvalue         = rr$pvalue,
  padj           = rr$padj,
  log2FoldChange = rr$log2FoldChange,
  log2FC_shrunk  = rr$log2FC_shrunk,
  p_de_post      = rr$p_de_post,
  p_de_only_post = rr$p_de_only_post,
  de_post_call   = rr$de_post_call,
  de_fdr_call    = rr$de_fdr_call,
  p_dv_post      = rr$p_dv_post,
  p_dv_only_post = rr$p_dv_only_post,
  dv_post_call   = rr$dv_post_call,
  dv_fdr_call    = rr$dv_fdr_call,
  p_dg_post      = rr$p_dg_post,
  p_dg_only_post = rr$p_dg_only_post,
  dg_post_call   = rr$dg_post_call,
  dg_fdr_call    = rr$dg_fdr_call,
  dv_pvalue      = rr$dv_pvalue,
  dv_padj        = rr$dv_padj,
  dv_pvalue_model = rr$dv_pvalue_model,
  dv_padj_model = rr$dv_padj_model,
  dv_pvalue_model_raw = rr$dv_pvalue_model_raw,
  dv_padj_model_raw = rr$dv_padj_model_raw,
  dv_pvalue_bf = rr$dv_pvalue_bf,
  dv_padj_bf = rr$dv_padj_bf,
  dv_pvalue_mix = rr$dv_pvalue_mix,
  dv_padj_mix = rr$dv_padj_mix,
  dv_log2_cv2    = rr$dv_log2_cv2_ratio,
  dg_gamma_pvalue = rr$dg_gamma_pvalue,
  dg_gamma_padj  = rr$dg_gamma_padj,
  pvalue_dd      = rr$pvalue_dd,
  padj_dd        = rr$padj_dd,
  SGCB_Score     = rr$SGCB_Score,
  SGCB_Score_p   = rr$SGCB_Score_p,
  SGCB_Score_padj = rr$SGCB_Score_padj,
  geodesic_dist  = rr$geodesic_dist
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  SGCB done: %.2f sec, %d genes\n", time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]
has_dv <- "dv_pvalue" %in% names(res_df)
has_dd <- "pvalue_dd" %in% names(res_df)
has_score <- "SGCB_Score_padj" %in% names(res_df)

metrics <- data.table(
  method       = "SGCB",
  n_genes      = nrow(res_df),
  time_sec     = time_sec,
  n_sig_fdr5   = sum(res_df$padj < 0.05),
  n_sig_fdr10  = sum(res_df$padj < 0.10),
  p_median     = median(res_df$pvalue),
  lfc_median   = median(abs(res_df$log2FoldChange)),
  n_post_de_0.5 = sum(res_df$p_de_post > 0.5, na.rm = TRUE),
  n_post_de_0.8 = sum(res_df$p_de_post > 0.8, na.rm = TRUE),
  n_post_deonly_0.8 = sum(res_df$p_de_only_post > 0.8, na.rm = TRUE),
  n_post_de_call = sum(res_df$de_post_call, na.rm = TRUE),
  n_de_fdr_call = sum(res_df$de_fdr_call, na.rm = TRUE),
  n_post_dv_0.8 = sum(res_df$p_dv_post > 0.8, na.rm = TRUE),
  n_post_dvonly_0.8 = sum(res_df$p_dv_only_post > 0.8, na.rm = TRUE),
  n_post_dv_call = sum(res_df$dv_post_call, na.rm = TRUE),
  n_dv_fdr_call = sum(res_df$dv_fdr_call, na.rm = TRUE),
  n_post_dg_0.8 = sum(res_df$p_dg_post > 0.8, na.rm = TRUE),
  n_post_dgonly_0.8 = sum(res_df$p_dg_only_post > 0.8, na.rm = TRUE),
  n_post_dg_call = sum(res_df$dg_post_call, na.rm = TRUE),
  n_dg_fdr_call = sum(res_df$dg_fdr_call, na.rm = TRUE),
  prior_de_used = as.numeric(attr(res, "prior_de_used")),
  prior_dv_used = as.numeric(attr(res, "prior_dv_used")),
  prior_dv_upper = as.numeric(attr(res, "prior_dv_upper")),
  prior_dg_used = as.numeric(attr(res, "prior_dg_used")),
  post_de_threshold = as.numeric(attr(res, "post_de_threshold")),
  post_dv_threshold = as.numeric(attr(res, "post_dv_threshold")),
  post_dg_threshold = as.numeric(attr(res, "post_dg_threshold")),
  post_de_realized_fdr = as.numeric(attr(res, "post_de_realized_fdr")),
  post_dv_realized_fdr = as.numeric(attr(res, "post_dv_realized_fdr")),
  post_dg_realized_fdr = as.numeric(attr(res, "post_dg_realized_fdr")),
  dataset      = DATASET,
  is_null      = IS_NULL,
  n_dv_fdr5    = ifelse(has_dv, sum(res_df$dv_padj < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_model_raw_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_model_raw < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_model_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_model < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_bf_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_bf < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_mix_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_mix < 0.05, na.rm = TRUE), NA_integer_),
  n_dd_fdr5    = ifelse(has_dd, sum(res_df$padj_dd < 0.05, na.rm = TRUE), NA_integer_),
  n_score_fdr5 = ifelse(has_score, sum(res_df$SGCB_Score_padj < 0.05, na.rm = TRUE), NA_integer_)
)

# Ground truth
has_truth_sim <- !is.null(sim_truth) && !IS_NULL
truth_cols <- data.table(n_true_de = NA_integer_, TP_de_fdr5 = NA_integer_,
                          FP_de_fdr5 = NA_integer_, prec_de = NA_real_,
                          recall_de = NA_real_, F1_de = NA_real_,
                          TP_de_post08 = NA_integer_, FP_de_post08 = NA_integer_,
                          prec_de_post08 = NA_real_, recall_de_post08 = NA_real_,
                          F1_de_post08 = NA_real_)
switch(as.character(has_truth_sim),
  "TRUE" = {
    matched <- sim_truth[match(res_df$gene, sim_truth$gene), ]
    td <- ifelse(is.na(matched$is_de), FALSE, matched$is_de)
    sig <- res_df$padj < 0.05
    sig_post <- res_df$p_de_only_post > 0.8
    tp <- sum(sig & td); fp <- sum(sig & !td); n_d <- sum(td)
    tp_p <- sum(sig_post & td); fp_p <- sum(sig_post & !td)
    pr <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rc <- ifelse(n_d == 0, 0, tp / n_d)
    f1 <- ifelse(pr + rc == 0, 0, 2 * pr * rc / (pr + rc))
    pr_p <- ifelse(tp_p + fp_p == 0, 0, tp_p / (tp_p + fp_p))
    rc_p <- ifelse(n_d == 0, 0, tp_p / n_d)
    f1_p <- ifelse(pr_p + rc_p == 0, 0, 2 * pr_p * rc_p / (pr_p + rc_p))
    truth_cols <- data.table(n_true_de = n_d, TP_de_fdr5 = tp, FP_de_fdr5 = fp,
                              prec_de = round(pr, 4), recall_de = round(rc, 4),
                              F1_de = round(f1, 4),
                              TP_de_post08 = tp_p, FP_de_post08 = fp_p,
                              prec_de_post08 = round(pr_p, 4),
                              recall_de_post08 = round(rc_p, 4),
                              F1_de_post08 = round(f1_p, 4))
    cat(sprintf("  GT DE: TP=%d FP=%d F1=%.4f\n", tp, fp, f1))
    NULL
  },
  "FALSE" = NULL
)
metrics <- cbind(metrics, truth_cols)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
