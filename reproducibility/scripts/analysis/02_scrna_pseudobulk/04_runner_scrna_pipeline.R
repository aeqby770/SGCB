# =============================================================================
#   method:  SGCB | DESeq2 | edgeR | limma | glmGamPoi | Seurat
#   dataset: kang18 | segerstolpe | xin | zilionis | jakel | muscat_sim
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = 1L)
data.table::setDTthreads(1L)

args    <- commandArgs(trailingOnly = TRUE)
METHOD  <- args[1]
DATASET <- args[2]
NULL_SEED <- as.integer(args[3])

BASE     <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk")
PB_DIR   <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/05_pseudobulk")
OUT_DIR  <- file.path(BASE, "results")
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

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== scRNA PB DE: %s / %s %s=====\n",
            Sys.time(), METHOD, DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

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

# Ground truth (muscat_sim only)
sim_truth <- ds$truth  # NULL for real datasets

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

t0 <- proc.time()

res_df <- switch(METHOD,

  "SGCB" = {
    library(SGCB)
    res <- do.call(sgcbDE, c(list(counts = counts, group = group), SGCB_ARGS))
    rr  <- as.data.frame(res)
    data.table(
      gene           = rr$gene_id,
      pvalue         = rr$pvalue,
      padj           = rr$padj,
      log2FoldChange = rr$log2FoldChange,
      log2FC_shrunk  = rr$log2FC_shrunk,
      p_de_post      = rr$p_de_post,
      p_de_only_post = rr$p_de_only_post,
      de_fdr_call    = rr$de_fdr_call,
      p_dv_post      = rr$p_dv_post,
      p_dv_only_post = rr$p_dv_only_post,
      dv_fdr_call    = rr$dv_fdr_call,
      p_dg_post      = rr$p_dg_post,
      p_dg_only_post = rr$p_dg_only_post,
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
  },

  # ----- DESeq2 -----
  "DESeq2" = {
    library(DESeq2)
    coldata <- data.frame(group = group, row.names = colnames(counts))
    dds <- DESeqDataSetFromMatrix(counts, coldata, ~group)
    dds <- DESeq(dds, quiet = TRUE)
    rr  <- results(dds)
    data.table(gene = rownames(rr), pvalue = rr$pvalue, padj = rr$padj,
               log2FoldChange = rr$log2FoldChange)
  },

  # ----- edgeR -----
  "edgeR" = {
    library(edgeR)
    design <- model.matrix(~group)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    tt  <- glmQLFTest(fit, coef = 2)
    rr  <- topTags(tt, n = Inf, sort.by = "none")$table
    data.table(gene = rownames(rr), pvalue = rr$PValue, padj = rr$FDR,
               log2FoldChange = rr$logFC)
  },

  # ----- limma-voom -----
  "limma" = {
    library(limma)
    library(edgeR)
    design <- model.matrix(~group)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    rr  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  },

  # ----- glmGamPoi -----
  "glmGamPoi" = {
    library(glmGamPoi)
    fit <- glm_gp(counts, design = ~group,
                   col_data = data.frame(group = group,
                                         row.names = colnames(counts)))
    rr  <- test_de(fit, contrast = colnames(fit$Beta)[2])
    data.table(gene = rr$name, pvalue = rr$pval, padj = rr$adj_pval,
               log2FoldChange = rr$lfc)
  },

  # ----- Seurat (Wilcoxon on pseudobulk) -----
  "Seurat" = {
    library(Seurat)
    obj <- CreateSeuratObject(counts = counts)
    obj$group <- group
    Idents(obj) <- "group"
    obj <- NormalizeData(obj, verbose = FALSE)
    markers <- FindMarkers(obj, ident.1 = levels(group)[2],
                           ident.2 = levels(group)[1],
                           test.use = "wilcox",
                           logfc.threshold = 0, min.pct = 0)
    data.table(gene = rownames(markers), pvalue = markers$p_val,
               padj = markers$p_val_adj, log2FoldChange = markers$avg_log2FC)
  }
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)

cat(sprintf("  %s done: %.2f sec, %d genes\n", METHOD, time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]
has_dv <- "dv_pvalue" %in% names(res_df)
has_dd <- "pvalue_dd" %in% names(res_df)
has_score <- "SGCB_Score_padj" %in% names(res_df)
has_post <- "p_de_post" %in% names(res_df)

metrics <- data.table(
  method       = METHOD,
  n_genes      = nrow(res_df),
  time_sec     = time_sec,
  n_sig_fdr1   = sum(res_df$padj < 0.01),
  n_sig_fdr5   = sum(res_df$padj < 0.05),
  n_sig_fdr10  = sum(res_df$padj < 0.10),
  n_sig_lfc0.5 = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5),
  n_sig_lfc1   = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0),
  p_median     = median(res_df$pvalue),
  p_mean       = mean(res_df$pvalue),
  p_lt_0.05    = mean(res_df$pvalue < 0.05),
  lfc_median   = median(abs(res_df$log2FoldChange)),
  n_post_de_0.5 = ifelse(has_post, sum(res_df$p_de_post > 0.5, na.rm = TRUE), NA_integer_),
  n_post_de_0.8 = ifelse(has_post, sum(res_df$p_de_post > 0.8, na.rm = TRUE), NA_integer_),
  n_post_deonly_0.8 = ifelse(has_post, sum(res_df$p_de_only_post > 0.8, na.rm = TRUE), NA_integer_),
  n_de_fdr_call = ifelse("de_fdr_call" %in% names(res_df), sum(res_df$de_fdr_call, na.rm = TRUE), NA_integer_),
  n_post_dv_0.8 = ifelse(has_post, sum(res_df$p_dv_post > 0.8, na.rm = TRUE), NA_integer_),
  n_post_dvonly_0.8 = ifelse(has_post, sum(res_df$p_dv_only_post > 0.8, na.rm = TRUE), NA_integer_),
  n_dv_fdr_call = ifelse("dv_fdr_call" %in% names(res_df), sum(res_df$dv_fdr_call, na.rm = TRUE), NA_integer_),
  n_post_dg_0.8 = ifelse(has_post, sum(res_df$p_dg_post > 0.8, na.rm = TRUE), NA_integer_),
  n_post_dgonly_0.8 = ifelse(has_post, sum(res_df$p_dg_only_post > 0.8, na.rm = TRUE), NA_integer_),
  n_dg_fdr_call = ifelse("dg_fdr_call" %in% names(res_df), sum(res_df$dg_fdr_call, na.rm = TRUE), NA_integer_),
  prior_de_used = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_de_used")), NA_real_),
  prior_dv_used = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_dv_used")), NA_real_),
  prior_dv_upper = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_dv_upper")), NA_real_),
  prior_dg_used = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_dg_used")), NA_real_),
  dataset      = paste0("PB_", DATASET),
  is_null      = IS_NULL,
  n_dv_fdr5    = ifelse(has_dv, sum(res_df$dv_padj < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_model_raw_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_model_raw < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_model_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_model < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_bf_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_bf < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_mix_fdr5 = ifelse(has_dv, sum(res_df$dv_padj_mix < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_fdr10   = ifelse(has_dv, sum(res_df$dv_padj < 0.10, na.rm = TRUE), NA_integer_),
  dv_p_median  = ifelse(has_dv, median(res_df$dv_pvalue, na.rm = TRUE), NA_real_),
  dv_p_lt_0.05 = ifelse(has_dv, mean(res_df$dv_pvalue < 0.05, na.rm = TRUE), NA_real_),
  n_dd_fdr5    = ifelse(has_dd, sum(res_df$padj_dd < 0.05, na.rm = TRUE), NA_integer_),
  n_dd_fdr10   = ifelse(has_dd, sum(res_df$padj_dd < 0.10, na.rm = TRUE), NA_integer_),
  n_score_fdr5 = ifelse(has_score, sum(res_df$SGCB_Score_padj < 0.05, na.rm = TRUE), NA_integer_),
  n_score_fdr10 = ifelse(has_score, sum(res_df$SGCB_Score_padj < 0.10, na.rm = TRUE), NA_integer_),
  score_median = ifelse(has_score, median(res_df$SGCB_Score, na.rm = TRUE), NA_real_)
)

has_truth_sim <- !is.null(sim_truth) && !IS_NULL

truth_cols <- data.table(
  n_true_de = NA_integer_, n_true_dv = NA_integer_,
  TP_de_fdr5 = NA_integer_, FP_de_fdr5 = NA_integer_,
  prec_de = NA_real_, recall_de = NA_real_, F1_de = NA_real_,
  TP_de_post08 = NA_integer_, FP_de_post08 = NA_integer_,
  prec_de_post08 = NA_real_, recall_de_post08 = NA_real_, F1_de_post08 = NA_real_,
  TP_dv_fdr5 = NA_integer_, FP_dv_fdr5 = NA_integer_,
  prec_dv = NA_real_, recall_dv = NA_real_, F1_dv = NA_real_
)

switch(as.character(has_truth_sim),
  "TRUE" = {
    matched_truth <- sim_truth[match(res_df$gene, sim_truth$gene), ]
    truth_de <- ifelse(is.na(matched_truth$is_de), FALSE, matched_truth$is_de)
    truth_dv <- ifelse(is.na(matched_truth$is_dv), FALSE, matched_truth$is_dv)
    sig_de <- res_df$padj < 0.05
    sig_de_post <- ifelse(has_post, res_df$p_de_only_post > 0.8, FALSE)
    tp_de <- sum(sig_de & truth_de)
    fp_de <- sum(sig_de & !truth_de)
    n_de  <- sum(truth_de)
    prec_de <- ifelse(tp_de + fp_de == 0, 0, tp_de / (tp_de + fp_de))
    rec_de  <- ifelse(n_de == 0, 0, tp_de / n_de)
    f1_de   <- ifelse(prec_de + rec_de == 0, 0, 2 * prec_de * rec_de / (prec_de + rec_de))

    tp_de_post <- sum(sig_de_post & truth_de, na.rm = TRUE)
    fp_de_post <- sum(sig_de_post & !truth_de, na.rm = TRUE)
    prec_de_post <- ifelse(tp_de_post + fp_de_post == 0, 0, tp_de_post / (tp_de_post + fp_de_post))
    rec_de_post  <- ifelse(n_de == 0, 0, tp_de_post / n_de)
    f1_de_post   <- ifelse(prec_de_post + rec_de_post == 0, 0, 2 * prec_de_post * rec_de_post / (prec_de_post + rec_de_post))

    n_dv <- sum(truth_dv)
    tp_dv <- NA_integer_; fp_dv <- NA_integer_
    prec_dv_v <- NA_real_; rec_dv_v <- NA_real_; f1_dv_v <- NA_real_
    switch(as.character(has_dv),
      "TRUE" = {
        sig_dv <- res_df$dv_padj < 0.05
        tp_dv  <- sum(sig_dv & truth_dv, na.rm = TRUE)
        fp_dv  <- sum(sig_dv & !truth_dv, na.rm = TRUE)
        prec_dv_v <- ifelse(tp_dv + fp_dv == 0, 0, tp_dv / (tp_dv + fp_dv))
        rec_dv_v  <- ifelse(n_dv == 0, 0, tp_dv / n_dv)
        f1_dv_v   <- ifelse(prec_dv_v + rec_dv_v == 0, 0, 2 * prec_dv_v * rec_dv_v / (prec_dv_v + rec_dv_v))
        NULL
      },
      "FALSE" = NULL
    )

    truth_cols <- data.table(
      n_true_de = n_de, n_true_dv = n_dv,
      TP_de_fdr5 = tp_de, FP_de_fdr5 = fp_de,
      prec_de = round(prec_de, 4), recall_de = round(rec_de, 4), F1_de = round(f1_de, 4),
      TP_de_post08 = tp_de_post, FP_de_post08 = fp_de_post,
      prec_de_post08 = round(prec_de_post, 4), recall_de_post08 = round(rec_de_post, 4), F1_de_post08 = round(f1_de_post, 4),
      TP_dv_fdr5 = tp_dv, FP_dv_fdr5 = fp_dv,
      prec_dv = round(prec_dv_v, 4), recall_dv = round(rec_dv_v, 4), F1_dv = round(f1_dv_v, 4)
    )
    cat(sprintf("  Ground truth DE: TP=%d FP=%d F1=%.4f\n", tp_de, fp_de, f1_de))
    NULL
  },
  "FALSE" = NULL
)

metrics <- cbind(metrics, truth_cols)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
