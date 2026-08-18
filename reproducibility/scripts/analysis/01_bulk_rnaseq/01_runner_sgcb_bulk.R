# =============================================================================
#   test_type: classic | simulation | gtex_pair | tcga | pseudobulk | null_fdr |
#              seqc | dv
#   dataset:   bottomly | sim_n10_pDE10_ef2.0 | liver_vs_kidney | brca | kang18 |
#              1 | seqc | cerebellum_vs_cortex | cerehemi_vs_cortex | hippo_vs_hypo
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

args       <- commandArgs(trailingOnly = TRUE)
TEST_TYPE  <- args[1]
DATASET    <- args[2]

BASE       <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq")
DATA_ROOT  <- file.path(BASE, "data/bulk")
OUT_DIR    <- file.path(BASE, "sgcb/results")
N_SUB      <- 50L
SGCB_ARGS  <- list(
  prior_de = "auto",
  prior_dv = "auto",
  prior_dg = "auto",
  prior_sd_de = 1.0,
  prior_sd_dv = 0.5,
  prior_sd_dg_alpha = 0.5,
  prior_sd_dg_gamma = 0.5
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_CSV  <- file.path(OUT_DIR, paste0(TEST_TYPE, "__", DATASET, ".csv"))
OUT_RDS  <- file.path(OUT_DIR, paste0(TEST_TYPE, "__", DATASET, "_full.rds"))

library(SGCB)
library(data.table)

cat(sprintf("[%s] ===== SGCB bulk: %s / %s =====\n", Sys.time(), TEST_TYPE, DATASET))

loaded <- switch(TEST_TYPE,

  "classic" = {
    ds <- readRDS(file.path(DATA_ROOT, "01_classic", paste0("real_", DATASET, ".rds")))
    list(counts = as.matrix(ds$counts),
         group  = factor(ds$sample_info$group),
         truth  = NULL)
  },

  "simulation" = {
    ds <- readRDS(file.path(DATA_ROOT, "02_simulation", paste0(DATASET, ".rds")))
    list(counts = as.matrix(ds[["counts"]]),
         group  = factor(ds[["sample_info"]][["group"]]),
         truth  = ds[["truth"]])
  },

  "gtex_pair" = {
    tt <- strsplit(DATASET, "_vs_")[[1]]
    d1 <- readRDS(file.path(DATA_ROOT, "03_gtex", paste0("gtex_", tt[1], ".rds")))
    d2 <- readRDS(file.path(DATA_ROOT, "03_gtex", paste0("gtex_", tt[2], ".rds")))
    cg <- intersect(rownames(d1$counts), rownames(d2$counts))
    set.seed(12345)
    i1 <- sample.int(ncol(d1$counts), min(N_SUB, ncol(d1$counts)))
    i2 <- sample.int(ncol(d2$counts), min(N_SUB, ncol(d2$counts)))
    cm <- cbind(d1$counts[cg, i1], d2$counts[cg, i2])
    gm <- factor(c(rep(tt[1], length(i1)), rep(tt[2], length(i2))))
    list(counts = as.matrix(cm), group = gm, truth = NULL)
  },

  "tcga" = {
    ds <- readRDS(file.path(DATA_ROOT, "04_tcga",
                            paste0("tcga_", DATASET, "_tumor_vs_normal.rds")))
    list(counts = as.matrix(ds$counts),
         group  = factor(ds$group),
         truth  = NULL)
  },

  "pseudobulk" = {
    flist <- list.files(file.path(DATA_ROOT, "05_pseudobulk"),
                        pattern = DATASET, full.names = TRUE)
    ds <- readRDS(flist[1])
    cm <- as.matrix(ds$counts)
    rownames(cm) <- make.unique(rownames(cm))
    list(counts = cm, group = factor(ds$group), truth = NULL)
  },

  "null_fdr" = {
    ds <- readRDS(file.path(DATA_ROOT, "03_gtex", "gtex_liver.rds"))
    seed_val <- as.integer(DATASET)
    set.seed(seed_val)
    idx <- sample.int(ncol(ds$counts), min(2L * N_SUB, ncol(ds$counts)))
    cs  <- ds$counts[, idx]
    nh  <- length(idx) %/% 2L
    gp  <- factor(sample(c(rep("A", nh), rep("B", length(idx) - nh))))
    list(counts = as.matrix(cs), group = gp, truth = NULL)
  },

  "seqc" = {
    seqc_file <- file.path(DATA_ROOT, "06_seqc", "GSE49712_HTSeq.txt.gz")
    seqc_raw  <- read.table(gzfile(seqc_file), header = TRUE, row.names = 1,
                            check.names = FALSE)
    cts <- as.matrix(seqc_raw)
    storage.mode(cts) <- "integer"
    gp  <- factor(ifelse(grepl("^A_", colnames(cts)), "A", "B"))

    library(DESeq2); library(edgeR); library(limma)
    keep_gs <- rowSums(cts > 5) >= 2
    cts_gs  <- cts[keep_gs, ]
    des     <- model.matrix(~gp)

    # DESeq2
    dds <- DESeqDataSetFromMatrix(cts_gs, data.frame(gp = gp, row.names = colnames(cts_gs)), ~gp)
    dds <- DESeq(dds, quiet = TRUE)
    rd  <- results(dds)
    de_d <- rownames(rd)[which(rd$padj < 0.001 & abs(rd$log2FoldChange) > 1)]

    # edgeR
    dge <- DGEList(counts = cts_gs, group = gp)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, des)
    fe  <- glmQLFit(dge, des)
    re  <- topTags(glmQLFTest(fe, coef = 2), n = Inf, sort.by = "none")$table
    de_e <- rownames(re)[which(re$FDR < 0.001 & abs(re$logFC) > 1)]

    # limma-voom
    v   <- voom(dge, des, plot = FALSE)
    fl  <- eBayes(lmFit(v, des))
    rl  <- topTable(fl, coef = 2, number = Inf, sort.by = "none")
    de_l <- rownames(rl)[which(rl$adj.P.Val < 0.001 & abs(rl$logFC) > 1)]

    gold <- Reduce(intersect, list(de_d, de_e, de_l))
    cat(sprintf("  SEQC gold standard: %d genes (D=%d, E=%d, L=%d)\n",
                length(gold), length(de_d), length(de_e), length(de_l)))

    truth_df <- data.frame(gene = rownames(cts), is_de = rownames(cts) %in% gold,
                           stringsAsFactors = FALSE)
    list(counts = cts, group = gp, truth = truth_df)
  },

  "dv" = {
    dv_file <- file.path(DATA_ROOT, "07_gtex_brain_dv",
                          paste0("gtex_brain_", DATASET, ".rds"))
    ds <- readRDS(dv_file)
    set.seed(12345)
    lvs  <- levels(ds$group)
    idx1 <- which(ds$group == lvs[1])
    idx2 <- which(ds$group == lvs[2])
    i1   <- idx1[sample.int(length(idx1), min(N_SUB, length(idx1)))]
    i2   <- idx2[sample.int(length(idx2), min(N_SUB, length(idx2)))]
    cm   <- as.matrix(ds$counts[, c(i1, i2)])
    gp   <- factor(c(rep(lvs[1], length(i1)), rep(lvs[2], length(i2))))
    list(counts = cm, group = gp, truth = NULL)
  }
)

counts <- loaded$counts
group  <- loaded$group
truth  <- loaded$truth

cat(sprintf("  Data: %d genes x %d samples, groups: %s\n",
            nrow(counts), ncol(counts),
            paste(paste0(levels(group), "=", table(group)), collapse=", ")))

t0 <- proc.time()
res <- do.call(sgcbDE, c(list(counts = counts, group = group), SGCB_ARGS))
time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)

rdf  <- as.data.frame(res)
cat(sprintf("  Done: %.2f sec, %d genes tested\n", time_sec, nrow(rdf)))

padj <- fifelse(is.na(rdf$padj), 1, rdf$padj)
pval <- fifelse(is.na(rdf$pvalue), 1, rdf$pvalue)
lfc  <- fifelse(is.na(rdf$log2FoldChange), 0, rdf$log2FoldChange)
p_de_post <- fifelse(is.na(rdf$p_de_post), 0, rdf$p_de_post)
p_de_only_post <- fifelse(is.na(rdf$p_de_only_post), 0, rdf$p_de_only_post)
p_dv_post <- fifelse(is.na(rdf$p_dv_post), 0, rdf$p_dv_post)
p_dv_only_post <- fifelse(is.na(rdf$p_dv_only_post), 0, rdf$p_dv_only_post)
p_dg_post <- fifelse(is.na(rdf$p_dg_post), 0, rdf$p_dg_post)
p_dg_only_post <- fifelse(is.na(rdf$p_dg_only_post), 0, rdf$p_dg_only_post)
de_post_call <- fifelse(is.na(rdf$de_post_call), FALSE, rdf$de_post_call)
dv_post_call <- fifelse(is.na(rdf$dv_post_call), FALSE, rdf$dv_post_call)
dg_post_call <- fifelse(is.na(rdf$dg_post_call), FALSE, rdf$dg_post_call)
de_fdr_call <- fifelse(is.na(rdf$de_fdr_call), FALSE, rdf$de_fdr_call)
dv_fdr_call <- fifelse(is.na(rdf$dv_fdr_call), FALSE, rdf$dv_fdr_call)
dg_fdr_call <- fifelse(is.na(rdf$dg_fdr_call), FALSE, rdf$dg_fdr_call)

DS_PREFIX <- c(gtex_pair = "GTEx_", tcga = "TCGA_", pseudobulk = "PB_")
prefix    <- DS_PREFIX[TEST_TYPE]
prefix[is.na(prefix)] <- ""
ds_label  <- paste0(prefix, DATASET)

metrics <- data.table(
  method       = "SGCB",
  n_genes      = nrow(rdf),
  time_sec     = time_sec,
  n_sig_fdr1   = sum(padj < 0.01),
  n_sig_fdr5   = sum(padj < 0.05),
  n_sig_fdr10  = sum(padj < 0.10),
  n_sig_lfc0.5 = sum(padj < 0.05 & abs(lfc) > 0.5),
  n_sig_lfc1   = sum(padj < 0.05 & abs(lfc) > 1.0),
  n_sig_lfc2   = sum(padj < 0.05 & abs(lfc) > 2.0),
  p_median     = median(pval),
  p_mean       = mean(pval),
  p_lt_0.05    = mean(pval < 0.05),
  lfc_median   = median(abs(lfc)),
  lfc_iqr      = IQR(abs(lfc)),
  n_post_de_0.5 = sum(p_de_post > 0.5),
  n_post_de_0.8 = sum(p_de_post > 0.8),
  n_post_de_0.9 = sum(p_de_post > 0.9),
  n_post_deonly_0.5 = sum(p_de_only_post > 0.5),
  n_post_deonly_0.8 = sum(p_de_only_post > 0.8),
  n_post_deonly_0.9 = sum(p_de_only_post > 0.9),
  n_post_de_call = sum(de_post_call),
  n_de_fdr_call = sum(de_fdr_call),
  n_post_dv_0.8 = sum(p_dv_post > 0.8),
  n_post_dvonly_0.8 = sum(p_dv_only_post > 0.8),
  n_post_dv_call = sum(dv_post_call),
  n_dv_fdr_call = sum(dv_fdr_call),
  n_post_dg_0.8 = sum(p_dg_post > 0.8),
  n_post_dgonly_0.8 = sum(p_dg_only_post > 0.8),
  n_post_dg_call = sum(dg_post_call),
  n_dg_fdr_call = sum(dg_fdr_call),
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
  dataset      = ds_label
)

has_truth <- !is.null(truth) && "is_de" %in% names(truth)

truth_cols <- data.table(
  TP_fdr1 = NA_integer_, FP_fdr1 = NA_integer_,
  precision_fdr1 = NA_real_, recall_fdr1 = NA_real_,
  F1_fdr1 = NA_real_, actual_FDR_fdr1 = NA_real_,
  TP_fdr5 = NA_integer_, FP_fdr5 = NA_integer_,
  precision_fdr5 = NA_real_, recall_fdr5 = NA_real_,
  F1_fdr5 = NA_real_, actual_FDR_fdr5 = NA_real_,
  TP_fdr10 = NA_integer_, FP_fdr10 = NA_integer_,
  precision_fdr10 = NA_real_, recall_fdr10 = NA_real_,
  F1_fdr10 = NA_real_, actual_FDR_fdr10 = NA_real_,
  TP_post_deonly_08 = NA_integer_, FP_post_deonly_08 = NA_integer_,
  precision_post_deonly_08 = NA_real_, recall_post_deonly_08 = NA_real_,
  F1_post_deonly_08 = NA_real_, actual_FDR_post_deonly_08 = NA_real_,
  AUC_ROC = NA_real_, AUC_PR = NA_real_
)

switch(as.character(has_truth),
  "TRUE" = {
    gene_ids  <- rdf$gene_id
    truth_dt  <- as.data.table(truth)
    is_de_vec <- truth_dt[match(gene_ids, truth_dt$gene), is_de]
    is_de_vec[is.na(is_de_vec)] <- FALSE
    # Fix: use filtered truth size as denominator so recall/F1/AUC are
    # comparable with other_methods (which score only over the filtered gene
    # list). Previously this used sum(truth_dt$is_de) [unfiltered], which
    # systematically under-estimated SGCB's recall when SGCB's internal
    # min_count filter dropped low-count true-DE genes.
    n_true_de <- sum(is_de_vec)

    sig01  <- padj < 0.01
    sig05  <- padj < 0.05
    sig10  <- padj < 0.10

    tp01 <- sum(sig01 & is_de_vec);  fp01 <- sum(sig01 & !is_de_vec)
    tp05 <- sum(sig05 & is_de_vec);  fp05 <- sum(sig05 & !is_de_vec)
    tp10 <- sum(sig10 & is_de_vec);  fp10 <- sum(sig10 & !is_de_vec)

    safe_div <- \(a, b) fifelse(b == 0, 0, a / b)

    prec01 <- safe_div(tp01, tp01 + fp01)
    prec05 <- safe_div(tp05, tp05 + fp05)
    prec10 <- safe_div(tp10, tp10 + fp10)
    rec01  <- safe_div(tp01, n_true_de)
    rec05  <- safe_div(tp05, n_true_de)
    rec10  <- safe_div(tp10, n_true_de)
    f1_01  <- safe_div(2 * prec01 * rec01, prec01 + rec01)
    f1_05  <- safe_div(2 * prec05 * rec05, prec05 + rec05)
    f1_10  <- safe_div(2 * prec10 * rec10, prec10 + rec10)
    fdr01  <- safe_div(fp01, tp01 + fp01)
    fdr05  <- safe_div(fp05, tp05 + fp05)
    fdr10  <- safe_div(fp10, tp10 + fp10)

    sig_post08 <- p_de_only_post > 0.8
    tp_p08 <- sum(sig_post08 & is_de_vec)
    fp_p08 <- sum(sig_post08 & !is_de_vec)
    prec_p08 <- safe_div(tp_p08, tp_p08 + fp_p08)
    rec_p08 <- safe_div(tp_p08, n_true_de)
    f1_p08 <- safe_div(2 * prec_p08 * rec_p08, prec_p08 + rec_p08)
    fdr_p08 <- safe_div(fp_p08, tp_p08 + fp_p08)

    ord        <- order(pval)
    is_de_ord  <- is_de_vec[ord]
    n_neg      <- sum(!is_de_vec)
    tpr_curve  <- cumsum(is_de_ord) / n_true_de
    fpr_curve  <- cumsum(!is_de_ord) / n_neg
    auc_roc    <- sum(diff(fpr_curve) * (head(tpr_curve, -1) + tail(tpr_curve, -1)) / 2)

    prec_curve <- cumsum(is_de_ord) / seq_along(is_de_ord)
    rec_curve  <- cumsum(is_de_ord) / n_true_de
    auc_pr     <- sum(diff(rec_curve) * (head(prec_curve, -1) + tail(prec_curve, -1)) / 2)

    truth_cols <- data.table(
      TP_fdr1 = tp01, FP_fdr1 = fp01, precision_fdr1 = round(prec01, 4),
      recall_fdr1 = round(rec01, 4), F1_fdr1 = round(f1_01, 4),
      actual_FDR_fdr1 = round(fdr01, 4),
      TP_fdr5 = tp05, FP_fdr5 = fp05, precision_fdr5 = round(prec05, 4),
      recall_fdr5 = round(rec05, 4), F1_fdr5 = round(f1_05, 4),
      actual_FDR_fdr5 = round(fdr05, 4),
      TP_fdr10 = tp10, FP_fdr10 = fp10, precision_fdr10 = round(prec10, 4),
      recall_fdr10 = round(rec10, 4), F1_fdr10 = round(f1_10, 4),
      actual_FDR_fdr10 = round(fdr10, 4),
      TP_post_deonly_08 = tp_p08, FP_post_deonly_08 = fp_p08,
      precision_post_deonly_08 = round(prec_p08, 4),
      recall_post_deonly_08 = round(rec_p08, 4),
      F1_post_deonly_08 = round(f1_p08, 4),
      actual_FDR_post_deonly_08 = round(fdr_p08, 4),
      AUC_ROC = round(auc_roc, 4), AUC_PR = round(auc_pr, 4)
    )

    cat(sprintf("  Simulation: TP@5%%=%d, FP@5%%=%d, F1@5%%=%.4f, AUC=%.4f\n",
                tp05, fp05, f1_05, auc_roc))
    NULL
  },
  "FALSE" = NULL
)

out <- cbind(metrics, truth_cols)
fwrite(out, OUT_CSV)
saveRDS(rdf, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
