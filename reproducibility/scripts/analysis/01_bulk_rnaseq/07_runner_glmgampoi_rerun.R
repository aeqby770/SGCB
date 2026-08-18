# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(glmGamPoi)
library(pROC)
CONFIG <- list(
  SEED        = 12345L,
  N_CORES     = parallel::detectCores(),
  DATA_DIR    = paste0(PROJECT_ROOT, "/benchmark/data"),
  BULK_DIR    = paste0(PROJECT_ROOT, "/benchmark/data/bulk"),
  RAW_DIR     = paste0(PROJECT_ROOT, "/benchmark/output"),
  PAPER_DIR   = paste0(PROJECT_ROOT, "/benchmark/output/paper_figures"),
  SEQC_FILE   = paste0(PROJECT_ROOT, "/benchmark/data/bulk/GSE49712_HTSeq.txt.gz")
)

cat(strrep("=", 70), "\n")
cat(" glmGamPoi-only (gene name fix)\n")
cat(" :", CONFIG$N_CORES, "\n")
cat(strrep("=", 70), "\n")

# =============================================================================
# =============================================================================

sim_files <- list.files(CONFIG$DATA_DIR, pattern = "^sim_.*\\.rds$", full.names = TRUE)
cat("\n", length(sim_files), " simulation \n")

sim_results <- lapply(sim_files, function(sim_file) {
  cat("  ", basename(sim_file), "... ")
  sim_name <- gsub("\\.rds$", "", basename(sim_file))
  sim <- readRDS(sim_file)

  counts <- sim$counts
  group  <- factor(sim$sample_info$group)

  truth <- character(0)
  is_truth <- "truth" %in% names(sim) && is.data.frame(sim$truth)
  is_de    <- "de_genes" %in% names(sim)
  truth <- data.table::fifelse(
    rep(is_truth, 1L), list(sim$truth$gene[sim$truth$is_de == TRUE]),
    data.table::fifelse(rep(is_de, 1L), list(sim$de_genes), list(character(0)))
  )[[1]]

  keep   <- rowSums(counts > 5) >= 2
  counts <- counts[keep, ]
  truth  <- intersect(truth, rownames(counts))

  t0 <- Sys.time()
  col_data <- data.frame(condition = group, row.names = colnames(counts))
  fit <- glm_gp(counts, design = ~condition, col_data = col_data,
                 size_factors = "normed_sum")
  contrast_str <- paste0("condition", levels(group)[2])
  res <- test_de(fit, contrast = contrast_str)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  dt <- data.table(
    gene     = res$name,
    pvalue   = res$pval,
    padj     = res$adj_pval,
    log2FC   = res$lfc,
    lfcSE    = NA_real_,
    baseMean = rowMeans(counts),
    stat     = res$f_statistic,
    method   = "glmGamPoi",
    time     = elapsed
  )

  dt_clean <- dt[!is.na(pvalue) & !is.na(padj)]

  metrics <- data.table(
    method  = "glmGamPoi",
    n_genes = nrow(dt_clean),
    time_sec = round(elapsed, 2)
  )

  fdr_levels <- c(0.01, 0.05, 0.1)
  lfc_cuts   <- c(0.5, 1, 2)

  metrics[, paste0("n_sig_fdr", fdr_levels * 100) :=
            lapply(fdr_levels, function(fdr) sum(dt_clean$padj < fdr, na.rm = TRUE))]
  metrics[, paste0("n_sig_lfc", lfc_cuts) :=
            lapply(lfc_cuts, function(lfc) sum(dt_clean$padj < 0.05 & abs(dt_clean$log2FC) > lfc, na.rm = TRUE))]

  metrics[, `:=`(
    p_median  = median(dt_clean$pvalue, na.rm = TRUE),
    p_mean    = mean(dt_clean$pvalue, na.rm = TRUE),
    p_lt_0.05 = mean(dt_clean$pvalue < 0.05, na.rm = TRUE),
    lfc_median = median(abs(dt_clean$log2FC), na.rm = TRUE),
    lfc_iqr    = IQR(dt_clean$log2FC, na.rm = TRUE)
  )]

  vapply(fdr_levels, function(fdr) {
    sig_genes <- dt_clean$gene[dt_clean$padj < fdr]
    TP  <- sum(sig_genes %in% truth)
    FP  <- sum(!(sig_genes %in% truth))
    FN  <- sum(!(truth %in% sig_genes))
    prec <- TP / max(1, TP + FP)
    rec  <- TP / max(1, length(truth))
    f1   <- 2 * prec * rec / max(0.001, prec + rec)
    fdr_a <- FP / max(1, TP + FP)

    sfx <- paste0("_fdr", fdr * 100)
    metrics[, paste0("TP",          sfx) := TP]
    metrics[, paste0("FP",          sfx) := FP]
    metrics[, paste0("precision",   sfx) := round(prec, 4)]
    metrics[, paste0("recall",      sfx) := round(rec, 4)]
    metrics[, paste0("F1",          sfx) := round(f1, 4)]
    metrics[, paste0("actual_FDR",  sfx) := round(fdr_a, 4)]
    0L
  }, integer(1))

  # AUC
  is_de_vec <- as.integer(dt_clean$gene %in% truth)
  n_pos <- sum(is_de_vec); n_neg <- length(is_de_vec) - n_pos
  auc_val <- NA_real_
  tryCatch({
    roc_obj <- pROC::roc(is_de_vec, 1 - dt_clean$pvalue, quiet = TRUE)
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 4)
  }, error = function(e) NULL)
  metrics[, AUC_ROC := auc_val]
  metrics[, AUC_PR  := auc_val]

  metrics[, dataset := sim_name]
  cat("TP=", metrics$TP_fdr5, " FP=", metrics$FP_fdr5, " F1=", metrics$F1_fdr5, "\n")
  metrics
})

sim_dt <- rbindlist(sim_results, fill = TRUE)
cat("\n===== Simulation =====\n")
print(sim_dt[, .(dataset, n_genes, TP_fdr5, FP_fdr5, F1_fdr5, actual_FDR_fdr5, AUC_ROC, time_sec)])

raw_sim_path <- file.path(CONFIG$RAW_DIR, "bulk_simulation_metrics.csv")
old_sim <- fread(raw_sim_path)
updated_sim <- rbind(old_sim[method != "glmGamPoi"], sim_dt, fill = TRUE)
fwrite(updated_sim, raw_sim_path)
cat("\n:", raw_sim_path, "\n")

combined_sim_path <- file.path(CONFIG$PAPER_DIR, "combined_simulation.csv")
old_combined <- fread(combined_sim_path)
updated_combined <- rbind(old_combined[method != "glmGamPoi"], sim_dt, fill = TRUE)
fwrite(updated_combined, combined_sim_path)
cat(":", combined_sim_path, "\n")

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" Part 2: SEQC glmGamPoi \n")
cat(strrep("=", 70), "\n")

library(DESeq2)
library(edgeR)
library(limma)

seqc_raw <- read.table(gzfile(CONFIG$SEQC_FILE), header = TRUE,
                        row.names = 1, check.names = FALSE)
counts_seqc <- as.matrix(seqc_raw)
storage.mode(counts_seqc) <- "integer"
group_seqc <- factor(ifelse(grepl("^A_", colnames(counts_seqc)), "A", "B"))
keep_gs <- rowSums(counts_seqc > 5) >= 2
counts_gs <- counts_seqc[keep_gs, ]
cat(" SEQC :", nrow(counts_gs), " :", ncol(counts_gs), "\n")

cd <- data.frame(condition = group_seqc, row.names = colnames(counts_gs))

dds_gs <- DESeqDataSetFromMatrix(counts_gs, cd, ~condition)
dds_gs <- estimateSizeFactors(dds_gs)
dds_gs <- estimateDispersions(dds_gs, quiet = TRUE)
dds_gs <- nbinomWaldTest(dds_gs)
res_gs <- results(dds_gs)
de_d <- rownames(res_gs)[which(res_gs$padj < 0.001 & abs(res_gs$log2FoldChange) > 1)]

dge_gs <- DGEList(counts = counts_gs, group = group_seqc)
dge_gs <- calcNormFactors(dge_gs)
des_gs <- model.matrix(~group_seqc)
dge_gs <- estimateDisp(dge_gs, des_gs)
fit_gs <- glmQLFit(dge_gs, des_gs)
test_gs <- glmQLFTest(fit_gs, coef = 2)
res_gs_e <- topTags(test_gs, n = Inf, sort.by = "none")$table
de_e <- rownames(res_gs_e)[which(res_gs_e$FDR < 0.001 & abs(res_gs_e$logFC) > 1)]

v_gs <- voom(dge_gs, des_gs, plot = FALSE)
fit_gs_l <- lmFit(v_gs, des_gs)
fit_gs_l <- eBayes(fit_gs_l)
res_gs_l <- topTable(fit_gs_l, coef = 2, number = Inf, sort.by = "none")
de_l <- rownames(res_gs_l)[which(res_gs_l$adj.P.Val < 0.001 & abs(res_gs_l$logFC) > 1)]

gold_genes <- Reduce(intersect, list(de_d, de_e, de_l))
cat("  Gold standard:", length(gold_genes), "genes\n")

cat("  Running glmGamPoi... ")
t0 <- Sys.time()
col_data_s <- data.frame(condition = group_seqc, row.names = colnames(counts_gs))
fit_s <- glm_gp(counts_gs, design = ~condition, col_data = col_data_s,
                 size_factors = "normed_sum")
contrast_s <- paste0("condition", levels(group_seqc)[2])
res_s <- test_de(fit_s, contrast = contrast_s)
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

dt_s <- data.table(
  gene     = res_s$name,
  pvalue   = res_s$pval,
  padj     = res_s$adj_pval,
  log2FC   = res_s$lfc,
  baseMean = rowMeans(counts_gs),
  method   = "glmGamPoi",
  time     = elapsed_s
)

dt_sc <- dt_s[!is.na(padj)]
sig <- dt_sc$gene[dt_sc$padj < 0.05]
TP_s  <- sum(sig %in% gold_genes)
FP_s  <- sum(!(sig %in% gold_genes))
prec_s <- TP_s / max(1, TP_s + FP_s)
rec_s  <- TP_s / max(1, length(gold_genes))
f1_s   <- 2 * prec_s * rec_s / max(0.001, prec_s + rec_s)
fdr_s  <- FP_s / max(1, TP_s + FP_s)

is_de_s <- as.integer(dt_sc$gene %in% gold_genes)
pv_s <- dt_sc$pvalue; pv_s[is.na(pv_s)] <- 1
roc_s <- roc(is_de_s, 1 - pv_s, quiet = TRUE)
auc_s <- round(as.numeric(auc(roc_s)), 4)

# LFC concordance
ref_lfc <- log2(rowMeans(counts_gs[, group_seqc == "A"]) + 1) -
           log2(rowMeans(counts_gs[, group_seqc == "B"]) + 1)
common <- intersect(names(ref_lfc), dt_sc$gene)
lfc_r <- round(cor(ref_lfc[common], setNames(dt_sc$log2FC, dt_sc$gene)[common],
                   use = "complete.obs"), 4)

seqc_new <- data.table(
  TP = TP_s, FP = FP_s,
  precision = round(prec_s, 4), recall = round(rec_s, 4),
  F1 = round(f1_s, 4), actual_FDR = round(fdr_s, 4),
  n_sig = length(sig), method = "glmGamPoi",
  AUC = auc_s, LFC_r = lfc_r
)
cat("OK  AUC=", auc_s, " F1=", round(f1_s, 4), " TP=", TP_s, " FP=", FP_s, "\n")

seqc_path <- file.path(CONFIG$PAPER_DIR, "test4_seqc_metrics.csv")
old_seqc <- fread(seqc_path)
updated_seqc <- rbind(old_seqc[method != "glmGamPoi"], seqc_new, fill = TRUE)
fwrite(updated_seqc, seqc_path)
cat(":", seqc_path, "\n")

seqc8_path <- file.path(CONFIG$PAPER_DIR, "test4_seqc_metrics_8methods.csv")
old_seqc8 <- fread(seqc8_path)
updated_seqc8 <- rbind(old_seqc8[method != "glmGamPoi"], seqc_new, fill = TRUE)
fwrite(updated_seqc8, seqc8_path)
cat(":", seqc8_path, "\n")

roc_df <- data.table(
  sensitivity = roc_s$sensitivities,
  specificity = 1 - roc_s$specificities,
  method = "glmGamPoi"
)
roc8_path <- file.path(CONFIG$PAPER_DIR, "test4_seqc_roc_8methods.csv.gz")
old_roc8 <- fread(roc8_path)
updated_roc8 <- rbind(old_roc8[method != "glmGamPoi"], roc_df, fill = TRUE)
fwrite(updated_roc8, roc8_path)
cat(":", roc8_path, "\n")

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" glmGamPoi \n")
cat(" :\n")
cat("    - ", raw_sim_path, "\n")
cat("    - ", combined_sim_path, "\n")
cat("    - ", seqc_path, "\n")
cat("    - ", seqc8_path, "\n")
cat("    - ", roc8_path, "\n")
cat(strrep("=", 70), "\n")
