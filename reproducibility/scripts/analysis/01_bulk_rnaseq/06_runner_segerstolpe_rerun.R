PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  BULK_DIR = paste0(PROJECT_ROOT, "/benchmark/data/bulk"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/rerun_large_missing"),
  DATASET = "PB_segerstolpe"
)

dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(edgeR)
library(limma)
library(NOISeq)
library(glmGamPoi)

pb <- readRDS(file.path(CONFIG$BULK_DIR, "segerstolpe_t2d_pseudobulk.rds"))
counts <- as.matrix(pb$counts)
group <- pb$group
keep <- rowSums(counts > 5) >= 2
counts <- counts[keep, ]

counts_limma <- counts
group_limma <- group
t0_limma <- Sys.time()
design_limma <- model.matrix(~group_limma)
dge_limma <- DGEList(counts = counts_limma, group = group_limma)
dge_limma <- calcNormFactors(dge_limma)
v_limma <- voom(dge_limma, design_limma, plot = FALSE)
fit_limma <- lmFit(v_limma, design_limma)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = 2, number = Inf, sort.by = "none")
time_sec_limma <- as.numeric(difftime(Sys.time(), t0_limma, units = "secs"))
dt_limma <- data.table(
  gene = rownames(res_limma),
  pvalue = res_limma$P.Value,
  padj = res_limma$adj.P.Val,
  log2FC = res_limma$logFC,
  lfcSE = sqrt(fit_limma$s2.post) * fit_limma$stdev.unscaled[, 2],
  baseMean = rowMeans(counts_limma),
  stat = res_limma$t,
  method = "limma",
  time = time_sec_limma
)
dt_limma <- dt_limma[!is.na(pvalue) & !is.na(padj)]
metrics_limma <- data.table(
  method = "limma",
  n_genes = nrow(dt_limma),
  time_sec = round(time_sec_limma, 2),
  n_sig_fdr1 = sum(dt_limma$padj < 0.01, na.rm = TRUE),
  n_sig_fdr5 = sum(dt_limma$padj < 0.05, na.rm = TRUE),
  n_sig_fdr10 = sum(dt_limma$padj < 0.10, na.rm = TRUE),
  n_sig_lfc0.5 = sum(dt_limma$padj < 0.05 & abs(dt_limma$log2FC) > 0.5, na.rm = TRUE),
  n_sig_lfc1 = sum(dt_limma$padj < 0.05 & abs(dt_limma$log2FC) > 1, na.rm = TRUE),
  n_sig_lfc2 = sum(dt_limma$padj < 0.05 & abs(dt_limma$log2FC) > 2, na.rm = TRUE),
  p_median = median(dt_limma$pvalue, na.rm = TRUE),
  p_mean = mean(dt_limma$pvalue, na.rm = TRUE),
  p_lt_0.05 = mean(dt_limma$pvalue < 0.05, na.rm = TRUE),
  lfc_median = median(abs(dt_limma$log2FC), na.rm = TRUE),
  lfc_iqr = IQR(dt_limma$log2FC, na.rm = TRUE),
  dataset = CONFIG$DATASET
)

counts_noiseq <- counts
group_noiseq <- group
t0_noiseq <- Sys.time()
rownames(counts_noiseq) <- make.unique(as.character(rownames(counts_noiseq)))
factors_noiseq <- data.frame(condition = group_noiseq, row.names = colnames(counts_noiseq))
nsdata_noiseq <- NOISeq::readData(data = counts_noiseq, factors = factors_noiseq)
res_noiseq <- NOISeq::noiseqbio(nsdata_noiseq, k = 0.5, norm = "tmm", nclust = 15, factor = "condition", conditions = levels(factor(group_noiseq)), r = 20, filter = 0)
de_table_noiseq <- res_noiseq@results[[1]]
prob_noiseq <- de_table_noiseq$prob
pval_noiseq <- 1 - prob_noiseq
time_sec_noiseq <- as.numeric(difftime(Sys.time(), t0_noiseq, units = "secs"))
dt_noiseq <- data.table(
  gene = rownames(de_table_noiseq),
  pvalue = pval_noiseq,
  padj = p.adjust(pval_noiseq, method = "BH"),
  log2FC = de_table_noiseq$log2FC,
  lfcSE = NA_real_,
  baseMean = rowMeans(counts_noiseq),
  stat = de_table_noiseq$theta,
  method = "NOISeq",
  time = time_sec_noiseq
)
dt_noiseq <- dt_noiseq[!is.na(pvalue) & !is.na(padj)]
metrics_noiseq <- data.table(
  method = "NOISeq",
  n_genes = nrow(dt_noiseq),
  time_sec = round(time_sec_noiseq, 2),
  n_sig_fdr1 = sum(dt_noiseq$padj < 0.01, na.rm = TRUE),
  n_sig_fdr5 = sum(dt_noiseq$padj < 0.05, na.rm = TRUE),
  n_sig_fdr10 = sum(dt_noiseq$padj < 0.10, na.rm = TRUE),
  n_sig_lfc0.5 = sum(dt_noiseq$padj < 0.05 & abs(dt_noiseq$log2FC) > 0.5, na.rm = TRUE),
  n_sig_lfc1 = sum(dt_noiseq$padj < 0.05 & abs(dt_noiseq$log2FC) > 1, na.rm = TRUE),
  n_sig_lfc2 = sum(dt_noiseq$padj < 0.05 & abs(dt_noiseq$log2FC) > 2, na.rm = TRUE),
  p_median = median(dt_noiseq$pvalue, na.rm = TRUE),
  p_mean = mean(dt_noiseq$pvalue, na.rm = TRUE),
  p_lt_0.05 = mean(dt_noiseq$pvalue < 0.05, na.rm = TRUE),
  lfc_median = median(abs(dt_noiseq$log2FC), na.rm = TRUE),
  lfc_iqr = IQR(dt_noiseq$log2FC, na.rm = TRUE),
  dataset = CONFIG$DATASET
)

counts_glm <- counts
group_glm <- factor(make.names(as.character(group)))
t0_glm <- Sys.time()
rownames(counts_glm) <- make.unique(as.character(rownames(counts_glm)))
col_data_glm <- data.frame(condition = group_glm, row.names = colnames(counts_glm))
fit_glm <- glm_gp(counts_glm, design = ~condition, col_data = col_data_glm, size_factors = "normed_sum")
res_glm <- test_de(fit_glm, contrast = c(0, 1))
time_sec_glm <- as.numeric(difftime(Sys.time(), t0_glm, units = "secs"))
dt_glm <- data.table(
  gene = rownames(res_glm),
  pvalue = res_glm$pval,
  padj = res_glm$adj_pval,
  log2FC = res_glm$lfc,
  lfcSE = NA_real_,
  baseMean = rowMeans(counts_glm),
  stat = res_glm$f_statistic,
  method = "glmGamPoi",
  time = time_sec_glm
)
dt_glm <- dt_glm[!is.na(pvalue) & !is.na(padj)]
metrics_glm <- data.table(
  method = "glmGamPoi",
  n_genes = nrow(dt_glm),
  time_sec = round(time_sec_glm, 2),
  n_sig_fdr1 = sum(dt_glm$padj < 0.01, na.rm = TRUE),
  n_sig_fdr5 = sum(dt_glm$padj < 0.05, na.rm = TRUE),
  n_sig_fdr10 = sum(dt_glm$padj < 0.10, na.rm = TRUE),
  n_sig_lfc0.5 = sum(dt_glm$padj < 0.05 & abs(dt_glm$log2FC) > 0.5, na.rm = TRUE),
  n_sig_lfc1 = sum(dt_glm$padj < 0.05 & abs(dt_glm$log2FC) > 1, na.rm = TRUE),
  n_sig_lfc2 = sum(dt_glm$padj < 0.05 & abs(dt_glm$log2FC) > 2, na.rm = TRUE),
  p_median = median(dt_glm$pvalue, na.rm = TRUE),
  p_mean = mean(dt_glm$pvalue, na.rm = TRUE),
  p_lt_0.05 = mean(dt_glm$pvalue < 0.05, na.rm = TRUE),
  lfc_median = median(abs(dt_glm$log2FC), na.rm = TRUE),
  lfc_iqr = IQR(dt_glm$log2FC, na.rm = TRUE),
  dataset = CONFIG$DATASET
)

metrics_all <- rbindlist(list(metrics_limma, metrics_noiseq, metrics_glm), fill = TRUE)
fwrite(metrics_all, file.path(CONFIG$OUT_DIR, "pb_segerstolpe_missing_methods.csv"))
cat("Saved:", file.path(CONFIG$OUT_DIR, "pb_segerstolpe_missing_methods.csv"), "\n")
