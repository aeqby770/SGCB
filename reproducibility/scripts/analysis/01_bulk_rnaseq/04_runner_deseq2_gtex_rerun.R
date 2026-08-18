PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  BULK_DIR = paste0(PROJECT_ROOT, "/benchmark/data/bulk"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/rerun_large_missing"),
  DATASET = "GTEx_blood_vs_skin",
  METHOD = "DESeq2",
  N_SAMPLES = 50L
)

dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(DESeq2)

d1 <- readRDS(file.path(CONFIG$BULK_DIR, "gtex_blood.rds"))
d2 <- readRDS(file.path(CONFIG$BULK_DIR, "gtex_skin.rds"))
common_genes <- intersect(rownames(d1$counts), rownames(d2$counts))
idx1 <- sample(ncol(d1$counts), min(CONFIG$N_SAMPLES, ncol(d1$counts)))
idx2 <- sample(ncol(d2$counts), min(CONFIG$N_SAMPLES, ncol(d2$counts)))
counts <- cbind(d1$counts[common_genes, idx1], d2$counts[common_genes, idx2])
keep <- rowSums(counts > 10) >= 10
counts <- counts[keep, ]
group <- factor(c(rep("blood", length(idx1)), rep("skin", length(idx2))))

t0 <- Sys.time()
counts <- as.matrix(counts)
counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
counts[!is.finite(counts)] <- 0
counts <- round(counts)
counts <- pmax(counts, 0)
counts <- pmin(counts, .Machine$integer.max)
storage.mode(counts) <- "integer"
group <- factor(make.names(as.character(group)))
col_data <- data.frame(condition = group, row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

dt <- data.table(
  gene = rownames(res),
  pvalue = res$pvalue,
  padj = res$padj,
  log2FC = res$log2FoldChange,
  lfcSE = res$lfcSE,
  baseMean = res$baseMean,
  stat = res$stat,
  method = CONFIG$METHOD,
  time = time_sec
)
dt <- dt[!is.na(pvalue) & !is.na(padj)]

metrics <- data.table(
  method = CONFIG$METHOD,
  n_genes = nrow(dt),
  time_sec = round(time_sec, 2),
  n_sig_fdr1 = sum(dt$padj < 0.01, na.rm = TRUE),
  n_sig_fdr5 = sum(dt$padj < 0.05, na.rm = TRUE),
  n_sig_fdr10 = sum(dt$padj < 0.10, na.rm = TRUE),
  n_sig_lfc0.5 = sum(dt$padj < 0.05 & abs(dt$log2FC) > 0.5, na.rm = TRUE),
  n_sig_lfc1 = sum(dt$padj < 0.05 & abs(dt$log2FC) > 1, na.rm = TRUE),
  n_sig_lfc2 = sum(dt$padj < 0.05 & abs(dt$log2FC) > 2, na.rm = TRUE),
  p_median = median(dt$pvalue, na.rm = TRUE),
  p_mean = mean(dt$pvalue, na.rm = TRUE),
  p_lt_0.05 = mean(dt$pvalue < 0.05, na.rm = TRUE),
  lfc_median = median(abs(dt$log2FC), na.rm = TRUE),
  lfc_iqr = IQR(dt$log2FC, na.rm = TRUE),
  dataset = CONFIG$DATASET
)

fwrite(metrics, file.path(CONFIG$OUT_DIR, "gtex_blood_vs_skin_deseq2.csv"))
cat("Saved:", file.path(CONFIG$OUT_DIR, "gtex_blood_vs_skin_deseq2.csv"), "\n")
