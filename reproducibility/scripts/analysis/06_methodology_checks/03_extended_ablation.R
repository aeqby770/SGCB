PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(pROC)
library(DESeq2)
library(edgeR)
library(limma)
devtools::load_all(paste0(PROJECT_ROOT, "/SGCB"))

CONFIG <- list(
  SEED = 12345L,
  DATA_DIR = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/02_simulation"),
  BULK_DIR = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/03_gtex"),
  SEQC_FILE = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/06_seqc/GSE49712_HTSeq.txt.gz"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/methodology_checks"),
  FDR = 0.05,
  MIN_COUNT = 10L,
  MIN_SAMPLES = 2L,
  GS_FDR = 0.001,
  GS_LFC = 1,
  GTEX_NULL_N = 10L
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

variant_dt <- data.table(
  variant = c("SGCB_full", "no_gg_variance", "no_hierarchical", "no_firth", "no_gamma_submodel", "no_natural_grad"),
  use_hierarchical = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  use_firth = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  use_gamma_submodel = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  use_gg_variance = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  use_natural_grad = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)

sim_files <- list.files(CONFIG$DATA_DIR, pattern = "^sim_.*\\.rds$", full.names = TRUE)

sim_results <- rbindlist(lapply(sim_files, function(sim_file) {
  sim_obj <- readRDS(sim_file)
  counts <- sim_obj$counts
  group <- factor(sim_obj$sample_info$group)
  truth <- sim_obj$truth$gene[sim_obj$truth$is_de]
  keep <- rowSums(counts > 5) >= CONFIG$MIN_SAMPLES
  counts <- counts[keep, , drop = FALSE]
  truth <- intersect(truth, rownames(counts))

  rbindlist(lapply(seq_len(nrow(variant_dt)), function(i) {
    v <- variant_dt[i]
    t0 <- Sys.time()
    fit <- sgcbDE(
      counts,
      group,
      use_hierarchical = v$use_hierarchical,
      use_firth = v$use_firth,
      use_gamma_submodel = v$use_gamma_submodel,
      use_gg_variance = v$use_gg_variance,
      use_natural_grad = v$use_natural_grad
    )
    time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    sig <- rownames(fit)[which(!is.na(fit$padj) & fit$padj < CONFIG$FDR)]
    tp <- sum(sig %in% truth)
    fp <- sum(!(sig %in% truth))
    fn <- sum(!(truth %in% sig))
    precision <- tp / max(1, tp + fp)
    recall <- tp / max(1, length(truth))
    f1 <- 2 * precision * recall / max(1e-8, precision + recall)
    data.table(
      dataset = sub("\\.rds$", "", basename(sim_file)),
      variant = v$variant,
      TP = tp,
      FP = fp,
      FN = fn,
      precision = precision,
      recall = recall,
      F1 = f1,
      actual_FDR = fp / max(1, tp + fp),
      n_sig = length(sig),
      time_sec = time_sec
    )
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

sim_summary <- sim_results[, .(
  mean_precision = mean(precision, na.rm = TRUE),
  mean_recall = mean(recall, na.rm = TRUE),
  mean_F1 = mean(F1, na.rm = TRUE),
  mean_actual_FDR = mean(actual_FDR, na.rm = TRUE),
  mean_n_sig = mean(n_sig, na.rm = TRUE),
  mean_time_sec = mean(time_sec, na.rm = TRUE)
), by = variant]

seqc_raw <- read.table(gzfile(CONFIG$SEQC_FILE), header = TRUE, row.names = 1, check.names = FALSE)
counts_seqc <- as.matrix(seqc_raw)
storage.mode(counts_seqc) <- "integer"
group_seqc <- factor(ifelse(grepl("^A_", colnames(counts_seqc)), "A", "B"))
keep_seqc_gs <- rowSums(counts_seqc > 5) >= CONFIG$MIN_SAMPLES
counts_seqc_gs <- counts_seqc[keep_seqc_gs, , drop = FALSE]

cd <- data.frame(condition = group_seqc, row.names = colnames(counts_seqc_gs))
dds <- DESeqDataSetFromMatrix(counts_seqc_gs, cd, ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, quiet = TRUE)
dds <- nbinomWaldTest(dds)
res_d <- results(dds)
de_d <- rownames(res_d)[which(res_d$padj < CONFIG$GS_FDR & abs(res_d$log2FoldChange) > CONFIG$GS_LFC)]

dge <- DGEList(counts = counts_seqc_gs, group = group_seqc)
dge <- calcNormFactors(dge)
design_gs <- model.matrix(~group_seqc)
dge <- estimateDisp(dge, design_gs)
fit_e <- glmQLFit(dge, design_gs)
test_e <- glmQLFTest(fit_e, coef = 2)
res_e <- topTags(test_e, n = Inf, sort.by = "none")$table
de_e <- rownames(res_e)[which(res_e$FDR < CONFIG$GS_FDR & abs(res_e$logFC) > CONFIG$GS_LFC)]

v_voom <- voom(dge, design_gs, plot = FALSE)
fit_l <- lmFit(v_voom, design_gs)
fit_l <- eBayes(fit_l)
res_l <- topTable(fit_l, coef = 2, number = Inf, sort.by = "none")
de_l <- rownames(res_l)[which(res_l$adj.P.Val < CONFIG$GS_FDR & abs(res_l$logFC) > CONFIG$GS_LFC)]

gold_genes <- Reduce(intersect, list(de_d, de_e, de_l))
keep_seqc <- rowSums(counts_seqc > 5) >= CONFIG$MIN_SAMPLES
counts_seqc_fit <- counts_seqc[keep_seqc, , drop = FALSE]
truth_seqc <- intersect(gold_genes, rownames(counts_seqc_fit))

seqc_results <- rbindlist(lapply(seq_len(nrow(variant_dt)), function(i) {
  v <- variant_dt[i]
  t0 <- Sys.time()
  fit <- sgcbDE(
    counts_seqc_fit,
    group_seqc,
    use_hierarchical = v$use_hierarchical,
    use_firth = v$use_firth,
    use_gamma_submodel = v$use_gamma_submodel,
    use_gg_variance = v$use_gg_variance,
    use_natural_grad = v$use_natural_grad
  )
  time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sig <- rownames(fit)[which(!is.na(fit$padj) & fit$padj < CONFIG$FDR)]
  tp <- sum(sig %in% truth_seqc)
  fp <- sum(!(sig %in% truth_seqc))
  precision <- tp / max(1, tp + fp)
  recall <- tp / max(1, length(truth_seqc))
  f1 <- 2 * precision * recall / max(1e-8, precision + recall)
  score <- -log10(pmax(fit$pvalue, 1e-300))
  label <- as.integer(rownames(fit) %in% truth_seqc)
  valid <- is.finite(score)
  auc <- if (sum(label[valid]) > 0 && sum(label[valid]) < sum(valid)) as.numeric(pROC::auc(pROC::roc(label[valid], score[valid], quiet = TRUE, direction = "<"))) else NA_real_
  data.table(
    variant = v$variant,
    AUC = auc,
    precision = precision,
    recall = recall,
    F1 = f1,
    actual_FDR = fp / max(1, tp + fp),
    n_sig = length(sig),
    time_sec = time_sec
  )
}), use.names = TRUE, fill = TRUE)

gtex_obj <- readRDS(file.path(CONFIG$BULK_DIR, "gtex_liver.rds"))
set.seed(CONFIG$SEED)
idx_null <- sample(seq_len(ncol(gtex_obj$counts)), 2L * CONFIG$GTEX_NULL_N)
counts_null <- round(gtex_obj$counts[, idx_null, drop = FALSE])
group_null <- factor(rep(c("A", "B"), each = CONFIG$GTEX_NULL_N))
keep_null <- rowSums(counts_null > 5) >= CONFIG$MIN_SAMPLES
counts_null <- counts_null[keep_null, , drop = FALSE]

null_results <- rbindlist(lapply(seq_len(nrow(variant_dt)), function(i) {
  v <- variant_dt[i]
  fit <- sgcbDE(
    counts_null,
    group_null,
    use_hierarchical = v$use_hierarchical,
    use_firth = v$use_firth,
    use_gamma_submodel = v$use_gamma_submodel,
    use_gg_variance = v$use_gg_variance,
    use_natural_grad = v$use_natural_grad
  )
  pv <- fit$pvalue[!is.na(fit$pvalue)]
  data.table(
    variant = v$variant,
    ks_p = ks.test(pv, "punif")$p.value,
    null_fdr_05 = mean(!is.na(fit$padj) & fit$padj < CONFIG$FDR),
    mean_pvalue = mean(pv),
    prop_p_lt_0.01 = mean(pv < 0.01)
  )
}), use.names = TRUE, fill = TRUE)

fwrite(variant_dt, file.path(CONFIG$OUT_DIR, "ablation_variant_settings.csv"))
fwrite(sim_results, file.path(CONFIG$OUT_DIR, "ablation_simulation_results.csv"))
fwrite(sim_summary, file.path(CONFIG$OUT_DIR, "ablation_simulation_summary.csv"))
fwrite(seqc_results, file.path(CONFIG$OUT_DIR, "ablation_seqc_results.csv"))
fwrite(null_results, file.path(CONFIG$OUT_DIR, "ablation_null_results.csv"))
