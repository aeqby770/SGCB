ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
INPUT <- file.path(file.path(ROOT, "data", "original_inputs"), "benchmark_v2", "01_bulk_rnaseq", "data", "bulk",
                    "05_pseudobulk", "kang18_pseudobulk.rds")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
.libPaths(c(Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library")),.libPaths()))
library(SGCB)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")
ds <- readRDS(INPUT)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.5.7")
PREFIX <- paste0("sgcb_",RECORD_ID,"_pseudobulk_kang18")
res <- readRDS(file.path(ROOT,"results",paste0(PREFIX,".rds")))
verify_study_record(file.path(ROOT,"results",paste0(PREFIX,".rds")),attr(res,"native_source_id"),RECORD_ID)
stopifnot(
    identical(colnames(ds$counts),attr(res,"samples")$sample),
    identical(as.character(ds$group),as.character(attr(res,"samples")$group)))
counts <- ds$counts[res$gene_id, , drop = FALSE]
sf <- attr(res,"size_factors")
logx <- log2(sweep(counts, 2, sf, "/") + 0.5)
map <- data.table(sample = colnames(counts), column = seq_len(ncol(counts)),
    donor = sub("_.*$", "", colnames(counts)), condition = as.character(ds$group))
paired <- dcast(map, donor ~ condition, value.var = "column")
difference <- logx[, paired$stim, drop = FALSE] - logx[, paired$ctrl, drop = FALSE]
effect <- rowMeans(difference)
se <- matrixStats::rowSds(difference) / sqrt(nrow(paired))
p <- 2 * pt(-abs(effect / pmax(se, 1e-12)), nrow(paired) - 1)
out <- data.table(gene_id = res$gene_id, paired_log2FC = effect, paired_p = p,
    paired_q = p.adjust(p, "BH"), two_group_log2FC = res$log2FoldChange,
    two_group_p = res$pvalue, two_group_q = res$padj)
fwrite(out, file.path(ROOT, "results", paste0(PREFIX,"_paired_genes.csv")))
fwrite(map, file.path(ROOT, "results", paste0(PREFIX,"_biological_units.csv")))
summary <- out[, .(n_genes = .N, n_donors = nrow(paired), n_profiles = nrow(map),
    paired_calls = sum(paired_q < 0.05), two_group_calls = sum(two_group_q < 0.05),
    overlap = sum(paired_q < 0.05 & two_group_q < 0.05),
    effect_max_difference = max(abs(paired_log2FC - two_group_log2FC)),
    p_rank_correlation = cor(paired_p, two_group_p, method = "spearman"))]
fwrite(summary, file.path(ROOT, "results", paste0(PREFIX,"_paired_summary.csv")))
print(summary)
