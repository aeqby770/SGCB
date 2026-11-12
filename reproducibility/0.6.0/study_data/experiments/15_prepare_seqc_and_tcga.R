ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
PROJECT <- file.path(ROOT, "data", "original_inputs")
DATA <- file.path(ROOT, "data")
SEQ_INPUT <- file.path(PROJECT, "benchmark_v2", "01_bulk_rnaseq", "data", "bulk", "06_seqc", "GSE49712_HTSeq.txt.gz")
SEQ_ARCHIVE <- file.path(PROJECT, "benchmark_v2", "01_bulk_rnaseq", "other_methods", "results", "_archive", "seqc_multimethod", "seqc_multimethod_results.rds")
THREADS <- parallel::detectCores()
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
seqc <- as.matrix(read.table(gzfile(SEQ_INPUT), header = TRUE, row.names = 1,
                            check.names = FALSE, comment.char = ""))
seqc_group <- factor(ifelse(grepl("^A_", colnames(seqc)), "A", "B"))
stopifnot(identical(as.integer(table(seqc_group)), c(5L, 5L)))
saveRDS(list(counts = seqc, group = seqc_group), file.path(DATA, "seqc_input.rds"))
archived <- readRDS(SEQ_ARCHIVE)
fwrite(as.data.table(archived$all_results)[method != "SGCB"],
       file.path(DATA, "seqc_archived_comparator_scores.csv"))
fwrite(as.data.table(archived$metrics), file.path(DATA, "seqc_archived_consensus_metrics.csv"))

kirc <- readRDS(file.path(DATA, "tcga_kirc_cache", "tcga_kirc_tumor_vs_normal.rds"))
pid <- substr(colnames(kirc$counts), 1, 12)
overlap <- intersect(pid[kirc$group == "tumor"], pid[kirc$group == "normal"])
# Keep all discovery tumors and one observation per biological patient;
# the 13 matching normals are retained in the original-input sensitivity file.
independent <- !(kirc$group == "normal" & pid %in% overlap)
stopifnot(!anyDuplicated(pid[independent]))
saveRDS(list(counts = kirc$counts[, independent], group = droplevels(kirc$group[independent])),
        file.path(DATA, "kirc_independent_discovery.rds"))
fwrite(data.table(sample = colnames(kirc$counts), patient = pid, group = kirc$group,
                  primary_discovery = independent, exclude_from_survival_validation = TRUE),
       file.path(ROOT, "results", "kirc_discovery_sample_membership.csv"))
print(data.table(dataset = c("SEQC", "KIRC independent discovery"),
                  n_features = c(nrow(seqc), nrow(kirc$counts)),
                  group1 = c(sum(seqc_group == "A"), sum(kirc$group[independent] == "normal")),
                  group2 = c(sum(seqc_group == "B"), sum(kirc$group[independent] == "tumor"))))
