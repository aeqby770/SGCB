ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
LABEL <- commandArgs(trailingOnly = TRUE)[1]
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor"))))
OUTPUT <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_nb_",LABEL))
INPUT <- file.path(file.path(DATA_ROOT, "data", "original_inputs"), "benchmark_v2", "01_bulk_rnaseq", "data", "bulk",
                    "02_simulation", paste0(LABEL, ".rds"))
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
.libPaths(c(Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor")), .libPaths()))
library(SGCB)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")
ds <- readRDS(INPUT)
group <- factor(ds$sample_info$group)
stopifnot(identical(colnames(ds$counts), ds$sample_info$sample_id))
elapsed <- system.time(res <- sgcbDE(ds$counts, group, alpha = 0.05, threads=THREADS))
stopifnot(identical(attr(res,"native_source_id"),SOURCE_ID))
attr(res,"api_elapsed_seconds") <- elapsed[["elapsed"]]
attr(res,"api_cpu_seconds") <- unname(elapsed[["user.self"]]+elapsed[["sys.self"]])
process_memory <- ps::ps_memory_info(ps::ps_handle())
attr(res,"process_peak_working_set_bytes") <- unname(process_memory[["peak_wset"]])
attr(res,"process_peak_pagefile_bytes") <- unname(process_memory[["peak_pagefile"]])
attr(res,"input_file") <- normalizePath(INPUT,winslash="/")
attr(res,"input_id") <- unname(tools::md5sum(INPUT))
attr(res,"samples") <- data.frame(sample=colnames(ds$counts),group=as.character(group))
out <- merge(as.data.table(ds$truth), as.data.table(res)[, .(gene = gene_id, pvalue, padj,
                  log2FoldChange, lfcSE, de_inference_available, fit_converged, dv_inference_available, DE)], by = "gene", all.x = TRUE)
out[, called := !is.na(DE) & DE]
summary <- out[, .(dataset = LABEL,run_id=attr(res,"run_id"),
    native_source_id=attr(res,"native_source_id"),n_genes=nrow(res),
    n_ctrl=attr(res,"n_ctrl"),n_treat=attr(res,"n_treat"),
    n_input = .N, n_tested = sum(is.finite(pvalue)),
    n_group = sum(group == levels(group)[1]), n_true = sum(is_de),
    TP = sum(called & is_de), FP = sum(called & !is_de), FN = sum(!called & is_de),
    TN = sum(!called & !is_de), elapsed_seconds = elapsed[["elapsed"]],version=VERSION,
    cpu_seconds=attr(res,"api_cpu_seconds"),
    process_peak_working_set_bytes=attr(res,"process_peak_working_set_bytes"),
    process_peak_pagefile_bytes=attr(res,"process_peak_pagefile_bytes"))]
summary[, `:=`(FDP = FP / pmax(TP + FP, 1), sensitivity = TP / pmax(TP + FN, 1),
               F1 = 2 * TP / pmax(2 * TP + FP + FN, 1))]
summary[,c("mixture_weight_KKT_gap","mixture_active_states"):=list(
    max(attr(res,"mean_model")$convergence$kkt_gap),
    nrow(attr(res,"mean_model")$pairs))]
saveRDS(res, paste0(OUTPUT,".rds"))
fwrite(out, paste0(OUTPUT,"_genes.csv"))
fwrite(summary, paste0(OUTPUT,"_summary.csv"))
print(summary)

