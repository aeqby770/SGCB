ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
DATASET <- ARGS[1]
NORMALIZATION <- "ratio_mode"
RNG_INPUT <- Sys.getenv("SGCB_RNG_INPUT","")
RNG_OUTPUT <- Sys.getenv("SGCB_RNG_OUTPUT","")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS"))
ALPHA <- .05
DATA_DIR <- file.path(DATA_ROOT,"data","original_inputs","benchmark_v2","03_proteomics","data")
LOADER <- file.path(DATA_ROOT,"data","original_inputs","CSBJ_submission","reproducibility",
    "scripts","analysis","03_proteomics","05_load_proteomics_data.R")
PREFIX <- file.path(DATA_ROOT,"results",paste0("protein_",DATASET))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
library(data.table)
data.table::setDTthreads(THREADS)
source(LOADER,encoding="UTF-8")

saved <- readRDS(paste0(PREFIX,"_preprocessing_diagnosis_",NORMALIZATION,".rds"))
logx <- loaded$prot_log2
group <- loaded$group
complete <- rowSums(is.finite(logx))==ncol(logx)
common <- unique(fread(paste0(PREFIX,"_common_genes.csv"))$gene)
offsets <- saved$offsets$normalization_log2_offset[match(colnames(logx),saved$offsets$sample)]
corrected_logx <- sweep(logx[complete,,drop=FALSE],2L,offsets,"-")
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")
archived <- readRDS(paste0(PREFIX,"_matched_preprocessing.rds"))
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library"))))
OUTPUT <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_protein_",DATASET))
x <- 2^corrected_logx
switch(as.character(nzchar(RNG_INPUT)),"TRUE"={
    rng_initial <- readRDS(RNG_INPUT)
    assign(".Random.seed",rng_initial$rng_state,envir=.GlobalEnv)
})
elapsed <- system.time(fit <- SGCB::sgcbDE(x,group,alpha=ALPHA,
    normalization="none",offset=0,min_count=0,min_samples=1,threads=THREADS))
attr(fit,"api_elapsed_seconds") <- elapsed[["elapsed"]]
attr(fit,"api_cpu_seconds") <- sum(elapsed[c("user.self","sys.self")])
process_memory <- ps::ps_memory_info(ps::ps_handle())
attr(fit,"process_peak_working_set_bytes") <- unname(process_memory[["peak_wset"]])
attr(fit,"process_peak_pagefile_bytes") <- unname(process_memory[["peak_pagefile"]])
attr(fit,"samples") <- data.frame(sample=colnames(x),group=as.character(group))
switch(as.character(nzchar(RNG_INPUT)),"TRUE"={
    attr(fit,"initial_rng_state") <- rng_initial$rng_state
    attr(fit,"initial_rng_file") <- normalizePath(RNG_INPUT,winslash="/")
})
attr(fit,"matched_normalization_offsets") <- saved$offsets
stopifnot(identical(fit$gene_id,rownames(x)),identical(archived$offsets,saved$offsets),
    identical(fit$padj,p.adjust(fit$pvalue,"BH",n=nrow(fit))),
    all(fit$pvalue>=0 & fit$pvalue<=1))
scores <- rbindlist(list(
    as.data.table(fit)[,.(method="SGCB",gene=gene_id,effect=log2FoldChange,p=pvalue)],
    archived$scores[method!="SGCB",.(method,gene,effect,p)]
))
scores <- merge(scores[gene%in%common],
    data.table(gene=rownames(logx),is_spike=loaded$is_spike),by="gene")
scores[,q:=p.adjust(p,"BH",n=length(common)),by=method]
scores[,called:=is.finite(q)&q<ALPHA]
old_comparators <- archived$scores[method!="SGCB",.(method,gene,p_old=p,q_old=q)]
readback <- merge(scores[method!="SGCB"],old_comparators,by=c("method","gene"))
stopifnot(nrow(readback)==nrow(old_comparators),
    identical(readback$p,readback$p_old),max(abs(readback$q-readback$q_old))<1e-14)
metrics <- scores[,.(features=.N,spike_proteins=sum(is_spike),
    null_proteins=sum(!is_spike),TP=sum(called&is_spike),
    FP=sum(called&!is_spike),null_negative=sum(called&!is_spike&effect<0),
    null_positive=sum(called&!is_spike&effect>0),
    median_null_effect=median(effect[!is_spike]),
    raw_null_rejection=mean(p[!is_spike]<ALPHA)),by=method]
metrics[,c("FDP","power","dataset","normalization"):=list(
    FP/pmax(TP+FP,1),TP/spike_proteins,DATASET,NORMALIZATION)]
metrics[,fit_seconds:=archived$metrics$fit_seconds[match(method,archived$metrics$method)]]
metrics[method=="SGCB",fit_seconds:=elapsed[["elapsed"]]]
metrics[,fit_cpu_seconds:=NA_real_]
metrics[method=="SGCB",fit_cpu_seconds:=attr(fit,"api_cpu_seconds")]
metrics[,c("process_peak_working_set_bytes","process_peak_pagefile_bytes"):=list(NA_real_,NA_real_)]
metrics[method=="SGCB",c("process_peak_working_set_bytes","process_peak_pagefile_bytes"):=list(
    attr(fit,"process_peak_working_set_bytes"),attr(fit,"process_peak_pagefile_bytes"))]
metrics[,score_source:=ifelse(method=="SGCB","current arithmetic-mean public API",
    "unchanged archived matched-input fit")]
stopifnot(all(metrics$features==length(common)),setequal(rownames(x),saved$fit$gene_id))
saveRDS(fit,paste0(OUTPUT,".rds"))
fwrite(as.data.table(fit),paste0(OUTPUT,"_genes.csv"))
fwrite(scores,paste0(OUTPUT,"_scores.csv"))
fwrite(metrics,paste0(OUTPUT,"_metrics.csv"))
summary <- data.table(dataset=basename(OUTPUT),version=VERSION,run_id=attr(fit,"run_id"),
    native_source_id=attr(fit,"native_source_id"),n_genes=nrow(fit),
    n_ctrl=attr(fit,"n_ctrl"),n_treat=attr(fit,"n_treat"),
    elapsed_seconds=elapsed[["elapsed"]],cpu_seconds=attr(fit,"api_cpu_seconds"),
    process_peak_working_set_bytes=attr(fit,"process_peak_working_set_bytes"),
    process_peak_pagefile_bytes=attr(fit,"process_peak_pagefile_bytes"),
    DE=sum(fit$DE),DV=sum(fit$DV),DV_available=sum(fit$dv_inference_available))
fwrite(summary,paste0(OUTPUT,"_summary.csv"))
print(metrics)
print(attr(fit,"mean_timing"))

switch(as.character(nzchar(RNG_OUTPUT)),"TRUE"={
    saveRDS(list(rng_state=get(".Random.seed",envir=.GlobalEnv),source=basename(OUTPUT)),RNG_OUTPUT)
})
