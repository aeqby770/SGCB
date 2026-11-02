ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
INPUT <- ARGS[1]
LABEL <- ARGS[2]
TYPE <- ARGS[3]
RNG_INPUT <- Sys.getenv("SGCB_RNG_INPUT","")
RNG_OUTPUT <- Sys.getenv("SGCB_RNG_OUTPUT","")
PARENT_OUTPUT <- Sys.getenv("SGCB_MEAN_PARENT_OUTPUT","")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
.libPaths(c(Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor")),.libPaths()))
library(SGCB)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")
switch(as.character(nzchar(RNG_INPUT)),"TRUE"={
    rng_initial <- readRDS(RNG_INPUT)
    assign(".Random.seed",rng_initial$rng_state,envir=.GlobalEnv)
})
ds <- readRDS(INPUT)
samples <- switch(TYPE,liver={
    original <- readRDS(file.path(ROOT,"data","gtex_liver_original_null_split_1.rds"))
    switch(ARGS[4],"1"=data.table(sample=original$sample_ids,group=as.character(original$group)),
        fread(file.path(ROOT,"results",paste0("null_revision_",ARGS[4],"_samples.csv"))))
},NULL)
x <- switch(TYPE,counts=as.matrix(ds$counts),abundance=,abundance_zero=ds$abundance,normalized=as.matrix(ds$abundance),
    liver=as.matrix(ds$counts[,match(samples$sample,colnames(ds$counts))]))
group <- factor(switch(TYPE,liver=samples$group,ds$group))
index <- switch(TYPE,counts=seq_len(ncol(x)),normalized=seq_len(ncol(x)),liver=seq_len(ncol(x)),abundance=,abundance_zero=c(
    which(group==levels(group)[1])[seq_len(as.integer(ARGS[4]))],
    which(group==levels(group)[2])[seq_len(as.integer(ARGS[5]))]))
elapsed <- system.time(res <- sgcbDE(x[,index,drop=FALSE],group[index],alpha=.05,
    normalization=switch(TYPE,counts="tmm",liver="tmm",abundance=,abundance_zero="none",normalized="none"),
    offset=switch(TYPE,counts=.5,liver=.5,abundance=0,abundance_zero=.5,normalized=.5),
    min_count=switch(TYPE,counts=10,liver=10,abundance=,abundance_zero=0,normalized=1),
    min_samples=switch(TYPE,counts=2,liver=2,abundance=,abundance_zero=1,normalized=2),threads=THREADS))
switch(as.character(nzchar(PARENT_OUTPUT)),"TRUE"={
    attr(res,"mean_parent_file") <- normalizePath(PARENT_OUTPUT,winslash="/",mustWork=TRUE)
})
attr(res,"api_elapsed_seconds") <- elapsed[["elapsed"]]
attr(res,"api_cpu_seconds") <- unname(elapsed[["user.self"]]+elapsed[["sys.self"]])
process_memory <- ps::ps_memory_info(ps::ps_handle())
attr(res,"process_peak_working_set_bytes") <- unname(process_memory[["peak_wset"]])
attr(res,"process_peak_pagefile_bytes") <- unname(process_memory[["peak_pagefile"]])
switch(as.character(nzchar(RNG_INPUT)),"TRUE"={
    attr(res,"initial_rng_state") <- rng_initial$rng_state
    attr(res,"initial_rng_file") <- normalizePath(RNG_INPUT,winslash="/")
})
attr(res,"input_file") <- normalizePath(INPUT,winslash="/")
attr(res,"input_id") <- unname(tools::md5sum(INPUT))
attr(res,"samples") <- data.frame(sample=colnames(x)[index],group=as.character(group[index]))
saveRDS(res,file.path(ROOT,"results",paste0(LABEL,".rds")))
fwrite(as.data.table(res),file.path(ROOT,"results",paste0(LABEL,"_genes.csv")))
native <- as.data.table(attr(res,"distribution_fits"))
native[,`:=`(gene_id=res$gene_id,zero_fraction=rowMeans(x[res$gene_id,index,drop=FALSE]==0))]
fwrite(native,file.path(ROOT,"results",paste0(LABEL,"_fits.csv")))
fit <- attr(res,"fit_diagnostics")
summary <- data.table(dataset=LABEL,version=as.character(packageVersion("SGCB")),
    run_id=attr(res,"run_id"),fit_started=attr(res,"fit_started"),
    native_source_id=attr(res,"native_source_id"),input_id=attr(res,"input_id"),
    n_genes=nrow(res),n_ctrl=attr(res,"n_ctrl"),n_treat=attr(res,"n_treat"),
    elapsed_seconds=elapsed[["elapsed"]],cpu_seconds=attr(res,"api_cpu_seconds"),
    process_peak_working_set_bytes=attr(res,"process_peak_working_set_bytes"),
    process_peak_pagefile_bytes=attr(res,"process_peak_pagefile_bytes"),
    converged_ctrl=mean(fit$ctrl$converged),
    converged_treat=mean(fit$treat$converged),converged_both=mean(res$fit_converged),
    DV_available=sum(res$dv_inference_available),DV_moments_finite=sum(res$dv_moments_finite),
    DE=sum(res$DE),DV=sum(res$DV),DE_only=sum(res$channel=="DE_only"),
    DV_only=sum(res$channel=="DV_only"),DE_and_DV=sum(res$channel=="DE_and_DV"),
    shape_exploratory=sum(res$shape_Q_padj<.05,na.rm=TRUE),
    logspread=sum(res$logspread_padj<.05,na.rm=TRUE),
    joint_DE=sum(res$DE_joint_q<.05),joint_DV=sum(res$DV_joint_q<.05))
fwrite(summary,file.path(ROOT,"results",paste0(LABEL,"_summary.csv")))
fwrite(data.table(component=names(attr(res,"elapsed_components")),
    elapsed_seconds=unname(attr(res,"elapsed_components"))),
    file.path(ROOT,"results",paste0(LABEL,"_timing.csv")))
print(summary)
switch(as.character(nzchar(RNG_OUTPUT)),"TRUE"={
    saveRDS(list(rng_state=get(".Random.seed",envir=.GlobalEnv),source=LABEL),RNG_OUTPUT)
})
