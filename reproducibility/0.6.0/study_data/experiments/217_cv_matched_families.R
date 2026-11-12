ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/")
ARGS <- commandArgs(trailingOnly=TRUE)
CONFIG <- list(panel=ARGS[1],design=ARGS[2],threads=as.integer(Sys.getenv("SGCB_THREADS")),
    alpha=.05,paths=list(results=file.path(ROOT,"results")))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,
    MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID",as.character(packageVersion("SGCB")))
LABEL <- paste0(CONFIG$panel,"_n",CONFIG$design)
PREFIX <- file.path(CONFIG$paths$results,paste0("cv_matched_families_",RECORD_ID,"_",LABEL))
reference_path <- file.path(CONFIG$paths$results,paste0("sgcb_",RECORD_ID,"_",LABEL,".rds"))
reference <- readRDS(reference_path)
verify_study_record(reference_path,attr(reference,"native_source_id"),RECORD_ID)

INPUT <- file.path(ROOT,"data",switch(CONFIG$panel,channel="channel_simulation.rds",zero="zero_calibration.rds"))
stopifnot(identical(attr(reference,"input_id"),unname(tools::md5sum(INPUT))))
dataset <- readRDS(INPUT)
samples <- attr(reference,"samples")
x <- sweep(dataset$abundance[reference$gene_id,samples$sample,drop=FALSE],
    2L,attr(reference,"size_factors"),"/")
group <- as.integer(factor(samples$group))-1L
evaluation <- readRDS(file.path(CONFIG$paths$results,paste0("sgcb_",RECORD_ID,"_",LABEL,"_evaluation.rds")))
truth <- copy(evaluation$genes[endpoint=="DV",.(gene_id,scenario,truth,truth_finite,truth_effect)])
stopifnot(identical(truth$gene_id,reference$gene_id))
gamma_time <- system.time(gamma <- SGCB:::.populationGG(x,group,CONFIG$threads,
    CONFIG$alpha,family="Gamma",truth=truth$truth_effect))
lognormal_time <- system.time(lognormal <- SGCB:::.populationGG(x,group,CONFIG$threads,
    CONFIG$alpha,family="lognormal",truth=truth$truth_effect))
gg <- copy(evaluation$genes[endpoint=="DV"])
genes <- rbindlist(list(data.table(method="GG",gene_id=gg$gene_id,
        effect=gg$effect,se=gg$se,pvalue=gg$pvalue,
        available=gg$available,covered=as.logical(gg$native_covered)),
    data.table(method="Gamma",gene_id=reference$gene_id,effect=gamma$genes$effect,
        se=gamma$genes$se,pvalue=gamma$genes$p,available=gamma$genes$available==1,
        covered=as.logical(gamma$genes$covered)),
    data.table(method="lognormal",gene_id=reference$gene_id,effect=lognormal$genes$effect,
        se=lognormal$genes$se,pvalue=lognormal$genes$p,available=lognormal$genes$available==1,
        covered=as.logical(lognormal$genes$covered))))
genes <- merge(genes,truth,by="gene_id",sort=FALSE)
genes[available==FALSE,pvalue:=NA_real_]
genes[,adjusted:=p.adjust(pvalue,"BH",n=.N),by=method]
genes[,called:=available & !is.na(adjusted) & adjusted<CONFIG$alpha]
metrics <- genes[,.(features=.N,available=sum(available),alternatives=sum(truth & truth_finite),
    TP=sum(called & truth & truth_finite),FP=sum(called & !truth & truth_finite),
    undefined_calls=sum(called & !truth_finite),
    nonnull_coverage=mean(covered[available & truth & truth_finite],na.rm=TRUE),
    nonnull_RMSE=sqrt(mean((effect[available & truth & truth_finite]-
        truth_effect[available & truth & truth_finite])^2))),by=method]
metrics[,c("power","FDP"):=list(TP/alternatives,FP/pmax(TP+FP,1))]
prior_time <- attr(reference,"population_timing")
metrics[,seconds:=c(unname(prior_time["total"]),gamma_time[["elapsed"]],lognormal_time[["elapsed"]])]
metrics[,timing_source:=c("stored GG API DV component","current Gamma DV refit","current lognormal DV refit")]
saveRDS(list(genes=genes,metrics=metrics,gamma=gamma,lognormal=lognormal,
    gene_ids=reference$gene_id,samples=samples,size_factors=attr(reference,"size_factors"),
    source_id=SOURCE_ID,gg_source_id=attr(reference,"native_source_id"),
    gg_input=file.path(CONFIG$paths$results,paste0("sgcb_",RECORD_ID,"_",LABEL,".rds")),config=CONFIG),paste0(PREFIX,".rds"))
fwrite(genes,paste0(PREFIX,"_genes.csv"))
fwrite(metrics,paste0(PREFIX,"_metrics.csv"))
print(metrics)
print(gamma$convergence)
print(lognormal$convergence)
