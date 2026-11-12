root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
args <- commandArgs(trailingOnly=TRUE)
threads <- as.integer(Sys.getenv("SGCB_THREADS","4"))
edition <- args[1L]; ids <- as.integer(args[-1L])
stopifnot(edition=="final")
lib <- file.path(root,"library")
library(SGCB,lib.loc=lib)
stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
source(file.path(root,"experiments/revision_calibration.R"))
d <- readRDS(file.path(root,"experiments/design.rds"))
template <- readRDS(file.path(root,"experiments/control_template.rds"))
out <- file.path(root,"results","final"); dir.create(out,showWarnings=FALSE)
run_one <- function(index) {
    destination <- file.path(out,sprintf("matrix_%03d.rds",index))
    if (file.exists(destination)) return(invisible(index))
    job <- d$jobs[index,]; sim <- revision_matrix(job,d,template)
    if (edition=="final") {
        paired <- readRDS(file.path(root,"results","reference_scoring",sprintf("matrix_%03d.rds",index)))
        stopifnot(identical(sim$truth,paired$truth))
        baseline <- paired$moments; baseline_time <- unname(paired$seconds["moments"])
    } else baseline_time <- system.time(baseline <- moment_workflow(sim$x,sim$group))[["elapsed"]]
    time <- system.time(fit <- sgcbDE(sim$x,sim$group,threads=threads,min_count=0,min_samples=1,normalization="none"))[["elapsed"]]
    stopifnot(identical(fit$gene_id,sim$truth$id),all(is.finite(fit$padj)),all(is.finite(sim$truth$dv_effect)))
    sg <- data.frame(mean_p=fit$pvalue,mean_q=fit$padj,dv_p=fit$dv_pvalue,dv_q=fit$dv_padj,
        estimate=fit$dv_log2_cv2_ratio,lower=fit$dv_ci_lower,upper=fit$dv_ci_upper,
        mean_available=fit$de_inference_available,dv_available=fit$dv_inference_available)
    answer <- list(job=job,truth=sim$truth,SGCB=sg,moments=baseline,seconds=c(SGCB=time,moments=baseline_time),
        version=as.character(packageVersion("SGCB",lib.loc=lib)),nthreads=threads,
        mean_score_regime=attr(fit,"mean_model")$control$mean_score_regime)
    saveRDS(answer,paste0(destination,".partial"))
    stopifnot(file.rename(paste0(destination,".partial"),destination))
    cat(edition,index,job$panel,job$size,round(time,2),"seconds",format(Sys.time()),"\n"); flush.console()
    invisible(index)
}
invisible(lapply(ids,run_one))