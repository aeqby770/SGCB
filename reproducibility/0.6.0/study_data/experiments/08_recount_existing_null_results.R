ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
RESULTS <- file.path(ROOT,"results")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS"))
ALPHA <- .05
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
library(data.table)
setDTthreads(THREADS)
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.6.0")
LABELS <- paste0("sgcb_",RECORD_ID,"_readcount_null_",1:10)
LABELS[3] <- paste0(LABELS[3],ifelse(RECORD_ID=="0.5.9","_streamed",""))
PATHS <- file.path(RESULTS,paste0(LABELS,".rds"))
OUTPUT <- file.path(RESULTS,paste0("sgcb_",RECORD_ID,"_readcount_null"))

# Read the complete current API results; split 3 retains its actual streamed run.
runs <- list(
    readRDS(PATHS[1]),readRDS(PATHS[2]),readRDS(PATHS[3]),readRDS(PATHS[4]),
    readRDS(PATHS[5]),readRDS(PATHS[6]),readRDS(PATHS[7]),readRDS(PATHS[8]),
    readRDS(PATHS[9]),readRDS(PATHS[10]))
fits <- rbindlist(runs,idcol="split")
mean_post <- rbindlist(lapply(runs,attr,"mean_posterior"),idcol="split")
cv_post <- rbindlist(lapply(runs,attr,"population_posterior"),idcol="split")
stopifnot(identical(fits$split,mean_post$split),identical(fits$split,cv_post$split))
fits[,mean_covered:=mean_post$lower<=0 & mean_post$upper>=0]
fits[,dv_covered:=dv_ci_lower<=0 & dv_ci_upper>=0]
metadata <- data.table(split=1:10,result_file=basename(PATHS),
    run_id=unlist(lapply(runs,attr,"run_id")),
    input_id=unlist(lapply(runs,attr,"input_id")),
    n_A=unlist(lapply(runs,attr,"n_ctrl")),n_B=unlist(lapply(runs,attr,"n_treat")),
    api_seconds=unlist(lapply(runs,attr,"api_elapsed_seconds")),
    api_peak_working_set_bytes=unlist(lapply(runs,attr,"process_peak_working_set_bytes")))
resources <- rbindlist(lapply(file.path(RESULTS,paste0(LABELS,"_process_runtime.csv")),fread),
    idcol="split",fill=TRUE)
stopifnot(all(resources$execution_outcome=="completed"),all(resources$exit_code==0))
metadata[,c("process_seconds","process_peak_working_set_bytes"):=list(
    resources$elapsed_seconds,resources$peak_working_set_bytes)]
de <- fits[,.(method="SGCB",target="Mean expression",n_common_genes=.N,
    n_valid_p=sum(de_inference_available),calls=sum(DE),any_discovery=as.integer(any(DE)),
    raw_p005_full=sum(pvalue<ALPHA,na.rm=TRUE)/.N,
    raw_p005_available=mean(pvalue[de_inference_available]<ALPHA,na.rm=TRUE),
    null_coverage=mean(mean_covered[de_inference_available],na.rm=TRUE)),by=split]
dv <- fits[,.(method="SGCB",target="CV2",n_common_genes=.N,
    n_valid_p=sum(dv_inference_available),calls=sum(DV),any_discovery=as.integer(any(DV)),
    raw_p005_full=sum(dv_pvalue<ALPHA,na.rm=TRUE)/.N,
    raw_p005_available=mean(dv_pvalue[dv_inference_available]<ALPHA,na.rm=TRUE),
    null_coverage=mean(dv_covered[dv_inference_available],na.rm=TRUE)),by=split]
metrics <- merge(rbindlist(list(de,dv)),metadata,by="split",sort=FALSE)
metrics[,positive_gene_fraction:=calls/n_common_genes]
summary <- metrics[,.(partitions=.N,with_discovery=sum(any_discovery),
    mean_calls=mean(calls),min_calls=min(calls),max_calls=max(calls),
    partition_fraction_with_calls=mean(any_discovery),
    mean_raw_p005=mean(raw_p005_full),mean_conditional_coverage=mean(null_coverage),
    mean_api_seconds=mean(api_seconds),mean_process_seconds=mean(process_seconds),
    maximum_peak_bytes=max(process_peak_working_set_bytes)),
    by=.(method,target)]
# These fixed partitions overlap; do not attach an independent-binomial interval.
fwrite(metrics,paste0(OUTPUT,"_metrics.csv"),bom=FALSE,eol="\n")
fwrite(summary,paste0(OUTPUT,"_summary.csv"),bom=FALSE,eol="\n")
print(metrics[,.(split,target,n_common_genes,n_valid_p,calls,null_coverage,api_seconds)])
print(summary)
