ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
PANEL <- ARGS[1]
N0 <- as.integer(ARGS[2])
N1 <- as.integer(ARGS[3])
THREADS <- as.integer(Sys.getenv("SGCB_THREADS"))
ALPHA <- .05
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library"))))
PREFIX <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_",PANEL,"_n",N0,"_",N1))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")

res <- readRDS(paste0(PREFIX,".rds"))
mean_timing <- attr(res,"mean_timing")
fwrite(data.table(component=names(mean_timing),seconds=unname(mean_timing)),
    paste0(PREFIX,"_mean_timing.csv"))
print(mean_timing)
ds <- readRDS(file.path(ROOT,"data",switch(PANEL,
    channel="channel_simulation.rds",zero="zero_calibration.rds")))
group_all <- factor(ds$group)
indices <- c(which(group_all==levels(group_all)[1])[seq_len(N0)],
    which(group_all==levels(group_all)[2])[seq_len(N1)])
x <- ds$abundance[,indices,drop=FALSE]
group <- as.integer(droplevels(group_all[indices]))-1L
truth <- as.data.table(ds$truth)[match(res$gene_id,gene_id)]
de_truth <- switch(PANEL,channel=truth$de_arithmetic_truth,
    zero=abs(log2(truth$m1/truth$m0))>1e-8)
de_effect <- switch(PANEL,channel=log2(truth$mean1/truth$mean0),
    zero=log2(truth$m1/truth$m0))
dv_effect <- log2(truth$cv21/truth$cv20)
dv_finite <- is.finite(dv_effect)
dv_truth <- dv_finite & abs(dv_effect)>1e-8
dv_effect[dv_finite & !dv_truth] <- 0

# Use the intervals already computed by the current sgcbDE entry point.
reference <- as.data.frame(attr(res,"population_posterior"))
ci <- data.frame(lower=res$dv_ci_lower,upper=res$dv_ci_upper,
    width=res$dv_ci_upper-res$dv_ci_lower,
    covered=as.integer(dv_effect>=res$dv_ci_lower & dv_effect<=res$dv_ci_upper))
mean_posterior <- as.data.table(attr(res,"mean_posterior"))
mean_covered <- de_effect>=mean_posterior$lower & de_effect<=mean_posterior$upper
de <- data.table(gene_id=res$gene_id,scenario=truth$scenario,endpoint="DE",
    available=res$de_inference_available,pvalue=res$pvalue,qvalue=res$padj,
    called=res$DE,truth_finite=TRUE,truth=de_truth,truth_effect=de_effect,
    effect=res$log2FoldChange,se=res$lfcSE,
    lower=mean_posterior$lower,upper=mean_posterior$upper,
    width=mean_posterior$upper-mean_posterior$lower,covered=mean_covered,
    native_covered=mean_covered)
dv <- data.table(gene_id=res$gene_id,scenario=truth$scenario,endpoint="DV",
    available=res$dv_inference_available,pvalue=res$dv_pvalue,qvalue=res$dv_padj,
    called=res$DV,truth_finite=dv_finite,truth=dv_truth,truth_effect=dv_effect,
    effect=res$dv_log2_cv2_ratio,se=res$dv_se_log2_cv2,
    lower=ci$lower,upper=ci$upper,width=ci$width,
    covered=ci$covered==1,native_covered=ci$covered==1)
genes <- rbindlist(list(de,dv))
genes[,c("panel","design","n0","n1"):=list(PANEL,paste0(N0,":",N1),N0,N1)]
genes[,c("error","interval_available"):=list(effect-truth_effect,!is.na(lower)&!is.na(upper))]
metrics <- genes[,.(features=.N,available=sum(available),finite_truth=sum(truth_finite),
    nulls=sum(truth_finite & !truth),alternatives=sum(truth_finite & truth),
    TP=sum(called & truth_finite & truth),FP=sum(called & truth_finite & !truth),
    undefined_calls=sum(called & !truth_finite),
    raw_null_rate=sum(pvalue<ALPHA & truth_finite & !truth,na.rm=TRUE)/sum(truth_finite & !truth),
    raw_null_rate_available=sum(pvalue<ALPHA & available & truth_finite & !truth,na.rm=TRUE)/
        sum(available & truth_finite & !truth),
    RMSE=sqrt(mean(error[available & truth_finite]^2)),bias=mean(error[available & truth_finite]),
    median_se=median(se[available & truth_finite],na.rm=TRUE),
    interval_available=sum(interval_available & truth_finite),
    coverage_all=sum(covered & available & truth_finite,na.rm=TRUE)/sum(truth_finite),
    coverage_available=sum(covered & available & truth_finite,na.rm=TRUE)/
        sum(available & truth_finite),
    native_coverage_available=sum(native_covered & available & truth_finite,na.rm=TRUE)/
        sum(available & truth_finite),
    median_width=median(width[available & truth_finite]),
    mean_width=mean(width[available & truth_finite])),
    by=.(panel,design,endpoint)]
metrics[,c("power","FDP","FDP_with_undefined","api_seconds","api_cpu_seconds",
    "process_peak_working_set_bytes","process_peak_pagefile_bytes"):=list(
    TP/alternatives,FP/pmax(1,TP+FP),(FP+undefined_calls)/pmax(1,TP+FP+undefined_calls),
    attr(res,"api_elapsed_seconds"),attr(res,"api_cpu_seconds"),
    attr(res,"process_peak_working_set_bytes"),attr(res,"process_peak_pagefile_bytes"))]
mean_posterior <- as.data.table(attr(res,"mean_posterior"))
metrics[endpoint=="DE",c("mixture_weight_KKT_gap","mixture_active_states",
    "called_zero_se"):=list(
    max(attr(res,"mean_model")$convergence$kkt_gap),
    nrow(attr(res,"mean_model")$pairs),
    sum(res$DE & res$lfcSE==0))]
scenario <- genes[,.(features=.N,available=sum(available),calls=sum(called),
    truth_finite=sum(truth_finite),alternatives=sum(truth & truth_finite),
    TP=sum(called & truth_finite & truth),FP=sum(called & truth_finite & !truth),
    undefined_calls=sum(called & !truth_finite),raw_p005=sum(pvalue<ALPHA,na.rm=TRUE)/.N,
    RMSE=sqrt(mean(error[available & truth_finite]^2)),
    bias=mean(error[available & truth_finite]),
    coverage_available=sum(covered & available & truth_finite,na.rm=TRUE)/
        sum(available & truth_finite),
    native_coverage_available=sum(native_covered & available & truth_finite,na.rm=TRUE)/
        sum(available & truth_finite),
    median_width=median(width[available & truth_finite]),
    mean_width=mean(width[available & truth_finite])),
    by=.(panel,design,endpoint,scenario)]
checks <- list(runtime_scope="complete sgcbDE entry point",
    reference=attr(res,"dv_score_type"),input_file=attr(res,"input_file"),
    api_seconds=attr(res,"api_elapsed_seconds"),
    source_id=attr(res,"native_source_id"))
saveRDS(list(genes=genes,metrics=metrics,scenario=scenario,checks=checks),
    paste0(PREFIX,"_evaluation.rds"))
fwrite(genes,paste0(PREFIX,"_evaluation_genes.csv"))
fwrite(metrics,paste0(PREFIX,"_metrics.csv"))
fwrite(scenario,paste0(PREFIX,"_scenario_metrics.csv"))
print(metrics)
print(scenario[endpoint=="DE"&grepl("undefined|near_second|equal_cv",scenario)])
