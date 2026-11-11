root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
library(data.table); setDTthreads(1)
out <- file.path(root,"results"); d <- readRDS(file.path(root,"experiments/design.rds"))
files <- file.path(out,"final",sprintf("matrix_%03d.rds",d$jobs$id))
stopifnot(all(file.exists(files)))
replicates <- lapply(files,readRDS)
ratio <- function(x,y) ifelse(y>0,x/y,NA_real_)
mean_ok <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
summarize_method <- function(run,method,z,endpoint,keep=rep(TRUE,nrow(z)),state="all",cutoff=.05) {
    tr <- run$truth
    p <- z[[paste0(endpoint,"_p")]][keep]; q <- z[[paste0(endpoint,"_q")]][keep]
    e <- tr[[paste0(endpoint,"_effect")]][keep]
    inf <- z[[paste0(endpoint,"_available")]][keep]
    alt <- abs(e)>1e-9; call <- is.finite(q)&q<cutoff
    n <- length(e); false <- sum(call & !alt); true <- sum(call & alt)
    ci <- rep(NA_real_,4)
    if(endpoint=="dv"){
        lo <- z$lower[keep]; hi <- z$upper[keep]; valid <- is.finite(lo)&is.finite(hi)
        in_ci <- valid & lo<=e & hi>=e
        score <- hi-lo+40*pmax(lo-e,0)+40*pmax(e-hi,0)
        ci <- c(mean(in_ci),ratio(sum(in_ci),sum(valid)),mean_ok((hi-lo)[valid]),mean_ok(score[valid]))
    }
    data.table(id=run$job$id,replicate=run$job$replicate,panel=run$job$panel,size=run$job$size,
        state=state,method=method,endpoint=endpoint,genes=n,nulls=sum(!alt),alternatives=sum(alt),available=sum(inf),
        tp=true,fp=false,calls=sum(call),FDP=false/max(1,sum(call)),any_call=as.numeric(any(call)),
        power=ratio(true,sum(alt)),type1=ratio(sum(is.finite(p)&p<.05 & !alt),sum(!alt)),
        availability=mean(inf),covered_all=ci[1],coverage=ci[2],width=ci[3],interval_score=ci[4],
        seconds=unname(run$seconds[if(method=="Moments") "moments" else "SGCB"]))
}
assemble <- function(run) {
    z <- run$SGCB; m <- run$moments; check <- z
    check$mean_p <- pmax(ifelse(is.finite(z$mean_p),z$mean_p,1),ifelse(is.finite(m$mean_p),m$mean_p,1))
    check$mean_q <- p.adjust(check$mean_p,"BH",n=nrow(check)); check$mean_available <- z$mean_available & m$mean_available
    methods <- list(SGCB=z,Moments=m)
    overall <- rbindlist(c(lapply(names(methods),function(nm) rbindlist(lapply(c("mean","dv"),function(ep)
        summarize_method(run,nm,methods[[nm]],ep)))),list(summarize_method(run,"Corroborated",check,"mean"))))
    labels <- unique(run$truth$state)
    strata <- rbindlist(lapply(labels,function(st) rbindlist(lapply(names(methods),function(nm)
        summarize_method(run,nm,methods[[nm]],"dv",run$truth$state==st,st)))))
    rbind(overall,strata)
}
metrics <- rbindlist(lapply(replicates,assemble))
fwrite(metrics,file.path(out,"replicate_metrics.tsv"),sep="\t",na="")
# Matrix-level intervals retain between-experiment variability and feature dependence.
measurements <- c("FDP","any_call","power","type1","availability","covered_all","coverage","width","interval_score","seconds")
rate_interval <- function(x,metric,panel,ep){
    x <- x[is.finite(x)]; n <- length(x)
    if(n==0L) return(list(value=NA_real_,lower=NA_real_,upper=NA_real_,mcse=NA_real_,matrices=0L,interval_method="not_defined"))
    av <- mean(x); se <- sd(x)/sqrt(n)
    bounded <- !metric%in%c("width","interval_score","seconds")
    zero_truth <- panel%in%c("null","equal_depth") || (panel=="depth" & ep=="mean")
    exact <- metric=="any_call" || (metric=="FDP" && zero_truth)
    boundary <- bounded && (all(x==0)|all(x==1))
    if(exact || boundary){
        lim <- binom.test(sum(x>0),n)$conf.int
        label <- if(exact) "exact_binomial_over_matrices" else "boundary_bound_via_nonzero_matrices"
    } else {
        lim <- av+c(-1,1)*qt(.975,n-1L)*se
        lim <- pmax(0,lim); if(bounded) lim <- pmin(1,lim)
        label <- "replicate_t_MC_interval"
    }
    list(value=av,lower=lim[1],upper=lim[2],mcse=se,matrices=n,interval_method=label)
}
pooled <- melt(metrics,id.vars=c("id","replicate","panel","size","state","method","endpoint"),measure.vars=measurements,variable.name="metric",variable.factor=FALSE)
summary <- pooled[,rate_interval(value,metric[1],panel[1],endpoint[1]),by=.(panel,size,state,method,endpoint,metric)]
fwrite(summary,file.path(out,"independent_operating_characteristics.tsv"),sep="\t",na="")
# Two reporting cutoffs applied to the same fitted, full-family arithmetic-mean scores.
primary_means <- copy(metrics[state=="all" & endpoint=="mean" & method=="SGCB"])
primary_means[,cutoff:=.05]
stricter_means <- rbindlist(lapply(replicates,function(run)
    summarize_method(run,"SGCB",run$SGCB,"mean",cutoff=.01)))
stricter_means[,cutoff:=.01]
fields <- c("panel","size","replicate","alternatives","tp","fp","FDP","power","any_call","cutoff")
mean_screens <- rbind(primary_means[,..fields],stricter_means[,..fields])[,study:="primary_matrices"]

# Use the same two cutoffs in the completed boundary and genuine-shift studies.
# The original technical-null matrices are already in the primary panel.
boundary_files <- list.files(file.path(out,"boundary"),pattern="[.]rds$",full.names=TRUE)
boundary_screens <- rbindlist(lapply(boundary_files,function(f) {
    stage <- sub("_[0-9]+[.]rds$","",basename(f))
    replicate_id <- as.integer(sub(".*_([0-9]+)[.]rds$","\\1",basename(f)))
    assays <- readRDS(f)
    if(stage=="technical") assays <- assays["depth_mixed"]
    rbindlist(lapply(names(assays),function(arm) {
        assay <- assays[[arm]]
        q <- p.adjust(assay$p$ratio_limit,"BH")
        cutoffs <- c(.05,.01)
        calls <- outer(q,cutoffs,"<") & is.finite(q)
        alt <- assay$truth
        tp <- colSums(calls[alt,,drop=FALSE]); fp <- colSums(calls[!alt,,drop=FALSE])
        name <- if(arm=="depth_null") "depth" else arm
        design <- if(stage=="positive" && replicate_id>8L) "50_50" else "10_10"
        data.table(panel=name,size=design,replicate=replicate_id,alternatives=sum(alt),tp=tp,fp=fp,
            FDP=fp/pmax(tp+fp,1),power=if(any(alt)) tp/sum(alt) else NA_real_,
            any_call=as.integer(tp+fp>0),cutoff=cutoffs,study=stage)
    }))
}))
mean_screens <- rbind(mean_screens,boundary_screens,use.names=TRUE)
thresholds <- mean_screens[,{
    fdr <- rate_interval(FDP,"FDP",panel[1L],"mean")
    sensitivity <- rate_interval(power,"power",panel[1L],"mean")
    family <- rate_interval(any_call,"any_call",panel[1L],"mean")
    list(matrices=.N,alternatives=sum(alternatives),true_calls=sum(tp),false_calls=sum(fp),
        mean_FDP=fdr$value,FDP_low=fdr$lower,FDP_high=fdr$upper,
        mean_power=sensitivity$value,power_low=sensitivity$lower,power_high=sensitivity$upper,
        matrices_with_calls=sum(any_call),any_probability=family$value,any_low=family$lower,any_high=family$upper)
},by=.(study,panel,size,cutoff)]
setorder(thresholds,study,panel,size,-cutoff)
fwrite(thresholds,file.path(out,"mean_screen_thresholds.tsv"),sep="\t",na="")
held <- metrics[state=="all" & endpoint=="mean" & (replicate>fifelse(panel%in%c("depth","equal_depth"),8L,16L))]
hs <- melt(held,id.vars=c("id","replicate","panel","size","state","method","endpoint"),measure.vars=c("FDP","power","any_call","type1"),variable.factor=FALSE)
holdsummary <- hs[,rate_interval(value,variable[1],panel[1],endpoint[1]),by=.(panel,size,method,endpoint,variable)]
wide <- dcast(hs,id+replicate+panel+size+variable~method,value.var="value")
wide[,diff:=Corroborated-SGCB]
paired <- wide[is.finite(diff),.(matrices=.N,delta=mean(diff),mcse=sd(diff)/sqrt(.N),
    lower=mean(diff)-qt(.975,.N-1)*sd(diff)/sqrt(.N),upper=mean(diff)+qt(.975,.N-1)*sd(diff)/sqrt(.N)),by=.(panel,size,variable)]
fwrite(holdsummary,file.path(out,"mean_candidate_holdout.tsv"),sep="\t",na="")
fwrite(paired,file.path(out,"mean_candidate_paired_differences.tsv"),sep="\t",na="")
# Compare fresh candidate and reference on identical inputs and all retained numeric endpoints.
ref_paths <- list.files(file.path(out,"reference"),pattern="\\.rds$",full.names=TRUE)
checks <- rbindlist(lapply(ref_paths,function(f) {
    ref <- readRDS(f); current <- replicates[[ref$job$id]]
    stopifnot(identical(ref$truth,current$truth))
    rbindlist(lapply(names(ref$SGCB),function(nm){
        a <- as.numeric(ref$SGCB[[nm]]); b <- as.numeric(current$SGCB[[nm]])
        data.table(id=ref$job$id,field=nm,missing_disagreements=sum(xor(is.na(a),is.na(b))),
            max_absolute_delta=max(abs(a[is.finite(a)&is.finite(b)]-b[is.finite(a)&is.finite(b)]),0),
            differing_calls=if(nm%in%c("mean_q","dv_q"))sum(xor(is.finite(a)&a<.05,is.finite(b)&b<.05)) else NA_integer_)
    }))
}))
fwrite(checks,file.path(out,"new_input_reference_agreement.tsv"),sep="\t",na="")
cat("Completed",length(replicates),"independently generated matrices and",length(ref_paths),"separate reference fits.\n")