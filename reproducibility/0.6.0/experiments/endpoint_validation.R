root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); lib <- file.path(root,"library")
library(SGCB,lib.loc=lib);library(data.table);setDTthreads(1)
source(file.path(root,"experiments/revision_calibration.R"))
d <- readRDS(file.path(root,"experiments/design.rds"));t <- readRDS(file.path(root,"experiments/control_template.rds"))
threads <- as.integer(Sys.getenv("SGCB_THREADS","3"))
cases <- c(1L,65L,129L,193L,257L,321L)
check_one <- function(i) {
    id <- cases[i]
    job <- d$jobs[match(id,d$jobs$id),,drop=FALSE]
    sim <- revision_matrix(job,d,t); labels <- c("both","mean","variability")
    order <- labels[(seq_along(labels)+i-2L)%%3L+1L]
    run <- lapply(order,function(e) {
        invisible(gc(FALSE))
        wall <- system.time(z <- sgcbDE(sim$x,sim$group,min_count=0,min_samples=1,normalization="none",threads=threads,endpoint=e))[["elapsed"]]
        list(z=z,seconds=wall)
    }); names(run) <- order
    rbindlist(lapply(labels,function(e){
        z <- run[[e]]$z; full <- run$both$z
        keys <- switch(e,both=c("padj","dv_padj"),mean=c("log2FoldChange","lfcSE","pvalue","padj","DE"),
            variability=c("dv_log2_cv2_ratio","dv_se_log2_cv2","dv_pvalue","dv_padj","dv_ci_lower","dv_ci_upper","DV"))
        diffs <- vapply(keys,function(k) max(c(0,abs(as.numeric(z[[k]])-as.numeric(full[[k]]))),na.rm=TRUE),numeric(1))
        miss <- sum(vapply(keys,function(k) sum(xor(is.na(z[[k]]),is.na(full[[k]]))),numeric(1)))
        calls <- intersect(keys,c("DE","DV")); changed <- sum(vapply(calls,function(k) sum(z[[k]]!=full[[k]],na.rm=TRUE),numeric(1)))
        data.table(id=id,design=job$panel,size=job$size,mode=e,run_order=match(e,order),threads=threads,
            seconds=run[[e]]$seconds,speedup=run$both$seconds/run[[e]]$seconds,call_changes=changed,missingness_changes=miss,
            max_score_delta=max(diffs),genes=nrow(z))
    }))
}
ans <- rbindlist(lapply(seq_along(cases),check_one))
fwrite(ans,file.path(root,"results","endpoint_runtime_equivalence.tsv"),sep="\t")
cat("Paired comparisons",nrow(ans),"calls changed",sum(ans$call_changes),"\n")