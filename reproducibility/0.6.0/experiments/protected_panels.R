root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); lib <- file.path(root,"library")
library(SGCB,lib.loc=lib); library(data.table)
stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
ref <- file.path(root,"reference_0.5.9"); dest <- file.path(root,"results","protected")
dir.create(dest,showWarnings=FALSE); threads <- as.integer(Sys.getenv("SGCB_THREADS","3"))
plan <- data.frame(panel=rep(c("zero","channel"),3),n0=c(10L,10L,10L,10L,50L,50L),n1=c(10L,10L,50L,50L,50L,50L))
run_panel <- function(j) {
    p <- plan[j,]; id <- sprintf("%s_n%i_%i",p$panel,p$n0,p$n1)
    out <- file.path(dest,paste0(id,".tsv")); if(file.exists(out)) return(invisible(j))
    source <- if(p$panel=="zero") "zero_calibration.rds" else "channel_simulation.rds"
    d <- readRDS(file.path(ref,"data",source)); g <- factor(d$group)
    index <- c(head(which(g==levels(g)[1]),p$n0),head(which(g==levels(g)[2]),p$n1))
    old <- readRDS(file.path(ref,"results",paste0("sgcb_0.5.9_",id,".rds")))
    if(!is.null(attr(old,"initial_rng_state"))) assign(".Random.seed",attr(old,"initial_rng_state"),.GlobalEnv)
    z <- sgcbDE(d$abundance[,index,drop=FALSE],g[index],min_count=0,min_samples=1,normalization="none",
        offset=if(p$panel=="zero") .5 else 0,threads=threads)
    stopifnot(identical(z$gene_id,old$gene_id))
    columns <- c("pvalue","padj","log2FoldChange","dv_pvalue","dv_padj","dv_log2_cv2_ratio","dv_ci_lower","dv_ci_upper",
        "de_inference_available","dv_inference_available","DE","DV")
    tab <- rbindlist(lapply(columns,function(col) {
        a <- as.numeric(old[[col]]); b <- as.numeric(z[[col]])
        same <- is.finite(a)&is.finite(b)
        data.table(panel=id,genes=length(a),field=col,finite_both=sum(same),missing_disagreements=sum(xor(is.na(a),is.na(b))),
            max_absolute_delta=if(any(same)) max(abs(a[same]-b[same])) else NA_real_,
            differing_entries=sum(a[same]!=b[same]),old_sum=if(is.logical(old[[col]]))sum(a,na.rm=TRUE) else NA_real_,
            new_sum=if(is.logical(z[[col]]))sum(b,na.rm=TRUE) else NA_real_)
    }))
    fwrite(tab,out,sep="\t")
    cat("protected",id,"call differences",sum(tab[field%in%c("DE","DV"),differing_entries]),"\n");flush.console()
}
invisible(lapply(seq_len(nrow(plan)),run_panel))