root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); args <- commandArgs(TRUE)
edition <- args[1]; lib <- if(edition=="reference") file.path(root,"reference_0.5.9/library") else file.path(root,"library")
threads <- as.integer(args[2]); panels <- args[-(1:2)]
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads)
data.table::setDTthreads(threads);library(SGCB,lib.loc=lib)
ref <- file.path(root,"reference_0.5.9");out <- file.path(root,"results","protected")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
match_one <- function(panel) {
    old <- readRDS(file.path(ref,"results",paste0("sgcb_0.5.9_",panel,".rds")))
    kind <- substr(panel,1,regexpr("_",panel)-1); size <- as.integer(strsplit(sub(".*_n","",panel),"_")[[1]])
    data <- readRDS(file.path(ref,"data",if(kind=="channel")"channel_simulation.rds" else "zero_calibration.rds"))
    g <- factor(data$group); idx <- c(head(which(g==levels(g)[1]),size[1]),head(which(g==levels(g)[2]),size[2]))
    z <- sgcbDE(data$abundance[,idx,drop=FALSE],g[idx],min_count=0,min_samples=1,normalization="none",offset=if(kind=="channel")0 else .5,threads=threads)
    stopifnot(identical(z$gene_id,old$gene_id))
    target <- c("gene_id","pvalue","padj","DE","dv_log2_cv2_ratio","dv_pvalue","dv_padj","dv_ci_lower","dv_ci_upper","dv_inference_available","DV")
    res <- data.frame(z[,target]); saveRDS(res,file.path(out,sprintf("%s_%s_t%02d.rds",panel,edition,threads)))
    data.table::fwrite(data.frame(edition=edition,threads=threads,panel=panel,features=nrow(z),mean_calls=sum(z$DE),dv_calls=sum(z$DV),
        mean_call_changes=sum(xor(z$DE,old$DE)),dv_call_changes=sum(xor(z$DV,old$DV)),max_DV_delta=max(abs(z$dv_padj-old$dv_padj),na.rm=TRUE)),
        file.path(out,sprintf("%s_%s_t%02d.tsv",panel,edition,threads)),sep="\t")
}
invisible(lapply(panels,match_one))