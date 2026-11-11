root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
library(SGCB,lib.loc=file.path(root,"library"));library(data.table);setDTthreads(1)
source(file.path(root,"experiments/revision_calibration.R")); design <- readRDS(file.path(root,"experiments/design.rds"))
template <- readRDS(file.path(root,"experiments/control_template.rds"));out <- file.path(root,"results","exposure_final");dir.create(out,showWarnings=FALSE)
threads <- as.integer(Sys.getenv("SGCB_THREADS","2")); args <- as.integer(commandArgs(TRUE))
source(file.path(root,"experiments/known_depth_inputs.R"))
one <- function(rep_id) {
    dest <- file.path(out,sprintf("exposure_%02d.tsv",rep_id));if(file.exists(dest))return(invisible(rep_id))
    inputs <- known_depth_inputs(rep_id,design,template)
    id <- inputs$id; raw <- inputs$raw; raw_alt <- inputs$raw_alt; g <- inputs$g; exposure <- inputs$exposure
    G <- inputs$G; mu <- inputs$mu; m1 <- inputs$m1; alt_stream <- inputs$alt_stream
    process <- function(condition,counts,means_differ,ref_stream) {
        assign(".Random.seed",ref_stream,.GlobalEnv)
        # Exact thinning of Poisson observations to a common known exposure.
        rate <- rep(min(exposure)/exposure,each=G)
        thin <- matrix(rbinom(length(counts),size=c(counts),prob=rate),nrow=G,dimnames=dimnames(counts))/min(exposure)
        org <- sweep(counts,2,exposure,"/")
        fit_mean <- function(x) {
            sgcbDE(x,g,min_count=0,min_samples=1,normalization="none",endpoint="mean",threads=threads)$padj
        }
        p <- list(supplied_depth=fit_mean(org),common_depth=fit_mean(thin))
        rbindlist(lapply(names(p),function(arm){
            call <- is.finite(p[[arm]]) & p[[arm]]<.05
            tp <- sum(call & means_differ);fp <- sum(call&!means_differ)
            data.table(replicate=rep_id,condition=condition,arm=arm,mean_alternatives=sum(means_differ),
                tp=tp,fp=fp,FDP=fp/max(1,sum(call)),any=as.integer(any(call)),power=if(any(means_differ))tp/sum(means_differ) else NA_real_)
        }))
    }
    # A synthetic-equivalence control uses exactly the same null data as the original depth experiment.
    null_result <- process("technical_only",raw,rep(FALSE,G),parallel::nextRNGSubStream(alt_stream))
    alt_result <- process("technical_and_biological",raw_alt,m1!=mu,parallel::nextRNGSubStream(parallel::nextRNGSubStream(alt_stream)))
    fwrite(rbind(null_result,alt_result),dest,sep="\t")
}
invisible(lapply(args,one))