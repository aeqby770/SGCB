root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
lib <- file.path(root,"library"); library(SGCB,lib.loc=lib)
stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
source(file.path(root,"experiments/revision_calibration.R"))
design <- readRDS(file.path(root,"experiments/design.rds")); template <- readRDS(file.path(root,"experiments/control_template.rds"))
threads <- as.integer(Sys.getenv("SGCB_THREADS","3")); jobs <- as.integer(commandArgs(TRUE))
out <- file.path(root,"results","outer"); dir.create(out,showWarnings=FALSE)
refit_one <- function(id) {
    dest <- file.path(out,sprintf("outer_%03d.rds",id))
    if(file.exists(dest)) return(invisible(id))
    sim <- revision_matrix(design$jobs[id,],design,template)
    # Fixed, balanced probes; simulated gene labels select evaluation strata only.
    probe <- sort(unlist(lapply(split(seq_len(nrow(sim$x)),sim$truth$state),head,12L),use.names=FALSE))
    background <- setdiff(seq_len(nrow(sim$x)),probe)
    assign(".Random.seed",parallel::nextRNGSubStream(design$streams[[id]]),envir=.GlobalEnv)
    B <- 16L
    backgrounds <- replicate(B,sample(background,length(background),replace=TRUE))
    fit_data <- function(indices) {
        x <- sim$x[indices,,drop=FALSE]; rownames(x) <- sprintf("unit%04d",seq_len(nrow(x)))
        z <- sgcbDE(x,sim$group,min_count=0,min_samples=1,normalization="none",
            endpoint="variability",threads=threads)[seq_along(probe),]
        z[,c("dv_log2_cv2_ratio","dv_se_log2_cv2","dv_ci_lower","dv_ci_upper","dv_pvalue")]
    }
    base <- fit_data(c(probe,background))
    boot <- lapply(asplit(backgrounds,2),function(indices) fit_data(c(probe,indices)))
    centres <- do.call(cbind,lapply(boot,"[[","dv_log2_cv2_ratio"))
    lowers <- do.call(cbind,lapply(boot,"[[","dv_ci_lower"))
    uppers <- do.call(cbind,lapply(boot,"[[","dv_ci_upper"))
    fits_available <- rowSums(is.finite(lowers)&is.finite(uppers))
    outer_sd <- matrixStats::rowSds(centres,na.rm=TRUE)
    outer_sd[fits_available<B] <- NA_real_
    # Sensitivity band: add between-refit scale in quadrature to each original half-width.
    # Retain the original endpoints when between-background variance is zero.
    low_distance <- abs(base$dv_log2_cv2_ratio-base$dv_ci_lower)
    high_distance <- abs(base$dv_ci_upper-base$dv_log2_cv2_ratio)
    extra <- stats::qnorm(.975)*outer_sd
    sensitivity_lower <- base$dv_ci_lower-(sqrt(low_distance^2+extra^2)-low_distance)
    sensitivity_upper <- base$dv_ci_upper+(sqrt(high_distance^2+extra^2)-high_distance)
    ans <- list(id=id,size=design$jobs$size[id],B=B,features=cbind(sim$truth[probe,],base,
        available_refits=fits_available,outer_sd=outer_sd,sensitivity_lower=sensitivity_lower,
        sensitivity_upper=sensitivity_upper),estimates=centres,lowers=lowers,uppers=uppers)
    saveRDS(ans,paste0(dest,".partial")); stopifnot(file.rename(paste0(dest,".partial"),dest))
    cat("outer",id,"completed",format(Sys.time()),"\n");flush.console()
}
invisible(lapply(jobs,refit_one))