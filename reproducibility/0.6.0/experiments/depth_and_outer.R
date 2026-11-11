root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); library(data.table); setDTthreads(1)
out <- file.path(root,"results"); design <- readRDS(file.path(root,"experiments/design.rds"))
args <- commandArgs(TRUE); stage <- match.arg(args[1],c("depth","outer"))
if(stage=="depth") {
    source(file.path(root,"experiments/revision_calibration.R"))
    template <- readRDS(file.path(root,"experiments/control_template.rds"))
    jobs <- design$jobs[design$jobs$panel%in%c("depth","equal_depth"),]
    calc <- function(i) {
        job <- jobs[i,]; sim <- revision_matrix(job,design,template)
        s <- c(1,if(job$panel=="depth") .25 else 1)
        block <- lapply(0:1,function(g){
            x <- sim$x[,sim$group==g,drop=FALSE];m <- rowMeans(x);v <- matrixStats::rowVars(x)
            cbind(mean=m,raw_cv=v/m^2,latent_cv=(v-m/s[g+1])/m^2)
        })
        good <- is.finite(block[[1]][,"raw_cv"]) & is.finite(block[[2]][,"raw_cv"])
        # Compare raw and noise-subtracted moment estimators on identical finite genes. Negative
        # unbiased-variance estimates are kept, not truncated or silently removed.
        observed_delta <- block[[2]][,"raw_cv"]-block[[1]][,"raw_cv"]
        latent_delta <- block[[2]][,"latent_cv"]-block[[1]][,"latent_cv"]
        data.table(id=job$id,panel=job$panel,replicate=job$replicate,genes=nrow(sim$x),
            finite=sum(good),observed_contrast_median=median(observed_delta[good]),
            corrected_contrast_median=median(latent_delta[good]),mean_contrast=mean(block[[2]][,"mean"]-block[[1]][,"mean"]),
            corrected_negative0=sum(block[[1]][good,"latent_cv"]<0),corrected_negative1=sum(block[[2]][good,"latent_cv"]<0),
            theoretical_log2_cv_ratio_median=median(sim$truth$dv_effect),theoretical_mean_contrast=0,theoretical_latent_contrast=0)
    }
    num <- rbindlist(lapply(seq_len(nrow(jobs)),calc)); fwrite(num,file.path(out,"observation_depth_moments.tsv"),sep="\t")
} else {
    ids <- c(2:13,258:263)
    paths <- file.path(out,"outer",sprintf("outer_%03d.rds",ids));stopifnot(all(file.exists(paths)))
    details <- rbindlist(lapply(lapply(paths,readRDS),function(z) {
        f <- as.data.table(z$features);f[,':='(matrix_id=z$id,size=z$size,B=z$B)]
        f[,refit_intervals_covering:=rowSums(z$lowers<=f$dv_effect & z$uppers>=f$dv_effect,na.rm=TRUE)]
        f
    }))
    fwrite(details,file.path(out,"outer_refit_probes.tsv"),sep="\t",na="")
    interval_sum <- function(lo,hi,e) {
        valid <- is.finite(lo)&is.finite(hi)
        list(n=length(lo),available=sum(valid),coverage=if(any(valid)) mean(lo[valid]<=e[valid] & hi[valid]>=e[valid]) else NA_real_,
            joint_coverage=mean(valid & lo<=e & hi>=e),width=if(any(valid))mean(hi[valid]-lo[valid]) else NA_real_,
            interval_score=if(any(valid))mean(hi[valid]-lo[valid]+40*pmax(lo[valid]-e[valid],0)+40*pmax(e[valid]-hi[valid],0)) else NA_real_)
    }
    # Matched probes for interval-quality comparisons; preserve unavailable counts in the feature table.
    details[,both:=is.finite(dv_ci_lower)&is.finite(dv_ci_upper)&is.finite(sensitivity_lower)&is.finite(sensitivity_upper)]
    retained <- details[both==TRUE]
    old <- retained[,interval_sum(dv_ci_lower,dv_ci_upper,dv_effect),by=.(matrix_id,size,state)][,interval:="conditional"]
    new <- retained[,interval_sum(sensitivity_lower,sensitivity_upper,dv_effect),by=.(matrix_id,size,state)][,interval:="outer_scale"]
    tabs <- rbind(old,new)
    fwrite(tabs,file.path(out,"outer_refit_matrices.tsv"),sep="\t",na="")
    wide <- dcast(tabs,matrix_id+size+state ~ interval,value.var=c("coverage","width","interval_score"))
    ans <- wide[,.(datasets=.N,conditional_coverage=mean(coverage_conditional),outer_coverage=mean(coverage_outer_scale),
        coverage_delta=mean(coverage_outer_scale-coverage_conditional),
        coverage_delta_se=sd(coverage_outer_scale-coverage_conditional)/sqrt(.N),
        conditional_width=mean(width_conditional),outer_width=mean(width_outer_scale),
        conditional_score=mean(interval_score_conditional),outer_score=mean(interval_score_outer_scale),
        score_delta_se=sd(interval_score_outer_scale-interval_score_conditional)/sqrt(.N)),by=.(size,state)]
    mechanism <- details[,.(probes=.N,matched=sum(both),conditional_available=sum(is.finite(dv_ci_lower)&is.finite(dv_ci_upper)),
        below_four=sum(count0<4 | count1<4),conditional_below_truth=sum(dv_ci_upper<dv_effect,na.rm=TRUE),
        median_outer_sd=median(outer_sd,na.rm=TRUE),median_conditional_se=median(dv_se_log2_cv2,na.rm=TRUE),
        median_fitted_effect=median(dv_log2_cv2_ratio,na.rm=TRUE),mean_refits_covering=mean(refit_intervals_covering[both]/B[both])),
        by=.(size,state)]
    fwrite(merge(ans,mechanism,by=c("size","state"),all=TRUE),file.path(out,"outer_refit_summary.tsv"),sep="\t",na="")
}