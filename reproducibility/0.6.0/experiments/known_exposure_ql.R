# Same simulated integer counts and known exposures as the thinning comparison.
root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE);library(data.table);setDTthreads(1)
source(file.path(root,"experiments/revision_calibration.R"));source(file.path(root,"experiments/known_depth_inputs.R"))
d <- readRDS(file.path(root,"experiments/design.rds"));template <- readRDS(file.path(root,"experiments/control_template.rds"))
score <- function(counts,group,exposure,alt) {
    keep <- rowSums(counts)>0; design <- model.matrix(~group)
    y <- edgeR::DGEList(counts=counts[keep,,drop=FALSE],lib.size=1e6*exposure,norm.factors=rep(1,length(exposure)))
    y <- edgeR::estimateDisp(y,design,robust=TRUE)
    fit <- edgeR::glmQLFit(y,design,robust=TRUE)
    p <- rep(1,nrow(counts)); p[keep] <- edgeR::glmQLFTest(fit,coef=2)$table$PValue
    call <- p.adjust(p,"BH")<.05
    tp<-sum(call&alt); fp<-sum(call&!alt)
    list(mean_alternatives=sum(alt),tp=tp,fp=fp,FDP=fp/max(1,sum(call)),any=as.integer(any(call)),power=if(any(alt))tp/sum(alt) else NA_real_)
}
each <- function(i) {
    a<-known_depth_inputs(i,d,template)
    b <- rbindlist(list(technical_only=score(a$raw,a$g,a$exposure,rep(FALSE,a$G)),
        technical_and_biological=score(a$raw_alt,a$g,a$exposure,a$m1!=a$mu)),idcol="condition")
    b[,':='(replicate=i,arm="edgeR_known_offset",method_version=as.character(packageVersion("edgeR")))]
    b
}
x<-rbindlist(lapply(1:24,each));fwrite(x,file.path(root,"results/known_exposure_ql_matrices.tsv"),sep="\t",na="")
s<-x[,.(replicates=.N,any_datasets=sum(any),any_probability=mean(any),any_lower=binom.test(sum(any),.N)$conf.int[1],
    any_upper=binom.test(sum(any),.N)$conf.int[2],mean_FDP=mean(FDP),mean_FDP_se=sd(FDP)/sqrt(.N),
    mean_power=if(all(is.na(power)))NA_real_ else mean(power),power_se=sd(power)/sqrt(.N),false_calls=sum(fp),true_calls=sum(tp)),
    by=.(condition,arm)]
fwrite(s,file.path(root,"results/known_exposure_ql_summary.tsv"),sep="\t",na="");print(s)