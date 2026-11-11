root<-normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE);p<-file.path(root,"results");library(data.table)
dir<-file.path(p,"boundary");files<-list.files(dir,pattern="\\.rds$",full.names=TRUE);stopifnot(length(files)==72)
z<-rbindlist(lapply(files,function(file){
    results<-readRDS(file);t<-rbindlist(lapply(results,"[[","stat"),idcol="arm")
    t[,':='(stage=sub("_[0-9]+\\.rds$","",basename(file)),replicate=as.integer(sub(".*_([0-9]+)\\.rds$","\\1",basename(file))))]
    t
}))
fwrite(z,file.path(p,"boundary_matrices.tsv"),sep="\t",na="")
a<-z[,.(matrices=.N,boundary=sum(total_null==0),any=sum(any),family_probability=mean(any),
    family_lower=binom.test(sum(any),.N)$conf.int[1],family_upper=binom.test(sum(any),.N)$conf.int[2],
    true_calls=sum(TP),false_calls=sum(FP),mean_FDP=mean(FDP),FDP_mcse=sd(FDP)/sqrt(.N),
    power=if(all(is.na(power)))NA_real_ else mean(power),power_mcse=sd(power)/sqrt(.N)),by=.(stage,arm,candidate)]
fwrite(a,file.path(p,"boundary_operating_characteristics.tsv"),sep="\t",na="")
b<-dcast(z,stage+arm+replicate~candidate,value.var=c("FDP","power","TP","FP"))
paired<-b[,.(datasets=.N,delta_FDP=mean(FDP_ratio_limit-FDP_released),
    se_FDP=sd(FDP_ratio_limit-FDP_released)/sqrt(.N),delta_power=mean(power_ratio_limit-power_released),
    se_power=sd(power_ratio_limit-power_released)/sqrt(.N),delta_TP=sum(TP_ratio_limit-TP_released),delta_FP=sum(FP_ratio_limit-FP_released)),
    by=.(stage,arm)]
fwrite(paired,file.path(p,"boundary_paired_differences.tsv"),sep="\t",na="")
# Direct calls through the installed release confirm the independent formula implementation.
source(file.path(root,"experiments/revision_calibration.R"));source(file.path(root,"experiments/known_depth_inputs.R"))
library(SGCB,lib.loc=file.path(root,"library"));design<-readRDS(file.path(root,"experiments/design.rds"))
template<-readRDS(file.path(root,"experiments/control_template.rds"));tr<-as.integer(Sys.getenv("SGCB_THREADS","3"))
verify<-function(stage,i,arm){
    if(stage=="positive"){
        sim<-revision_matrix(design$jobs[64L+i,],design,template);x<-sim$x;g<-sim$group;x[,g==1]<-x[,g==1]*2
    } else {
        input<-design
        if(stage=="fresh_depth"){id<-384L+i;input$jobs[id,]<-design$jobs[385,];input$jobs$id[id]<-id;input$jobs$replicate[id]<-i
            input$streams[[id]]<-Reduce(function(s,j)parallel::nextRNGStream(s),seq_len(i),init=tail(design$streams,1)[[1]])}
        a<-known_depth_inputs(i,input,template);g<-a$g;x<-sweep(if(arm=="depth_null") a$raw else a$raw_alt,2,a$exposure,"/")
    }
    fit<-sgcbDE(x,g,min_count=0,min_samples=1,normalization="none",endpoint="mean",threads=tr)
    original<-readRDS(file.path(dir,sprintf("%s_%03d.rds",stage,i)))[[arm]]$p$ratio_limit
    data.table(stage=stage,rep=i,arm=arm,genes=nrow(fit),max_score_difference=max(abs(fit$pvalue-original)),
        different_BH_calls=sum(xor(fit$DE,p.adjust(original,"BH")<.05)))
}
fresh<-z[stage=="fresh_depth"&arm=="depth_null"&total_null==0,unique(replicate)]
verified<-rbindlist(list(verify("technical",1L,"depth_null"),verify("positive",1L,"all_double"),verify("fresh_depth",fresh[1],"depth_null")))
stopifnot(all(verified$different_BH_calls==0),all(verified$max_score_difference<1e-12))
fwrite(verified,file.path(p,"installed_boundary_verification.tsv"),sep="\t")
