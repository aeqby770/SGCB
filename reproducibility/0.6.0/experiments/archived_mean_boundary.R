root<-normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE);library(data.table)
r<-file.path(root,"study_data/results"); files<-list.files(r,pattern="^sgcb_0[.]5[.][79]_.*[.]rds$",full.names=TRUE)
read_one<-function(file) {
    z<-readRDS(file)
    if(!is.data.frame(z) || !all(c("gene_id","pvalue","padj")%in%names(z)))return(NULL)
    post<-attr(z,"mean_posterior");model<-attr(z,"mean_model")
    available<-if("de_inference_available"%in%names(z)) z$de_inference_available else is.finite(z$pvalue)
    prior<-model$control$prior_null
    total<-if(!is.null(post$posterior_null))sum(post$posterior_null) else NA_real_
    data.table(source=basename(file),genes=nrow(z),available=sum(available),calls=sum(z$padj<.05,na.rm=TRUE),
        nonzero_score_available=sum(z$pvalue[available]>0,na.rm=TRUE),total_posterior_null_mass=total,
        fitted_prior_null=if(is.null(prior))NA_real_ else prior,
        zero_mass_boundary=if(is.finite(total)) total==0 else if(any(z$pvalue[available]>0,na.rm=TRUE)) FALSE else NA)
}
a<-rbindlist(lapply(files,read_one));fwrite(a,file.path(root,"results/archived_mean_boundary.tsv"),sep="\t",na="")
print(a[,.(fits=.N,zero_mass=sum(zero_mass_boundary==TRUE,na.rm=TRUE),indeterminate=sum(is.na(zero_mass_boundary)))])

# Check whether the current constant-positive eligibility correction affects retained study inputs.
library(data.table);data.table::setDTthreads(3)
path<-file.path(root,"study_data/results")
files<-list.files(path,pattern="^sgcb_0[.]5[.]9_(channel_n[0-9_]+|zero_n[0-9_]+|readcount_null_[0-9]+(_streamed)?|kirc|sle_classical_monocytes_2021|masld_liver_2023)[.]rds$",full.names=TRUE)
stopifnot(length(files)==19L)
check<-function(file) {
    fit<-readRDS(file);input<-attr(fit,"input_file")
    relative<-sub("^.*?/data/","data/",input,perl=TRUE)
    d<-readRDS(file.path(root,"study_data",relative))
    x<-switch(as.character("counts"%in%names(d)),"TRUE"=as.matrix(d$counts),"FALSE"=as.matrix(d$abundance))
    samples<-attr(fit,"samples");sf<-attr(fit,"size_factors")
    x<-x[fit$gene_id,as.character(samples$sample),drop=FALSE]
    stopifnot(ncol(x)==length(sf),all(is.finite(sf)),all(sf>0))
    x<-sweep(x,2,sf,"/");g<-factor(samples$group)
    constant<-function(y) {
        p<-rowSums(y>0)
        lo<-matrixStats::rowMins(replace(y,y<=0,Inf))
        hi<-matrixStats::rowMaxs(replace(y,y<=0,-Inf))
        p>=4 & lo==hi
    }
    a<-constant(x[,g==levels(g)[1],drop=FALSE]);b<-constant(x[,g==levels(g)[2],drop=FALSE])
    affected<-(a|b)&fit$dv_inference_available
    data.table(record=basename(file),genes=nrow(fit),constant_positive_control=sum(a),constant_positive_treatment=sum(b),
        previously_available_constant_groups=sum(affected),affected_gene_ids=paste(fit$gene_id[affected],collapse=";"))
}
answer<-rbindlist(lapply(files,check))
fwrite(answer,file.path(root,"results/retained_input_eligibility.tsv"),sep="\t",na="")
print(answer[,.(record,previously_available_constant_groups)])
cat("Total retained DV records requiring eligibility updates:",sum(answer$previously_available_constant_groups),"\n")
m<-fread(file.path(root,"results/archived_mean_boundary.tsv"))
stopifnot(nrow(m)==31L,!any(is.na(m$zero_mass_boundary)),!any(m$zero_mass_boundary))
cat("All 31 saved mean-fit records remain outside the zero-null-mass boundary.\n")

