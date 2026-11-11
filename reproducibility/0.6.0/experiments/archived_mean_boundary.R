root<-normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE);library(data.table)
r<-file.path(root,"reference_0.5.9/results"); files<-list.files(r,pattern="^sgcb_0[.]5[.][79]_.*[.]rds$",full.names=TRUE)
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
