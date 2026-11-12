# Compare the analytic one-sided limit with the original zero-mass convention.
root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
lib<-file.path(root,"library");library(SGCB,lib.loc=lib);library(data.table);setDTthreads(1)
stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
source(file.path(root,"experiments/revision_calibration.R"));source(file.path(root,"experiments/known_depth_inputs.R"))
design <- readRDS(file.path(root,"experiments/design.rds"));template <- readRDS(file.path(root,"experiments/control_template.rds"))
# Use the current installed mean model; the reference convention is a score calculation only.
args<-commandArgs(TRUE); stage<-args[1]; ids<-as.integer(args[-1]); dir<-file.path(root,"results/boundary");dir.create(dir,showWarnings=FALSE)
threads<-as.integer(Sys.getenv("SGCB_THREADS","2"))
run_one<-function(i) {
    dest<-file.path(dir,sprintf("%s_%03d.rds",stage,i))
    if(file.exists(dest)) return(invisible(i))
    if(stage%in%c("technical","fresh_depth")){
        input_design<-design
        if(stage=="fresh_depth") {
            stream<-Reduce(function(s,j)parallel::nextRNGStream(s),seq_len(i),init=tail(design$streams,1)[[1]])
            id<-384L+i; input_design$jobs[id,]<-design$jobs[385,];input_design$jobs$id[id]<-id
            input_design$jobs$replicate[id]<-i;input_design$streams[[id]]<-stream
        }
        a<-known_depth_inputs(i,input_design,template);g<-a$g
        dat<-list(depth_null=list(x=sweep(a$raw,2,a$exposure,"/"),alt=rep(FALSE,a$G)),
            depth_mixed=list(x=sweep(a$raw_alt,2,a$exposure,"/"),alt=a$mu!=a$m1))
    } else {
        idx<-if(i<=8) 64L+i else 320L+i-8L; s<-revision_matrix(design$jobs[idx,],design,template);g<-s$group
        make<-function(mult){x<-s$x;x[,g==1]<-x[,g==1,drop=FALSE]*mult;list(x=x,alt=mult!=1)}
        G<-nrow(s$x); dat<-list(all_double=make(rep(2,G)),all_quadruple=make(rep(4,G)),balanced_fourfold=make(rep(c(.25,4),length.out=G)))
    }
    one<-function(z) {
        m<-SGCB:::.populationMean(z$x,as.integer(g),threads)
        av<-m$genes$available==1
        mass<-m$genes$posterior_null; ranked<-ifelse(av,mass,Inf)
        total<-sum(mass)
        old<-cumsum(mass[order(ranked)])[data.table::frank(ranked,ties.method="max")]/ifelse(total>0,total,1)
        old[!av]<-1
        candidates<-list(released=old,ratio_limit=m$genes$p,boundary_Welch=old)
        if(total==0 && any(av)){
            welch<-moment_workflow(z$x,g)$mean_p
            candidates$boundary_Welch<-ifelse(av & is.finite(welch),welch,1)
        }
        stat<-rbindlist(lapply(names(candidates),function(k){
            p<-candidates[[k]];calls<-p.adjust(p,"BH")<.05
            data.table(candidate=k,prior_null=m$control$prior_null,total_null=sum(m$genes$posterior_null),
                eligible=sum(av),TP=sum(calls&z$alt),FP=sum(calls&!z$alt),alternatives=sum(z$alt),nulls=sum(!z$alt),
                FDP=sum(calls&!z$alt)/max(1,sum(calls)),any=as.integer(any(calls)),
                power=if(any(z$alt))sum(calls&z$alt)/sum(z$alt) else NA_real_)
        }))
        list(stat=stat,p=candidates,available=av,truth=z$alt,ratio=NULL)
    }
    ans<-lapply(dat,one);saveRDS(ans,paste0(dest,".partial"));stopifnot(file.rename(paste0(dest,".partial"),dest))
    cat("Boundary",stage,i,format(Sys.time()),"\n");flush.console()
}
invisible(lapply(ids,run_one))