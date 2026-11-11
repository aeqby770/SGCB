# Independent positive-null-weight limit using the R1 likelihood.
root<-normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE);lib<-file.path(root,"reference_0.5.9/library")
library(SGCB,lib.loc=lib); env<-new.env(parent=asNamespace("SGCB"))
txt<-paste(readLines(file.path(root,"reference_0.5.9/software/SGCB/R/population.R"),encoding="UTF-8"),collapse="\n")
marker<-"    moments <- summary\n    genes <- data.frame(effect=direction*moments"
addition<-'    eta <- 1e-8
    nuisance <- rowSums(matrix(joint_fit$weights,ncol=2L))
    z <- anchored_location_summary_cpp(joint_kernel$dictionary,c(eta*nuisance,(1-eta)*nuisance),.05,threads,control$probability_tolerance)[,"posterior_null"]
    observed <- ifelse(available,z,Inf)
    control$positive_weight_reference <- cumsum(z[order(observed)])[data.table::frank(observed,ties.method="max")]/sum(z)
    control$positive_weight_reference[!available] <- 1
'
expr<-parse(text=sub(marker,paste0(addition,marker),txt,fixed=TRUE))
i<-which(vapply(expr,function(e) is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name(".populationMean")),logical(1)))
eval(expr[[i]],env)
x<-matrix(exp(sin(seq_len(240))/2+seq_len(240)/1000),nrow=12);x[1,1:7]<-0
g<-rep(0:1,each=10); forward<-env$.populationMean(x,g,1L);reverse<-env$.populationMean(x,1L-g,1L)
stopifnot(forward$control$prior_null==0,reverse$control$prior_null==0)
saveRDS(list(forward=forward$control$positive_weight_reference,reverse=reverse$control$positive_weight_reference),
    file.path(root,"software/SGCB/tests/testthat/fixtures/mean_boundary_directions.rds"))
