# Run from the extracted research directory with Rscript --vanilla experiments/reproduce.R.
file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
root<-normalizePath(file.path(dirname(sub("^--file=","",file_arg)),".."),winslash="/",mustWork=TRUE)
setwd(root);Sys.setenv(SGCB_R2_ROOT=root,SGCB_THREADS=3,OMP_NUM_THREADS=3,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
dir.create(file.path(root,"logs"),showWarnings=FALSE)
dir.create(file.path(root,"results"),showWarnings=FALSE)
workspace_tmp<-file.path(root,"tmp");dir.create(workspace_tmp,showWarnings=FALSE)
# Keep the system R temporary directory; workspace paths may contain spaces.
rscript<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
run<-function(script,args=character()) {
    log<-file.path(root,"logs",paste0("reproduce_",tools::file_path_sans_ext(basename(script)),".log"))
    message(basename(script)," ",if(length(args))paste(head(args,2),collapse=" ") else "")
    status<-system2(rscript,c("--vanilla","--encoding=UTF-8",shQuote(file.path(root,script)),args),stdout=log,stderr=log)
    if(status!=0L)stop("Stage failed. See ",log)
}
e<-function(name,args=character())run(file.path("experiments",name),args)
# The saved original-scoring records are comparison data; all fits use the current package.
e("prepare_inputs.R");e("replicate.R",c("final",1:432))
e("matched_reference.R",c("final","16","channel_n10_10","channel_n10_50","channel_n50_50","zero_n10_10","zero_n10_50","zero_n50_50"))
e("outer_refit.R",c(2:13,258:263));e("depth_and_outer.R","depth");e("depth_and_outer.R","outer")
e("exposure_sensitivity.R",1:24);e("known_exposure_ql.R");e("summarize_exposure.R")
e("null_boundary.R",c("technical",1:24));e("null_boundary.R",c("positive",1:16));e("null_boundary.R",c("fresh_depth",1:32))
e("boundary_summary.R");e("retained_evidence.R");e("mean_corroboration.R");e("archived_mean_boundary.R")
e("summarize_revision.R");e("final_agreement.R");e("provenance.R");e("package_tests.R");e("endpoint_validation.R")
run("visualization/figure_S9.R")
message("Numerical evidence is in ",file.path(root,"results"),". Follow REPRODUCE.md for typesetting and author review.")
