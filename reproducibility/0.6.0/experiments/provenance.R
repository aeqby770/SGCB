root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
r1 <- file.path(root,"reference_0.5.9"); old <- file.path(r1,"software/SGCB"); nw <- file.path(root,"software/SGCB")
library(data.table)
numerical <- c(list.files(file.path(old,"src"),full.names=FALSE),list.files(file.path(old,"R"),full.names=FALSE))
paths <- c(file.path("src",list.files(file.path(old,"src"))),file.path("R",list.files(file.path(old,"R"))))
same <- function(f) identical(readLines(file.path(old,f),warn=FALSE,encoding="UTF-8"),readLines(file.path(nw,f),warn=FALSE,encoding="UTF-8"))
compare <- data.table(file=paths,identical_source=vapply(paths,same,logical(1)))
env_old <- new.env(); env_new <- new.env()
sys.source(file.path(old,"R/differential.R"),env_old);sys.source(file.path(nw,"R/differential.R"),env_new)
b0 <- as.list(body(env_old$sgcbDE)); b1 <- as.list(body(env_new$sgcbDE))
# In single-target dispatch only, remove endpoint selection and early returns. The unsupplied
# cache is an empty environment and may be allocated before descriptive work without changing the fit.
string <- function(x) paste(deparse(x,width.cutoff=500L),collapse=" ")
key <- function(x) vapply(x,string,character(1))
a <- key(b0); b <- key(b1)
a <- a[!startsWith(a,"kernel_cache <- new.env")]
b <- b[!startsWith(b,"kernel_cache <- new.env") & !startsWith(b,"endpoint <- match.arg(endpoint)") & !grepl('^if \\(endpoint == "(mean|variability)"\\)',b)]
agree <- identical(a,b);stopifnot(agree,setequal(compare[identical_source==FALSE,file],c("R/differential.R","R/population.R","src/population_gg.cpp")))
get_text<-function(path)readLines(path,warn=FALSE,encoding="UTF-8")
outside_score<-function(x){
    first<-which(trimws(x)=="posterior_null <- joint_fit$posterior_null")
    last<-which(trimws(x)=="moments <- summary")
    stopifnot(length(first)==1L,length(last)==1L,first<last)
    x[-seq.int(first+2L,last-1L)]
}
old_outside <- paste(outside_score(get_text(file.path(old,"R/population.R"))),collapse="\n")
new_outside <- paste(outside_score(get_text(file.path(nw,"R/population.R"))),collapse="\n")
scale_block <- paste(c(
    "    log_scale <- .5*log(cbind(matrixStats::rowVars(positive_log[,group==0,drop=FALSE],na.rm=TRUE),",
    "        matrixStats::rowVars(positive_log[,group==1,drop=FALSE],na.rm=TRUE)))"),collapse="\n")
new_head <- paste(c(scale_block,
    "    positive_usable <- positive_counts>=4L & is.finite(log_scale)",
    "    control$preserve_nested <- any(x==0)",
    "    # No two-group CV2 posterior is available without four variable positives per group.",
    "    switch(as.character(any(rowSums(positive_usable)==2L)),"),collapse="\n")
old_head <- paste(c("    control$preserve_nested <- any(x==0)",
    "    # No two-group CV2 posterior is available without four positives per group.",
    "    switch(as.character(any(rowSums(positive_counts>=4L)==2L)),"),collapse="\n")
new_range <- "    scale_range <- range(log_scale[positive_usable],finite=TRUE)"
old_range <- paste(scale_block,"    scale_range <- range(log_scale[positive_counts>=4],finite=TRUE)",sep="\n")
stopifnot(sum(gregexpr(new_head,new_outside,fixed=TRUE)[[1]]>0L)==1L,
    sum(gregexpr(new_range,new_outside,fixed=TRUE)[[1]]>0L)==1L)
new_outside <- sub(new_head,old_head,new_outside,fixed=TRUE)
new_outside <- sub(new_range,old_range,new_outside,fixed=TRUE)
stopifnot(identical(old_outside,new_outside))
fwrite(compare,file.path(root,"results","source_dependency_comparison.tsv"),sep="\t")
v <- data.table(
    analysis=c("Fixed positive and observed-zero panels; Fig. 2 and 3","qPCR Fig. 4 and the six negative-binomial mean benchmarks, protein spike-ins and Kang mean",
               "Kernel restrictions","GTEx liver partitions","KIRC, SLE and liver primary Fig. 5 and 6","Covariate-adjusted moments and held-out/external survival",
               "Repeated null, mixed and technical panels; Fig. S9","Welch and empirical log-CV-squared comparator","Paired qPCR AUC differences",
               "Background resampling and interval-scale sensitivity","Exposure-alignment sensitivity","Mean corroboration exploratory comparison",
               "Endpoint timing and package tests","Fig. 1"),
    fitted_by=c("0.5.9 + 0.6.0","0.5.7","0.5.9","0.5.9","0.5.9","R1 downstream models","0.6.0","two-moment workflow","0.5.7 SGCB + native external scores","0.6.0","0.6.0 mean endpoint","0.6.0 + Welch","0.6.0","schematic"),
    disposition=c("refit with R1 thread and input settings","retained SGCB arithmetic-mean estimates; competitor tables retained under original names",
       "matched GG/Gamma/lognormal results reused","existing overlapping partitions, descriptive","R1 primary fits retained; both-channel numerical path compared",
       "separate same-cohort sensitivity and survival models; retained","fresh generation and independent population refit for each matrix",
       "same newly generated inputs; two endpoints BH-adjusted separately","paired differences newly calculated; existing predictions reused",
       "fixed probe measurements; background draws refitted","same raw counts and fixed known exposures; separately identified simulated alternatives",
       "raw-score corroboration vs unchanged full endpoint","fresh execution","conceptual illustration of the same primary targets (separate submission attachment)"),
    source=c("results/protected/*_final_t16.rds","reference_0.5.9/results/sgcb_0.5.7_*","reference_0.5.9/results/cv_matched_families_0.5.9_*",
       "reference_0.5.9/results/sgcb_0.5.9_readcount_null_*","reference_0.5.9/results/sgcb_0.5.9_{kirc,sle_classical_monocytes_2021,masld_liver_2023}*",
       "reference_0.5.9/results/{*covariate*,*survival*}","experiments/design.csv; results/final/; results/mean_screen_thresholds.tsv","experiments/revision_calibration.R",
       "results/paired_qpcr_auc.tsv","results/outer/","results/exposure_final/","results/mean_candidate_holdout.tsv","results/endpoint_runtime_equivalence.tsv; results/package_tests.tsv","Figure_1.pdf"))
v<-rbind(v,data.table(analysis=c("Zero-null scoring and global-shift counterexamples","Known exposure count-likelihood reference","Retained R1 mean score boundary identity"),
    fitted_by=c("0.5.9 + 0.6.0","edgeR 4.8.2, robust QL with known exposure","original 0.5.7 and 0.5.9 objects"),
    disposition=c("fixed fitted effects and nuisance distribution; analytic boundary limit checked against installed public entry point",
        "new same-input count comparator; the method label remains edgeR","31 original fits checked against exact boundary eligibility"),
    source=c("results/boundary_operating_characteristics.tsv","results/known_exposure_ql_summary.tsv","results/archived_mean_boundary.tsv")))
switch(as.character(file.exists(file.path(root,"results/eligibility_refit_comparison.tsv"))),"TRUE"={
    v[analysis=="Repeated null, mixed and technical panels; Fig. S9",disposition:=
        "R2 records retained except six original depth matrices refitted after constant-positive eligibility correction; original random streams and three threads reused"]
    v[analysis=="Endpoint timing and package tests",disposition:=
        "original R2 benchmark runtime on 5700X3D retained; 15 constant-input boundary assertions passed after correction; six refit runtimes on 5900HX retained separately"]
    v[analysis=="Endpoint timing and package tests",source:=
        "results/endpoint_runtime_equivalence.tsv; results/package_tests.tsv; software/SGCB/tests/testthat/test-eligibility.R; results/eligibility_refit_comparison.tsv"]
})
fwrite(v,file.path(root,"results","analysis_versions.tsv"),sep="\t")
cat("Shared full-analysis expressions",agree,"; unchanged source files",sum(compare$identical_source),"of",nrow(compare),"\n")
