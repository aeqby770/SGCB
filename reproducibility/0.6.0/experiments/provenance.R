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
agree <- identical(a,b);stopifnot(agree,setequal(compare[identical_source==FALSE,file],c("R/differential.R","R/population.R")))
get_text<-function(path)readLines(path,warn=FALSE,encoding="UTF-8")
outside_score<-function(x){
    first<-which(trimws(x)=="posterior_null <- joint_fit$posterior_null")
    last<-which(trimws(x)=="moments <- summary")
    stopifnot(length(first)==1L,length(last)==1L,first<last)
    x[-seq.int(first+2L,last-1L)]
}
stopifnot(identical(outside_score(get_text(file.path(old,"R/population.R"))),
    outside_score(get_text(file.path(nw,"R/population.R")))))
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
fwrite(v,file.path(root,"results","analysis_versions.tsv"),sep="\t")
cat("Shared full-analysis expressions",agree,"; unchanged source files",sum(compare$identical_source),"of",nrow(compare),"\n")