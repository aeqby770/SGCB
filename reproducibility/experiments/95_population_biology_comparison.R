ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
LABEL <- ARGS[1]
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor"))))
TYPE <- switch(LABEL,long_covid_2025="counts",masld_liver_2023="normalized",
    alzheimers_blood_2026="normalized")
OUTPUT <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_",LABEL))
ALPHA <- .05
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")

stored <- readRDS(file.path(DATA_ROOT,"results",paste0("current_",LABEL,"_analysis_data.rds")))
primary <- readRDS(paste0(OUTPUT,".rds"))
current <- as.data.table(primary)
samples <- attr(primary,"samples")
stopifnot(setequal(rownames(stored$normalized_abundance),current$gene_id),
    identical(colnames(stored$normalized_abundance),samples$sample),
    identical(as.character(stored$group),samples$group))
group <- as.integer(stored$group)-1L
N0 <- sum(group==0)
N1 <- sum(group==1)
original <- fread(file.path(DATA_ROOT,"results",paste0("current_",LABEL,"_biological_results.csv")))
baseline <- fread(file.path(DATA_ROOT,"results",paste0("current_",LABEL,"_genes.csv")))
stopifnot(!anyDuplicated(current$gene_id),!anyDuplicated(baseline$gene_id),
    setequal(current$gene_id,baseline$gene_id),setequal(original$gene_id,current$gene_id))
baseline <- baseline[match(current$gene_id,gene_id)]
current[,c("gene_wise_DV","gene_wise_effect","gene_wise_available"):=list(
    baseline$DV,baseline$dv_log2_cv2_ratio,baseline$dv_inference_available)]

original <- merge(current,original[,c("gene_id",setdiff(names(original),names(current))),with=FALSE],
    by="gene_id",sort=FALSE)
candidate <- current[,.(gene_id,population_effect=dv_log2_cv2_ratio,
    population_se=dv_se_log2_cv2,population_lower=dv_ci_lower,population_upper=dv_ci_upper,
    population_pvalue=dv_pvalue,
    population_BH=dv_padj,population_DV=DV,population_available=dv_inference_available)]
stopifnot(!anyDuplicated(original$gene_id),!anyDuplicated(candidate$gene_id),
    setequal(original$gene_id,candidate$gene_id),nrow(original)==nrow(stored$normalized_abundance))
genes <- merge(original,candidate,by="gene_id",sort=FALSE)
genes[,`:=`(population_without_mean_DE=population_DV & !DE & !any_conventional_DE,
    population_adjusted_retained=population_DV & adjusted_q<ALPHA &
        sign(population_effect)==sign(adjusted_effect))]
genes[,population_nonmean_adjusted_retained:=population_without_mean_DE &
    population_adjusted_retained]
genes[,shared_effect_available:=population_available & gene_wise_available &
    is.finite(population_effect) & is.finite(gene_wise_effect)]

# The full-family mean and covariate fits use the same abundance-CV2 target
# and input, so their saved estimates can be compared without fitting again.
moments <- fread(file.path(DATA_ROOT,"results",paste0("current_",LABEL,"_covariate_genes.csv")))
moments <- merge(moments,candidate,by="gene_id",sort=FALSE)
covariates <- moments[,.(features=.N,available=sum(is.finite(pvalue)),
    moment_calls=sum(qvalue<ALPHA,na.rm=TRUE),population_selected=sum(population_DV),
    selected_available=sum(population_DV & is.finite(pvalue)),
    selected_direction_retained=sum(population_DV &
        sign(log2_cv2_ratio)==sign(population_effect),na.rm=TRUE),
    selected_same_direction_significant=sum(population_DV & qvalue<ALPHA &
        sign(log2_cv2_ratio)==sign(population_effect),na.rm=TRUE)),by=model]
covariates[,dataset:=LABEL]
summary <- genes[,.(dataset=LABEL,native_source_id=attr(primary,"native_source_id"),
    mean_inference_source_id=SOURCE_ID,n_control=N0,n_case=N1,features=.N,
    available=sum(population_available),original_DV=sum(gene_wise_DV),
    population_DV=sum(population_DV),shared_DV=sum(gene_wise_DV & population_DV),
    population_only=sum(population_DV & !gene_wise_DV),original_only=sum(gene_wise_DV & !population_DV),
    population_without_mean_DE=sum(population_without_mean_DE),
    population_higher=sum(population_DV & population_effect>0),
    population_lower=sum(population_DV & population_effect<0),
    comparable_effects=sum(shared_effect_available),
    effect_correlation_spearman=cor(population_effect[shared_effect_available],
        gene_wise_effect[shared_effect_available],method="spearman"),
    effect_direction_agreement=mean(sign(population_effect[shared_effect_available])==
        sign(gene_wise_effect[shared_effect_available])),
    adjusted_retained=sum(population_adjusted_retained,na.rm=TRUE),
    nonmean_adjusted_retained=sum(population_nonmean_adjusted_retained,na.rm=TRUE))]
fwrite(genes,paste0(OUTPUT,"_biological_results.csv"))
fwrite(covariates,paste0(OUTPUT,"_covariate_summary.csv"))
fwrite(summary,paste0(OUTPUT,"_comparison_summary.csv"))
print(summary)
print(covariates)
