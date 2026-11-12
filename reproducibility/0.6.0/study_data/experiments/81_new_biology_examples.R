ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
LABEL <- ARGS[1]
STATE_SUFFIX <- c(ARGS[-1L],"_refined_continued")[1L]
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.6.0")
CURRENT_PREFIX <- file.path(DATA_ROOT,"results",paste0("current_",LABEL))
analysis <- readRDS(paste0(CURRENT_PREFIX,"_analysis_data.rds"))
group <- as.integer(analysis$group)-1L
PREFIX <- file.path(ROOT,"results",paste0("sgcb_",RECORD_ID,"_",LABEL))
EXAMPLES <- switch(LABEL,long_covid_2025=c("MTIF3","OMA1","C1QBP","COX20"),
    masld_liver_2023=c("SIAH2","ABCG8","EPHX1"),alzheimers_blood_2026=c("NSMAF","TNRC18"))
THREADS <- as.integer(Sys.getenv("SGCB_THREADS"))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)

results <- fread(paste0(PREFIX,"_biological_results.csv"))
metadata <- fread(paste0(CURRENT_PREFIX,"_covariate_samples.csv"))
moment <- fread(paste0(CURRENT_PREFIX,"_covariate_genes.csv"))
example_indices <- unique(c(match(EXAMPLES,results$gene_name),
    which(results$population_without_mean_DE & results$population_adjusted_retained)))
examples <- results[example_indices]
stopifnot(identical(examples$gene_name[seq_along(EXAMPLES)],EXAMPLES))
observation <- as.data.table(analysis$normalized_abundance[examples$gene_id,,drop=FALSE],
    keep.rownames="gene_id")
observation <- melt(observation,id.vars="gene_id",variable.name="sample",value.name="abundance")
observation <- merge(observation,examples[,.(gene_id,gene_name)],by="gene_id",sort=FALSE)
observation <- merge(observation,metadata,by="sample",sort=FALSE)
observation[,log2_abundance:=log2(abundance+.5)]
empirical <- observation[,.(n=.N,nonzero=sum(abundance>0),mean=mean(abundance),
    median=median(abundance),variance=mean((abundance-mean(abundance))^2),
    CV2=mean((abundance-mean(abundance))^2)/mean(abundance)^2,
    maximum_variance_fraction=max((abundance-mean(abundance))^2)/sum((abundance-mean(abundance))^2)),
    by=.(gene_id,gene_name,group)]
gg <- examples[,.(gene_id,gene_name,model="SGCB",effect=population_effect,
    se=population_se,df=NA_real_,pvalue=population_pvalue,qvalue=population_BH,
    lower=population_lower,upper=population_upper,
    interval_type="primary_95_percent_reference_interval")]
conditional <- merge(moment[gene_id %in% examples$gene_id],examples[,.(gene_id,gene_name)],by="gene_id")
conditional <- conditional[,.(gene_id,gene_name,model,effect=log2_cv2_ratio,se=se_log2_cv2,
    df=uniqueN(metadata$donor)-1,pvalue,qvalue)]
conditional[,c("lower","upper","interval_type"):=list(
    effect-qt(.975,df)*se,effect+qt(.975,df)*se,"frequentist_95_percent_CI")]
effects <- rbindlist(list(gg,conditional))
fwrite(examples,paste0(PREFIX,"_example_summary.csv"))
fwrite(observation,paste0(PREFIX,"_example_observations.csv"))
fwrite(empirical,paste0(PREFIX,"_example_empirical.csv"))
fwrite(effects,paste0(PREFIX,"_example_effects.csv"))
print(examples[,.(gene_name,log2FoldChange,padj,limma_q,edgeR_q,population_effect,population_BH,
    population_available,ctrl_nonzero_fraction,treat_nonzero_fraction,adjusted_effect,adjusted_q)])
