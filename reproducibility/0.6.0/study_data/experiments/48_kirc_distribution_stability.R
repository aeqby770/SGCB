ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.6.0")
PREFIX <- paste0("sgcb_",RECORD_ID,"_kirc")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
.libPaths(c(Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library")),.libPaths()))
library(SGCB)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

res <- readRDS(file.path(ROOT,"results",paste0(PREFIX,".rds")))
verify_study_record(file.path(ROOT,"results",paste0(PREFIX,".rds")),attr(res,"native_source_id"),RECORD_ID)
# Reuse the covariate estimates only for the same input and normalization.
reference <- readRDS(file.path(DATA_ROOT,"results","population_kirc.rds"))
covariate_design <- readRDS(file.path(DATA_ROOT,"results","covariate_CV2_kirc_discovery_design.rds"))
stopifnot(identical(attr(res,"input_id"),attr(reference,"input_id")),
    identical(attr(res,"samples"),attr(reference,"samples")),
    identical(unname(attr(res,"size_factors")),unname(attr(reference,"size_factors"))),
    identical(covariate_design$complete_samples,attr(res,"samples")$sample),
    identical(unname(covariate_design$normalization_factors),unname(attr(res,"size_factors"))))
covariates <- fread(file.path(DATA_ROOT,"results","population_kirc_covariate_genes.csv"))
stopifnot(!anyDuplicated(covariates$gene_id),setequal(covariates$gene_id,res$gene_id))
covariates <- covariates[match(res$gene_id,gene_id)]
covariates[,c("primary_GG_effect","primary_GG_q","channel","population_available","population_DV"):=
    list(res$dv_log2_cv2_ratio,res$dv_padj,res$channel,res$dv_inference_available,res$DV)]
fwrite(covariates,file.path(ROOT,"results",paste0(PREFIX,"_covariate_genes.csv")))
covariate_summary <- fread(file.path(DATA_ROOT,"results","population_kirc_covariate_summary.csv"))
selected <- covariates[population_DV & is.finite(primary_GG_effect) &
    is.finite(adjusted_log2_cv2_ratio) & is.finite(unadjusted_log2_cv2_ratio)]
covariate_summary[,c("population_calls","selected_available","selected_sign_retained",
    "selected_adjusted_significant","selected_median_effect_change","selection_method"):=list(
    sum(res$DV),nrow(selected),
    mean(sign(selected$adjusted_log2_cv2_ratio)==sign(selected$primary_GG_effect)),
    sum(selected$adjusted_qvalue<.05,na.rm=TRUE),
    median(abs(selected$adjusted_log2_cv2_ratio-selected$unadjusted_log2_cv2_ratio)),
    paste0("SGCB ",RECORD_ID,"; ",c(attr(res,"runtime_scope"),"full public API")[1]))]
fwrite(covariate_summary,file.path(ROOT,"results",paste0(PREFIX,"_covariate_summary.csv")))
input <- readRDS(file.path(DATA_ROOT,"data","kirc_independent_discovery.rds"))
ds <- readRDS(file.path(DATA_ROOT,"data","tcga_kirc_cache","tcga_kirc.rds"))
membership <- fread(file.path(DATA_ROOT,"results","kirc_discovery_sample_membership.csv"))
genes <- intersect(res$gene_id,rownames(ds$expression))
discovery_index <- match(substr(colnames(input$counts),1,15),substr(colnames(ds$expression),1,15))
stopifnot(!anyNA(discovery_index),!anyDuplicated(discovery_index))
discovery_samples <- colnames(ds$expression)[discovery_index]
heldout_samples <- colnames(ds$expression)[substr(colnames(ds$expression),14,15)=="01" &
    !substr(colnames(ds$expression),1,12) %in% membership$patient]
samples <- c(discovery_samples,heldout_samples)
# The full TCGA cache uses a different gene annotation release from discovery.
# Compare matched genes and patients within this single cache, using CPM on all rows.
x <- sweep(ds$expression[,samples,drop=FALSE],2,
           colSums(ds$expression[,samples,drop=FALSE]),"/")*1e6
x <- x[genes,,drop=FALSE]
discovery_groups <- as.character(input$group)
normal <- which(c(discovery_groups,rep("heldout",length(heldout_samples)))=="normal")
tumor <- which(c(discovery_groups,rep("heldout",length(heldout_samples)))=="tumor")
heldout <- match(heldout_samples,samples)
cv0 <- matrixStats::rowVars(x[,normal,drop=FALSE])/rowMeans(x[,normal,drop=FALSE])^2
cv1 <- matrixStats::rowVars(x[,tumor,drop=FALSE])/rowMeans(x[,tumor,drop=FALSE])^2
cvh <- matrixStats::rowVars(x[,heldout,drop=FALSE])/rowMeans(x[,heldout,drop=FALSE])^2
out <- as.data.table(res)[match(genes,gene_id)]
out[,`:=`(empirical_cv2_normal=cv0,empirical_cv2_discovery=cv1,empirical_cv2_heldout=cvh,
    empirical_log2_cv2_ratio_discovery=log2(cv1/cv0),
    empirical_log2_cv2_ratio_heldout=log2(cvh/cv0))]
out[,effect_available:=is.finite(dv_log2_cv2_ratio) & is.finite(empirical_log2_cv2_ratio_heldout)]
out[,direction_agreement:=sign(dv_log2_cv2_ratio)==sign(empirical_log2_cv2_ratio_heldout)]
summary <- out[,.(n_genes=.N,n_effects=sum(effect_available),
    rho_tumor_cv2=cor(log2(empirical_cv2_discovery),log2(empirical_cv2_heldout),
                     method="spearman",use="pairwise.complete.obs"),
    rho_effect=cor(dv_log2_cv2_ratio[effect_available],
                   empirical_log2_cv2_ratio_heldout[effect_available],method="spearman"),
    rho_empirical_effect=cor(empirical_log2_cv2_ratio_discovery,
                            empirical_log2_cv2_ratio_heldout,method="spearman",use="pairwise.complete.obs"),
    direction_agreement=mean(direction_agreement[effect_available])),by=channel]
summary[,`:=`(normal_samples=length(normal),discovery_tumors=length(tumor),
              heldout_tumors=length(heldout),shared_normal_reference=TRUE)]
fwrite(out,file.path(ROOT,"results",paste0(PREFIX,"_distribution_stability.csv")))
fwrite(summary,file.path(ROOT,"results",paste0(PREFIX,"_distribution_stability_summary.csv")))
coverage <- data.table(original_discovery_genes=nrow(res),matched_genes=length(genes),
    original_discovery_normal=sum(input$group=="normal"),original_discovery_tumor=sum(input$group=="tumor"),
    matched_discovery_normal=length(normal),matched_discovery_tumor=length(tumor),
    heldout_tumors=length(heldout),different_aliquot_barcode=sum(discovery_samples!=colnames(input$counts)),
    normal_reference="shared",normalization="CPM on full cache gene set")
fwrite(coverage,file.path(ROOT,"results",paste0(PREFIX,"_stability_coverage.csv")))
fwrite(data.table(discovery_barcode=colnames(input$counts),cache_barcode=discovery_samples,
    patient=substr(discovery_samples,1,12),group=discovery_groups),
    file.path(ROOT,"results",paste0(PREFIX,"_stability_sample_mapping.csv")))

# Retain the three examples selected for the original biological analysis.
# Their population-model channels are reported even if classification changes.
examples <- out[gene_id %in% c("SF3B1","FTL","GPX3")]
values <- as.data.table(t(x[examples$gene_id,,drop=FALSE]),keep.rownames="sample")
values <- melt(values,id.vars="sample",variable.name="gene_id",value.name="abundance")
sample_groups <- data.table(sample=samples,
    group=c(discovery_groups,rep("heldout_tumor",length(heldout_samples))))
values <- merge(values,sample_groups,by="sample")
values <- merge(values,examples[,.(gene_id,channel,padj,dv_padj,dv_log2_cv2_ratio)],by="gene_id")
fwrite(values,file.path(ROOT,"results",paste0(PREFIX,"_distribution_examples.csv")))
fwrite(examples,file.path(ROOT,"results",paste0(PREFIX,"_example_features.csv")))
print(summary)
print(examples[,.(gene_id,channel,baseMean,log2FoldChange,dv_log2_cv2_ratio,
                  empirical_log2_cv2_ratio_heldout,padj,dv_padj)])
