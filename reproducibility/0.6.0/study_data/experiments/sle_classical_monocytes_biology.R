ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
LABEL <- "sle_classical_monocytes_2021"
TYPE <- "counts"
INPUT <- file.path(DATA_ROOT,"data","new_biology",LABEL,"input.rds")
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.5.9")
PREFIX <- file.path(ROOT,"results",paste0("sgcb_",RECORD_ID,"_",LABEL))
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
ALPHA <- .05
MIN_SET <- 10L
MAX_SET <- 500L
KEYTYPE <- "SYMBOL"
EXAMPLES <- c("TNFSF13B","IL1RN","IL1R2","MTCO1P40","HBA1","HBA2","HBB")
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")
Rcpp::sourceCpp(file.path(ROOT,"development","57_cv_moment_regression.cpp"),
    rebuild=FALSE,cacheDir=file.path(ROOT,"tmp","cv_moment_cpp"))

ds <- readRDS(INPUT)
primary <- readRDS(paste0(PREFIX,".rds"))
stopifnot(
    identical(attr(primary,"native_source_id"),SOURCE_ID))
input <- switch(TYPE,counts=as.matrix(ds$counts),normalized=as.matrix(ds$abundance))
input <- input[primary$gene_id,,drop=FALSE]
group <- factor(ds$group)
sf <- attr(primary,"size_factors")
x <- sweep(input,2,sf,"/")
logx <- log2(x+.5)
annotation <- as.data.table(ds$annotation)[match(primary$gene_id,gene_id)]
meta <- as.data.table(ds$metadata)
stopifnot(identical(meta$sample,colnames(x)),!anyDuplicated(meta$donor),
    identical(annotation$gene_id,primary$gene_id))
meta[,group:=as.integer(group)-1L]

stopifnot(identical(normalizePath(attr(primary,"input_file"),winslash="/"),INPUT),
    identical(colnames(ds$counts),attr(primary,"samples")$sample),
    identical(as.character(group),attr(primary,"samples")$group),
    identical(levels(group),c("control","sle")),identical(ds$covariates,c("age","sex")))
# Sorted cMo: age/sex sensitivity without bulk blood-composition proxies.
complete <- complete.cases(meta[,ds$covariates,with=FALSE])
m <- meta[complete]
clinical <- scale(model.matrix(reformulate(ds$covariates),data=m)[,-1,drop=FALSE])
baseline_design <- cbind("(Intercept)"=1,group=m$group)
clinical_design <- cbind(baseline_design,clinical)
stopifnot(qr(clinical_design)$rank==ncol(clinical_design))
cluster <- as.integer(factor(m$donor))-1L
moment0 <- as.data.table(cv_moment_regression_cpp(x[,complete,drop=FALSE],
    baseline_design,cluster,THREADS))
moment1 <- as.data.table(cv_moment_regression_cpp(x[,complete,drop=FALSE],
    clinical_design,cluster,THREADS))
moment <- rbindlist(list(moment0[,model:="unadjusted"],moment1[,model:="clinical"]))
moment[,gene_id:=rep(primary$gene_id,2L)]
moment[,qvalue:=p.adjust(pvalue,"BH",n=.N),by=model]
moment <- merge(moment,as.data.table(primary)[,.(gene_id,SGCB_DV=DV,
    SGCB_effect=dv_log2_cv2_ratio,SGCB_q=dv_padj)],by="gene_id",sort=FALSE)
moment_summary <- moment[,.(n_genes=.N,n_samples=nrow(m),n_group0=sum(m$group==0),
    n_group1=sum(m$group==1),available=sum(is.finite(pvalue)),
    calls=sum(qvalue<ALPHA,na.rm=TRUE),SGCB_selected=sum(SGCB_DV),
    selected_available=sum(SGCB_DV & is.finite(pvalue)),
    selected_sign_retained=sum(SGCB_DV & sign(log2_cv2_ratio)==sign(SGCB_effect),na.rm=TRUE),
    selected_same_direction_significant=sum(SGCB_DV & qvalue<ALPHA &
        sign(log2_cv2_ratio)==sign(SGCB_effect),na.rm=TRUE)),by=model]

# The mean-effect comparison uses SGCB.
results <- merge(as.data.table(primary),annotation,by="gene_id",sort=FALSE)
results[,c("population_effect","population_se","population_lower","population_upper","population_pvalue","population_BH","population_DV","population_available") :=
    list(dv_log2_cv2_ratio,dv_se_log2_cv2,dv_ci_lower,dv_ci_upper,dv_pvalue,dv_padj,DV,dv_inference_available)]
adjusted <- moment[model=="clinical",.(gene_id,adjusted_effect=log2_cv2_ratio,adjusted_q=qvalue)]
results <- merge(results,adjusted,by="gene_id",sort=FALSE)
results[,population_without_mean_DE:=population_DV & !DE]
results[,population_adjusted_retained:=population_DV & is.finite(adjusted_q) & adjusted_q<ALPHA &
    sign(population_effect)==sign(adjusted_effect)]
results[,population_nonmean_adjusted_retained:=population_without_mean_DE & population_adjusted_retained]
results[,retained_adjusted_DV:=population_adjusted_retained]
results[,c("mean_comparison","conventional_mean_available"):=list("SGCB_only",FALSE)]
results[,DE_direction:=fcase(DE & log2FoldChange>0,"increase",DE & log2FoldChange<0,"decrease",default="not_called")]
results[,DV_direction:=fcase(DV & population_effect>0,"increase",DV & population_effect<0,"decrease",default="not_called")]
channel_genes <- results[channel!="neither",.(gene_id,gene_name,ensembl,channel,DE_direction,DV_direction,
    log2FoldChange,padj,dv_log2_cv2_ratio,dv_padj,dv_cv2_ctrl,dv_cv2_treat,adjusted_effect,adjusted_q)]
channel_summary <- channel_genes[,.(genes=.N),by=.(channel,DE_direction,DV_direction)]
summary <- results[,.(dataset=LABEL,features=.N,n_genes=.N,n_control=sum(group==levels(group)[1]),
    n_case=sum(group==levels(group)[2]),SGCB_DE=sum(DE),SGCB_DV=sum(DV),population_DV=sum(DV),
    DE_only=sum(channel=="DE_only"),DV_only=sum(channel=="DV_only"),DE_and_DV=sum(channel=="DE_and_DV"),
    neither=sum(channel=="neither"),population_without_mean_DE=sum(population_without_mean_DE),
    nonmean_adjusted_retained=sum(population_nonmean_adjusted_retained),
    covariate_samples=nrow(m),covariate_missing=sum(!complete),mean_comparison="SGCB_only",
    conventional_mean_available=FALSE)]
fwrite(results,paste0(PREFIX,"_biological_results.csv"))
fwrite(channel_genes,paste0(PREFIX,"_channel_genes.csv"))
fwrite(channel_summary,paste0(PREFIX,"_channel_summary.csv"))
fwrite(summary,paste0(PREFIX,"_comparison_summary.csv"))
fwrite(moment,paste0(PREFIX,"_covariate_genes.csv"))
fwrite(moment_summary,paste0(PREFIX,"_covariate_summary.csv"))
fwrite(m,paste0(PREFIX,"_covariate_samples.csv"))
analysis <- list(normalized_abundance=x,group=group,annotation=annotation,complete=complete,
    baseline_design=baseline_design,clinical_design=clinical_design,samples=meta,input=INPUT)
saveRDS(analysis,paste0(PREFIX,"_analysis_data.rds"))
go <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","go_term_maps_human.rds"))
kegg <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","kegg_term_maps_hsa.rds"))
keys <- results[,.(gene_id,mapping_key=gene_id)]
mapping <- as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys=intersect(keys$mapping_key,AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype=KEYTYPE)),
    columns="ENTREZID",keytype=KEYTYPE))
setnames(mapping,KEYTYPE,"mapping_key")
mapping <- unique(merge(keys,mapping[!is.na(ENTREZID)],by="mapping_key"))
mapping[,feature_multiplicity:=length(unique(ENTREZID)),by=gene_id]
mapping[,entrez_multiplicity:=length(unique(gene_id)),by=ENTREZID]
mapping <- mapping[feature_multiplicity==1 & entrez_multiplicity==1]

bp <- as.data.table(AnnotationDbi::toTable(GO.db::GOBPANCESTOR))
mf <- as.data.table(AnnotationDbi::toTable(GO.db::GOMFANCESTOR))
cc <- as.data.table(AnnotationDbi::toTable(GO.db::GOCCANCESTOR))
setnames(bp,c("GO","ancestor"))
setnames(mf,c("GO","ancestor"))
setnames(cc,c("GO","ancestor"))
ancestors <- rbindlist(list(bp[,ONTOLOGY:="BP"],mf[,ONTOLOGY:="MF"],cc[,ONTOLOGY:="CC"]))
direct <- as.data.table(go$term2gene)[ENTREZID %in% mapping$ENTREZID]
inherited <- merge(direct,ancestors,by=c("GO","ONTOLOGY"),allow.cartesian=TRUE)
go_links <- unique(rbindlist(list(direct[,.(ENTREZID,term=GO,collection=paste0("GO:",ONTOLOGY))],
    inherited[!is.na(ancestor),.(ENTREZID,term=ancestor,collection=paste0("GO:",ONTOLOGY))])))
go_names <- as.data.table(AnnotationDbi::select(GO.db::GO.db,keys=unique(go_links$term),
    columns="TERM",keytype="GOID"))[,.(term=GOID,description=TERM)]
reactome_links <- as.data.table(AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID))
reactome_names <- as.data.table(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))
terms <- unique(rbindlist(list(go_links,
    as.data.table(kegg$term2gene)[,.(ENTREZID,term=pathway_id,collection="KEGG")],
    reactome_links[grepl("^R-HSA-",DB_ID),.(ENTREZID=gene_id,term=DB_ID,collection="Reactome")])))
term_names <- unique(rbindlist(list(go_names,
    as.data.table(kegg$term2name)[,.(term=pathway_id,description=pathway_name)],
    reactome_names[grepl("^R-HSA-",DB_ID),.(term=DB_ID,
        description=sub("^Homo sapiens: ","",path_name))])))
background <- rbindlist(list(
    results[is.finite(pvalue),.(gene_id,gene_name,category="DE_up",selected=DE & log2FoldChange>0)],
    results[is.finite(pvalue),.(gene_id,gene_name,category="DE_down",selected=DE & log2FoldChange<0)],
    results[is.finite(pvalue),.(gene_id,gene_name,category="DE_only",selected=channel=="DE_only")],
    results[is.finite(pvalue),.(gene_id,gene_name,category="DE_and_DV",selected=channel=="DE_and_DV")],
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_increase",
        selected=population_DV & population_effect>0)],
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_decrease",
        selected=population_DV & population_effect<0)],
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_without_mean_DE",
        selected=population_DV & !DE)],
    results[population_available==TRUE & is.finite(adjusted_q),
        .(gene_id,gene_name,category="DV_without_mean_DE_adjusted",
            selected=retained_adjusted_DV & !DE)]))
mapped <- merge(background,mapping[,.(gene_id,ENTREZID)],by="gene_id")
membership <- merge(mapped,terms,by="ENTREZID",allow.cartesian=TRUE)
sizes <- unique(membership[,.(category,collection,ENTREZID,selected)])[,.(N=.N,n=sum(selected)),
    by=.(category,collection)]
enrichment <- membership[,.(M=length(unique(ENTREZID)),k=sum(selected),
    genes=paste(sort(unique(gene_name[selected])),collapse="/")),by=.(category,collection,term)]
enrichment <- merge(enrichment[M>=MIN_SET & M<=MAX_SET],sizes,by=c("category","collection"))
enrichment[,`:=`(pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),
    fold_enrichment=fifelse(n>0,(k/n)/(M/N),NA_real_),gene_ratio=k/n,background_ratio=M/N)]
enrichment[,qvalue:=p.adjust(pvalue,"BH"),by=.(category,collection)]
enrichment <- merge(enrichment,term_names,by="term",all.x=TRUE)
setorder(enrichment,category,collection,qvalue,pvalue,term)
summary <- enrichment[,.(tested_terms=.N,significant_terms=sum(qvalue<ALPHA),
    background_genes=unique(N),selected_annotated=unique(n),minimum_q=min(qvalue)),
    by=.(category,collection)]
mapping_summary <- merge(background[,.(tested_genes=.N,selected_genes=sum(selected)),by=category],
    mapped[,.(mapped_background=.N,mapped_selected=sum(selected)),by=category],by="category")
fwrite(enrichment,paste0(PREFIX,"_enrichment.csv"))
fwrite(summary,paste0(PREFIX,"_enrichment_summary.csv"))
fwrite(mapped,paste0(PREFIX,"_enrichment_gene_sets.csv"))
fwrite(mapping_summary,paste0(PREFIX,"_enrichment_mapping.csv"))
fwrite(results,paste0(PREFIX,"_biological_results.csv"))
fwrite(enrichment[grepl("interferon",description,ignore.case=TRUE)],paste0(PREFIX,"_interferon_terms.csv"))
metadata <- meta
examples <- results[gene_id %in% EXAMPLES | population_DV]

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
print(channel_summary)
