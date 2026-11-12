ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.6.0")
PREFIX <- paste0("sgcb_",RECORD_ID,"_kirc")
ANNOTATION_FILE <- file.path(DATA_ROOT,"data","biological_story_annotation_cache","kirc_enrichment_annotation.rds")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
MIN_SET <- 10L
MAX_SET <- 500L
ALPHA <- .05
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

primary <- readRDS(file.path(ROOT,"results",paste0(PREFIX,".rds")))
verify_study_record(file.path(ROOT,"results",paste0(PREFIX,".rds")),attr(primary,"native_source_id"),RECORD_ID)
results <- as.data.table(primary)
# Resolved study annotations preserve the original gene universe and ontology.
annotation <- readRDS(ANNOTATION_FILE)
stopifnot(identical(attr(primary,"input_id"),annotation$input_id),
    setequal(results$gene_id,annotation$gene_ids))
mapping <- as.data.table(annotation$mapping)
terms <- as.data.table(annotation$terms)
term_names <- as.data.table(annotation$term_names)

# Each background contains all genes eligible for the relevant primary test.
background <- rbindlist(list(
    results[is.finite(pvalue),.(gene_id,category="DE_up",selected=DE & log2FoldChange>0)],
    results[is.finite(pvalue),.(gene_id,category="DE_down",selected=DE & log2FoldChange<0)],
    results[is.finite(pvalue) & dv_inference_available,
        .(gene_id,category="DV_only",selected=channel=="DV_only")]))
mapped <- merge(background,mapping[,.(gene_id,ENTREZID)],by="gene_id")
membership <- merge(mapped,terms,by="ENTREZID",allow.cartesian=TRUE)
sizes <- unique(membership[,.(category,collection,ENTREZID,selected)])[,.(N=.N,n=sum(selected)),
    by=.(category,collection)]
enrichment <- membership[,.(M=uniqueN(ENTREZID),k=sum(selected),
    genes=paste(sort(unique(gene_id[selected])),collapse="/")),by=.(category,collection,term)]
enrichment <- merge(enrichment[M>=MIN_SET & M<=MAX_SET],sizes,by=c("category","collection"))
enrichment[,`:=`(pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),
    fold_enrichment=(k/n)/(M/N),gene_ratio=k/n,background_ratio=M/N)]
enrichment[,qvalue:=p.adjust(pvalue,"BH"),by=.(category,collection)]
enrichment <- merge(enrichment,term_names,by="term",all.x=TRUE)
setorder(enrichment,category,collection,qvalue,pvalue,term)
summary <- enrichment[,.(tested_terms=.N,significant_terms=sum(qvalue<ALPHA),
    background_genes=unique(N),selected_annotated=unique(n),minimum_q=min(qvalue)),
    by=.(category,collection)]
mapping_summary <- merge(background[,.(tested_genes=.N,selected_genes=sum(selected)),by=category],
    mapped[,.(mapped_background=.N,mapped_selected=sum(selected)),by=category],by="category")
fwrite(enrichment,file.path(ROOT,"results",paste0(PREFIX,"_enrichment.csv")))
fwrite(summary,file.path(ROOT,"results",paste0(PREFIX,"_enrichment_summary.csv")))
fwrite(mapped,file.path(ROOT,"results",paste0(PREFIX,"_enrichment_gene_sets.csv")))
fwrite(mapping_summary,file.path(ROOT,"results",paste0(PREFIX,"_enrichment_mapping.csv")))
annotation_versions <- as.data.table(annotation$versions)
fwrite(annotation_versions,file.path(ROOT,"results",paste0(PREFIX,"_annotation_versions.csv")))
print(mapping_summary)
print(summary)
print(enrichment[qvalue<ALPHA,head(.SD,4),by=.(category,collection)][,
    .(category,collection,description,k,n,fold_enrichment,qvalue)])
