ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
ARGS <- commandArgs(trailingOnly=TRUE)
LABEL <- ARGS[1]
STATE_SUFFIX <- c(ARGS[-1L],"_refined_continued")[1L]
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor"))))
CURRENT_PREFIX <- file.path(DATA_ROOT,"results",paste0("current_",LABEL))
analysis <- readRDS(paste0(CURRENT_PREFIX,"_analysis_data.rds"))
group <- as.integer(analysis$group)-1L
PREFIX <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_",LABEL))
KEYTYPE <- switch(LABEL,masld_liver_2023="SYMBOL","ENSEMBL")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS"))
MIN_SET <- 10L
MAX_SET <- 500L
ALPHA <- .05
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)

results <- fread(paste0(PREFIX,"_biological_results.csv"))
results[,retained_adjusted_DV:=population_adjusted_retained]
go <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","go_term_maps_human.rds"))
kegg <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","kegg_term_maps_hsa.rds"))
keys <- results[,.(gene_id,mapping_key=switch(KEYTYPE,ENSEMBL=ensembl,SYMBOL=gene_name))]
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
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_increase",
        selected=population_DV & population_effect>0)],
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_decrease",
        selected=population_DV & population_effect<0)],
    results[population_available==TRUE,.(gene_id,gene_name,category="DV_without_mean_DE",
        selected=population_DV & !DE & !any_conventional_DE)],
    results[population_available==TRUE & is.finite(adjusted_q),
        .(gene_id,gene_name,category="DV_without_mean_DE_adjusted",
            selected=retained_adjusted_DV & !DE & !any_conventional_DE)]))
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
print(mapping_summary)
print(enrichment[category %in% c("DV_without_mean_DE","DV_without_mean_DE_adjusted") &
    qvalue<ALPHA,head(.SD,3),by=.(category,collection)][,
    .(category,collection,description,k,n,qvalue)])
