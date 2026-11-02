ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
VERSION <- as.character(utils::packageVersion("SGCB",lib.loc=Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library"))))
PREFIX <- paste0("sgcb_",VERSION,"_kirc")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
MIN_SET <- 10L
MAX_SET <- 500L
ALPHA <- .05
Sys.setenv(OMP_NUM_THREADS=THREADS,OPENBLAS_NUM_THREADS=THREADS,MKL_NUM_THREADS=THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

primary <- readRDS(file.path(ROOT,"results",paste0(PREFIX,".rds")))
stopifnot(identical(c(attr(primary,"mean_inference_source_id"),
    attr(primary,"native_source_id"))[1],SOURCE_ID))
results <- as.data.table(primary)
go <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","go_term_maps_human.rds"))
kegg <- readRDS(file.path(DATA_ROOT,"data","biological_story_annotation_cache","kegg_term_maps_hsa.rds"))
mapping <- as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys=intersect(results$gene_id,AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype="SYMBOL")),
    columns="ENTREZID",keytype="SYMBOL"))
setnames(mapping,"SYMBOL","gene_id")
mapping <- unique(mapping[!is.na(ENTREZID)])
mapping[,symbol_multiplicity:=uniqueN(ENTREZID),by=gene_id]
mapping[,entrez_multiplicity:=uniqueN(gene_id),by=ENTREZID]
mapping <- mapping[symbol_multiplicity==1 & entrez_multiplicity==1]

# Propagate the cached direct GO annotations through the installed GO graph.
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
print(mapping_summary)
print(summary)
print(enrichment[qvalue<ALPHA,head(.SD,4),by=.(category,collection)][,
    .(category,collection,description,k,n,fold_enrichment,qvalue)])
