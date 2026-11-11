root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
library(data.table); setDTthreads(1)
source(file.path(root,"experiments/revision_calibration.R"))
ref <- file.path(root,"reference_0.5.9"); result <- file.path(ref,"results"); out <- file.path(root,"results")
indices <- CJ(panel=c("channel","zero"),design=c("10_10","10_50","50_50"))
fixed <- function(panel,design) {
    id <- paste0(panel,"_n",design)
    raw <- readRDS(file.path(ref,"data",if (panel=="channel") "channel_simulation.rds" else "zero_calibration.rds"))
    old <- readRDS(file.path(result,paste0("sgcb_0.5.9_",id,".rds")))
    truth <- fread(file.path(result,paste0("sgcb_0.5.9_",id,"_evaluation_genes.csv")))[endpoint=="DE"]
    g<-as.integer(factor(raw$group))-1L; n<-as.integer(strsplit(design,"_")[[1]])
    ind<-c(head(which(g==0),n[1]),head(which(g==1),n[2]))
    additional <- moment_workflow(raw$abundance[,ind,drop=FALSE],g[ind])$mean_p
    secondary <- p.adjust(pmax(old$pvalue,ifelse(is.finite(additional),additional,1)),"BH")
    sg <- old$DE; guard <- secondary<.05; stopifnot(identical(old$gene_id,truth$gene_id))
    data.table(panel=panel,size=design,SGCB_TP=sum(sg & truth$truth),proposed_TP=sum(guard&truth$truth),
        SGCB_FP=sum(sg&!truth$truth),proposed_FP=sum(guard&!truth$truth),new_calls=sum(guard&!sg),lost_calls=sum(sg&!guard))
}
fwrite(rbindlist(Map(fixed,indices$panel,indices$design)),file.path(out,"mean_candidate_protected_fixed.tsv"),sep="\t")
z <- readRDS(file.path(result,"sgcb_0.5.7_seqc.rds"))
seqc <- readRDS(file.path(ref,"data","seqc_input.rds")); x <- as.matrix(seqc$counts)
stopifnot(all(z$gene_id%in%rownames(x))); x <- x[z$gene_id,,drop=FALSE]
sf <- attr(z,"size_factors");stopifnot(length(sf)==ncol(x))
m <- moment_workflow(sweep(x,2,sf,"/"),as.integer(factor(seqc$group))-1L)
pq <- pmax(z$pvalue,ifelse(is.finite(m$mean_p),m$mean_p,1)); sq <- p.adjust(pq,"BH",n=nrow(z))
raw <- fread(file.path(result,"sgcb_0.5.7_seqc_qpcr_gene_scores.csv"))[method=="SGCB"]
matches <- match(raw$gene,z$gene_id); stopifnot(nrow(raw)==683L,all(!is.na(matches)))
y <- abs(raw$log2FC_qpcr_B_over_A)>=1
score_auc <- function(p) {
    scores <- -log10(pmax(p,1e-300)); delta <- outer(scores[y],scores[!y],"-")
    mean((delta>0)+.5*(delta==0))
}
fwrite(data.table(method=c("SGCB","Corroborated"),n=length(y),positive=sum(y),
    auc=c(score_auc(z$pvalue[matches]),score_auc(pq[matches])),
    TP=c(sum(z$padj[matches]<.05 & y),sum(sq[matches]<.05 & y)),
    FP=c(sum(z$padj[matches]<.05 & !y),sum(sq[matches]<.05 & !y))),
    file.path(out,"mean_candidate_qpcr.tsv"),sep="\t")