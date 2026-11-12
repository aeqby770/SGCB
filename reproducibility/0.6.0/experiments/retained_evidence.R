root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE)
library(data.table); r <- file.path(root,"study_data","results"); out <- file.path(root,"results")
scores <- fread(file.path(r,"sgcb_0.5.7_seqc_qpcr_gene_scores.csv"))
methods <- c("SGCB","DESeq2","edgeR")
wide <- dcast(scores[method%in%methods],gene + log2FC_qpcr_B_over_A+log2FC_qpcr_arithmetic_B_over_A ~ method,value.var="pvalue")
stopifnot(nrow(wide)==683L,all(is.finite(as.matrix(wide[,..methods]))))
auc_rows <- function(target,cutoff) {
    observed <- abs(wide[[target]])>=cutoff
    scores <- -log10(pmax(as.matrix(wide[,..methods]),1e-300))
    place <- lapply(seq_len(ncol(scores)),function(m) {
        pos <- scores[observed,m]; neg <- scores[!observed,m]
        delta <- outer(pos,neg,"-"); ties <- (delta>0)+.5*(delta==0)
        list(pos=rowMeans(ties),neg=colMeans(ties))
    })
    rbindlist(lapply(2:3,function(m) {
        a <- place[[1]]; b <- place[[m]]; d <- mean(a$pos)-mean(b$pos)
        se <- sqrt(var(a$pos-b$pos)/sum(observed)+var(a$neg-b$neg)/sum(!observed))
        data.table(reference=target,threshold=cutoff,comparison=methods[m],genes=nrow(wide),
            positive=sum(observed),negative=sum(!observed),sgcb_auc=mean(a$pos),comparator_auc=mean(b$pos),difference=d,
            paired_se=se,lower=d-qnorm(.975)*se,upper=d+qnorm(.975)*se,p=2*pnorm(-abs(d/se)))
    }))
}
contrasts <- CJ(target=c("log2FC_qpcr_B_over_A","log2FC_qpcr_arithmetic_B_over_A"),cutoff=c(.5,1,1.5))
aucs <- rbindlist(Map(auc_rows,contrasts$target,contrasts$cutoff))
aucs[,p_holm_two:=p.adjust(p,"holm"),by=.(reference,threshold)]
fwrite(aucs,file.path(out,"paired_qpcr_auc.tsv"),sep="\t")
params <- fread(file.path(r,"zero_calibration_parameters.csv"))
pop <- unique(params[,.(scenario,group,pi)])
p <- dcast(pop,scenario~group,value.var="pi"); setnames(p,c("0","1"),c("p0","p1"))
files <- list.files(r,pattern="^sgcb_0[.]5[.]9_zero_n.*_evaluation_genes[.]csv$",full.names=TRUE)
old <- rbindlist(lapply(files,fread)); old <- old[endpoint=="DV" & truth_finite==TRUE]
old[,ci:=is.finite(lower)&is.finite(upper)]
g <- old[,.(N=.N,alternatives=sum(as.logical(truth)),available=sum(ci),calls=sum(called==TRUE),
    covered=sum(ci & lower<=truth_effect & upper>=truth_effect),covered_conditional=mean((lower<=truth_effect & upper>=truth_effect)[ci]),
    power=ifelse(any(as.logical(truth)),sum(called==TRUE & as.logical(truth))/sum(as.logical(truth)),NA_real_),
    raw_rejection_all=sum(pvalue<.05,na.rm=TRUE)/.N,median_effect=median(effect,na.rm=TRUE)),by=.(n0,n1,scenario)]
g <- merge(g,p,by="scenario"); g[,probability_at_least_four_each:=pbinom(3,n0,p0,lower.tail=FALSE)*pbinom(3,n1,p1,lower.tail=FALSE)]
fwrite(g,file.path(out,"fixed_sparse_guidance.tsv"),sep="\t")
time_paths <- list.files(r,pattern="^sgcb_0[.]5[.]9_.*process_runtime[.]csv$",full.names=TRUE)
resources <- rbindlist(lapply(time_paths,fread),fill=TRUE)
resources[,peak_GiB:=peak_working_set_bytes/2^30];fwrite(resources,file.path(out,"retained_resources.tsv"),sep="\t")
hrs <- fread(file.path(r,"sgcb_0.5.9_kirc_heldout_survival.csv"))[term=="DV_only"]
ext <- fread(file.path(r,"sgcb_0.5.9_kirc_external_survival.csv"))[term=="DV_only"]
fwrite(rbindlist(list(held_out=hrs,external=ext),idcol="cohort",fill=TRUE),file.path(out,"retained_primary_survival.tsv"),sep="\t")
# Retain actual generating versions, with a concrete filename for every released numerical file.
raw <- list.files(r,pattern="[.](rds|csv|tsv)$",full.names=FALSE)
map <- data.table(file=file.path("study_data/results",raw),author_build_id=NA_character_,role="published companion table or model; consult the cited source input")
map[grepl("/sgcb_0[.]5[.]7_",file),':='(author_build_id="0.5.7",role="retained mean analysis or its downstream evaluation; comparator rows retain original method identity")]
map[grepl("/sgcb_0[.]5[.]9_",file),':='(author_build_id="0.5.9",role="R1 full or downstream result; external method rows retain their own generating method")]
map[grepl("/cv_matched_families_0[.]5[.]9_",file),':='(author_build_id="0.5.9",role="matched-kernel restriction of SGCB")]
map[,bundle:="SGCB 0.6.0 study data"]; map[,software_for_refit:="SGCB 0.6.0"]; fwrite(map,file.path(out,"r1_file_index.tsv"),sep="\t",na="")
print(aucs[threshold==1]); print(resources[grepl("kirc|masld|sle",stage),.(stage,elapsed_seconds,peak_GiB)])