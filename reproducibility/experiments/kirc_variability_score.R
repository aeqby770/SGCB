CONFIG <- list(root=Sys.getenv("SGCB_REVISION_ROOT"),
    threads=as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores()))))
ROOT <- normalizePath(CONFIG$root,winslash="/",mustWork=TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,
    MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")
VERSION <- as.character(packageVersion("SGCB"))
PREFIX <- file.path(ROOT,"results",paste0("sgcb_",VERSION,"_kirc_variability"))
res <- readRDS(file.path(ROOT,"results",paste0("sgcb_",VERSION,"_kirc.rds")))
stopifnot(identical(attr(res,"native_source_id"),SOURCE_ID))
ds <- readRDS(file.path(DATA_ROOT,"data/kirc_independent_discovery.rds"))
# The predefined score uses DV genes regardless of their mean-effect status.
library(survival)
raw_tcga <- readRDS(file.path(DATA_ROOT,"data/tcga_kirc_cache/tcga_kirc.rds"))
clinical <- as.data.table(readRDS(file.path(DATA_ROOT,"data/tcga_kirc_cache/tcga_kirc_clinical.rds")))
expr <- raw_tcga$expression[,substr(colnames(raw_tcga$expression),14,15)=="01",drop=FALSE]
pid <- substr(colnames(expr),1,12)
discovery <- substr(colnames(ds$counts),1,12)
train <- pid %in% discovery[ds$group=="tumor"]
validation <- !pid %in% discovery
logx <- log2(sweep(expr,2,colSums(expr),"/")*1e6+.5)
training_mean <- rowMeans(logx[,train,drop=FALSE])
training_sd <- matrixStats::rowSds(logx[,train,drop=FALSE])
selected <- res$DV & is.finite(res$dv_log2_cv2_ratio) & is.finite(res$dv_se_log2_cv2)
genes <- intersect(res$gene_id[selected],rownames(logx)[training_sd>0])
stopifnot(length(genes)>1L)
parameters <- res[match(genes,res$gene_id),]
weight <- parameters$dv_log2_cv2_ratio/sqrt(parameters$dv_log2_cv2_ratio^2+parameters$dv_se_log2_cv2^2)
weight <- weight/sum(abs(weight))
z <- sweep(sweep(logx[genes,,drop=FALSE],1,training_mean[genes],"-"),1,training_sd[genes],"/")
scores <- cbind(mean_Z=colMeans(z),distribution_deviation=drop(crossprod(weight,z^2-1)))
scores <- sweep(scores,2,matrixStats::colSds(scores[train,,drop=FALSE]),"/")
immune_genes <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","LST1","TYROBP")
immune_z <- sweep(sweep(logx[immune_genes,,drop=FALSE],1,training_mean[immune_genes],"-"),
    1,training_sd[immune_genes],"/")
immune_proxy <- colMeans(immune_z)
immune_proxy <- immune_proxy/stats::sd(immune_proxy[train])
published <- fread(file.path(DATA_ROOT,"data/kirc_published_purity.csv"))
published[,patient:=substr(`Sample ID`,1,12)]
purity <- published[,.(CPE=median(CPE,na.rm=TRUE)),by=patient]
clinical[,`:=`(patient=submitter_id,
    time=fifelse(vital_status=="Dead",as.numeric(days_to_death),as.numeric(days_to_last_follow_up)),
    event=as.integer(vital_status=="Dead"),age=as.numeric(age_at_diagnosis)/365.25,
    stage=factor(ajcc_pathologic_stage,levels=c("Stage I","Stage II","Stage III","Stage IV")),
    grade=factor(fifelse(tumor_grade %in% c("G1","G2"),"G1/G2",
        fifelse(tumor_grade %in% c("G3","G4"),"G3/G4",NA_character_))))]
cohort <- merge(data.table(patient=pid,validation=validation,scores,immune_proxy=immune_proxy),
    clinical[,.(patient,time,event,age,stage,grade)],by="patient")
cohort <- merge(cohort,purity,by="patient",all.x=TRUE)
fwrite(cohort,paste0(PREFIX,"_tcga_patient_scores.csv"))
cohort <- cohort[validation & is.finite(time) & time>0 & complete.cases(cohort[,.(event,age,stage,grade)])]
mean_fit <- coxph(Surv(time,event)~mean_Z+age+strata(stage)+grade,data=cohort)
deviation_fit <- coxph(Surv(time,event)~distribution_deviation+age+strata(stage)+grade,data=cohort)
heldout <- rbindlist(list(
    as.data.table(broom::tidy(mean_fit,exponentiate=TRUE,conf.int=TRUE))[term=="mean_Z"],
    as.data.table(broom::tidy(deviation_fit,exponentiate=TRUE,conf.int=TRUE))[term=="distribution_deviation"]))
heldout[,`:=`(cohort="TCGA held-out",n=nrow(cohort),events=sum(cohort$event),genes=length(genes),
    qvalue=p.adjust(p.value,"BH"),model="clinical")]
composition <- cohort[is.finite(CPE) & is.finite(immune_proxy)]
mean_composition <- coxph(Surv(time,event)~mean_Z+age+strata(stage)+grade+CPE+immune_proxy,data=composition)
deviation_composition <- coxph(Surv(time,event)~distribution_deviation+age+strata(stage)+grade+CPE+immune_proxy,data=composition)
composition_result <- rbindlist(list(
    as.data.table(broom::tidy(mean_composition,exponentiate=TRUE,conf.int=TRUE))[term=="mean_Z"],
    as.data.table(broom::tidy(deviation_composition,exponentiate=TRUE,conf.int=TRUE))[term=="distribution_deviation"]))
composition_result[,`:=`(cohort="TCGA held-out",n=nrow(composition),events=sum(composition$event),
    genes=length(genes),qvalue=p.adjust(p.value,"BH"),model="clinical_plus_composition")]
heldout <- rbindlist(list(heldout,composition_result))
fwrite(heldout,paste0(PREFIX,"_heldout_survival.csv"))
fwrite(data.table(gene_id=genes,weight=weight,training_mean=training_mean[genes],
    training_sd=training_sd[genes]),paste0(PREFIX,"_discovery_score_parameters.csv"))

raw <- fread(file.path(DATA_ROOT,"data/ccRCC_exp_log_quantile_normalized.txt"))
samples <- grep("^ccRCC-",names(raw),value=TRUE)
raw[,refseq:=sub("\\..*$","",SystematicName)]
annotation <- as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys=unique(raw$refseq),keytype="REFSEQ",columns="SYMBOL"))
annotation <- unique(annotation[!is.na(SYMBOL) & SYMBOL!=""])
annotation <- annotation[,.(SYMBOL=SYMBOL[1],n_symbols=uniqueN(SYMBOL)),by=REFSEQ][n_symbols==1]
mapped <- merge(raw[,c("REF","refseq",samples),with=FALSE],annotation,by.x="refseq",by.y="REFSEQ")
external <- rowsum(as.matrix(mapped[,..samples]),mapped$SYMBOL,reorder=TRUE)
external <- external/as.numeric(table(mapped$SYMBOL)[rownames(external)])
external_sd <- matrixStats::rowSds(external)
external_genes <- intersect(genes,rownames(external)[external_sd>0])
external_weight <- weight[match(external_genes,genes)]
external_weight <- external_weight/sum(abs(external_weight))
z <- sweep(sweep(external[external_genes,,drop=FALSE],1,
    rowMeans(external[external_genes,,drop=FALSE]),"-"),1,external_sd[external_genes],"/")
scores <- cbind(mean_Z=colMeans(z),distribution_deviation=drop(crossprod(external_weight,z^2-1)))
scores <- sweep(scores,2,matrixStats::colSds(scores),"/")
clinical <- fread(file.path(DATA_ROOT,"data/sato_external_clinical.csv"))
clinical[,`:=`(patient=`sample ID`,time=as.numeric(`observation period (month)`),
    event=as.integer(outcome=="dead"),age=as.numeric(Age),
    stage_high=as.integer(as.integer(sub("^.*T([1-4]).*$","\\1",`Stage at diagnosis`))>=3),
    grade_high=as.integer(as.numeric(`Fuhrman grade`)>=3))]
cohort <- merge(data.table(patient=samples,scores),clinical[,.(patient,time,event,age,stage_high,grade_high)],by="patient")
cohort <- cohort[is.finite(time) & time>0 & complete.cases(cohort[,.(event,age,stage_high,grade_high)])]
fwrite(cohort,paste0(PREFIX,"_external_patient_scores.csv"))
mean_fit <- coxph(Surv(time,event)~mean_Z+age+stage_high+grade_high,data=cohort)
deviation_fit <- coxph(Surv(time,event)~distribution_deviation+age+stage_high+grade_high,data=cohort)
external_result <- rbindlist(list(
    as.data.table(broom::tidy(mean_fit,exponentiate=TRUE,conf.int=TRUE))[term=="mean_Z"],
    as.data.table(broom::tidy(deviation_fit,exponentiate=TRUE,conf.int=TRUE))[term=="distribution_deviation"]))
external_result[,`:=`(cohort="E-MTAB-1980",n=nrow(cohort),events=sum(cohort$event),genes=length(external_genes),
    qvalue=p.adjust(p.value,"BH"),model="clinical")]
fwrite(external_result,paste0(PREFIX,"_external_survival.csv"))
print(rbindlist(list(heldout,external_result)))
