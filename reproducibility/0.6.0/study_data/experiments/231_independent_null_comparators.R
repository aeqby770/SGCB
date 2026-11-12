ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"),winslash="/")
ARGS <- commandArgs(trailingOnly=TRUE)
CONFIG <- list(split=as.integer(ARGS[1]),method=ARGS[2],alpha=.05,
    threads=as.integer(Sys.getenv("SGCB_THREADS")),
    paths=list(data_root=normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/"),
        output=file.path(ROOT,"results")))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,
    MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
.libPaths(c(file.path(CONFIG$paths$data_root,"library_comparators"),file.path(CONFIG$paths$data_root,"library"),.libPaths()))
Sys.setenv(R_LIBS_USER=paste(.libPaths(),collapse=";"))
source(file.path(ROOT,"development/00_current_setup.R"),encoding="UTF-8")
INPUT <- file.path(CONFIG$paths$data_root,"data/gtex_readcount_independent_null",paste0("split_",CONFIG$split,".rds"))
PREFIX <- file.path(CONFIG$paths$output,paste0("readcount_independent_null_comparator_",CONFIG$method,"_",CONFIG$split))
shared <- readRDS(INPUT)
stopifnot(identical(shared$unit,"AUC-scaled integer read counts"))
counts <- shared$counts
x <- shared$normalized
sf <- shared$size_factors
group <- shared$group
n0 <- sum(group==levels(group)[1])
n1 <- sum(group==levels(group)[2])
design <- model.matrix(~group)
logx <- log(x+.5)
stopifnot(identical(rownames(counts),shared$gene_id),identical(colnames(counts),shared$samples$sample),
    !anyDuplicated(shared$samples$donor))
fit <- NULL
elapsed <- system.time(result <- switch(CONFIG$method,
    edgeR={
        dge <- edgeR::DGEList(counts)
        dge$samples$norm.factors <- sf/(dge$samples$lib.size/exp(mean(log(dge$samples$lib.size))))
        dge <- edgeR::estimateDisp(dge,design,robust=TRUE)
        fit <- edgeR::glmQLFTest(edgeR::glmQLFit(dge,design,robust=TRUE),coef=2)$table
        data.table(gene_id=shared$gene_id,pvalue=fit$PValue,effect=fit$logFC,target="Mean expression")
    },
    limma={
        dge <- edgeR::DGEList(counts)
        dge$samples$norm.factors <- sf/(dge$samples$lib.size/exp(mean(log(dge$samples$lib.size))))
        voom <- limma::voom(dge,design,plot=FALSE)
        fit <- limma::eBayes(limma::lmFit(voom,design),robust=TRUE)
        data.table(gene_id=shared$gene_id,pvalue=fit$p.value[,2],
            effect=fit$coefficients[,2],target="Mean expression")
    },
    DESeq2={
        stopifnot(all(counts==round(counts)),max(counts)<.Machine$integer.max)
        dds <- DESeq2::DESeqDataSetFromMatrix(counts,
            data.frame(group=group,row.names=colnames(counts)),~group)
        DESeq2::sizeFactors(dds) <- sf
        dds <- DESeq2::DESeq(dds,quiet=TRUE,parallel=FALSE)
        fit <- DESeq2::results(dds,independentFiltering=FALSE)
        data.table(gene_id=rownames(fit),pvalue=fit$pvalue,
            effect=fit$log2FoldChange,target="Mean expression")
    },
    diffVar={
        source(file.path(CONFIG$paths$data_root,"data/external_methods/DiffVar/DiffVar.R"),encoding="UTF-8")
        dge <- edgeR::DGEList(counts)
        dge$samples$norm.factors <- sf/(dge$samples$lib.size/exp(mean(log(dge$samples$lib.size))))
        fit <- varFit(dge,design=design,coef=c(1,2),type="AD")
        data.table(gene_id=rownames(logx),pvalue=fit$p.value[,2],
            effect=fit$coefficients[,2],target="Log spread")
    },
    clrDV={
        library(sn)
        source(file.path(CONFIG$paths$data_root,"data/external_methods/clrDV/clrDV.R"),encoding="UTF-8")
        clr_logx <- log(x+.5)
        clr <- sweep(clr_logx,2,colMeans(clr_logx),"-")
        fit <- clrDV(clr,group)
        data.table(gene_id=rownames(fit),pvalue=fit$pval,effect=fit$DV,
            original_qvalue=fit$adj_pval,target="CLR dispersion")
    },
    MDSeq={
        library(MDSeq)
        parallel:::setDefaultClusterOptions(rscript_args="--vanilla")
        md_design <- MDSeq::get.model.matrix(group)
        mean_covariate <- md_design$mean[,colnames(md_design$mean)!="(Intercept)",drop=FALSE]
        dispersion_covariate <- md_design$dispersion[,colnames(md_design$dispersion)!="(Intercept)",drop=FALSE]
        fit <- MDSeq::MDSeq(counts,X=mean_covariate,U=dispersion_covariate,
            offsets=sf,mc.cores=CONFIG$threads,verbose=FALSE)
        native <- as.data.table(fit$Dat,keep.rownames="gene_id")
        result <- data.table(gene_id=shared$gene_id)[native,on="gene_id"]
        result[,pvalue:=ifelse(is.finite(se.gamma.group1) & se.gamma.group1>0,
            2*pnorm(-abs(gamma.group1/se.gamma.group1)),NA_real_)]
        result[,.(gene_id,pvalue,effect=-2*gamma.group1/log(2),target="Count dispersion")]
    },
    GAMLSS={
        library(gamlss)
        dge <- edgeR::DGEList(counts)
        dge$samples$norm.factors <- sf/(dge$samples$lib.size/exp(mean(log(dge$samples$lib.size))))
        log_cpm <- edgeR::cpm(dge,log=TRUE,prior.count=1)
        v0 <- matrixStats::rowVars(log_cpm[,group==levels(group)[1],drop=FALSE])
        v1 <- matrixStats::rowVars(log_cpm[,group==levels(group)[2],drop=FALSE])
        mad0 <- matrixStats::rowMads(log_cpm[,group==levels(group)[1],drop=FALSE])
        mad1 <- matrixStats::rowMads(log_cpm[,group==levels(group)[2],drop=FALSE])
        estimable <- is.finite(v0) & is.finite(v1) & is.finite(mad0) & is.finite(mad1) &
            v0>0 & v1>0 & mad0>0 & mad1>0
        long <- as.data.table(log_cpm)
        long[,gene_id:=shared$gene_id]
        long <- melt(long,id.vars="gene_id",
            variable.name="sample",value.name="y")
        long[,g:=as.integer(shared$group[match(as.character(sample),colnames(counts))])-1L]
        fit <- long[gene_id %chin% shared$gene_id[estimable],{
            native_data <- data.frame(y=y,g=g)
            h1 <- do.call(gamlss::gamlss,list(formula=y~g,sigma.formula=~g,family=gamlss.dist::NO(),
                data=native_data,trace=FALSE,
                control=gamlss::gamlss.control(n.cyc=50,trace=FALSE)))
            h0 <- do.call(gamlss::gamlss,list(formula=y~g,sigma.formula=~1,family=gamlss.dist::NO(),
                data=native_data,trace=FALSE,
                control=gamlss::gamlss.control(n.cyc=50,trace=FALSE)))
            lrt <- max(as.numeric(2*(logLik(h1)-logLik(h0))),0)
            .(pvalue=ifelse(h1$converged & h0$converged,
                pchisq(lrt,1,lower.tail=FALSE),NA_real_),
              effect=2*coef(h1,what="sigma")[["g"]]/log(2),
              alternative_deviance=h1$G.deviance,null_deviance=h0$G.deviance,
              alternative_converged=h1$converged,null_converged=h0$converged)
        },by=gene_id]
        result <- merge(data.table(gene_id=shared$gene_id),fit,by="gene_id",all.x=TRUE,sort=FALSE)
        result[,target:="Log spread"]
        result
    },

    Brown_Forsythe={
        left <- logx[,group==levels(group)[1],drop=FALSE]
        right <- logx[,group==levels(group)[2],drop=FALSE]
        d0 <- abs(left-matrixStats::rowMedians(left))
        d1 <- abs(right-matrixStats::rowMedians(right))
        mean0 <- rowMeans(d0)
        mean1 <- rowMeans(d1)
        within <- rowSums((d0-mean0)^2)+rowSums((d1-mean1)^2)
        statistic <- (n0*n1/(n0+n1))*(mean1-mean0)^2/(within/(n0+n1-2))
        data.table(gene_id=rownames(logx),pvalue=pf(statistic,1,n0+n1-2,lower.tail=FALSE),
            effect=mean1-mean0,target="Log spread")
    },
    Gamma={
        fit <- SGCB:::.populationGG(x,as.integer(group)-1L,CONFIG$threads,
            CONFIG$alpha,family="Gamma",truth=rep(0,nrow(x)))
        data.table(gene_id=shared$gene_id,pvalue=replace(fit$genes$p,fit$genes$available==0,NA_real_),
            effect=fit$genes$effect,covered=as.logical(fit$genes$covered),target="CV2")
    },
    lognormal={
        fit <- SGCB:::.populationGG(x,as.integer(group)-1L,CONFIG$threads,
            CONFIG$alpha,family="lognormal",truth=rep(0,nrow(x)))
        data.table(gene_id=shared$gene_id,pvalue=replace(fit$genes$p,fit$genes$available==0,NA_real_),
            effect=fit$genes$effect,covered=as.logical(fit$genes$covered),target="CV2")
    }))
stopifnot(nrow(result)==length(shared$gene_id),setequal(result$gene_id,shared$gene_id))
result <- result[match(shared$gene_id,gene_id)]
result[,available:=is.finite(pvalue)]
result[,qvalue:=p.adjust(pvalue,"BH",n=nrow(counts))]
result[,called:=!is.na(qvalue) & qvalue<CONFIG$alpha]
result[,c("split","method"):=list(CONFIG$split,CONFIG$method)]
metrics <- result[,.(split=CONFIG$split,method=CONFIG$method,target=target[1],
    genes=.N,available=sum(available),calls=sum(called),has_call=as.integer(any(called)),
    raw_p005_full=sum(pvalue<CONFIG$alpha,na.rm=TRUE)/.N,
    raw_p005_available=mean(pvalue[available]<CONFIG$alpha),
    seconds=elapsed[["elapsed"]],n_A=n0,n_B=n1)]
saveRDS(list(genes=result,metrics=metrics,fit=fit,input=INPUT,samples=shared$samples,
    size_factors=sf,source_id=SOURCE_ID,version=as.character(packageVersion("SGCB")),
    config=CONFIG,session_info=sessionInfo()),paste0(PREFIX,".rds"),compress=FALSE)
fwrite(result,paste0(PREFIX,"_genes.csv"),bom=FALSE,eol="\n")
fwrite(metrics,paste0(PREFIX,"_metrics.csv"),bom=FALSE,eol="\n")
print(metrics)
