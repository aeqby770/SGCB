# =============================================================================
# Bulk RNA-seq Benchmark - Other Methods Only
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(ggplot2)
library(patchwork)

DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark/data")
BULK_DIR <- file.path(DATA_DIR, "bulk")
RESULTS_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# =============================================================================
#
# (SGCB removed - tested separately)
# 7. samr      - SAM/SAMseq
#
#
# =============================================================================

# =============================================================================
# =============================================================================
#
#
# - sim_n3_pDE10_ef1.5.rds  : n=3, 10% DE, effect=1.5
# - sim_n3_pDE10_ef2.0.rds  : n=3, 10% DE, effect=2.0
# - sim_n5_pDE10_ef1.5.rds  : n=5, 10% DE, effect=1.5
# - sim_n5_pDE10_ef2.0.rds  : n=5, 10% DE, effect=2.0
# - sim_n10_pDE10_ef1.5.rds : n=10, 10% DE, effect=1.5
# - sim_n10_pDE10_ef2.0.rds : n=10, 10% DE, effect=2.0
#
# - GTEx: adipose, blood, brain, colon, heart, kidney, liver, lung, muscle, skin
# - TCGA: brca, coad, kirc, lihc, luad, lusc, prad
#
# =============================================================================

cat("\n========== ==========\n")
library(DESeq2)
library(edgeR)
library(limma)
library(NOISeq)
library(EBSeq)
library(samr)
library(glmGamPoi)

EXTERNAL_METHODS_DIR <- paste0(PROJECT_ROOT, "/external_methods")
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/MDSeq.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/check.input.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/glm.ZIMD.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/class.ZIMD.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/get.model.matrix.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/formula.NB.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/formula.ZINB.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/initial.value.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "MDSeq/R/inequality.test.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "clrDV/R/clr.SN.fit.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "clrDV/R/clrSeq.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "ComBat_Seq/src/ComBat_seq.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "ComBat_Seq/src/helper_seq.R"), local = TRUE)
source(file.path(EXTERNAL_METHODS_DIR, "ComBat-ref/Combat_ref.R"), local = TRUE)

library(sva)
library(RUVSeq)
library(gamlss)
library(DiPhiSeq)
library(sn)
# missMethyl removed - use limma::varFit directly for DV

# =============================================================================
# =============================================================================

run_deseq2 <- function(counts, group) {
    t0 <- Sys.time()
    col_data <- data.frame(condition = group, row.names = colnames(counts))
    dds <- DESeqDataSetFromMatrix(counts, col_data, ~condition)
    dds <- tryCatch(
        estimateSizeFactors(dds, type = "poscounts"), 
        error = function(e) estimateSizeFactors(dds)
    )
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res), 
        pvalue = res$pvalue, 
        padj = res$padj,
        log2FC = res$log2FoldChange, 
        lfcSE = res$lfcSE,
        baseMean = res$baseMean,
        stat = res$stat,
        method = "DESeq2", 
        time = time
    )
}

run_edger <- function(counts, group) {
    t0 <- Sys.time()
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    res <- topTags(qlf, n = Inf, sort.by = "none")$table
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res), 
        pvalue = res$PValue, 
        padj = res$FDR,
        log2FC = res$logFC, 
        lfcSE = NA_real_,
        baseMean = res$logCPM,
        stat = res$F,
        method = "edgeR", 
        time = time
    )
}

run_limma <- function(counts, group) {
    t0 <- Sys.time()
    design <- model.matrix(~group)
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res), 
        pvalue = res$P.Value, 
        padj = res$adj.P.Val,
        log2FC = res$logFC, 
        lfcSE = sqrt(fit$s2.post) * fit$stdev.unscaled[rownames(res), 2],
        baseMean = res$AveExpr,
        stat = res$t,
        method = "limma", 
        time = time
    )
}

run_noiseq <- function(counts, group) {
    t0 <- Sys.time()
    factors <- data.frame(condition = group, row.names = colnames(counts))
    noiseq_data <- NOISeq::readData(data = as.matrix(counts), factors = factors)
    res <- NOISeq::noiseqbio(noiseq_data, k = 0.5, norm = "tmm", factor = "condition", 
                              conditions = levels(group), r = 20)
    res_df <- res@results[[1]]
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res_df),
        pvalue = 1 - res_df$prob,
        padj = p.adjust(1 - res_df$prob, method = "BH"),
        log2FC = res_df$log2FC,
        lfcSE = NA_real_,
        baseMean = rowMeans(counts),
        stat = res_df$prob,
        method = "NOISeq",
        time = time
    )
}

run_ebseq <- function(counts, group) {
    t0 <- Sys.time()
    sizes <- MedianNorm(counts)
    res <- EBTest(Data = counts, Conditions = as.character(group), sizeFactors = sizes, maxround = 5)
    ppmat <- GetPPMat(res)
    fc <- PostFC(res)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(ppmat),
        pvalue = 1 - ppmat[, "PPDE"],
        padj = p.adjust(1 - ppmat[, "PPDE"], method = "BH"),
        log2FC = log2(fc$PostFC),
        lfcSE = NA_real_,
        baseMean = rowMeans(counts),
        stat = ppmat[, "PPDE"],
        method = "EBSeq",
        time = time
    )
}

run_samr <- function(counts, group) {
    t0 <- Sys.time()
    y <- as.numeric(factor(group))
    x_mat <- as.matrix(log2(counts + 1))
    samr_data <- list(x = x_mat, y = y, geneid = rownames(counts), genenames = rownames(counts), logged2 = TRUE)
    samr_obj <- samr(samr_data, resp.type = "Two class unpaired", nperms = 100)
    delta_table <- samr.compute.delta.table(samr_obj)
    if (nrow(delta_table) > 0 && "FDR" %in% colnames(delta_table)) {
        fdr_col <- delta_table[, "FDR"]
        valid_idx <- which(!is.na(fdr_col) & is.finite(fdr_col))
        if (length(valid_idx) > 0) {
            delta_opt <- delta_table[valid_idx[which.min(abs(fdr_col[valid_idx] - 0.05))], "Delta"]
        } else {
            delta_opt <- delta_table[1, "Delta"]
        }
    } else {
        delta_opt <- 0.5
    }
    siggenes <- samr.compute.siggenes.table(samr_obj, delta_opt, samr_data, delta.table = delta_table)
    sig_genes_up <- if (!is.null(siggenes$genes.up) && nrow(siggenes$genes.up) > 0) siggenes$genes.up[, "Gene ID"] else character(0)
    sig_genes_down <- if (!is.null(siggenes$genes.lo) && nrow(siggenes$genes.lo) > 0) siggenes$genes.lo[, "Gene ID"] else character(0)
    sig_genes <- c(sig_genes_up, sig_genes_down)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    padj_val <- rep(1, nrow(counts))
    padj_val[rownames(counts) %in% sig_genes] <- 0.01
    lfc <- rowMeans(log2(counts[, group == levels(group)[2], drop=FALSE] + 1)) - 
           rowMeans(log2(counts[, group == levels(group)[1], drop=FALSE] + 1))
    data.table(
        gene = rownames(counts),
        pvalue = padj_val,
        padj = padj_val,
        log2FC = lfc,
        lfcSE = NA_real_,
        baseMean = rowMeans(counts),
        stat = NA_real_,
        method = "samr",
        time = time
    )
}

run_glmgampoi <- function(counts, group) {
    t0 <- Sys.time()
    col_data <- data.frame(condition = group, row.names = colnames(counts))
    fit <- glm_gp(counts, design = ~condition, col_data = col_data, size_factors = "normed_sum")
    contrast_str <- paste0("condition", levels(group)[2])
    res <- test_de(fit, contrast = contrast_str)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = res$name,
        pvalue = res$pval,
        padj = res$adj_pval,
        log2FC = res$lfc,
        lfcSE = NA_real_,
        baseMean = rowMeans(counts),
        stat = res$f_statistic,
        method = "glmGamPoi",
        time = time
    )
}

METHODS <- list(
    DESeq2 = run_deseq2,
    edgeR = run_edger,
    limma = run_limma,
    NOISeq = run_noiseq,
    EBSeq = run_ebseq,
    samr = run_samr,
    glmGamPoi = run_glmgampoi
)

# =============================================================================
# =============================================================================

run_mdseq_dv <- function(counts, group) {
    t0 <- Sys.time()
    counts <- round(counts)
    storage.mode(counts) <- "integer"
    group <- factor(group)
    design <- get.model.matrix(group)
    fit <- MDSeq(
        counts,
        X = NULL,
        U = NULL,
        contrast = design,
        offsets = rep(1, ncol(counts)),
        verbose = FALSE,
        mc.cores = 1
    )
    res <- extract.ZIMD(
        fit,
        compare = list(A = levels(group)[2], B = levels(group)[1]),
        p.adj = "BH"
    )
    mean_lfc_col <- grep("mean\\.log2FC", colnames(res), value = TRUE)[1]
    disp_lfc_col <- grep("dispersion\\.log2FC", colnames(res), value = TRUE)[1]
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res),
        pvalue_dv = res[, "Pvalue.dispersion"],
        padj_dv = res[, "FDR.dispersion"],
        log2VarRatio = res[, disp_lfc_col],
        pvalue_de = res[, "Pvalue.mean"],
        padj_de = res[, "FDR.mean"],
        log2FC = res[, mean_lfc_col],
        method = "MDSeq",
        time = time
    )
}

run_clrdv <- function(counts, group) {
    t0 <- Sys.time()
    group <- factor(group)
    clr_counts <- log2(counts + 1) - rowMeans(log2(counts + 1))
    group_vec <- as.numeric(group) - 1
    res <- clrSeq(clr_counts, group_vec)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    pvals_sigma1 <- res$p.sigma1
    pvals_sigma2 <- res$p.sigma2
    pvals <- pmin(pvals_sigma1, pvals_sigma2, na.rm = TRUE) * 2
    pvals[is.na(pvals)] <- 1
    pvals <- pmin(pvals, 1)
    data.table(
        gene = rownames(counts),
        pvalue_dv = pvals,
        padj_dv = p.adjust(pvals, method = "BH"),
        log2VarRatio = log2((res$sigma2 + 1e-8) / (res$sigma1 + 1e-8)),
        method = "clrDV",
        time = time
    )
}

run_diffvar <- function(counts, group) {
    t0 <- Sys.time()
    design <- model.matrix(~group)
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    residuals <- abs(residuals(fit, v))
    fit_var <- lmFit(residuals, design)
    fit_var <- eBayes(fit_var)
    res <- topTable(fit_var, coef = 2, number = Inf, sort.by = "none")
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(res),
        pvalue_dv = res$P.Value,
        padj_dv = res$adj.P.Val,
        log2VarRatio = res$logFC,
        method = "diffVar",
        time = time
    )
}

run_gamlss_dv <- function(counts, group) {
    t0 <- Sys.time()
    group <- factor(group)
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
    ctrl_idx <- which(group == levels(group)[1])
    treat_idx <- which(group == levels(group)[2])
    g <- as.numeric(group) - 1
    pvals <- vapply(seq_len(nrow(log_cpm)), function(i) {
        y <- log_cpm[i, ]
        p <- NA_real_
        tryCatch({
            fit_h1 <- gamlss::gamlss(y ~ g, sigma.formula = ~ g, family = gamlss.dist::NO(),
                                      trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50))
            fit_h0 <- gamlss::gamlss(y ~ g, sigma.formula = ~ 1, family = gamlss.dist::NO(),
                                      trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50))
            lr <- -2 * (logLik(fit_h0) - logLik(fit_h1))
            lr <- max(as.numeric(lr), 0)
            p <- pchisq(lr, df = 1, lower.tail = FALSE)
        }, error = function(e) {})
        p
    }, numeric(1))
    var_ctrl <- matrixStats::rowVars(log_cpm[, ctrl_idx])
    var_treat <- matrixStats::rowVars(log_cpm[, treat_idx])
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(counts),
        pvalue_dv = pvals,
        padj_dv = p.adjust(pvals, method = "BH"),
        log2VarRatio = log2(var_treat / var_ctrl),
        method = "GAMLSS",
        time = time
    )
}

run_diphiseq_dv <- function(counts, group) {
    t0 <- Sys.time()
    group <- factor(group)
    counts <- round(counts)
    storage.mode(counts) <- "integer"
    classlab <- as.integer(group == levels(group)[2]) + 1L
    res <- DiPhiSeq::diphiseq(countmat = counts, classlab = classlab)
    tab <- res$tab
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data.table(
        gene = rownames(counts),
        pvalue_dv = tab[, "p.value.phi"],
        padj_dv = tab[, "fdr.phi"],
        log2VarRatio = log2((tab[, "phi2"] + 1e-8) / (tab[, "phi1"] + 1e-8)),
        method = "DiPhiSeq",
        time = time
    )
}

DV_METHODS <- list(
    MDSeq = run_mdseq_dv,
    clrDV = run_clrdv,
    diffVar = run_diffvar,
    GAMLSS = run_gamlss_dv,
    DiPhiSeq = run_diphiseq_dv
)

# =============================================================================
# =============================================================================

run_combat_seq <- function(counts, group, batch) {
    t0 <- Sys.time()
    adj_counts <- ComBat_seq(counts, batch = batch, group = group)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(counts = adj_counts, time = time, method = "ComBat-seq")
}

run_ruvseq <- function(counts, group, k = 1) {
    t0 <- Sys.time()
    set <- newSeqExpressionSet(as.matrix(counts), phenoData = data.frame(group = group, row.names = colnames(counts)))
    set <- betweenLaneNormalization(set, which = "upper")
    design <- model.matrix(~group)
    y <- DGEList(counts = counts(set), group = group)
    y <- calcNormFactors(y, method = "upper")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    ctrl_genes <- rownames(topTags(lrt, n = 5000, sort.by = "none")$table)[topTags(lrt, n = 5000, sort.by = "none")$table$PValue > 0.5]
    ctrl_genes <- head(ctrl_genes, min(1000, length(ctrl_genes)))
    set_ruv <- RUVg(set, ctrl_genes, k = k)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(counts = normCounts(set_ruv), W = pData(set_ruv)$W_1, time = time, method = "RUVSeq")
}

run_sva_correct <- function(counts, group) {
    t0 <- Sys.time()
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
    mod <- model.matrix(~group)
    mod0 <- model.matrix(~1, data = data.frame(row.names = colnames(counts)))
    sv_obj <- sva(log_cpm, mod, mod0, n.sv = NULL)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(sv = sv_obj$sv, n.sv = sv_obj$n.sv, time = time, method = "SVA")
}

run_combat_ref <- function(counts, group, batch) {
    t0 <- Sys.time()
    adj_counts <- ComBat_ref(counts = counts, batch = batch, group = group)
    time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(counts = adj_counts, time = time, method = "ComBat-ref")
}

BATCH_METHODS <- list(
    ComBat_seq = run_combat_seq,
    RUVSeq = run_ruvseq,
    SVA = run_sva_correct,
    ComBat_ref = run_combat_ref
)

# =============================================================================
# =============================================================================

compute_metrics <- function(dt, truth = NULL, fdr_levels = c(0.01, 0.05, 0.1)) {
    dt <- dt[!is.na(pvalue) & !is.na(padj)]
    
    metrics <- data.table(
        method = dt$method[1],
        n_genes = nrow(dt),
        time_sec = round(dt$time[1], 2)
    )
    
    for (fdr in fdr_levels) {
        col_name <- paste0("n_sig_fdr", fdr * 100)
        metrics[[col_name]] <- sum(dt$padj < fdr, na.rm = TRUE)
    }
    
    for (lfc in c(0.5, 1, 2)) {
        col_name <- paste0("n_sig_lfc", lfc)
        metrics[[col_name]] <- sum(dt$padj < 0.05 & abs(dt$log2FC) > lfc, na.rm = TRUE)
    }
    
    metrics$p_median <- median(dt$pvalue, na.rm = TRUE)
    metrics$p_mean <- mean(dt$pvalue, na.rm = TRUE)
    metrics$p_lt_0.05 <- mean(dt$pvalue < 0.05, na.rm = TRUE)
    
    metrics$lfc_median <- median(abs(dt$log2FC), na.rm = TRUE)
    metrics$lfc_iqr <- IQR(dt$log2FC, na.rm = TRUE)
    
    if (!is.null(truth) && length(truth) > 0) {
        for (fdr in fdr_levels) {
            sig_genes <- dt$gene[dt$padj < fdr]
            TP <- sum(sig_genes %in% truth)
            FP <- sum(!(sig_genes %in% truth))
            FN <- sum(!(truth %in% sig_genes))
            TN <- nrow(dt) - TP - FP - FN
            
            precision <- TP / max(1, TP + FP)
            recall <- TP / max(1, length(truth))
            f1 <- 2 * precision * recall / max(0.001, precision + recall)
            fdr_actual <- FP / max(1, TP + FP)
            
            suffix <- paste0("_fdr", fdr * 100)
            metrics[[paste0("TP", suffix)]] <- TP
            metrics[[paste0("FP", suffix)]] <- FP
            metrics[[paste0("precision", suffix)]] <- round(precision, 4)
            metrics[[paste0("recall", suffix)]] <- round(recall, 4)
            metrics[[paste0("F1", suffix)]] <- round(f1, 4)
            metrics[[paste0("actual_FDR", suffix)]] <- round(fdr_actual, 4)
        }
        
        if (requireNamespace("pROC", quietly = TRUE)) {
            is_de <- as.integer(dt$gene %in% truth)
            if (sum(is_de) > 0 && sum(is_de) < length(is_de)) {
                roc_obj <- pROC::roc(is_de, 1 - dt$pvalue, quiet = TRUE)
                metrics$AUC_ROC <- round(as.numeric(pROC::auc(roc_obj)), 4)
                
                pr_obj <- pROC::roc(is_de, 1 - dt$pvalue, quiet = TRUE)
                metrics$AUC_PR <- round(as.numeric(pROC::auc(pr_obj)), 4)
            }
        }
    }
    
    metrics
}

compute_concordance <- function(results_list, fdr = 0.05, top_n = NULL) {
    methods <- names(results_list)
    n_methods <- length(methods)
    
    sig_genes <- lapply(results_list, function(dt) {
        if (is.null(top_n)) {
            dt$gene[dt$padj < fdr]
        } else {
            dt$gene[order(dt$pvalue)][1:min(top_n, nrow(dt))]
        }
    })
    
    concordance_mat <- matrix(0, n_methods, n_methods, dimnames = list(methods, methods))
    for (i in seq_len(n_methods)) {
        for (j in seq_len(n_methods)) {
            g1 <- sig_genes[[i]]
            g2 <- sig_genes[[j]]
            concordance_mat[i, j] <- length(intersect(g1, g2)) / length(union(g1, g2))
        }
    }
    
    concordance_mat
}

compute_lfc_correlation <- function(results_list) {
    methods <- names(results_list)
    n_methods <- length(methods)
    
    if (n_methods < 2) {
        return(list(pearson = matrix(1), spearman = matrix(1), n_genes = 0))
    }
    
    all_genes <- Reduce(union, lapply(results_list, function(x) x$gene))
    
    if (length(all_genes) < 10) {
        return(list(pearson = matrix(NA, n_methods, n_methods),
                    spearman = matrix(NA, n_methods, n_methods), n_genes = length(all_genes)))
    }
    
    lfc_mat <- sapply(results_list, function(dt) {
        setNames(dt$log2FC, dt$gene)[all_genes]
    })
    n_pair <- crossprod(!is.na(lfc_mat))
    
    cor_pearson <- cor(lfc_mat, use = "pairwise.complete.obs", method = "pearson")
    cor_spearman <- cor(lfc_mat, use = "pairwise.complete.obs", method = "spearman")
    
    list(pearson = cor_pearson, spearman = cor_spearman, n_genes = n_pair)
}

compute_dv_metrics <- function(dt, truth_dv = NULL, fdr = 0.05) {
    dt <- dt[!is.na(pvalue_dv) & !is.na(padj_dv)]
    
    metrics <- data.table(
        method = dt$method[1],
        n_genes = nrow(dt),
        n_sig_dv = sum(dt$padj_dv < fdr, na.rm = TRUE),
        time_sec = round(dt$time[1], 2)
    )
    
    if (!is.null(truth_dv) && length(truth_dv) > 0) {
        sig_genes <- dt$gene[dt$padj_dv < fdr]
        TP <- sum(sig_genes %in% truth_dv)
        FP <- sum(!(sig_genes %in% truth_dv))
        FN <- sum(!(truth_dv %in% sig_genes))
        precision <- TP / max(1, TP + FP)
        recall <- TP / max(1, length(truth_dv))
        f1 <- 2 * precision * recall / max(0.001, precision + recall)
        actual_fdr <- FP / max(1, TP + FP)
        
        metrics$TP <- TP
        metrics$FP <- FP
        metrics$FN <- FN
        metrics$precision <- round(precision, 4)
        metrics$recall <- round(recall, 4)
        metrics$F1 <- round(f1, 4)
        metrics$actual_FDR <- round(actual_fdr, 4)
    }
    
    metrics
}

compute_batch_metrics <- function(counts_orig, counts_adj, batch, group, truth_de = NULL) {
    pca_orig <- prcomp(t(log2(counts_orig + 1)), scale. = TRUE)
    pca_adj <- prcomp(t(log2(counts_adj + 1)), scale. = TRUE)
    
    batch_var_orig <- summary(aov(pca_orig$x[, 1] ~ batch))[[1]][1, 2] / sum(summary(aov(pca_orig$x[, 1] ~ batch))[[1]][, 2])
    batch_var_adj <- summary(aov(pca_adj$x[, 1] ~ batch))[[1]][1, 2] / sum(summary(aov(pca_adj$x[, 1] ~ batch))[[1]][, 2])
    
    group_var_orig <- summary(aov(pca_orig$x[, 1] ~ group))[[1]][1, 2] / sum(summary(aov(pca_orig$x[, 1] ~ group))[[1]][, 2])
    group_var_adj <- summary(aov(pca_adj$x[, 1] ~ group))[[1]][1, 2] / sum(summary(aov(pca_adj$x[, 1] ~ group))[[1]][, 2])
    
    de_metrics <- NULL
    if (!is.null(truth_de) && length(truth_de) > 0) {
        col_data <- data.frame(condition = group, row.names = colnames(counts_adj))
        dds <- DESeqDataSetFromMatrix(round(counts_adj), col_data, ~condition)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds)
        sig_genes <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
        TP <- sum(sig_genes %in% truth_de)
        FP <- sum(!(sig_genes %in% truth_de))
        de_metrics <- data.table(TP = TP, FP = FP, 
                                  precision = round(TP / max(1, TP + FP), 4),
                                  recall = round(TP / max(1, length(truth_de)), 4))
    }
    
    list(
        batch_var_orig = round(batch_var_orig, 4),
        batch_var_adj = round(batch_var_adj, 4),
        batch_reduction = round(1 - batch_var_adj / max(0.001, batch_var_orig), 4),
        group_var_orig = round(group_var_orig, 4),
        group_var_adj = round(group_var_adj, 4),
        signal_retention = round(group_var_adj / max(0.001, group_var_orig), 4),
        de_metrics = de_metrics
    )
}

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 1. Bulk RNA-seq \n")
cat(strrep("=", 70), "\n")

classic_datasets <- list(
    bottomly = "real_bottomly.rds",
    hammer = "real_hammer.rds",
    pasilla = "real_pasilla.rds",
    airway = "real_airway.rds",
    sultan = "real_sultan.rds",
    wang = "real_wang.rds"
)

classic_results <- list()
classic_metrics <- list()

for (ds_name in names(classic_datasets)) {
    ds_file <- file.path(DATA_DIR, classic_datasets[[ds_name]])
    
    cat("\n========== ", ds_name, " ==========\n", sep = "")
    
    tryCatch({
        ds <- readRDS(ds_file)
        
        if ("counts" %in% names(ds) && "sample_info" %in% names(ds)) {
            counts <- ds$counts
            group <- ds$sample_info$group
        } else if ("counts" %in% names(ds) && "group" %in% names(ds)) {
            counts <- ds$counts
            group <- ds$group
        } else if (inherits(ds, "SummarizedExperiment")) {
            counts <- SummarizedExperiment::assay(ds, "counts")
            group <- SummarizedExperiment::colData(ds)$condition
        } else if (inherits(ds, "DESeqDataSet")) {
            counts <- DESeq2::counts(ds)
            group <- SummarizedExperiment::colData(ds)$condition
        } else {
 cat(" \n")
            next
        }
        
        counts <- as.matrix(counts)
        group <- factor(group)
        
        keep <- rowSums(counts > 5) >= 2
        counts <- counts[keep, ]
        
 cat(" : ", ncol(counts), " (", table(group)[1], " vs ", table(group)[2], ")\n", sep = "")
 cat(" : ", nrow(counts), "\n", sep = "")
        
        ds_results <- list()
        for (m_name in names(METHODS)) {
            cat("    Running ", m_name, "... ", sep = "")
            ds_results[[m_name]] <- tryCatch({
                res <- METHODS[[m_name]](counts, group)
                cat("OK (", round(res$time[1], 1), "s)\n", sep = "")
                res
            }, error = function(e) {
                cat("FAILED: ", conditionMessage(e), "\n", sep = "")
                NULL
            })
        }
        
        ds_results <- ds_results[!sapply(ds_results, is.null)]
        classic_results[[ds_name]] <- ds_results
        
        ds_metrics <- rbindlist(lapply(ds_results, compute_metrics), fill = TRUE)
        ds_metrics[, dataset := ds_name]
        classic_metrics[[ds_name]] <- ds_metrics
        
        concord <- compute_concordance(ds_results)
        lfc_cor <- compute_lfc_correlation(ds_results)
        
 cat("\n (FDR < 0.05)\n")
        print(ds_metrics[, .(method, n_sig_fdr5, n_sig_lfc1, time_sec)])
        
 cat("\n LFC Spearman \n")
        print(round(lfc_cor$spearman, 3))
        
    }, error = function(e) {
 cat(" : ", conditionMessage(e), "\n", sep = "")
    })
}

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 2. (FDR/Power )\n")
cat(strrep("=", 70), "\n")

sim_files <- list.files(DATA_DIR, pattern = "^sim_.*\\.rds$", full.names = TRUE)
sim_metrics <- list()

for (sim_file in sim_files) {
    sim_name <- gsub("\\.rds$", "", basename(sim_file))
    
    cat("\n========== ", sim_name, " ==========\n", sep = "")
    
    tryCatch({
        sim <- readRDS(sim_file)
        
        counts <- sim$counts
        group <- factor(sim$sample_info$group)
        
        if ("truth" %in% names(sim) && is.data.frame(sim$truth)) {
            truth <- sim$truth$gene[sim$truth$is_de == TRUE]
        } else if ("de_genes" %in% names(sim)) {
            truth <- sim$de_genes
        } else {
            truth <- character(0)
        }
        
        keep <- rowSums(counts > 5) >= 2
        counts <- counts[keep, ]
        truth <- intersect(truth, rownames(counts))
        
 cat(" : ", ncol(counts), "\n", sep = "")
 cat(" : ", nrow(counts), " ( DE: ", length(truth), ")\n", sep = "")
        
        sim_results <- list()
        for (m_name in names(METHODS)) {
            cat("    Running ", m_name, "... ", sep = "")
            sim_results[[m_name]] <- tryCatch({
                res <- METHODS[[m_name]](counts, group)
                cat("OK\n")
                res
            }, error = function(e) {
                cat("FAILED\n")
                NULL
            })
        }
        
        sim_results <- sim_results[!sapply(sim_results, is.null)]
        
        sim_ds_metrics <- rbindlist(lapply(sim_results, function(x) compute_metrics(x, truth)), fill = TRUE)
        sim_ds_metrics[, dataset := sim_name]
        sim_metrics[[sim_name]] <- sim_ds_metrics
        
 cat("\n FDR Power (FDR < 0.05)\n")
        print(sim_ds_metrics[, .(method, TP_fdr5, FP_fdr5, actual_FDR_fdr5, recall_fdr5, F1_fdr5)])
        
    }, error = function(e) {
 cat(" : ", conditionMessage(e), "\n", sep = "")
    })
}

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 3. (GTEx/TCGA)\n")
cat(strrep("=", 70), "\n")

large_metrics <- list()

gtex_pairs <- list(
    c("liver", "kidney"),
    c("brain", "heart"),
    c("lung", "colon"),
    c("adipose_tissue", "muscle"),
    c("blood", "skin")
)

for (pair in gtex_pairs) {
    t1 <- pair[1]
    t2 <- pair[2]
    
    f1 <- file.path(BULK_DIR, paste0("gtex_", t1, ".rds"))
    f2 <- file.path(BULK_DIR, paste0("gtex_", t2, ".rds"))
    
    if (!file.exists(f1) || !file.exists(f2)) next
    
    cat("\n========== GTEx: ", t1, " vs ", t2, " ==========\n", sep = "")
    
    tryCatch({
        d1 <- readRDS(f1)
        d2 <- readRDS(f2)
        
        common_genes <- intersect(rownames(d1$counts), rownames(d2$counts))
        n_samples <- 50
        
        idx1 <- sample(ncol(d1$counts), min(n_samples, ncol(d1$counts)))
        idx2 <- sample(ncol(d2$counts), min(n_samples, ncol(d2$counts)))
        
        counts <- cbind(d1$counts[common_genes, idx1], d2$counts[common_genes, idx2])
        keep <- rowSums(counts > 10) >= 10
        counts <- counts[keep, ]
        
        group <- factor(c(rep(t1, length(idx1)), rep(t2, length(idx2))))
        
 cat(" : ", ncol(counts), " : ", nrow(counts), "\n", sep = "")
        
        gtex_results <- list()
        for (m_name in names(METHODS)) {
            cat("    Running ", m_name, "... ", sep = "")
            gtex_results[[m_name]] <- tryCatch({
                res <- METHODS[[m_name]](counts, group)
                cat("OK (", round(res$time[1], 1), "s)\n", sep = "")
                res
            }, error = function(e) {
                cat("FAILED\n")
                NULL
            })
        }
        
        gtex_results <- gtex_results[!sapply(gtex_results, is.null)]
        
        gtex_metrics <- rbindlist(lapply(gtex_results, compute_metrics), fill = TRUE)
        ds_name <- paste0("GTEx_", t1, "_vs_", t2)
        gtex_metrics[, dataset := ds_name]
        large_metrics[[ds_name]] <- gtex_metrics
        
 cat("\n \n")
        print(gtex_metrics[, .(method, n_sig_fdr5, n_sig_lfc1, time_sec)])
        
    }, error = function(e) {
 cat(" : ", conditionMessage(e), "\n", sep = "")
    })
}

tcga_tvn_cancers <- c("brca", "luad", "kirc", "lihc", "coad", "hnsc", "stad")

for (cancer in tcga_tvn_cancers) {
    tcga_file <- file.path(BULK_DIR, paste0("tcga_", cancer, "_tumor_vs_normal.rds"))
    if (!file.exists(tcga_file)) next

    cat("\n========== TCGA: ", toupper(cancer), " (tumor vs normal) ==========\n", sep = "")

    tryCatch({
        tcga <- readRDS(tcga_file)
        tcga_counts <- as.matrix(tcga$counts)
        tcga_group <- tcga$group

        keep <- rowSums(tcga_counts > 10) >= 5
        tcga_counts <- tcga_counts[keep, ]

        cat("  ", sum(tcga_group == "tumor"), " tumor vs ",
            sum(tcga_group == "normal"), " normal | Genes: ", nrow(tcga_counts), "\n", sep = "")

        tcga_results <- list()
        for (m_name in names(METHODS)) {
            cat("    Running ", m_name, "... ", sep = "")
            tcga_results[[m_name]] <- tryCatch({
                res <- METHODS[[m_name]](tcga_counts, tcga_group)
                cat("OK (", round(res$time[1], 1), "s)\n", sep = "")
                res
            }, error = function(e) {
                cat("FAILED\n")
                NULL
            })
        }

        tcga_results <- tcga_results[!sapply(tcga_results, is.null)]

        tcga_ds_metrics <- rbindlist(lapply(tcga_results, compute_metrics), fill = TRUE)
        ds_name <- paste0("TCGA_", toupper(cancer))
        tcga_ds_metrics[, dataset := ds_name]
        large_metrics[[ds_name]] <- tcga_ds_metrics

 cat("\n \n")
        print(tcga_ds_metrics[, .(method, n_sig_fdr5, n_sig_lfc1, time_sec)])

    }, error = function(e) {
 cat(" : ", conditionMessage(e), "\n", sep = "")
    })
}

pseudo_bulk_datasets <- list(
    kang18 = list(file = "kang18_pseudobulk.rds", label = "Kang18 PBMC stim vs ctrl"),
    segerstolpe = list(file = "segerstolpe_t2d_pseudobulk.rds", label = "Segerstolpe T2D vs Normal"),
    xin = list(file = "xin_t2d_pseudobulk.rds", label = "Xin T2D vs Healthy"),
    zilionis = list(file = "zilionis_lung_pseudobulk.rds", label = "Zilionis Lung tumor vs blood")
)

for (pb_name in names(pseudo_bulk_datasets)) {
    pb_info <- pseudo_bulk_datasets[[pb_name]]
    pb_file <- file.path(BULK_DIR, pb_info$file)
    if (!file.exists(pb_file)) next

    cat("\n========== Pseudo-bulk: ", pb_info$label, " ==========\n", sep = "")

    tryCatch({
        pb <- readRDS(pb_file)
        pb_counts <- as.matrix(pb$counts)
        pb_group <- pb$group

        keep <- rowSums(pb_counts > 5) >= 2
        pb_counts <- pb_counts[keep, ]

        cat("  Groups: ", paste(table(pb_group), collapse = " vs "),
            " | Genes: ", nrow(pb_counts), "\n", sep = "")

        pb_results <- list()
        for (m_name in names(METHODS)) {
            cat("    Running ", m_name, "... ", sep = "")
            pb_results[[m_name]] <- tryCatch({
                res <- METHODS[[m_name]](pb_counts, pb_group)
                cat("OK (", round(res$time[1], 1), "s)\n", sep = "")
                res
            }, error = function(e) {
                cat("FAILED\n")
                NULL
            })
        }

        pb_results <- pb_results[!sapply(pb_results, is.null)]

        pb_metrics <- rbindlist(lapply(pb_results, compute_metrics), fill = TRUE)
        ds_name <- paste0("PB_", pb_name)
        pb_metrics[, dataset := ds_name]
        large_metrics[[ds_name]] <- pb_metrics

 cat("\n \n")
        print(pb_metrics[, .(method, n_sig_fdr5, n_sig_lfc1, time_sec)])

    }, error = function(e) {
 cat(" : ", conditionMessage(e), "\n", sep = "")
    })
}

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 4. Null (, FDR )\n")
cat(strrep("=", 70), "\n")

null_results <- list()
n_null_sim <- 10

gtex_liver <- tryCatch(readRDS(file.path(BULK_DIR, "gtex_liver.rds")), error = function(e) NULL)

if (!is.null(gtex_liver)) {
 cat("\n GTEx Liver Null (", n_null_sim, " )\n", sep = "")
    
    for (i in seq_len(n_null_sim)) {
        cat("  Simulation ", i, "/", n_null_sim, "\n", sep = "")
        
        n_samples <- min(40, ncol(gtex_liver$counts))
        idx <- sample(ncol(gtex_liver$counts), n_samples)
        null_counts <- gtex_liver$counts[, idx]
        
        keep <- rowSums(null_counts > 10) >= 5
        null_counts <- null_counts[keep, ]
        
        half <- n_samples %/% 2
        null_group <- factor(c(rep("A", half), rep("B", n_samples - half)))
        
        for (m_name in names(METHODS)) {
            null_res <- tryCatch({
                res <- METHODS[[m_name]](null_counts, null_group)
                data.table(
                    method = m_name,
                    fdr01 = mean(res$padj < 0.01, na.rm = TRUE),
                    fdr05 = mean(res$padj < 0.05, na.rm = TRUE),
                    fdr10 = mean(res$padj < 0.10, na.rm = TRUE),
                    sim = i
                )
            }, error = function(e) NULL)
            
            if (!is.null(null_res)) {
                null_results[[length(null_results) + 1]] <- null_res
            }
        }
    }
}

null_metrics <- rbindlist(null_results, fill = TRUE)
null_summary <- null_metrics[, .(
    mean_FDR01 = round(mean(fdr01, na.rm = TRUE), 4),
    sd_FDR01 = round(sd(fdr01, na.rm = TRUE), 4),
    mean_FDR05 = round(mean(fdr05, na.rm = TRUE), 4),
    sd_FDR05 = round(sd(fdr05, na.rm = TRUE), 4),
    mean_FDR10 = round(mean(fdr10, na.rm = TRUE), 4),
    sd_FDR10 = round(sd(fdr10, na.rm = TRUE), 4)
), by = method]

cat("\nNull FDR (: 0.01/0.05/0.10)\n")
print(null_summary)

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 5. \n")
cat(strrep("=", 70), "\n")

all_classic <- rbindlist(classic_metrics, fill = TRUE)
all_sim <- rbindlist(sim_metrics, fill = TRUE)
all_large <- rbindlist(large_metrics, fill = TRUE)

saveRDS(list(
    classic = all_classic,
    simulation = all_sim,
    large_scale = all_large,
    null_fdr = null_summary
), file.path(RESULTS_DIR, "bulk_benchmark_results.rds"))

fwrite(all_classic, file.path(RESULTS_DIR, "bulk_classic_metrics.csv"))
fwrite(all_sim, file.path(RESULTS_DIR, "bulk_simulation_metrics.csv"))
fwrite(all_large, file.path(RESULTS_DIR, "bulk_large_scale_metrics.csv"))
fwrite(null_summary, file.path(RESULTS_DIR, "bulk_null_fdr.csv"))

# =============================================================================
# =============================================================================

cat("\n...\n")

if (nrow(all_sim) > 0 && "F1_fdr5" %in% names(all_sim)) {
    p_sim <- ggplot(all_sim, aes(x = dataset, y = F1_fdr5, fill = method)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        labs(title = "Simulation: F1 Score (FDR < 0.05)", x = "", y = "F1") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
    p_sim <- ggplot() + theme_void() + labs(title = "No simulation data")
}

if (nrow(all_classic) > 0) {
    p_classic <- ggplot(all_classic, aes(x = dataset, y = n_sig_fdr5, fill = method)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        labs(title = "Classic Datasets: DEGs (FDR < 0.05)", x = "", y = "DEGs") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
    p_classic <- ggplot() + theme_void() + labs(title = "No classic data")
}

time_data <- rbind(
    all_classic[, .(method, time_sec, type = "Classic")],
    all_large[, .(method, time_sec, type = "Large-scale")],
    fill = TRUE
)

if (nrow(time_data) > 0) {
    time_summary <- time_data[, .(mean_time = mean(time_sec, na.rm = TRUE)), by = .(method, type)]
    
    p_time <- ggplot(time_summary, aes(x = method, y = mean_time, fill = type)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        labs(title = "Running Time Comparison", x = "", y = "Time (seconds)") +
        theme_minimal()
} else {
    p_time <- ggplot() + theme_void() + labs(title = "No time data")
}

# 4. Null FDR
if (nrow(null_summary) > 0) {
    null_long <- melt(null_summary, id.vars = "method", 
                      measure.vars = c("mean_FDR01", "mean_FDR05", "mean_FDR10"),
                      variable.name = "threshold", value.name = "FDR")
    null_long[, target := c(0.01, 0.05, 0.10)[match(threshold, c("mean_FDR01", "mean_FDR05", "mean_FDR10"))]]
    
    p_null <- ggplot(null_long, aes(x = method, y = FDR, fill = threshold)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        geom_hline(yintercept = c(0.01, 0.05, 0.10), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(title = "Null FDR Control", x = "", y = "Actual FDR") +
        theme_minimal()
} else {
    p_null <- ggplot() + theme_void() + labs(title = "No null FDR data")
}

combined <- (p_sim | p_classic) / (p_time | p_null) +
    plot_annotation(
        title = "Bulk RNA-seq Benchmark Results",
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

ggsave(file.path(RESULTS_DIR, "bulk_benchmark_results.pdf"), combined, width = 14, height = 10)
ggsave(file.path(RESULTS_DIR, "bulk_benchmark_results.png"), combined, width = 14, height = 10, dpi = 300)

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" Bulk RNA-seq \n")
cat(strrep("=", 70), "\n")

cat("\n\n")
if (nrow(all_classic) > 0) {
    classic_summary <- all_classic[, .(
        mean_DEGs = round(mean(n_sig_fdr5, na.rm = TRUE)),
        mean_time = round(mean(time_sec, na.rm = TRUE), 2)
    ), by = method]
    print(classic_summary[order(-mean_DEGs)])
}

cat("\n (FDR/Power)\n")
if (nrow(all_sim) > 0 && "F1_fdr5" %in% names(all_sim)) {
    sim_summary <- all_sim[, .(
        mean_F1 = round(mean(F1_fdr5, na.rm = TRUE), 4),
        mean_recall = round(mean(recall_fdr5, na.rm = TRUE), 4),
        mean_FDR = round(mean(actual_FDR_fdr5, na.rm = TRUE), 4)
    ), by = method]
    print(sim_summary[order(-mean_F1)])
}

cat("\n\n")
if (nrow(all_large) > 0) {
    large_summary <- all_large[, .(
        mean_DEGs = round(mean(n_sig_fdr5, na.rm = TRUE)),
        mean_time = round(mean(time_sec, na.rm = TRUE), 2)
    ), by = method]
    print(large_summary[order(-mean_DEGs)])
}

cat("\nNull FDR \n")
print(null_summary)

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 5. (DV) \n")
cat(strrep("=", 70), "\n")

dv_metrics <- list()

generate_dv_data <- function(n_genes = 5000, n_samples = 10, pDV = 0.1, var_ratio = 3) {
    n_per_group <- n_samples / 2
    base_mean <- exp(rnorm(n_genes, mean = 4, sd = 2))
    base_mean <- pmax(base_mean, 1)
    base_var <- base_mean * runif(n_genes, 0.5, 2)
    
    counts <- matrix(0, n_genes, n_samples)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Sample", seq_len(n_samples))
    
    n_dv <- round(n_genes * pDV)
    dv_idx <- sample(n_genes, n_dv)
    dv_genes <- rownames(counts)[dv_idx]
    
    for (i in seq_len(n_genes)) {
        mu <- base_mean[i]
        var_ctrl <- base_var[i]
        size_ctrl <- mu^2 / (var_ctrl - mu)
        size_ctrl <- pmax(size_ctrl, 0.1)
        counts[i, 1:n_per_group] <- rnbinom(n_per_group, size = size_ctrl, mu = mu)
        
        if (i %in% dv_idx) {
            var_treat <- var_ctrl * var_ratio
        } else {
            var_treat <- var_ctrl
        }
        size_treat <- mu^2 / (var_treat - mu)
        size_treat <- pmax(size_treat, 0.1)
        counts[i, (n_per_group + 1):n_samples] <- rnbinom(n_per_group, size = size_treat, mu = mu)
    }
    
    group <- factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group)))
    list(counts = counts, group = group, dv_genes = dv_genes)
}

var_ratios <- c(2, 3, 5)

for (vr in var_ratios) {
 cat("\n--- DV : = ", vr, "x ---\n", sep = "")
    
    sim_dv <- generate_dv_data(n_genes = 5000, n_samples = 10, pDV = 0.1, var_ratio = vr)
    
    keep <- rowSums(sim_dv$counts > 5) >= 2
    counts_dv <- sim_dv$counts[keep, ]
    truth_dv <- intersect(sim_dv$dv_genes, rownames(counts_dv))
    
 cat(" : ", nrow(counts_dv), " | DV: ", length(truth_dv), "\n", sep = "")
    
    for (m_name in names(DV_METHODS)) {
        cat("    ", m_name, "... ", sep = "")
        res <- tryCatch({
            dt <- DV_METHODS[[m_name]](counts_dv, sim_dv$group)
            metrics <- compute_dv_metrics(dt, truth_dv)
            metrics[, var_ratio := vr]
            cat("OK\n")
            metrics
        }, error = function(e) {
            cat("FAILED: ", conditionMessage(e), "\n")
            NULL
        })
        if (!is.null(res)) dv_metrics[[length(dv_metrics) + 1]] <- res
    }
}

dv_dt <- rbindlist(dv_metrics, fill = TRUE)
cat("\nDV (F1 Score)\n")
print(dcast(dv_dt, method ~ var_ratio, value.var = "F1"))

fwrite(dv_dt, file.path(RESULTS_DIR, "dv_benchmark_results.csv"))

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 6. \n")
cat(strrep("=", 70), "\n")

batch_benchmark_results <- list()

generate_batch_benchmark_data <- function(n_genes = 5000, n_samples = 24, n_batch = 3,
                                           pDE = 0.1, effect_size = 2, batch_effect = 2) {
    n_per_group <- n_samples / 2
    base_mean <- exp(rnorm(n_genes, mean = 4, sd = 2))
    base_mean <- pmax(base_mean, 1)
    alpha <- runif(n_genes, 0.05, 0.3)
    
    counts <- matrix(0, n_genes, n_samples)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Sample", seq_len(n_samples))
    
    n_de <- round(n_genes * pDE)
    de_idx <- sample(n_genes, n_de)
    de_genes <- rownames(counts)[de_idx]
    
    batch <- factor(rep(1:n_batch, length.out = n_samples))
    condition <- factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group)))
    
    batch_effects <- matrix(1, n_genes, n_batch)
    for (b in 2:n_batch) {
        batch_effects[, b] <- exp(rnorm(n_genes, mean = 0, sd = log(batch_effect)))
    }
    
    for (j in seq_len(n_samples)) {
        b <- as.integer(batch[j])
        is_treat <- condition[j] == "treat"
        for (i in seq_len(n_genes)) {
            mu <- base_mean[i] * batch_effects[i, b]
            if (is_treat && (i %in% de_idx)) mu <- mu * effect_size
            size <- 1 / alpha[i]
            counts[i, j] <- rnbinom(1, size = size, mu = mu)
        }
    }
    
    list(counts = counts, condition = condition, batch = batch, de_genes = de_genes)
}

batch_strengths <- c(1.5, 2, 3)

for (be in batch_strengths) {
 cat("\n--- : ", be, "x ---\n", sep = "")
    
    sim_batch <- generate_batch_benchmark_data(
        n_genes = 5000, n_samples = 24, n_batch = 3,
        pDE = 0.1, effect_size = 2, batch_effect = be
    )
    
    keep <- rowSums(sim_batch$counts > 5) >= 2
    counts_batch <- sim_batch$counts[keep, ]
    truth_de <- intersect(sim_batch$de_genes, rownames(counts_batch))
    
 cat(" : ", nrow(counts_batch), " | DE: ", length(truth_de), "\n", sep = "")
    
    cat("    DESeq2 (no batch)... ")
    res_baseline <- tryCatch({
        col_data <- data.frame(condition = sim_batch$condition, row.names = colnames(counts_batch))
        dds <- DESeqDataSetFromMatrix(counts_batch, col_data, ~condition)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds)
        sig <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
        TP <- sum(sig %in% truth_de)
        FP <- sum(!(sig %in% truth_de))
        f1 <- 2 * TP / max(1, 2 * TP + FP + length(truth_de) - TP)
        cat("OK\n")
        data.table(method = "DESeq2_no_batch", batch_effect = be, TP = TP, FP = FP, F1 = round(f1, 4))
    }, error = function(e) { cat("FAILED\n"); NULL })
    if (!is.null(res_baseline)) batch_benchmark_results[[length(batch_benchmark_results) + 1]] <- res_baseline
    
    # DESeq2 with batch in design
    cat("    DESeq2 (with batch)... ")
    res_deseq2_batch <- tryCatch({
        col_data <- data.frame(condition = sim_batch$condition, batch = sim_batch$batch, row.names = colnames(counts_batch))
        dds <- DESeqDataSetFromMatrix(counts_batch, col_data, ~ batch + condition)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds)
        sig <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
        TP <- sum(sig %in% truth_de)
        FP <- sum(!(sig %in% truth_de))
        f1 <- 2 * TP / max(1, 2 * TP + FP + length(truth_de) - TP)
        cat("OK\n")
        data.table(method = "DESeq2_with_batch", batch_effect = be, TP = TP, FP = FP, F1 = round(f1, 4))
    }, error = function(e) { cat("FAILED\n"); NULL })
    if (!is.null(res_deseq2_batch)) batch_benchmark_results[[length(batch_benchmark_results) + 1]] <- res_deseq2_batch
    
    # ComBat-seq
    cat("    ComBat-seq... ")
    res_combat <- tryCatch({
        adj <- run_combat_seq(counts_batch, sim_batch$condition, sim_batch$batch)
        col_data <- data.frame(condition = sim_batch$condition, row.names = colnames(adj$counts))
        dds <- DESeqDataSetFromMatrix(round(adj$counts), col_data, ~condition)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds)
        sig <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
        TP <- sum(sig %in% truth_de)
        FP <- sum(!(sig %in% truth_de))
        f1 <- 2 * TP / max(1, 2 * TP + FP + length(truth_de) - TP)
        cat("OK\n")
        data.table(method = "ComBat-seq", batch_effect = be, TP = TP, FP = FP, F1 = round(f1, 4))
    }, error = function(e) { cat("FAILED: ", conditionMessage(e), "\n"); NULL })
    if (!is.null(res_combat)) batch_benchmark_results[[length(batch_benchmark_results) + 1]] <- res_combat
    
    # ComBat-ref
    cat("    ComBat-ref... ")
    res_combat_ref <- tryCatch({
        adj <- run_combat_ref(counts_batch, sim_batch$condition, sim_batch$batch)
        col_data <- data.frame(condition = sim_batch$condition, row.names = colnames(adj$counts))
        dds <- DESeqDataSetFromMatrix(round(pmax(adj$counts, 0)), col_data, ~condition)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds)
        sig <- rownames(res)[res$padj < 0.05 & !is.na(res$padj)]
        TP <- sum(sig %in% truth_de)
        FP <- sum(!(sig %in% truth_de))
        f1 <- 2 * TP / max(1, 2 * TP + FP + length(truth_de) - TP)
        cat("OK\n")
        data.table(method = "ComBat-ref", batch_effect = be, TP = TP, FP = FP, F1 = round(f1, 4))
    }, error = function(e) { cat("FAILED: ", conditionMessage(e), "\n"); NULL })
    if (!is.null(res_combat_ref)) batch_benchmark_results[[length(batch_benchmark_results) + 1]] <- res_combat_ref
    
    # (SGCB implicit batch correction removed - tested separately)
}

batch_benchmark_dt <- rbindlist(batch_benchmark_results, fill = TRUE)
cat("\n (F1 Score)\n")
print(dcast(batch_benchmark_dt, method ~ batch_effect, value.var = "F1"))

fwrite(batch_benchmark_dt, file.path(RESULTS_DIR, "batch_benchmark_results.csv"))

cat("\n: ", RESULTS_DIR, "\n", sep = "")
cat(strrep("=", 70), "\n")

# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 7. GG (Differential Distribution)\n")
cat(" : /, DE \n")
cat(strrep("=", 70), "\n")

dd_sim_metrics <- list()

generate_gg_shape_data <- function(n_genes = 5000, n_samples = 20, pDD = 0.1,
                                    alpha_shift = 0.5, gamma_shift = 0.3) {
    n_per_group <- n_samples / 2

    base_mean <- exp(rnorm(n_genes, mean = 4, sd = 2))
    base_mean <- pmax(base_mean, 1)

    alpha_base <- runif(n_genes, 1, 5)
    gamma_base <- runif(n_genes, 0.5, 2)
    beta_base <- base_mean / (gamma(alpha_base + 1 / gamma_base) / gamma(alpha_base))
    beta_base <- pmax(beta_base, 0.01)

    n_dd <- round(n_genes * pDD)
    dd_idx <- sample(n_genes, n_dd)
    dd_genes <- paste0("Gene", dd_idx)

    counts <- matrix(0L, n_genes, n_samples)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Sample", seq_len(n_samples))

    sample_gg <- function(n, alpha, beta, gamma_param) {
        y <- rgamma(n, shape = alpha, rate = 1)
        x <- beta * y^(1 / gamma_param)
        pmax(round(x), 0L)
    }

    sapply(seq_len(n_genes), function(i) {
        a <- alpha_base[i]; b <- beta_base[i]; g <- gamma_base[i]

        counts[i, 1:n_per_group] <<- sample_gg(n_per_group, a, b, g)

        a_treat <- a; g_treat <- g
        b_treat <- b

        is_dd <- i %in% dd_idx
        a_treat <- a + is_dd * alpha_shift
        g_treat <- g + is_dd * gamma_shift
        mean_orig <- b * gamma(a + 1/g) / gamma(a)
        mean_new <- gamma(a_treat + 1/g_treat) / gamma(a_treat)
        b_treat <- mean_orig / pmax(mean_new, 1e-6)

        counts[i, (n_per_group + 1):n_samples] <<- sample_gg(n_per_group, a_treat, b_treat, g_treat)
        NULL
    })

    group <- factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group)))
    list(counts = counts, group = group, dd_genes = dd_genes,
         alpha_shift = alpha_shift, gamma_shift = gamma_shift)
}

shape_shifts <- list(
    list(alpha = 0.3, gamma = 0.2),
    list(alpha = 0.5, gamma = 0.3),
    list(alpha = 1.0, gamma = 0.5)
)

cat("\n--- DE GG ---\n")

dd_de_results <- list()
dd_dv_results <- list()

for (shift in shape_shifts) {
    shift_label <- paste0("a", shift$alpha, "_g", shift$gamma)
 cat("\n--- : =", shift$alpha, ", =", shift$gamma, " ---\n", sep = "")

    sim_dd <- generate_gg_shape_data(
        n_genes = 5000, n_samples = 20, pDD = 0.1,
        alpha_shift = shift$alpha, gamma_shift = shift$gamma
    )

    keep <- rowSums(sim_dd$counts > 5) >= 2
    counts_dd <- sim_dd$counts[keep, ]
    truth_dd <- intersect(sim_dd$dd_genes, rownames(counts_dd))

 cat(" : ", nrow(counts_dd), " | DD: ", length(truth_dd), "\n", sep = "")

    for (m_name in names(METHODS)) {
        cat("    [DE] ", m_name, "... ", sep = "")
        res <- tryCatch({
            dt <- METHODS[[m_name]](counts_dd, sim_dd$group)
            sig <- dt$gene[dt$padj < 0.05 & !is.na(dt$padj)]
            TP <- sum(sig %in% truth_dd)
            FP <- sum(!(sig %in% truth_dd))
            recall <- TP / max(1, length(truth_dd))
            precision <- TP / max(1, TP + FP)
            f1 <- 2 * precision * recall / max(0.001, precision + recall)
            cat("recall=", round(recall, 3), "\n")
            data.table(method = m_name, type = "DE", shift = shift_label,
                       TP = TP, FP = FP, recall = round(recall, 4),
                       precision = round(precision, 4), F1 = round(f1, 4))
        }, error = function(e) { cat("FAILED\n"); NULL })
        dd_de_results[[length(dd_de_results) + 1]] <- res
    }

    for (m_name in names(DV_METHODS)) {
        cat("    [DV] ", m_name, "... ", sep = "")
        res <- tryCatch({
            dt <- DV_METHODS[[m_name]](counts_dd, sim_dd$group)
            sig <- dt$gene[dt$padj_dv < 0.05 & !is.na(dt$padj_dv)]
            TP <- sum(sig %in% truth_dd)
            FP <- sum(!(sig %in% truth_dd))
            recall <- TP / max(1, length(truth_dd))
            precision <- TP / max(1, TP + FP)
            f1 <- 2 * precision * recall / max(0.001, precision + recall)
            cat("recall=", round(recall, 3), "\n")
            data.table(method = m_name, type = "DV", shift = shift_label,
                       TP = TP, FP = FP, recall = round(recall, 4),
                       precision = round(precision, 4), F1 = round(f1, 4))
        }, error = function(e) { cat("FAILED\n"); NULL })
        dd_dv_results[[length(dd_dv_results) + 1]] <- res
    }
}

dd_all <- rbindlist(c(dd_de_results, dd_dv_results), fill = TRUE)

cat("\nGG : DE vs DV (F1 Score)\n")
print(dcast(dd_all, method + type ~ shift, value.var = "F1"))
cat("\nGG : Recall ()\n")
print(dcast(dd_all, method + type ~ shift, value.var = "recall"))

fwrite(dd_all, file.path(RESULTS_DIR, "gg_shape_dd_benchmark.csv"))

# =============================================================================
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 8. GTEx Brain DV ()\n")
cat(strrep("=", 70), "\n")

gtex_brain_file <- file.path(BULK_DIR, "gtex_brain_cerebellum_vs_cortex.rds")

gtex_brain_dv_results <- list()

tryCatch({
    gtex_cb_fc <- readRDS(gtex_brain_file)

    n_sub <- 50
    ctrl_idx <- sample(which(gtex_cb_fc$group == "Cerebellum"), min(n_sub, sum(gtex_cb_fc$group == "Cerebellum")))
    treat_idx <- sample(which(gtex_cb_fc$group == "FrontalCortex"), min(n_sub, sum(gtex_cb_fc$group == "FrontalCortex")))

    sel <- c(ctrl_idx, treat_idx)
    counts_brain <- gtex_cb_fc$counts[, sel]
    group_brain <- factor(c(rep("Cerebellum", length(ctrl_idx)), rep("FrontalCortex", length(treat_idx))))

    keep <- rowSums(counts_brain > 10) >= 10
    counts_brain <- counts_brain[keep, ]

    cat("  Cerebellum: ", length(ctrl_idx), " | FrontalCortex: ", length(treat_idx),
        " | Genes: ", nrow(counts_brain), "\n", sep = "")

    for (m_name in names(DV_METHODS)) {
        cat("    ", m_name, "... ", sep = "")
        res <- tryCatch({
            dt <- DV_METHODS[[m_name]](counts_brain, group_brain)
            n_sig <- sum(dt$padj_dv < 0.05, na.rm = TRUE)
            cat(n_sig, " DV genes\n")
            dt[, dataset := "GTEx_Brain_CB_vs_FC"]
            dt
        }, error = function(e) { cat("FAILED: ", conditionMessage(e), "\n"); NULL })
        gtex_brain_dv_results[[length(gtex_brain_dv_results) + 1]] <- res
    }

    gtex_brain_dv_dt <- rbindlist(gtex_brain_dv_results, fill = TRUE)

    dv_sig_by_method <- split(gtex_brain_dv_dt[padj_dv < 0.05]$gene, gtex_brain_dv_dt[padj_dv < 0.05]$method)
    n_m <- length(dv_sig_by_method)
    jaccard_mat <- matrix(0, n_m, n_m, dimnames = list(names(dv_sig_by_method), names(dv_sig_by_method)))
    sapply(seq_len(n_m), function(i) {
        sapply(seq_len(n_m), function(j) {
            g1 <- dv_sig_by_method[[i]]; g2 <- dv_sig_by_method[[j]]
            jaccard_mat[i, j] <<- length(intersect(g1, g2)) / max(1, length(union(g1, g2)))
            NULL
        })
        NULL
    })

 cat("\nGTEx Brain DV \n")
    dv_counts <- gtex_brain_dv_dt[, .(n_dv = sum(padj_dv < 0.05, na.rm = TRUE)), by = method]
    print(dv_counts)

 cat("\nDV Jaccard \n")
    print(round(jaccard_mat, 3))

    fwrite(dv_counts, file.path(RESULTS_DIR, "gtex_brain_dv_counts.csv"))

}, error = function(e) {
 cat(" GTEx Brain DV : ", conditionMessage(e), "\n", sep = "")
})

# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat(" 9. Tabula Muris Senis DV \n")
cat(" : , DV \n")
cat(strrep("=", 70), "\n")

tms_file <- file.path(BULK_DIR, "tms_bulk.rds")

tms_dv_results <- list()

tryCatch({
    tms <- readRDS(tms_file)

    cat("  TMS bulk: ", nrow(tms$counts), " genes x ", ncol(tms$counts), " samples\n", sep = "")

    age_col <- grep("age", colnames(tms$metadata), ignore.case = TRUE, value = TRUE)[1]
    organ_col <- grep("organ|tissue|source", colnames(tms$metadata), ignore.case = TRUE, value = TRUE)[1]

 cat(" :", age_col, " | :", organ_col, "\n")
 cat(" :\n")
    print(table(tms$metadata[[age_col]]))
 cat(" :\n")
    print(table(tms$metadata[[organ_col]]))

    organs <- names(sort(table(tms$metadata[[organ_col]]), decreasing = TRUE))
    target_organ <- organs[1]

    organ_mask <- tms$metadata[[organ_col]] == target_organ
    ages <- tms$metadata[[age_col]][organ_mask]
    unique_ages <- sort(unique(ages))

 cat("\n :", target_organ, "\n")
 cat(" :", paste(unique_ages, collapse=", "), "\n")

    young_age <- unique_ages[1]
    old_age <- unique_ages[length(unique_ages)]

    young_idx <- which(organ_mask & tms$metadata[[age_col]] == young_age)
    old_idx <- which(organ_mask & tms$metadata[[age_col]] == old_age)

    cat("  Young (", young_age, "): ", length(young_idx), " | Old (", old_age, "): ", length(old_idx), "\n", sep = "")

    sel_idx <- c(young_idx, old_idx)
    counts_tms <- as.matrix(tms$counts[, sel_idx])
    group_tms <- factor(c(rep("young", length(young_idx)), rep("old", length(old_idx))))

    keep <- rowSums(counts_tms > 5) >= 3
    counts_tms <- counts_tms[keep, ]

 cat(" : ", nrow(counts_tms), "\n", sep = "")

    for (m_name in names(DV_METHODS)) {
        cat("    ", m_name, "... ", sep = "")
        res <- tryCatch({
            dt <- DV_METHODS[[m_name]](counts_tms, group_tms)
            n_sig <- sum(dt$padj_dv < 0.05, na.rm = TRUE)
            cat(n_sig, " DV genes\n")
            dt[, dataset := paste0("TMS_", target_organ, "_", young_age, "_vs_", old_age)]
            dt
        }, error = function(e) { cat("FAILED: ", conditionMessage(e), "\n"); NULL })
        tms_dv_results[[length(tms_dv_results) + 1]] <- res
    }

    tms_dv_dt <- rbindlist(tms_dv_results, fill = TRUE)

 cat("\nTMS DV (", target_organ, ": ", young_age, " vs ", old_age, ")\n", sep = "")
    tms_dv_counts <- tms_dv_dt[, .(n_dv = sum(padj_dv < 0.05, na.rm = TRUE)), by = method]
    print(tms_dv_counts)

    fwrite(tms_dv_counts, file.path(RESULTS_DIR, "tms_aging_dv_counts.csv"))

}, error = function(e) {
 cat(" TMS DV : ", conditionMessage(e), "\n", sep = "")
})

cat("\n", strrep("=", 70), "\n")
cat(" Benchmark \n")
cat(strrep("=", 70), "\n")
