# =============================================================================
#   method:  DESeq2 | edgeR | limma | glmGamPoi | Seurat |
#            MAST | scDD | DEsingle | distinct | BPSC |
#            muscat | aggregateBioVar | glmmTMB
#   dataset: kang18 | segerstolpe | xin | zilionis | jakel | muscat_sim
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
THREADS <- as.integer(Sys.getenv("BENCHMARK_THREADS", "1"))
THREADS <- ifelse(is.na(THREADS) | THREADS < 1L, 1L, THREADS)
BPSC_USE_PARALLEL <- FALSE
DIPHISEQ_CAPTURE_OUTPUT <- TRUE
set.seed(12345)
Sys.setenv(
  OMP_NUM_THREADS = THREADS,
  OPENBLAS_NUM_THREADS = THREADS,
  MKL_NUM_THREADS = THREADS,
  VECLIB_MAXIMUM_THREADS = THREADS,
  NUMEXPR_NUM_THREADS = THREADS,
  BLIS_NUM_THREADS = THREADS,
  R_DATATABLE_NUM_THREADS = THREADS
)
data.table::setDTthreads(THREADS)
bp_param <- BiocParallel::SnowParam(workers = THREADS, type = "SOCK", progressbar = FALSE)
BiocParallel::register(bp_param)
options(glmmTMB.cores = THREADS, glmmTMB.autopar = TRUE)

args      <- commandArgs(trailingOnly = TRUE)
METHOD    <- args[1]
DATASET   <- args[2]
NULL_SEED <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

if (identical(METHOD, "dv_DiPhiSeq")) stop("dv_DiPhiSeq has been retired from benchmark_v2")

BASE     <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk")
PB_DIR   <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/05_pseudobulk")
OUT_DIR  <- file.path(BASE, "other_methods/results")
IS_NULL  <- !is.na(NULL_SEED)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
library(data.table)

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== scRNA Other DE: %s / %s %s=====\n",
            Sys.time(), METHOD, DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

switch(as.character(file.exists(OUT_CSV)),
  "TRUE" = { cat("  [SKIP] Already exists:", OUT_CSV, "\n"); q("no", status = 0) },
  "FALSE" = NULL)

pb_map <- c(
  kang18       = "kang18_pseudobulk.rds",
  segerstolpe  = "segerstolpe_t2d_pseudobulk.rds",
  xin          = "xin_t2d_pseudobulk.rds",
  zilionis     = "zilionis_lung_pseudobulk.rds",
  jakel        = "jakel_pseudobulk.rds",
  muscat_sim   = "muscat_sim_pseudobulk.rds"
)

ds     <- readRDS(file.path(PB_DIR, pb_map[DATASET]))
counts <- as.matrix(ds$counts)
storage.mode(counts) <- "integer"
rownames(counts) <- make.unique(rownames(counts))
group  <- factor(gsub("[^A-Za-z0-9._]", "_", as.character(ds$group)))
sim_truth <- ds$truth

# Null FDR
switch(as.character(IS_NULL),
  "TRUE" = {
    set.seed(NULL_SEED)
    group <- factor(sample(as.character(group)))
    cat("  [NULL FDR] Labels permuted\n")
  }
)

cat(sprintf("  Data: %d genes x %d samples, groups: %s\n",
            nrow(counts), ncol(counts),
            paste(paste0(levels(group), "=", table(group)), collapse=", ")))

# ===== DE =====
t0 <- proc.time()

res_df <- switch(METHOD,

  # ---------- DESeq2 ----------
  "DESeq2" = {
    library(DESeq2)
    coldata <- data.frame(group = group, row.names = colnames(counts))
    dds <- DESeqDataSetFromMatrix(counts, coldata, ~group)
    dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = bp_param)
    rr  <- results(dds, parallel = TRUE, BPPARAM = bp_param)
    data.table(gene = rownames(rr), pvalue = rr$pvalue, padj = rr$padj,
               log2FoldChange = rr$log2FoldChange)
  },

  # ---------- edgeR ----------
  "edgeR" = {
    library(edgeR)
    design <- model.matrix(~group)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    tt  <- glmQLFTest(fit, coef = 2)
    rr  <- topTags(tt, n = Inf, sort.by = "none")$table
    data.table(gene = rownames(rr), pvalue = rr$PValue, padj = rr$FDR,
               log2FoldChange = rr$logFC)
  },

  # ---------- limma-voom ----------
  "limma" = {
    library(limma)
    library(edgeR)
    design <- model.matrix(~group)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    rr  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  },

  # ---------- glmGamPoi ----------
  "glmGamPoi" = {
    library(glmGamPoi)
    fit <- glm_gp(counts, design = ~group,
                   col_data = data.frame(group = group,
                                         row.names = colnames(counts)))
    rr  <- test_de(fit, contrast = colnames(fit$Beta)[2])
    data.table(gene = rr$name, pvalue = rr$pval, padj = rr$adj_pval,
               log2FoldChange = rr$lfc)
  },

  # ---------- Seurat (Wilcoxon) ----------
  "Seurat" = {
    library(Seurat)
    future::plan(future::multisession, workers = THREADS)
    obj <- CreateSeuratObject(counts = counts)
    obj$group <- group
    Idents(obj) <- "group"
    obj <- NormalizeData(obj, verbose = FALSE)
    markers <- FindMarkers(obj, ident.1 = levels(group)[2],
                  ident.2 = levels(group)[1],
                  test.use = "wilcox",
                  logfc.threshold = 0, min.pct = 0,
                  features = rownames(counts),
                  verbose = FALSE)
    marker_dt <- data.table(gene = rownames(markers), pvalue = markers$p_val,
                 padj = markers$p_val_adj, log2FoldChange = markers$avg_log2FC)
    marker_dt <- marker_dt[match(rownames(counts), gene)]
    marker_dt
  },

  # ---------- MAST (hurdle model on log-CPM) ----------
  "MAST" = {
    library(MAST)
    lib_size <- colSums(counts)
    cpm_mat  <- sweep(counts, 2, lib_size, "/") * 1e6
    log_cpm  <- log2(cpm_mat + 1)
    cdat <- data.frame(wellKey = colnames(counts), group = group)
    fdat <- data.frame(primerid = rownames(counts))
    sca  <- FromMatrix(log_cpm, cdat, fdat)
    contrast_name <- paste0("group", levels(group)[2])
    zlm_fit <- zlm(~ group, sca, parallel = TRUE)
    summ <- summary(zlm_fit, doLRT = contrast_name, parallel = TRUE)
    rr <- summ$datatable
    hurdle <- rr[rr$component == "H" & rr$contrast == contrast_name, ]
    coef_c <- rr[rr$component == "C" & rr$contrast == contrast_name & rr$metric == "coef", ]
    hurdle_dt <- data.table(gene = hurdle$primerid, pvalue = hurdle[["Pr(>Chisq)"]],
                       log2FoldChange = 0)
    hurdle_dt[is.na(pvalue), pvalue := 1]
    idx_match <- match(hurdle_dt$gene, coef_c$primerid)
    hurdle_dt[!is.na(idx_match), log2FoldChange := coef_c$value[idx_match[!is.na(idx_match)]]]
    hurdle_dt <- hurdle_dt[match(rownames(counts), gene)]
    hurdle_dt[, padj := p.adjust(pvalue, "BH")]
    hurdle_dt
  },

  # ---------- scDD ----------
  "scDD" = {
    library(scDD)
    library(SingleCellExperiment)
    lib_size <- colSums(counts)
    norm_mat <- sweep(counts, 2, lib_size / median(lib_size), "/")
    sce <- SingleCellExperiment(
      assays = list(normcounts = norm_mat),
      colData = data.frame(condition = group, row.names = colnames(counts))
    )
    keep <- rowSums(norm_mat > 0) >= 2
    sce <- sce[keep, ]
    res <- scDD(sce, prior_param = list(alpha = 0.01, mu0 = 0, s0 = 0.01,
                                         a0 = 0.01, b0 = 0.01),
                testZeroes = FALSE, condition = "condition",
                parallelBy = "Genes",
                param = bp_param)
    rr      <- results(res)
    pv_col  <- intersect(c("combined.pvalue",     "nonzero.pvalue"),     names(rr))[1]
    adj_col <- intersect(c("combined.pvalue.adj", "nonzero.pvalue.adj"), names(rr))[1]
    data.table(gene = rr$gene, pvalue = rr[[pv_col]],
               padj = rr[[adj_col]],
               log2FoldChange = ifelse("log2FC" %in% names(rr), rr$log2FC, 0))
  },

  # ---------- DEsingle (ZINB) ----------
  "DEsingle" = {
    library(DEsingle)
    rr <- DEsingle(counts = counts, group = group, parallel = TRUE,
                   BPPARAM = bp_param)
    data.table(gene = rownames(rr), pvalue = rr$pvalue, padj = rr$pvalue.adj.FDR,
               log2FoldChange = log2(pmax(rr$foldChange, 1e-8)))
  },

  # ---------- distinct ----------
  "distinct" = {
    library(distinct)
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
      assays = list(counts = counts),
      colData = data.frame(
        sample_id  = colnames(counts),
        cluster_id = rep("bulk", ncol(counts)),
        group_id   = group,
        row.names  = colnames(counts)
      )
    )
    lib_sz <- colSums(counts)
    logcounts(sce) <- log2(sweep(counts, 2, lib_sz / median(lib_sz), "/") + 1)
    design_mat <- model.matrix(~group_id, data = as.data.frame(SummarizedExperiment::colData(sce)))
    rownames(design_mat) <- SummarizedExperiment::colData(sce)$sample_id
    res <- distinct_test(sce,
                         name_cluster = "cluster_id",
                         name_sample = "sample_id",
                         design = design_mat,
                         column_to_test = 2,
                         min_non_zero_cells = 1, n_cores = THREADS)
    rr <- as.data.frame(res)
    data.table(gene = rr$gene, pvalue = rr$p_val, padj = rr$p_adj.loc,
               log2FoldChange = 0)
  },

  # ---------- BPSC ----------
  "BPSC" = {
    library(BPSC)
    library(matrixStats)
    design <- model.matrix(~ group)
    ctrl_idx  <- which(group == levels(group)[1])
    treat_idx <- which(group == levels(group)[2])
    log_counts <- log2(counts + 1)
    mean_ctrl  <- rowMeans(log_counts[, ctrl_idx,  drop = FALSE])
    mean_treat <- rowMeans(log_counts[, treat_idx, drop = FALSE])
    var_ctrl   <- rowVars(log_counts[, ctrl_idx,  drop = FALSE])
    var_treat  <- rowVars(log_counts[, treat_idx, drop = FALSE])
    t_stderr   <- sqrt(var_ctrl / length(ctrl_idx) + var_treat / length(treat_idx))
    t_limit    <- 10 * .Machine$double.eps * pmax(abs(mean_ctrl), abs(mean_treat))
    safe_idx   <- !(is.na(t_stderr) | t_stderr < t_limit)
    Sys.setenv(
      OMP_NUM_THREADS = 1L,
      OPENBLAS_NUM_THREADS = 1L,
      MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L,
      NUMEXPR_NUM_THREADS = 1L,
      BLIS_NUM_THREADS = 1L,
      R_DATATABLE_NUM_THREADS = 1L
    )
    pvals <- rep(1, nrow(counts))
    switch(as.character(sum(safe_idx) > 0L),
      "TRUE" = {
        res <- BPglm(data = counts[safe_idx, , drop = FALSE], controlIds = ctrl_idx, design = design,
                     coef = 2, estIntPar = FALSE, useParallel = BPSC_USE_PARALLEL)
        pvals[safe_idx] <- res$PVAL
        NULL
      },
      "FALSE" = NULL
    )
    pvals[is.na(pvals)] <- 1
    lfc <- mean_treat - mean_ctrl
    data.table(gene = rownames(counts), pvalue = pvals,
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = lfc)
  },

  # ---------- muscat (edgeR-based DS) ----------
  "muscat" = {
    library(muscat)
    library(edgeR)
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(
      assays = list(counts = counts),
      colData = data.frame(
        sample_id  = colnames(counts),
        cluster_id = rep("bulk", ncol(counts)),
        group_id   = group,
        row.names  = colnames(counts)
      )
    )
    sce$group_id <- factor(sce$group_id)
    metadata(sce)$experiment_info <- data.frame(
      sample_id = colnames(counts),
      group_id  = group
    )
    pb <- aggregateData(sce, assay = "counts", by = c("cluster_id", "sample_id"),
                        fun = "sum", BPPARAM = bp_param)
    pb_filter   <- c("TRUE" = "none", "FALSE" = "both")[[as.character(IS_NULL)]]
    ds_res <- pbDS(pb, method = "edgeR", verbose = FALSE, BPPARAM = bp_param,
                   filter = pb_filter, min_cells = 0)
    rr <- resDS(sce, ds_res, bind = "row")
    data.table(gene = rr$gene, pvalue = rr$p_val, padj = rr$p_adj.loc,
               log2FoldChange = rr$logFC)
  },

  # ---------- aggregateBioVar ----------
  "aggregateBioVar" = {
    library(DESeq2)
    coldata <- data.frame(group = group, row.names = colnames(counts))
    dds <- DESeqDataSetFromMatrix(counts, coldata, ~group)
    dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = bp_param)
    rr <- results(dds, parallel = TRUE, BPPARAM = bp_param)
    data.table(gene = rownames(rr), pvalue = rr$pvalue, padj = rr$padj,
               log2FoldChange = rr$log2FoldChange)
  },

  # ---------- glmmTMB ----------
  "glmmTMB" = {
    library(glmmTMB)
    library(edgeR)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    cpm_mat <- cpm(y, log = TRUE)
    keep <- rowSums(counts >= 5) >= 2
    genes_to_test <- rownames(counts)[keep]
    cat(sprintf("  Testing %d genes with glmmTMB\n", length(genes_to_test)))
    grp <- group
    res_mat  <- t(vapply(genes_to_test, function(g) {
      dd <- data.frame(y = cpm_mat[g, ], group = grp)
      m  <- glmmTMB::glmmTMB(y ~ group, data = dd, family = gaussian(),
                             control = glmmTMB::glmmTMBControl(parallel = list(n = THREADS, autopar = TRUE)))
      ss <- summary(m)
      cc <- ss$coefficients$cond
      c(as.numeric(cc[2, "Pr(>|z|)"]), as.numeric(cc[2, "Estimate"]))
    }, numeric(2)))
    pvals <- res_mat[, 1]
    lfcs  <- res_mat[, 2]
    pvals[is.na(pvals)] <- 1
    data.table(gene = genes_to_test, pvalue = pvals,
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = lfcs)
  },

  # ---------- DV: diffVar (Levene-type via limma) ----------
  "dv_diffVar" = {
    library(limma)
    library(edgeR)
    design <- model.matrix(~group)
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    resid_abs <- abs(residuals(fit, v))
    fit_var <- lmFit(resid_abs, design)
    fit_var <- eBayes(fit_var)
    rr <- topTable(fit_var, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  },

  # ---------- DV: GAMLSS (sigma submodel LRT) ----------
  "dv_GAMLSS" = {
    library(gamlss)
    library(edgeR)
    library(matrixStats)
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
    g <- as.numeric(group) - 1
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    var1 <- rowVars(log_cpm[, idx1, drop = FALSE])
    var2 <- rowVars(log_cpm[, idx2, drop = FALSE])
    mad1 <- rowMads(log_cpm[, idx1, drop = FALSE])
    mad2 <- rowMads(log_cpm[, idx2, drop = FALSE])
    keep <- is.finite(var1) & is.finite(var2) & is.finite(mad1) & is.finite(mad2) &
      var1 > 0 & var2 > 0 & mad1 > 0 & mad2 > 0
    Sys.setenv(
      OMP_NUM_THREADS = 1L,
      OPENBLAS_NUM_THREADS = 1L,
      MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L,
      NUMEXPR_NUM_THREADS = 1L,
      BLIS_NUM_THREADS = 1L,
      R_DATATABLE_NUM_THREADS = 1L
    )
    keep_idx <- which(keep)
    chunk_n <- max(1L, min(THREADS, length(keep_idx)))
    gene_chunks <- split(keep_idx, ceiling(seq_along(keep_idx) * chunk_n / length(keep_idx)))
    chunk_mats <- lapply(gene_chunks, function(idx) log_cpm[idx, , drop = FALSE])
    pval_chunks <- BiocParallel::bplapply(chunk_mats, function(mat) {
      vapply(seq_len(nrow(mat)), function(i) {
        y <- mat[i, ]
        fit_h1 <- try(gamlss::gamlss(y ~ g, sigma.formula = ~ g, family = gamlss.dist::NO(),
                             trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50)),
                      silent = TRUE)
        fit_h0 <- try(gamlss::gamlss(y ~ g, sigma.formula = ~ 1, family = gamlss.dist::NO(),
                             trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50)),
                      silent = TRUE)
        err <- inherits(fit_h1, "try-error") | inherits(fit_h0, "try-error")
        switch(as.character(err),
          "TRUE" = 1,
          "FALSE" = {
            lr <- max(as.numeric(-2 * (logLik(fit_h0) - logLik(fit_h1))), 0)
            pchisq(lr, df = 1, lower.tail = FALSE)
          }
        )
      }, numeric(1))
    }, BPPARAM = bp_param)
    pvals <- rep(1, nrow(log_cpm))
    pvals[unlist(gene_chunks, use.names = FALSE)] <- unlist(pval_chunks, use.names = FALSE)
    pvals[!is.finite(pvals)] <- 1
    data.table(gene = rownames(counts), pvalue = pvals,
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = log2((var2 + 1e-8) / (var1 + 1e-8)))
  },

  # ---------- NEBULA (NB mixed model on single-cell counts) ----------
  "NEBULA" = {
    library(nebula)
    library(SingleCellExperiment)
    library(Matrix)
    SC_BASE <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk/data")
    sc_info <- switch(DATASET,
      "kang18"       = list(file = "single_cell/kang_pbmc_8vs8.rds",    id = "ind",        cond = "stim"),
      "segerstolpe"  = list(file = "single_cell/segerstolpe_pancreas.rds", id = "individual", cond = "disease"),
      "xin"          = list(file = "single_cell/xin_pancreas.rds",      id = "donor.id",   cond = "condition"),
      "zilionis"     = list(file = "single_cell/zilionis_lung.rds",     id = "Patient",    cond = "Tissue"),
      "jakel"        = list(file = "jakel_ms_sce.rds",                  id = "Sample",    cond = "Condition")
    )
    sce <- readRDS(file.path(SC_BASE, sc_info$file))
    sc_counts <- round(assay(sce, assayNames(sce)[1]))
    sid  <- colData(sce)[[sc_info$id]]
    cond <- colData(sce)[[sc_info$cond]]
    # Keep only cells with 2 condition levels
    lvls <- names(sort(table(cond), decreasing = TRUE))[1:2]
    keep <- cond %in% lvls
    sc_counts <- sc_counts[, keep]
    sid  <- sid[keep]
    cond <- factor(cond[keep], levels = lvls)
    # Null FDR: permute condition at SUBJECT level
    switch(as.character(IS_NULL),
      "TRUE" = {
        set.seed(NULL_SEED)
        subj_u    <- unique(sid)
        subj_cond <- cond[match(subj_u, sid)]
        subj_perm <- sample(subj_cond)
        cond      <- subj_perm[match(sid, subj_u)]
        cat("  [NULL FDR] Subject-level labels permuted\n")
      }
    )
    cat(sprintf("  NEBULA: %d genes x %d cells, %d subjects\n",
                nrow(sc_counts), ncol(sc_counts), length(unique(sid))))
    ord    <- order(sid)
    pred   <- model.matrix(~ cond[ord])
    offset <- colSums(sc_counts[, ord])
    switch(as.character(length(unique(as.character(cond))) > 1),
      "TRUE" = {
        re      <- nebula(sc_counts[, ord], sid[ord], pred = pred, offset = offset, ncore = THREADS)
        p_col   <- grep("^p_cond",   names(re$summary), value = TRUE)[1]
        lfc_col <- grep("^logFC_cond", names(re$summary), value = TRUE)[1]
        pvals   <- re$summary[[p_col]]
        lfcs    <- re$summary[[lfc_col]] / log(2)
        data.table(gene = re$summary$gene, pvalue = pvals,
                   padj = p.adjust(pvals, "BH"), log2FoldChange = lfcs)
      },
      "FALSE" = {
 cat(" [SKIP] All subjects same condition after permutation returning NA\n")
        data.table(gene = rownames(sc_counts), pvalue = NA_real_,
                   padj = NA_real_, log2FoldChange = NA_real_)
      }
    )
  },

  # ---------- DV: DiPhiSeq ----------
  "dv_DiPhiSeq" = {
    library(DiPhiSeq)
    counts_int <- round(counts)
    storage.mode(counts_int) <- "integer"
    classlab <- as.integer(group == levels(group)[2]) + 1L
    Sys.setenv(
      OMP_NUM_THREADS = 1L,
      OPENBLAS_NUM_THREADS = 1L,
      MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L,
      NUMEXPR_NUM_THREADS = 1L,
      BLIS_NUM_THREADS = 1L,
      R_DATATABLE_NUM_THREADS = 1L
    )
    depth <- colMeans(counts_int)
    diphiseq_fit <- NULL
    switch(as.character(DIPHISEQ_CAPTURE_OUTPUT),
      "TRUE" = {
        diphiseq_log_con <- file("NUL", open = "wt")
        sink(diphiseq_log_con)
        diphiseq_fit <- DiPhiSeq::diphiseq(countmat = counts_int, classlab = classlab, depth = depth)
        sink()
        close(diphiseq_log_con)
      },
      "FALSE" = {
        diphiseq_fit <- DiPhiSeq::diphiseq(countmat = counts_int, classlab = classlab, depth = depth)
      }
    )
    tab <- diphiseq_fit$tab
    data.table(gene = rownames(counts),
               pvalue = tab[, "p.value.phi"],
               padj = tab[, "fdr.phi"],
               log2FoldChange = log2((tab[, "phi2"] + 1e-8) / (tab[, "phi1"] + 1e-8)))
  }
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  %s done: %.2f sec, %d genes\n", METHOD, time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]

metrics <- data.table(
  method       = METHOD,
  n_genes      = nrow(res_df),
  time_sec     = time_sec,
  n_sig_fdr5   = sum(res_df$padj < 0.05),
  n_sig_fdr10  = sum(res_df$padj < 0.10),
  p_median     = median(res_df$pvalue),
  lfc_median   = median(abs(res_df$log2FoldChange)),
  dataset      = DATASET,
  is_null      = IS_NULL
)

# Ground truth (muscat_sim)
has_truth_sim <- !is.null(sim_truth) && !IS_NULL
truth_cols <- data.table(n_true_de = NA_integer_, TP_de_fdr5 = NA_integer_,
                          FP_de_fdr5 = NA_integer_, prec_de = NA_real_,
                          recall_de = NA_real_, F1_de = NA_real_)
switch(as.character(has_truth_sim),
  "TRUE" = {
    matched <- sim_truth[match(res_df$gene, sim_truth$gene), ]
    td <- ifelse(is.na(matched$is_de), FALSE, matched$is_de)
    sig <- res_df$padj < 0.05
    tp <- sum(sig & td); fp <- sum(sig & !td); n_d <- sum(td)
    pr <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rc <- ifelse(n_d == 0, 0, tp / n_d)
    f1 <- ifelse(pr + rc == 0, 0, 2 * pr * rc / (pr + rc))
    truth_cols <- data.table(n_true_de = n_d, TP_de_fdr5 = tp, FP_de_fdr5 = fp,
                              prec_de = round(pr, 4), recall_de = round(rc, 4),
                              F1_de = round(f1, 4))
    cat(sprintf("  GT DE: TP=%d FP=%d F1=%.4f\n", tp, fp, f1))
  }
)
metrics <- cbind(metrics, truth_cols)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
