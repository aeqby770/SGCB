# =============================================================================
#   method:    DESeq2 | edgeR | limma | glmGamPoi | EBSeq | NOISeq | samr
#   test_type: classic | simulation | gtex_pair | tcga | pseudobulk | null_fdr | seqc
#   dataset:   bottomly | sim_n10_pDE10_ef2.0 | liver_vs_kidney | brca | kang18 | 1 | seqc
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
THREADS <- as.integer(Sys.getenv("BENCHMARK_THREADS", "1"))
THREADS <- ifelse(is.na(THREADS) | THREADS < 1L, 1L, THREADS)
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

args      <- commandArgs(trailingOnly = TRUE)
METHOD    <- args[1]
TEST_TYPE <- args[2]
DATASET   <- args[3]

if (identical(METHOD, "dv_DiPhiSeq")) stop("dv_DiPhiSeq has been retired from benchmark_v2")

BASE      <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq")
DATA_ROOT <- file.path(BASE, "data/bulk")
OUT_DIR   <- file.path(BASE, "other_methods/results")
N_SUB     <- 50L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
library(data.table)

OUT_CSV <- file.path(OUT_DIR, paste0(METHOD, "__", TEST_TYPE, "__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(METHOD, "__", TEST_TYPE, "__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== Bulk Other DE: %s / %s / %s =====\n",
            Sys.time(), METHOD, TEST_TYPE, DATASET))

switch(as.character(file.exists(OUT_CSV)),
  "TRUE"  = { cat("  [SKIP] Already exists:", OUT_CSV, "\n"); q("no", status = 0) },
  "FALSE" = NULL)

loaded <- switch(TEST_TYPE,

  "classic" = {
    ds <- readRDS(file.path(DATA_ROOT, "01_classic", paste0("real_", DATASET, ".rds")))
    list(counts = as.matrix(ds$counts),
         group  = factor(ds$sample_info$group),
         truth  = NULL)
  },

  "simulation" = {
    ds <- readRDS(file.path(DATA_ROOT, "02_simulation", paste0(DATASET, ".rds")))
    list(counts = as.matrix(ds[["counts"]]),
         group  = factor(ds[["sample_info"]][["group"]]),
         truth  = ds[["truth"]])
  },

  "gtex_pair" = {
    tt <- strsplit(DATASET, "_vs_")[[1]]
    d1 <- readRDS(file.path(DATA_ROOT, "03_gtex", paste0("gtex_", tt[1], ".rds")))
    d2 <- readRDS(file.path(DATA_ROOT, "03_gtex", paste0("gtex_", tt[2], ".rds")))
    cg <- intersect(rownames(d1$counts), rownames(d2$counts))
    set.seed(12345)
    i1 <- sample.int(ncol(d1$counts), min(N_SUB, ncol(d1$counts)))
    i2 <- sample.int(ncol(d2$counts), min(N_SUB, ncol(d2$counts)))
    cm <- cbind(d1$counts[cg, i1], d2$counts[cg, i2])
    gm <- factor(c(rep(tt[1], length(i1)), rep(tt[2], length(i2))))
    list(counts = as.matrix(cm), group = gm, truth = NULL)
  },

  "tcga" = {
    ds <- readRDS(file.path(DATA_ROOT, "04_tcga",
                            paste0("tcga_", DATASET, "_tumor_vs_normal.rds")))
    list(counts = as.matrix(ds$counts),
         group  = factor(ds$group),
         truth  = NULL)
  },

  "pseudobulk" = {
    flist <- list.files(file.path(DATA_ROOT, "05_pseudobulk"),
                        pattern = DATASET, full.names = TRUE)
    ds <- readRDS(flist[1])
    cm <- as.matrix(ds$counts)
    rownames(cm) <- make.unique(rownames(cm))
    list(counts = cm, group = factor(ds$group), truth = NULL)
  },

  "null_fdr" = {
    ds <- readRDS(file.path(DATA_ROOT, "03_gtex", "gtex_liver.rds"))
    seed_val <- as.integer(DATASET)
    set.seed(seed_val)
    idx <- sample.int(ncol(ds$counts), min(2L * N_SUB, ncol(ds$counts)))
    cs  <- ds$counts[, idx]
    nh  <- length(idx) %/% 2L
    gp  <- factor(sample(c(rep("A", nh), rep("B", length(idx) - nh))))
    cs_mat <- as.matrix(cs)
    cs_mat[] <- pmin(round(cs_mat), .Machine$integer.max)
    storage.mode(cs_mat) <- "integer"
    list(counts = cs_mat, group = gp, truth = NULL)
  },

  "seqc" = {
    seqc_file <- file.path(DATA_ROOT, "06_seqc", "GSE49712_HTSeq.txt.gz")
    seqc_raw  <- read.table(gzfile(seqc_file), header = TRUE, row.names = 1,
                            check.names = FALSE)
    cts <- as.matrix(seqc_raw)
    storage.mode(cts) <- "integer"
    gp  <- factor(ifelse(grepl("^A_", colnames(cts)), "A", "B"))
    list(counts = cts, group = gp, truth = NULL)
  },

  "dv" = {
    dv_file <- file.path(DATA_ROOT, "07_gtex_brain_dv",
                          paste0("gtex_brain_", DATASET, ".rds"))
    ds <- readRDS(dv_file)
    set.seed(12345)
    lvs  <- levels(ds$group)
    idx1 <- which(ds$group == lvs[1])
    idx2 <- which(ds$group == lvs[2])
    i1   <- idx1[sample.int(length(idx1), min(N_SUB, length(idx1)))]
    i2   <- idx2[sample.int(length(idx2), min(N_SUB, length(idx2)))]
    cm   <- as.matrix(ds$counts[, c(i1, i2)])
    gp   <- factor(c(rep(lvs[1], length(i1)), rep(lvs[2], length(i2))))
    list(counts = cm, group = gp, truth = NULL)
  }
)

counts <- loaded$counts
group  <- loaded$group
truth  <- loaded$truth
rownames(counts) <- make.unique(rownames(counts))

cat(sprintf("  Data: %d genes x %d samples, groups: %s\n",
            nrow(counts), ncol(counts),
            paste(paste0(levels(group), "=", table(group)), collapse = ", ")))

# ===== DE =====
t0 <- proc.time()

res_df <- switch(METHOD,

  "DESeq2" = {
    library(DESeq2)
    coldata <- data.frame(group = group, row.names = colnames(counts))
    dds <- DESeqDataSetFromMatrix(counts, coldata, ~group)
    dds <- DESeq(dds, quiet = TRUE, parallel = TRUE, BPPARAM = bp_param)
    rr  <- results(dds, parallel = TRUE, BPPARAM = bp_param)
    data.table(gene = rownames(rr), pvalue = rr$pvalue, padj = rr$padj,
               log2FoldChange = rr$log2FoldChange)
  },

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

  "glmGamPoi" = {
    library(glmGamPoi)
    fit <- glm_gp(counts, design = ~group,
                   col_data = data.frame(group = group,
                                         row.names = colnames(counts)))
    rr  <- test_de(fit, contrast = colnames(fit$Beta)[2])
    data.table(gene = rr$name, pvalue = rr$pval, padj = rr$adj_pval,
               log2FoldChange = rr$lfc)
  },

  "EBSeq" = {
    library(EBSeq)
    sizes <- MedianNorm(counts)
    res   <- EBTest(Data = counts, Conditions = as.character(group),
                     sizeFactors = sizes, maxround = 5)
    ppmat <- GetPPMat(res)
    ppde  <- 1 - ppmat[, "PPEE"]
    nc <- GetNormalizedMat(counts, sizes)
    g1 <- which(group == levels(group)[1])
    g2 <- which(group == levels(group)[2])
    lfc <- log2((rowMeans(nc[, g2]) + 1) / (rowMeans(nc[, g1]) + 1))
    data.table(gene = rownames(counts), pvalue = ppde,
               padj = p.adjust(ppde, "BH"), log2FoldChange = lfc)
  },

  "NOISeq" = {
    library(NOISeq)
    factors <- data.frame(condition = group, row.names = colnames(counts))
    noiseq_data <- readData(data = as.matrix(counts), factors = factors)
    res <- noiseqbio(noiseq_data, k = 0.5, norm = "tmm", factor = "condition",
                      conditions = levels(group), r = 20)
    res_df <- res@results[[1]]
    prob   <- res_df$prob
    prob[is.na(prob)] <- 0
    pvals  <- 1 - prob
    data.table(gene = rownames(res_df), pvalue = pvals,
               padj = p.adjust(pvals, "BH"), log2FoldChange = res_df$log2FC)
  },

  "samr" = {
    library(samr)
    y_vec   <- as.numeric(factor(group))
    x_mat   <- as.matrix(log2(counts + 1))
    samr_data <- list(x = x_mat, y = y_vec, geneid = rownames(counts),
                       genenames = rownames(counts), logged2 = TRUE)
    samr_obj  <- samr(samr_data, resp.type = "Two class unpaired", nperms = 100)
    tt_stat   <- samr_obj$tt
    n1 <- sum(group == levels(group)[1])
    n2 <- sum(group == levels(group)[2])
    df_val <- n1 + n2 - 2
    pvals  <- 2 * pt(-abs(tt_stat), df = df_val)
    lfc    <- rowMeans(x_mat[, which(group == levels(group)[2])]) -
              rowMeans(x_mat[, which(group == levels(group)[1])])
    data.table(gene = rownames(counts), pvalue = pvals,
               padj = p.adjust(pvals, "BH"), log2FoldChange = lfc)
  },

  # ---------- DV: MDSeq ----------
  "dv_MDSeq" = {
    library(MDSeq)
    library(edgeR)
    cnt_safe  <- pmin(round(counts), .Machine$integer.max)
    keep_rows <- rowSums(cnt_safe) > 0
    cnts_norm <- normalize.counts(cnt_safe[keep_rows, , drop = FALSE])
    Xm  <- get.model.matrix(group)
    X_e <- Xm$mean[, colnames(Xm$mean) != "(Intercept)", drop = FALSE]
    U_e <- Xm$dispersion[, colnames(Xm$dispersion) != "(Intercept)", drop = FALSE]
    res <- MDSeq(cnts_norm, X = X_e, U = U_e, mc.cores = THREADS, verbose = FALSE)
    rr  <- res$Dat
    z   <- rr$gamma.group1 / rr$se.gamma.group1
    pv  <- 2 * pnorm(-abs(z))
    out <- data.table(gene = rownames(counts), pvalue = 1, padj = 1, log2FoldChange = 0)
    out[match(rownames(rr), gene), `:=`(pvalue = pv, padj = p.adjust(pv, "BH"),
                                         log2FoldChange = 2 * rr$gamma.group1 / log(2))]
    out
  },

  # ---------- DV: clrDV ----------
  "dv_clrDV" = {
    library(edgeR)
    library(matrixStats)
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
    clr <- log_cpm - rowMeans(log_cpm)
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    var1 <- rowVars(clr[, idx1])
    var2 <- rowVars(clr[, idx2])
    n1 <- length(idx1); n2 <- length(idx2)
    f_stat <- var2 / (var1 + 1e-10)
    pvals <- 2 * pmin(pf(f_stat, n2 - 1, n1 - 1, lower.tail = FALSE),
                       pf(f_stat, n2 - 1, n1 - 1, lower.tail = TRUE))
    data.table(gene = rownames(counts), pvalue = pvals,
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = log2((var2 + 1e-8) / (var1 + 1e-8)))
  },

  # ---------- DV: diffVar (Levene via limma) ----------
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
               padj   = tab[, "fdr.phi"],
               log2FoldChange = log2((tab[, "phi2"] + 1e-8) / (tab[, "phi1"] + 1e-8)))
  }
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  %s done: %.2f sec, %d genes\n", METHOD, time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]

metrics <- data.table(
  method    = METHOD,
  test_type = TEST_TYPE,
  dataset   = DATASET,
  n_genes   = nrow(res_df),
  time_sec  = time_sec,
  n_sig_fdr5  = sum(res_df$padj < 0.05),
  n_sig_fdr10 = sum(res_df$padj < 0.10),
  p_median    = median(res_df$pvalue)
)

has_truth <- !is.null(truth) && !is.null(truth$is_de)
switch(as.character(has_truth),
  "TRUE" = {
    matched <- truth[match(res_df$gene, truth$gene), ]
    td  <- ifelse(is.na(matched$is_de), FALSE, matched$is_de)
    sig <- res_df$padj < 0.05
    tp  <- sum(sig & td); fp <- sum(sig & !td); n_d <- sum(td)
    pr  <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rc  <- ifelse(n_d == 0, 0, tp / n_d)
    f1  <- ifelse(pr + rc == 0, 0, 2 * pr * rc / (pr + rc))
    metrics[, `:=`(n_true_de = n_d, TP_fdr5 = tp, FP_fdr5 = fp,
                   precision = round(pr, 4), recall = round(rc, 4),
                   F1 = round(f1, 4))]
    cat(sprintf("  GT: TP=%d FP=%d F1=%.4f\n", tp, fp, f1))
    NULL
  },
  "FALSE" = NULL
)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
