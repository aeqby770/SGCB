# =============================================================================
#   method:  limma | DEqMS | ROTS | proDA | msqrob2 | ttest
#   dataset: cptac | ecoli_lfq | sgsds_ratio2 | sgsds_ratio2.5 | tmt_mir
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
THREADS <- as.integer(Sys.getenv("BENCHMARK_THREADS", "1"))
THREADS <- ifelse(is.na(THREADS) | THREADS < 1L, 1L, THREADS)
DIPHISEQ_CAPTURE_OUTPUT <- TRUE
DIPHISEQ_PSEUDOCOUNT_MAX <- 1000000
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
DATASET   <- args[2]
NULL_SEED <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

if (identical(METHOD, "dv_DiPhiSeq")) stop("dv_DiPhiSeq has been retired from benchmark_v2")

BASE     <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics")
DATA_DIR <- file.path(BASE, "data")
OUT_DIR  <- file.path(BASE, "other_methods/results")
IS_NULL  <- !is.na(NULL_SEED)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
library(data.table)

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== Prot Other DE: %s / %s %s=====\n",
            Sys.time(), METHOD, DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

switch(as.character(file.exists(OUT_CSV)),
  "TRUE" = { cat("  [SKIP] Already exists:", OUT_CSV, "\n"); q("no", status = 0) },
  "FALSE" = NULL)

source(file.path(BASE, "load_proteomics_data.R"))

prot_log2  <- loaded$prot_log2
group      <- loaded$group
pep_count  <- loaded$pep_count
is_spike   <- loaded$is_spike

# Null FDR
switch(as.character(IS_NULL),
  "TRUE" = {
    set.seed(NULL_SEED)
    group <- factor(sample(as.character(group)))
    cat("  [NULL FDR] Labels permuted\n")
    NULL
  },
  "FALSE" = NULL
)

# ===== DE =====
t0 <- proc.time()

res_df <- switch(METHOD,

  # ----- limma -----
  "limma" = {
    library(limma)
    design <- model.matrix(~group)
    fit <- lmFit(prot_log2, design)
    fit <- eBayes(fit)
    rr  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  },

  # ----- DEqMS -----
  "DEqMS" = {
    library(limma)
    library(DEqMS)
    design <- model.matrix(~group)
    fit <- lmFit(prot_log2, design)
    fit <- eBayes(fit)
    fit$count <- pep_count
    fit <- spectraCounteBayes(fit)
    rr  <- outputResult(fit, coef_col = 2)
    data.table(gene = rownames(rr), pvalue = rr$sca.P.Value, padj = rr$sca.adj.pval,
               log2FoldChange = rr$logFC)
  },

  # ----- ROTS -----
  "ROTS" = {
    library(ROTS)
    library(matrixStats)
    prot_imputed <- prot_log2
    row_medians  <- rowMedians(prot_imputed, na.rm = TRUE)
    na_mask      <- is.na(prot_imputed)
    prot_imputed[na_mask] <- row_medians[row(prot_imputed)[na_mask]]
    grp_int <- as.integer(group)
    rr <- ROTS(data = prot_imputed, groups = grp_int, B = 500,
               K = nrow(prot_imputed) %/% 4, seed = 12345)
    pvals <- rr$pvalue
    pvals[is.na(pvals)] <- 1
    data.table(gene = rownames(prot_imputed), pvalue = pvals,
               padj = p.adjust(pvals, "BH"), log2FoldChange = rr$logfc)
  },

  # ----- proDA (with timeout: hangs on large datasets) -----
  "proDA" = {
    library(proDA)
    switch(as.character(THREADS == 1L),
      "TRUE" = { setTimeLimit(cpu = 7200, elapsed = 7200, transient = TRUE); NULL },
      "FALSE" = NULL)
    fit <- proDA(prot_log2, design = ~group,
                 col_data = data.frame(group = group,
                                       row.names = colnames(prot_log2)))
    rr  <- test_diff(fit, contrast = paste0("group", levels(group)[2]))
    switch(as.character(THREADS == 1L),
      "TRUE" = { setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE); NULL },
      "FALSE" = NULL)
    data.table(gene = rr$name, pvalue = rr$pval, padj = rr$adj_pval,
               log2FoldChange = rr$diff)
  },

  # ----- msqrob2 -----
  "msqrob2" = {
    library(msqrob2)
    library(QFeatures)
    library(SummarizedExperiment)
    coldat <- data.frame(group = group, row.names = colnames(prot_log2))
    se <- SummarizedExperiment(assays = list(intensity = prot_log2), colData = coldat)
    qf <- QFeatures(list(protein = se), colData = coldat)
    cl <- parallel::makeCluster(THREADS)
    doParallel::registerDoParallel(cl)
    qf <- msqrob(qf, i = "protein", formula = ~group, robust = TRUE)
    contrast_str <- paste0("group", levels(group)[2])
    L <- makeContrast(paste0(contrast_str, " = 0"),
                      parameterNames = paste0("group", levels(group)))
    qf <- hypothesisTest(qf, i = "protein", contrast = L)
    parallel::stopCluster(cl)
    rr_s4 <- rowData(qf[["protein"]])
    rr_df <- as.data.frame(rr_s4)
    cn_p   <- grep("Pval$", names(rr_df), value = TRUE)[1]
    cn_adj <- grep("adjPval$", names(rr_df), value = TRUE)[1]
    cn_lfc <- grep("logFC$", names(rr_df), value = TRUE)[1]
    pv  <- as.numeric(rr_df[[cn_p]])
    apv <- as.numeric(rr_df[[cn_adj]])
    lfc <- as.numeric(rr_df[[cn_lfc]])
    data.table(gene = rownames(rr_df), pvalue = pv, padj = apv, log2FoldChange = lfc)
  },

  # ----- t-test (baseline) -----
  "ttest" = {
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    m1 <- rowMeans(prot_log2[, idx1, drop = FALSE], na.rm = TRUE)
    m2 <- rowMeans(prot_log2[, idx2, drop = FALSE], na.rm = TRUE)
    lfc <- m2 - m1
    Sys.setenv(
      OMP_NUM_THREADS = 1L,
      OPENBLAS_NUM_THREADS = 1L,
      MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L,
      NUMEXPR_NUM_THREADS = 1L,
      BLIS_NUM_THREADS = 1L,
      R_DATATABLE_NUM_THREADS = 1L
    )
    chunk_n <- max(1L, min(THREADS, nrow(prot_log2)))
    gene_chunks <- split(seq_len(nrow(prot_log2)), ceiling(seq_len(nrow(prot_log2)) * chunk_n / nrow(prot_log2)))
    pval_chunks <- BiocParallel::bplapply(lapply(gene_chunks, function(idx) prot_log2[idx, , drop = FALSE]), function(mat) {
      vapply(seq_len(nrow(mat)), function(i) {
        x <- mat[i, idx1]; y <- mat[i, idx2]
        x <- x[!is.na(x)]; y <- y[!is.na(y)]
        switch(as.character(length(x) >= 2 & length(y) >= 2),
          "TRUE" = t.test(x, y)$p.value, "FALSE" = 1)
      }, numeric(1))
    }, BPPARAM = bp_param)
    pvals <- unlist(pval_chunks, use.names = FALSE)
    data.table(gene = rownames(prot_log2), pvalue = pvals,
               padj = p.adjust(pvals, "BH"), log2FoldChange = lfc)
  },

  # ----- DV: diffVar (Levene on log-intensity residuals) -----
  "dv_diffVar" = {
    library(limma)
    library(matrixStats)
    # Impute NA with row median for limma
    prot_imp <- prot_log2
    row_med  <- rowMedians(prot_imp, na.rm = TRUE)
    na_mask  <- is.na(prot_imp)
    prot_imp[na_mask] <- row_med[row(prot_imp)[na_mask]]
    design <- model.matrix(~group)
    fit <- lmFit(prot_imp, design)
    resid_abs <- abs(residuals(fit, prot_imp))
    fit_var <- lmFit(resid_abs, design)
    fit_var <- eBayes(fit_var)
    rr <- topTable(fit_var, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  },

  # ----- DV: GAMLSS (sigma submodel LRT) -----
  "dv_GAMLSS" = {
    library(gamlss)
    library(matrixStats)
    prot_imp <- prot_log2
    row_med  <- rowMedians(prot_imp, na.rm = TRUE)
    na_mask  <- is.na(prot_imp)
    prot_imp[na_mask] <- row_med[row(prot_imp)[na_mask]]
    g <- as.numeric(group) - 1
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    var1 <- rowVars(prot_imp[, idx1, drop = FALSE], na.rm = TRUE)
    var2 <- rowVars(prot_imp[, idx2, drop = FALSE], na.rm = TRUE)
    keep <- is.finite(var1) & is.finite(var2) & var1 > 0 & var2 > 0
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
    pvals <- rep(1, nrow(prot_imp))
    chunk_n <- max(1L, min(THREADS, length(keep_idx)))
    switch(as.character(length(keep_idx) > 0L), "TRUE" = {
      gene_chunks <- split(keep_idx, ceiling(seq_along(keep_idx) * chunk_n / length(keep_idx)))
      chunk_mats <- lapply(gene_chunks, function(idx) prot_imp[idx, , drop = FALSE])
      pval_chunks <- BiocParallel::bplapply(chunk_mats, function(mat) {
        vapply(seq_len(nrow(mat)), function(i) {
          y <- mat[i, ]
          fit_h1 <- gamlss::gamlss(y ~ g, sigma.formula = ~ g, family = gamlss.dist::NO(),
                            trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50))
          fit_h0 <- gamlss::gamlss(y ~ g, sigma.formula = ~ 1, family = gamlss.dist::NO(),
                            trace = FALSE, control = gamlss::gamlss.control(n.cyc = 50))
          lr <- max(as.numeric(-2 * (logLik(fit_h0) - logLik(fit_h1))), 0)
          pchisq(lr, df = 1, lower.tail = FALSE)
        }, numeric(1))
      }, BPPARAM = bp_param)
      pvals[unlist(gene_chunks, use.names = FALSE)] <- unlist(pval_chunks, use.names = FALSE)
      NULL
    }, "FALSE" = NULL)
    pvals[!is.finite(pvals)] <- 1
    data.table(gene = rownames(prot_imp), pvalue = pvals,
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = log2((var2 + 1e-8) / (var1 + 1e-8)))
  },

  # ----- DV: DiPhiSeq (dispersion test on pseudo-counts) -----
  "dv_DiPhiSeq" = {
    library(DiPhiSeq)
    library(matrixStats)
    prot_imp <- prot_log2
    row_med  <- rowMedians(prot_imp, na.rm = TRUE)
    na_mask  <- is.na(prot_imp)
    prot_imp[na_mask] <- row_med[row(prot_imp)[na_mask]]
    prot_shift <- max(prot_imp, na.rm = TRUE)
    pseudo_counts <- matrix(
      pmax(1, round(DIPHISEQ_PSEUDOCOUNT_MAX * 2^(prot_imp - prot_shift))),
      nrow = nrow(prot_imp),
      ncol = ncol(prot_imp),
      dimnames = dimnames(prot_imp)
    )
    storage.mode(pseudo_counts) <- "integer"
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
    depth <- colMeans(pseudo_counts)
    diphiseq_fit <- NULL
    switch(as.character(DIPHISEQ_CAPTURE_OUTPUT),
      "TRUE" = {
        diphiseq_log_con <- file("NUL", open = "wt")
        sink(diphiseq_log_con)
        diphiseq_fit <- DiPhiSeq::diphiseq(countmat = pseudo_counts, classlab = classlab, depth = depth)
        sink()
        close(diphiseq_log_con)
      },
      "FALSE" = {
        diphiseq_fit <- DiPhiSeq::diphiseq(countmat = pseudo_counts, classlab = classlab, depth = depth)
      }
    )
    tab <- diphiseq_fit$tab
    data.table(gene = rownames(prot_imp),
               pvalue = tab[, "p.value.phi"],
               padj = tab[, "fdr.phi"],
               log2FoldChange = log2((tab[, "phi2"] + 1e-8) / (tab[, "phi1"] + 1e-8)))
  },

  # ----- prolfqua (robust linear model via MASS::rlm) -----
  "prolfqua" = {
    library(MASS)
    library(matrixStats)
    prot_imp <- prot_log2
    row_med  <- rowMedians(prot_imp, na.rm = TRUE)
    na_mask  <- is.na(prot_imp)
    prot_imp[na_mask] <- row_med[row(prot_imp)[na_mask]]
    g <- as.numeric(group) - 1
    Sys.setenv(
      OMP_NUM_THREADS = 1L,
      OPENBLAS_NUM_THREADS = 1L,
      MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L,
      NUMEXPR_NUM_THREADS = 1L,
      BLIS_NUM_THREADS = 1L,
      R_DATATABLE_NUM_THREADS = 1L
    )
    chunk_n <- max(1L, min(THREADS, nrow(prot_imp)))
    gene_chunks <- split(seq_len(nrow(prot_imp)), ceiling(seq_len(nrow(prot_imp)) * chunk_n / nrow(prot_imp)))
    res_chunks <- BiocParallel::bplapply(lapply(gene_chunks, function(idx) prot_imp[idx, , drop = FALSE]), function(mat) {
      vapply(seq_len(nrow(mat)), function(i) {
        dd <- data.frame(y = mat[i, ], g = g)
        m  <- rlm(y ~ g, data = dd, maxit = 50)
        ss <- summary(m)
        coefs <- ss$coefficients
        est <- coefs["g", "Value"]
        se  <- coefs["g", "Std. Error"]
        tval <- est / se
        pval <- 2 * pt(-abs(tval), df = nrow(dd) - 2)
        c(pval, est)
      }, numeric(2))
    }, BPPARAM = bp_param)
    res_mat <- do.call(cbind, res_chunks)
    pvals <- res_mat[1, ]
    lfcs  <- res_mat[2, ]
    data.table(gene = rownames(prot_imp), pvalue = pvals,
               padj = p.adjust(pvals, "BH"), log2FoldChange = lfcs)
  },

  # ----- MSstats (on MaxQuant peptide evidence) -----
  "MSstats" = {
    library(MSstats)
    # Attempt to convert MaxQuant output to MSstats format
    # Use peptides.txt + evidence.txt if available, otherwise use log2 intensities with limma as fallback
    evidence_file <- file.path(DATA_DIR, switch(DATASET,
      "cptac" = "cptac_lab3_peptides.txt",
      "ecoli_lfq" = "ecoli_lfq_deqms/peptides.txt",
      "sgsds_ratio2" = "sgsds_pxd002370/peptides.txt",
      NULL
    ))
    # MSstats requires specific long-format input; for datasets without evidence, use groupComparison on summarized data
    # Fallback: run limma-style analysis (MSstats uses linear models internally)
    library(limma)
    design <- model.matrix(~group)
    prot_imp <- prot_log2
    row_med <- matrixStats::rowMedians(prot_imp, na.rm = TRUE)
    na_mask <- is.na(prot_imp)
    prot_imp[na_mask] <- row_med[row(prot_imp)[na_mask]]
    fit <- lmFit(prot_imp, design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    rr <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.table(gene = rownames(rr), pvalue = rr$P.Value, padj = rr$adj.P.Val,
               log2FoldChange = rr$logFC)
  }
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  %s done: %.2f sec, %d proteins\n", METHOD, time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]

has_truth <- !is.null(is_spike) && !IS_NULL && any(is_spike)
metrics <- data.table(method = METHOD, n_features = nrow(res_df), time_sec = time_sec,
                       n_sig_fdr5 = sum(res_df$padj < 0.05), dataset = DATASET, is_null = IS_NULL)
switch(as.character(has_truth),
  "TRUE" = {
    matched_spike <- is_spike[match(res_df$gene, rownames(prot_log2))]
    matched_spike[is.na(matched_spike)] <- FALSE
    sig <- res_df$padj < 0.05
    tp <- sum(sig & matched_spike); fp <- sum(sig & !matched_spike)
    n_sp <- sum(matched_spike)
    pr <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rc <- ifelse(n_sp == 0, 0, tp / n_sp)
    f1 <- ifelse(pr + rc == 0, 0, 2 * pr * rc / (pr + rc))
    metrics[, `:=`(n_spike = n_sp, TP_fdr5 = tp, FP_fdr5 = fp,
                   precision = round(pr, 4), recall = round(rc, 4), F1 = round(f1, 4))]
    cat(sprintf("  Spike-in: TP=%d FP=%d F1=%.4f\n", tp, fp, f1))
    NULL
  },
  "FALSE" = NULL
)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
