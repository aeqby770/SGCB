# =============================================================================
#   method:  SGCB | DEqMS | limma | ROTS | proDA
#   dataset: cptac | sim_proteomics
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = 1L)
data.table::setDTthreads(1L)

args      <- commandArgs(trailingOnly = TRUE)
METHOD    <- args[1]
DATASET   <- args[2]
NULL_SEED <- as.integer(args[3])

BASE      <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics")
DATA_DIR  <- file.path(BASE, "data")
OUT_DIR   <- file.path(BASE, "results")
IS_NULL   <- !is.na(NULL_SEED)
SGCB_ARGS <- list(
  prior_de = "auto",
  prior_dv = "auto",
  prior_sd_de = 1.0,
  prior_sd_dv = 0.5
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(data.table)

tag     <- ifelse(IS_NULL, paste0("null_", NULL_SEED, "__"), "")
OUT_CSV <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, ".csv"))
OUT_RDS <- file.path(OUT_DIR, paste0(tag, METHOD, "__", DATASET, "_full.rds"))

cat(sprintf("[%s] ===== Proteomics DE: %s / %s %s=====\n",
            Sys.time(), METHOD, DATASET,
            ifelse(IS_NULL, paste0("(null seed=", NULL_SEED, ") "), "")))

loaded <- switch(DATASET,

  # ----- CPTAC Study 6 (MaxQuant peptide output, 6A vs 6B) -----
  # 6A = 0.25 fmol UPS1 in yeast; 6B = 0.74 fmol UPS1 in yeast
  # UPS1 (HUMAN) proteins = true DE; yeast (YEAST) = non-DE
  "cptac" = {
    pep_dt <- fread(file.path(DATA_DIR, "cptac_lab3_peptides.txt"), integer64 = "numeric")

    pep_dt <- pep_dt[is.na(Reverse) | Reverse != "+"]
    pep_dt <- pep_dt[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]

    int_cols <- grep("^Intensity ", names(pep_dt), value = TRUE)
    int_mat  <- as.matrix(pep_dt[, ..int_cols])
    int_mat[int_mat == 0] <- NA  # MaxQuant 0 = missing

    prot_id <- pep_dt[["Leading razor protein"]]

    log2_pep <- log2(int_mat)
    pep_tbl  <- data.table(protein = prot_id, log2_pep)
    pep_tbl  <- pep_tbl[!is.na(protein) & protein != ""]

    pep_counts <- pep_tbl[, .(n_pep = .N), by = protein]

    prot_long <- melt(pep_tbl, id.vars = "protein", variable.name = "sample",
                      value.name = "log2int")
    prot_wide <- dcast(prot_long, protein ~ sample, value.var = "log2int",
                       fun.aggregate = median, na.rm = TRUE)
    prot_mat  <- as.matrix(prot_wide[, -1])
    rownames(prot_mat) <- prot_wide$protein
    prot_mat[is.nan(prot_mat)] <- NA

    grp_vec   <- fifelse(grepl("6A", colnames(prot_mat)), "6A", "6B")
    n_valid_A <- rowSums(!is.na(prot_mat[, grp_vec == "6A", drop = FALSE]))
    n_valid_B <- rowSums(!is.na(prot_mat[, grp_vec == "6B", drop = FALSE]))
    keep      <- n_valid_A >= 2 & n_valid_B >= 2
    prot_mat  <- prot_mat[keep, ]

    # Ground truth: UPS1 (HUMAN) = spike-in DE; YEAST = non-DE
    is_spike <- grepl("HUMAN|UPS", rownames(prot_mat))

    # Peptide counts for DEqMS (matched)
    pep_count_matched <- pep_counts[match(rownames(prot_mat), pep_counts$protein), n_pep]
    pep_count_matched[is.na(pep_count_matched)] <- 1L

    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("6A", "6B")),
         pep_count = pep_count_matched, is_spike = is_spike)
  },

  # ----- Simulated proteomics (spike-in DE + DV from CPTAC base) -----
  "sim_proteomics" = {
    sim_d <- readRDS(file.path(DATA_DIR, "sim_proteomics.rds"))
    list(prot_log2 = sim_d$prot_log2, group = sim_d$group,
         pep_count = sim_d$pep_count, is_spike = sim_d$is_spike,
         sim_truth = sim_d$truth)
  },

  # ----- PXD000279: E.coli spike-in LFQ (DEqMS core benchmark) -----
  # H = Heavy (higher E.coli), L = Light (lower E.coli)
  # ECOLI proteins = DE; HUMAN proteins = non-DE
  "ecoli_lfq" = {
    pg <- fread(file = file.path(DATA_DIR, "ecoli_lfq_deqms/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(Contaminant) | Contaminant != "+"]

    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)

    grp_vec <- fifelse(grepl("^LFQ intensity H", lfq_cols), "H", "L")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "H", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "L", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]

    # Species from Fasta headers
    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("ECOLI", fasta, ignore.case = TRUE)

    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L

    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("L", "H")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # D vs E, spike-in (UPS/HUMAN) = DE, yeast background = non-DE
  "sgsds_ratio2" = {
    pg <- fread(file = file.path(DATA_DIR, "sgsds_pxd002370/txt/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]

    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)

    grp_vec <- fifelse(grepl("^LFQ intensity D", lfq_cols), "D", "E")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "D", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "E", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]

    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("HUMAN|UPS", fasta, ignore.case = TRUE)

    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L

    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("D", "E")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- PXD002370 Ratio2.5: SGSDS spike-in (C vs D) -----
  "sgsds_ratio2.5" = {
    pg <- fread(file = file.path(DATA_DIR, "sgsds_pxd002370/ratio2.5/txt/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]

    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)

    grp_vec <- fifelse(grepl("^LFQ intensity C", lfq_cols), "C", "D")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "C", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "D", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]

    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("HUMAN|UPS", fasta, ignore.case = TRUE)

    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L

    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("C", "D")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- PXD004163: TMT miR proteomics (DEqMS TMT benchmark) -----
  # TMT 10-plex: 126-128C = ctrl (4 channels), 129N-131 = miR-372 (4 channels)
  # No spike-in ground truth; real miR knockdown experiment
  "tmt_mir" = {
    pg <- fread(file = file.path(DATA_DIR, "tmt_mir_deqms/Yan_miR_Protein_table.txt"), integer64 = "numeric")
    tmt_cols <- grep("tmt10plex_1[23]", names(pg), value = TRUE)
    tmt_cols <- tmt_cols[!grepl("quanted", tmt_cols)]
    psm_cols <- tmt_cols[grepl("quanted", grep("tmt10plex", names(pg), value = TRUE))]

    prot_mat <- as.matrix(pg[, ..tmt_cols])
    rownames(prot_mat) <- pg[["Protein accession"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)

    # DEqMS vignette: 126,127N,127C,128N = ctrl; 128C,129N,129C,130N = miR-372
    grp_vec <- fifelse(grepl("126|127N|127C|128N", tmt_cols), "ctrl", "miR372")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "ctrl", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "miR372", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]

    # PSM count for DEqMS
    pep_cnt <- pg[["miR FASP_# PSMs"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L

    # No spike-in ground truth for TMT miR (real experiment)
    is_spike <- rep(FALSE, sum(keep))

    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("ctrl", "miR372")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  }
)

prot_log2  <- loaded$prot_log2
group      <- loaded$group
pep_count  <- loaded$pep_count
is_spike   <- loaded$is_spike
sim_truth  <- loaded$sim_truth  # NULL for cptac

switch(as.character(IS_NULL),
  "TRUE" = {
    set.seed(NULL_SEED)
    group <- factor(sample(as.character(group)))
    cat("  [NULL FDR] Labels permuted\n")
    NULL
  },
  "FALSE" = NULL
)

cat(sprintf("  Data: %d proteins x %d samples, groups: %s\n",
            nrow(prot_log2), ncol(prot_log2),
            paste(paste0(levels(group), "=", table(group)), collapse=", ")))
cat(sprintf("  Spike-in proteins: %d / %d\n", sum(is_spike), length(is_spike)))

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

  # ----- DEqMS (limma + spectra count moderation) -----
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
               padj = p.adjust(pvals, "BH"),
               log2FoldChange = rr$logfc)
  },

  "proDA" = {
    library(proDA)
    fit <- proDA(prot_log2, design = ~group,
                 col_data = data.frame(group = group,
                                       row.names = colnames(prot_log2)))
    rr  <- test_diff(fit, contrast = paste0("group", levels(group)[2]))
    data.table(gene = rr$name, pvalue = rr$pval, padj = rr$adj_pval,
               log2FoldChange = rr$diff)
  },

  # ----- SGCB (pseudo-counts from intensity, via sgcbProtein) -----
  "SGCB" = {
    library(SGCB)
    med_log2 <- median(prot_log2, na.rm = TRUE)
    shifted  <- prot_log2 - med_log2 + log2(1000)
    pseudo_counts <- round(2^shifted)
    pseudo_counts[is.na(pseudo_counts) | pseudo_counts < 0] <- 0L
    storage.mode(pseudo_counts) <- "integer"
    res <- do.call(sgcbProtein, c(list(counts = pseudo_counts, group = group), SGCB_ARGS))
    rr  <- as.data.frame(res)
    dt_out <- data.table(
      gene           = rr$gene_id,
      pvalue         = rr$pvalue,
      padj           = rr$padj,
      log2FoldChange = rr$log2FoldChange,
      p_de_post      = rr$p_de_post,
      p_de_only_post = rr$p_de_only_post,
      de_fdr_call    = rr$de_fdr_call
    )
    # DV columns if present
    has_dv_col <- "pvalue_dv" %in% names(rr)
    switch(as.character(has_dv_col),
      "TRUE" = {
        dt_out[, dv_pvalue := rr$pvalue_dv]
        dt_out[, dv_padj := rr$padj_dv]
        dt_out[, dv_pvalue_model := rr$pvalue_dv]
        dt_out[, dv_padj_model := rr$padj_dv]
        dt_out[, dv_log2_cv2 := rr$dv_log2ratio]
        dt_out[, p_dv_post := rr$p_dv_post]
        dt_out[, p_dv_only_post := rr$p_dv_only_post]
        dt_out[, p_de_dv_post := rr$p_de_dv_post]
        dt_out[, dv_fdr_call := rr$dv_fdr_call]
        NULL
      },
      "FALSE" = NULL
    )
    dt_out
  }
)

time_sec <- round(as.numeric((proc.time() - t0)[3]), 2)
cat(sprintf("  %s done: %.2f sec, %d proteins\n", METHOD, time_sec, nrow(res_df)))

res_df[is.na(pvalue), pvalue := 1]
res_df[is.na(padj),   padj   := 1]
res_df[is.na(log2FoldChange), log2FoldChange := 0]
has_dv <- "dv_pvalue" %in% names(res_df)
has_dd <- "pvalue_dd" %in% names(res_df)
has_score <- "SGCB_Score_padj" %in% names(res_df)
has_post <- "p_de_post" %in% names(res_df)

has_truth <- !is.null(is_spike) && !IS_NULL

truth_cols <- data.table(
  n_spike = NA_integer_, TP_fdr5 = NA_integer_, FP_fdr5 = NA_integer_,
  precision_fdr5 = NA_real_, recall_fdr5 = NA_real_, F1_fdr5 = NA_real_,
  TP_post08 = NA_integer_, FP_post08 = NA_integer_,
  precision_post08 = NA_real_, recall_post08 = NA_real_, F1_post08 = NA_real_
)

switch(as.character(has_truth),
  "TRUE" = {
    matched_spike <- is_spike[match(res_df$gene, rownames(prot_log2))]
    matched_spike[is.na(matched_spike)] <- FALSE
    sig05 <- res_df$padj < 0.05
    sig_post <- ifelse(has_post, res_df$p_de_only_post > 0.8, FALSE)
    tp05  <- sum(sig05 & matched_spike)
    fp05  <- sum(sig05 & !matched_spike)
    n_sp  <- sum(matched_spike)
    tp_p  <- sum(sig_post & matched_spike)
    fp_p  <- sum(sig_post & !matched_spike)
    prec  <- fifelse(tp05 + fp05 == 0, 0, tp05 / (tp05 + fp05))
    rec   <- fifelse(n_sp == 0, 0, tp05 / n_sp)
    f1    <- fifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
    prec_p <- fifelse(tp_p + fp_p == 0, 0, tp_p / (tp_p + fp_p))
    rec_p  <- fifelse(n_sp == 0, 0, tp_p / n_sp)
    f1_p   <- fifelse(prec_p + rec_p == 0, 0, 2 * prec_p * rec_p / (prec_p + rec_p))

    truth_cols <- data.table(
      n_spike = n_sp, TP_fdr5 = tp05, FP_fdr5 = fp05,
      precision_fdr5 = round(prec, 4), recall_fdr5 = round(rec, 4),
      F1_fdr5 = round(f1, 4),
      TP_post08 = tp_p, FP_post08 = fp_p,
      precision_post08 = round(prec_p, 4), recall_post08 = round(rec_p, 4),
      F1_post08 = round(f1_p, 4)
    )
    cat(sprintf("  Spike-in: TP=%d, FP=%d, F1=%.4f\n", tp05, fp05, f1))
    NULL
  },
  "FALSE" = NULL
)

metrics <- cbind(data.table(
  method       = METHOD,
  n_features   = nrow(res_df),
  time_sec     = time_sec,
  n_sig_fdr1   = sum(res_df$padj < 0.01),
  n_sig_fdr5   = sum(res_df$padj < 0.05),
  n_sig_fdr10  = sum(res_df$padj < 0.10),
  n_sig_lfc0.5 = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5),
  n_sig_lfc1   = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0),
  p_median     = median(res_df$pvalue),
  p_mean       = mean(res_df$pvalue),
  p_lt_0.05    = mean(res_df$pvalue < 0.05),
  lfc_median   = median(abs(res_df$log2FoldChange)),
  n_post_de_0.8 = ifelse(has_post, sum(res_df$p_de_post > 0.8, na.rm = TRUE), NA_integer_),
  n_post_deonly_0.8 = ifelse(has_post, sum(res_df$p_de_only_post > 0.8, na.rm = TRUE), NA_integer_),
  n_de_fdr_call = ifelse("de_fdr_call" %in% names(res_df), sum(res_df$de_fdr_call, na.rm = TRUE), NA_integer_),
  n_post_dv_0.8 = ifelse("p_dv_post" %in% names(res_df), sum(res_df$p_dv_post > 0.8, na.rm = TRUE), NA_integer_),
  n_post_dvonly_0.8 = ifelse("p_dv_only_post" %in% names(res_df), sum(res_df$p_dv_only_post > 0.8, na.rm = TRUE), NA_integer_),
  n_dv_fdr_call = ifelse("dv_fdr_call" %in% names(res_df), sum(res_df$dv_fdr_call, na.rm = TRUE), NA_integer_),
  prior_de_used = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_de_used")), NA_real_),
  prior_dv_used = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_dv_used")), NA_real_),
  prior_dv_upper = ifelse(METHOD == "SGCB", as.numeric(attr(res, "prior_dv_upper")), NA_real_),
  dataset      = DATASET,
  is_null      = IS_NULL,
  n_dv_fdr5    = ifelse(has_dv, sum(res_df$dv_padj < 0.05, na.rm = TRUE), NA_integer_),
  n_dv_fdr10   = ifelse(has_dv, sum(res_df$dv_padj < 0.10, na.rm = TRUE), NA_integer_),
  dv_p_median  = ifelse(has_dv, median(res_df$dv_pvalue, na.rm = TRUE), NA_real_),
  n_dd_fdr5    = ifelse(has_dd, sum(res_df$padj_dd < 0.05, na.rm = TRUE), NA_integer_),
  n_dd_fdr10   = ifelse(has_dd, sum(res_df$padj_dd < 0.10, na.rm = TRUE), NA_integer_),
  n_score_fdr5 = ifelse(has_score, sum(res_df$SGCB_Score_padj < 0.05, na.rm = TRUE), NA_integer_),
  n_score_fdr10 = ifelse(has_score, sum(res_df$SGCB_Score_padj < 0.10, na.rm = TRUE), NA_integer_),
  score_median = ifelse(has_score, median(res_df$SGCB_Score, na.rm = TRUE), NA_real_)
), truth_cols)

# ===== Simulation ground truth (sim_proteomics) =====
has_sim_truth <- !is.null(sim_truth) && !IS_NULL

sim_truth_cols <- data.table(
  n_true_de = NA_integer_, n_true_dv = NA_integer_,
  sim_TP_de = NA_integer_, sim_FP_de = NA_integer_,
  sim_prec_de = NA_real_, sim_recall_de = NA_real_, sim_F1_de = NA_real_,
  sim_TP_de_post08 = NA_integer_, sim_FP_de_post08 = NA_integer_,
  sim_prec_de_post08 = NA_real_, sim_recall_de_post08 = NA_real_, sim_F1_de_post08 = NA_real_,
  sim_TP_dv = NA_integer_, sim_FP_dv = NA_integer_,
  sim_prec_dv = NA_real_, sim_recall_dv = NA_real_, sim_F1_dv = NA_real_
)

switch(as.character(has_sim_truth),
  "TRUE" = {
    mt <- sim_truth[match(res_df$gene, sim_truth$gene), ]
    td <- ifelse(is.na(mt$is_de), FALSE, mt$is_de)
    tv <- ifelse(is.na(mt$is_dv), FALSE, mt$is_dv)
    sig <- res_df$padj < 0.05
    sig_post <- ifelse(has_post, res_df$p_de_only_post > 0.8, FALSE)
    tp_d <- sum(sig & td); fp_d <- sum(sig & !td); n_d <- sum(td)
    tp_d_post <- sum(sig_post & td, na.rm = TRUE)
    fp_d_post <- sum(sig_post & !td, na.rm = TRUE)
    pr_d <- ifelse(tp_d + fp_d == 0, 0, tp_d / (tp_d + fp_d))
    rc_d <- ifelse(n_d == 0, 0, tp_d / n_d)
    f1_d <- ifelse(pr_d + rc_d == 0, 0, 2 * pr_d * rc_d / (pr_d + rc_d))
    pr_d_post <- ifelse(tp_d_post + fp_d_post == 0, 0, tp_d_post / (tp_d_post + fp_d_post))
    rc_d_post <- ifelse(n_d == 0, 0, tp_d_post / n_d)
    f1_d_post <- ifelse(pr_d_post + rc_d_post == 0, 0, 2 * pr_d_post * rc_d_post / (pr_d_post + rc_d_post))

    n_v <- sum(tv)
    tp_v <- NA_integer_; fp_v <- NA_integer_
    pr_v <- NA_real_; rc_v <- NA_real_; f1_v <- NA_real_
    switch(as.character(has_dv),
      "TRUE" = {
        sig_v <- res_df$dv_padj < 0.05
        tp_v <- sum(sig_v & tv, na.rm = TRUE)
        fp_v <- sum(sig_v & !tv, na.rm = TRUE)
        pr_v <- ifelse(tp_v + fp_v == 0, 0, tp_v / (tp_v + fp_v))
        rc_v <- ifelse(n_v == 0, 0, tp_v / n_v)
        f1_v <- ifelse(pr_v + rc_v == 0, 0, 2 * pr_v * rc_v / (pr_v + rc_v))
        NULL
      },
      "FALSE" = NULL
    )
    sim_truth_cols <- data.table(
      n_true_de = n_d, n_true_dv = n_v,
      sim_TP_de = tp_d, sim_FP_de = fp_d,
      sim_prec_de = round(pr_d, 4), sim_recall_de = round(rc_d, 4), sim_F1_de = round(f1_d, 4),
      sim_TP_de_post08 = tp_d_post, sim_FP_de_post08 = fp_d_post,
      sim_prec_de_post08 = round(pr_d_post, 4), sim_recall_de_post08 = round(rc_d_post, 4), sim_F1_de_post08 = round(f1_d_post, 4),
      sim_TP_dv = tp_v, sim_FP_dv = fp_v,
      sim_prec_dv = round(pr_v, 4), sim_recall_dv = round(rc_v, 4), sim_F1_dv = round(f1_v, 4)
    )
    cat(sprintf("  Sim truth DE: TP=%d FP=%d F1=%.4f\n", tp_d, fp_d, f1_d))
    NULL
  },
  "FALSE" = NULL
)

metrics <- cbind(metrics, sim_truth_cols)

fwrite(metrics, OUT_CSV)
saveRDS(res_df, OUT_RDS)
cat(sprintf("  Saved: %s\n\n", OUT_CSV))
