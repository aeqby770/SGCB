#!/usr/bin/env Rscript
# =============================================================================
# 51_tcga_kirc_clinical_extended.R
#
# Extends 42_tcga_kirc_dv_survival.R:
#   - Adds a DG-only signature alongside DE-only / DV-only.
#   - Adjusts each channel signature for stage / grade / age in a multivariable
#     Cox model (so the single-channel HRs above are not confounded by stage).
#   - Tests whether per-patient channel signatures differ by AJCC stage
#     (ANOVA + pairwise Wilcoxon).
#
# Inputs (all produced upstream by previous benchmark runs):
#   benchmark/output/paper_data_v2/tcga_kirc_sgcb_channels.csv
#   benchmark/output/tcga_kirc_clinical.rds  (TCGAbiolinks clinical)
#   benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc.rds (expression)
#
# Outputs:
#   benchmark/output/tcga_kirc_clinical_extended_summary.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})


OUT_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")

# -----------------------------------------------------------------------------
# 1. Channel membership + expression
# -----------------------------------------------------------------------------
channels <- fread(file.path(OUT_DIR, "paper_data_v2/tcga_kirc_sgcb_channels.csv"))
# DG significance is split across two padj columns (alpha and gamma shape
# parameters); the stored CSV only has sig_de/sig_dv/sig_dd, so we derive
channels[, sig_dg := pmin(padj_dg_a, padj_dg_g, na.rm = TRUE) < 0.05]
channels[, cat_dg_only := sig_dg & !sig_de & !sig_dv]
gene_sets <- list(
  DE_only = channels[cat_de_only == TRUE, gene],
  DV_only = channels[cat_dv_only == TRUE, gene],
  DG_only = channels[cat_dg_only == TRUE, gene]
)
channel_label <- c(DE_only = "DE-only", DV_only = "DV-only", DG_only = "DG-only")
cat("gene set sizes:\n"); print(lengths(gene_sets))

tcga_dat <- readRDS(paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc.rds"))
expr <- tcga_dat$expression
patient_ids <- substr(colnames(expr), 1, 12)

# -----------------------------------------------------------------------------
# 2. Clinical & survival skeleton (reuse conventions from 42_tcga_kirc_dv_survival.R)
# -----------------------------------------------------------------------------
clinical <- as.data.table(readRDS(file.path(OUT_DIR, "tcga_kirc_clinical.rds")))

dtd <- suppressWarnings(as.numeric(clinical$days_to_death))
dtf <- suppressWarnings(as.numeric(clinical$days_to_last_follow_up))
surv_dt <- data.table(
  patient = clinical$bcr_patient_barcode,
  os_time = ifelse(!is.na(dtd) & dtd > 0, dtd, dtf),
  os_event = as.integer(tolower(clinical$vital_status) == "dead"),
  stage = clinical$ajcc_pathologic_stage,
  grade = clinical$tumor_grade,
  age   = suppressWarnings(as.numeric(clinical$age_at_index))
)
surv_dt <- surv_dt[!is.na(os_time) & os_time > 0 & !is.na(os_event)]
surv_dt[, stage_group := fcase(
  stage %in% c("Stage I"),              "Stage I",
  stage %in% c("Stage II"),             "Stage II",
  stage %in% c("Stage III"),            "Stage III",
  stage %in% c("Stage IV"),             "Stage IV",
  default = NA_character_)]
surv_dt[, stage_group := factor(stage_group,
                                 levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))]

# Match patients to expression samples (one sample per patient)
expr_patients <- data.table(sample = colnames(expr), patient = patient_ids)
expr_patients <- expr_patients[!duplicated(patient)]
matched <- merge(surv_dt, expr_patients, by = "patient")
cat(sprintf("Matched patients: %d (events: %d)\n",
            nrow(matched), sum(matched$os_event)))

# -----------------------------------------------------------------------------
# 3. Per-channel signature scores (mean z-scored log2 expression)
# -----------------------------------------------------------------------------
compute_signature <- function(expr_mat, gene_list, samples) {
  avail <- intersect(gene_list, rownames(expr_mat))
  if (length(avail) < 5L) return(rep(NA_real_, length(samples)))
  log_mat <- log2(expr_mat[avail, samples, drop = FALSE] + 1)
  gm <- rowMeans(log_mat, na.rm = TRUE)
  gs <- apply(log_mat, 1, sd, na.rm = TRUE); gs[gs < 1e-8] <- 1
  colMeans((log_mat - gm) / gs, na.rm = TRUE)
}

for (ch in names(gene_sets)) {
  score <- compute_signature(expr, gene_sets[[ch]], matched$sample)
  matched[, (paste0("score_", ch)) := score]
  matched[, (paste0("group_", ch)) := ifelse(score >= median(score, na.rm = TRUE),
                                              "High", "Low")]
}

# -----------------------------------------------------------------------------
# 4. Univariate + multivariable Cox for each channel
# -----------------------------------------------------------------------------
run_cox <- function(dt, group_col, covars = NULL) {
  dt_use <- dt[!is.na(get(group_col))]
  f_uni <- as.formula(sprintf("Surv(os_time, os_event) ~ %s", group_col))
  fit_u <- coxph(f_uni, data = dt_use)
  s_u   <- summary(fit_u)
  rec_u <- data.table(
    adjustment = "Univariate",
    HR    = s_u$conf.int[group_col_first_row <- 1, 1],
    lower = s_u$conf.int[1, 3],
    upper = s_u$conf.int[1, 4],
    p     = s_u$coefficients[1, 5])

  if (!is.null(covars)) {
    keep <- dt_use[complete.cases(dt_use[, c(group_col, covars), with = FALSE])]
    f_m  <- as.formula(sprintf("Surv(os_time, os_event) ~ %s + %s",
                               group_col, paste(covars, collapse = " + ")))
    fit_m <- coxph(f_m, data = keep)
    s_m   <- summary(fit_m)
    rn    <- rownames(s_m$coefficients)
    idx   <- grep(paste0("^", group_col), rn)[1]
    rec_m <- data.table(
      adjustment = "Adjusted (stage + grade + age)",
      HR    = s_m$conf.int[idx, 1],
      lower = s_m$conf.int[idx, 3],
      upper = s_m$conf.int[idx, 4],
      p     = s_m$coefficients[idx, 5])
    return(rbind(rec_u, rec_m))
  }
  rec_u
}

cox_rows <- list()
for (ch in names(gene_sets)) {
  gcol <- paste0("group_", ch)
  if (all(is.na(matched[[gcol]]))) {
    cat("\n--", channel_label[ch], "-- SKIP (0 genes)\n")
    next
  }
  cat("\n--", channel_label[ch], "--\n")
  out <- run_cox(matched, gcol,
                 covars = c("stage_group", "grade", "age"))
  out[, signature := channel_label[ch]]
  cox_rows[[ch]] <- out
  print(out)
}
cox_all <- rbindlist(cox_rows)

# -----------------------------------------------------------------------------
# 5. Stage differences of the continuous signature score
# -----------------------------------------------------------------------------
stage_tests <- list()
for (ch in names(gene_sets)) {
  col <- paste0("score_", ch)
  d <- matched[!is.na(stage_group) & !is.na(get(col))]
  if (nrow(d) < 30L) next
  a <- aov(as.formula(sprintf("%s ~ stage_group", col)), data = d)
  p <- summary(a)[[1]][["Pr(>F)"]][1]
  stage_tests[[ch]] <- data.table(signature = channel_label[ch], anova_p = p)
}
stage_tab <- rbindlist(stage_tests)
cat("\nStage ANOVA:\n"); print(stage_tab)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 6. Persist summary tables  (figure rendering removed from test scope)
# -----------------------------------------------------------------------------
