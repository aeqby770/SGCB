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
#   - Produces one compact figure (forest + KM + stage boxplot) for §H.
#
# Inputs:
#   data/paper_data_v2/tcga_kirc_sgcb_channels.csv
#   data/tcga_kirc_cache/tcga_kirc_clinical.rds
#   data/tcga_kirc_cache/tcga_kirc.rds
#
# Outputs are written to SGCB_FIGURE_OUT.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(PROJECT_ROOT, "scripts", "figures", "37_theme_v3.R"),
       encoding = "UTF-8")

OUT_DIR <- Sys.getenv(
  "SGCB_FIGURE_OUT",
  unset = file.path(PROJECT_ROOT, "results", "regenerated_figures")
)
FIG_DIR <- OUT_DIR
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Channel membership + expression
# -----------------------------------------------------------------------------
channels <- fread(file.path(PROJECT_ROOT, "data", "paper_data_v2",
                            "tcga_kirc_sgcb_channels.csv"))
# DG significance is split across two padj columns (alpha and gamma shape
# parameters); the stored CSV only has sig_de/sig_dv/sig_dd, so we derive
# sig_dg here the same way §H uses it.
channels[, sig_dg := pmin(padj_dg_a, padj_dg_g, na.rm = TRUE) < 0.05]
channels[, cat_dg_only := sig_dg & !sig_de & !sig_dv]
gene_sets <- list(
  DE_only = channels[cat_de_only == TRUE, gene],
  DV_only = channels[cat_dv_only == TRUE, gene],
  DG_only = channels[cat_dg_only == TRUE, gene]
)
channel_label <- c(DE_only = "DE-only", DV_only = "DV-only", DG_only = "DG-only")
cat("gene set sizes:\n"); print(lengths(gene_sets))

tcga_dat <- readRDS(file.path(PROJECT_ROOT, "data", "tcga_kirc_cache",
                              "tcga_kirc.rds"))
expr <- tcga_dat$expression
patient_ids <- substr(colnames(expr), 1, 12)

# -----------------------------------------------------------------------------
# 2. Clinical & survival skeleton (reuse conventions from 42_tcga_kirc_dv_survival.R)
# -----------------------------------------------------------------------------
clinical <- as.data.table(readRDS(file.path(PROJECT_ROOT, "data", "tcga_kirc_cache",
                                           "tcga_kirc_clinical.rds")))

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
# 6. Figures (forest / KM / stage boxplot) and combined panel
# -----------------------------------------------------------------------------
cox_plot_dt <- copy(cox_all)
cox_plot_dt[, signature := factor(signature, levels = rev(unname(channel_label)))]
cox_plot_dt[, adjustment := factor(adjustment,
  levels = c("Univariate", "Adjusted (stage + grade + age)"))]
cox_plot_dt[, lab := sprintf("HR = %.2f (%.2f\u2013%.2f), p = %s",
                              HR, lower, upper,
                              formatC(p, format = "g", digits = 2))]

label_x <- max(cox_plot_dt$upper, na.rm = TRUE) * 1.15
x_hi    <- max(cox_plot_dt$upper, na.rm = TRUE) * 4.5

p_forest <- ggplot(cox_plot_dt,
                   aes(x = HR, y = signature, color = adjustment)) +
  geom_vline(xintercept = 1, linetype = "21", color = "grey55", linewidth = 0.35) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.6),
                  size = 0.6, linewidth = 0.8) +
  geom_text(aes(x = label_x, label = lab),
            position = position_dodge(width = 0.6),
            hjust = 0, size = 3.0, show.legend = FALSE) +
  scale_color_manual(values = c("Univariate" = "#3B4CC0",
                                 "Adjusted (stage + grade + age)" = "#D4A017"),
                     name = NULL) +
  scale_x_continuous(trans = "log2",
                     limits = c(0.5, x_hi),
                     breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(title = "**A.** Cox HR per SGCB channel (High vs Low signature)",
       x = "Hazard ratio (log\u2082 scale)", y = NULL) +
  THEME_JBHI() +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())

# KM curve for DV-only (the signal channel here)
km_fit <- survfit(Surv(os_time, os_event) ~ group_DV_only, data = matched)
km_plot <- ggsurvplot(km_fit, data = matched,
                      palette = c("#E07B5E", "#3B4CC0"),
                      legend.title = "DV-only signature", legend.labs = c("High", "Low"),
                      pval = TRUE, pval.size = 3.2, pval.coord = c(100, 0.08),
                      risk.table = FALSE, conf.int = TRUE, conf.int.alpha = 0.12,
                      censor.size = 2, size = 0.65,
                      xlab = "Time (days)", ylab = "Overall survival",
                      ggtheme = THEME_JBHI())
p_km <- km_plot$plot +
  labs(title = "**B.** Kaplan\u2013Meier: DV-only signature")

# Stage boxplot for DV-only
stage_box_dt <- matched[!is.na(stage_group)]
p_stage <- ggplot(stage_box_dt,
                  aes(x = stage_group, y = score_DV_only, fill = stage_group)) +
  geom_violin(alpha = 0.35, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.22, outlier.size = 0.6, outlier.alpha = 0.4,
               color = "grey25", linewidth = 0.3) +
  scale_fill_manual(values = c("Stage I" = "#4DBBD5", "Stage II" = "#7CC68E",
                                "Stage III" = "#E07B5E", "Stage IV" = "#B40426"),
                    guide = "none") +
  annotate("text", x = 0.6, y = max(stage_box_dt$score_DV_only, na.rm = TRUE),
           hjust = 0, vjust = 1, size = 3,
           label = sprintf("ANOVA p = %s",
                           formatC(stage_tab[signature == "DV-only", anova_p],
                                    format = "g", digits = 2))) +
  labs(title = "**C.** DV-only signature by AJCC stage",
       x = NULL, y = "DV-only signature score") +
  THEME_JBHI() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

fig_combined <- p_forest / (p_km | p_stage) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "TCGA-KIRC survival association of SGCB channel signatures",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

out_png <- file.path(FIG_DIR, "FigS5_tcga_kirc_clinical.png")
ggsave(out_png, fig_combined, width = 14, height = 9.5, dpi = 300, bg = "white")
cat("Saved: ", out_png, "\n")

# -----------------------------------------------------------------------------
# 7. Persist summary tables
# -----------------------------------------------------------------------------
fwrite(cox_all[, .(signature, adjustment, HR, lower, upper, p)],
       file.path(OUT_DIR, "tcga_kirc_clinical_extended_summary.csv"))
fwrite(stage_tab,
       file.path(OUT_DIR, "tcga_kirc_stage_anova.csv"))
cat("Done.\n")
