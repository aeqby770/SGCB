#!/usr/bin/env Rscript
# =============================================================================
# 42_tcga_kirc_dv_survival.R
#
# Strategy:
#   1. Load SGCB channel classification (DE-only / DV-only / Both)
#   2. Download TCGA-KIRC clinical data via TCGAbiolinks
#   3. For DV-only genes: compute per-patient CV-based signature score
#   5. Compare with DE-only gene signature as control
#
# Output:
#   - benchmark/output/tcga_kirc_dv_survival.csv  (per-patient scores + survival)
#   - benchmark/output/tcga_kirc_dv_survival_summary.csv  (KM/Cox results)
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(TCGAbiolinks)
})

OUT_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")
CLINICAL_CACHE <- file.path(OUT_DIR, "tcga_kirc_clinical.rds")

# =============================================================================
# 1. Load SGCB channel classification
# =============================================================================
channels <- fread(paste0(PROJECT_ROOT, "/benchmark/output/paper_data_v2/tcga_kirc_sgcb_channels.csv"))
cat(sprintf("Loaded %d genes, %d DE-only, %d DV-only, %d Both\n",
            nrow(channels),
            sum(channels$cat_de_only),
            sum(channels$cat_dv_only),
            sum(channels$cat_both)))

dv_only_genes <- channels[cat_dv_only == TRUE, gene]
de_only_genes <- channels[cat_de_only == TRUE, gene]

# =============================================================================
# 2. Load TCGA-KIRC expression (tumor samples)
# =============================================================================
tcga_dat <- readRDS(paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc.rds"))
expr <- tcga_dat$expression
cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# Extract patient barcodes (first 12 chars of TCGA barcode)
patient_ids <- substr(colnames(expr), 1, 12)

# =============================================================================
# 3. Download/load clinical data
# =============================================================================
if (file.exists(CLINICAL_CACHE)) {
  cat("Loading cached clinical data...\n")
  clinical <- readRDS(CLINICAL_CACHE)
} else {
  cat("Downloading TCGA-KIRC clinical data via GDCquery_clinic()...\n")
  clinical_raw <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")
  clinical <- as.data.table(clinical_raw)
  saveRDS(clinical, CLINICAL_CACHE)
  cat(sprintf("Clinical data: %d patients, columns: %s\n",
              nrow(clinical), paste(head(names(clinical), 15), collapse = ", ")))
}

# Extract survival info
# Try multiple column name formats (TCGA clinical data can vary)
os_time_col <- intersect(c("days_to_death", "days_to_last_followup",
                            "OS.time", "os_time"), names(clinical))
os_status_col <- intersect(c("vital_status", "OS", "os_status"), names(clinical))
patient_col <- intersect(c("bcr_patient_barcode", "patient", "submitter_id"), names(clinical))

cat("Clinical columns found:\n")
cat("  time:", paste(os_time_col, collapse = ", "), "\n")
cat("  status:", paste(os_status_col, collapse = ", "), "\n")
cat("  patient:", paste(patient_col, collapse = ", "), "\n")

# Build survival table
if (length(patient_col) == 0) {
  stop("Cannot find patient barcode column in clinical data. Columns: ",
       paste(head(names(clinical), 20), collapse = ", "))
}

surv_dt <- data.table(patient = clinical[[patient_col[1]]])

# OS time: days_to_death for dead, days_to_last_follow_up for alive
# GDCquery_clinic uses days_to_last_follow_up (with underscores)
cat("  Available time cols:", paste(grep("days_to", names(clinical), value=TRUE), collapse=", "), "\n")
dtd <- suppressWarnings(as.numeric(clinical[["days_to_death"]]))
followup_col <- intersect(c("days_to_last_follow_up", "days_to_last_followup",
                            "days_to_last_known_alive"), names(clinical))
if (length(followup_col) > 0) {
  dtf <- suppressWarnings(as.numeric(clinical[[followup_col[1]]]))
} else {
  dtf <- rep(NA_real_, nrow(clinical))
}
surv_dt[, os_time := ifelse(!is.na(dtd) & dtd > 0, dtd, dtf)]

# OS status: 1 = dead, 0 = alive
if ("vital_status" %in% names(clinical)) {
  vs <- tolower(clinical$vital_status)
  surv_dt[, os_event := as.integer(vs == "dead")]
} else if ("OS" %in% names(clinical)) {
  surv_dt[, os_event := as.integer(clinical[["OS"]])]
} else {
  surv_dt[, os_event := suppressWarnings(as.integer(clinical[[os_status_col[1]]]))]
}

surv_dt <- surv_dt[!is.na(os_time) & os_time > 0 & !is.na(os_event)]
cat(sprintf("Survival data: %d patients (events: %d)\n",
            nrow(surv_dt), sum(surv_dt$os_event)))

# =============================================================================
# 4. Compute per-patient signature scores
# =============================================================================
# Match patients between expression and clinical
expr_patients <- data.table(sample = colnames(expr), patient = patient_ids)
# Some patients may have multiple samples; keep first
expr_patients <- expr_patients[!duplicated(patient)]

matched <- merge(surv_dt, expr_patients, by = "patient")
cat(sprintf("Matched patients: %d\n", nrow(matched)))

if (nrow(matched) < 30) {
  stop("Too few matched patients (", nrow(matched), "). Check barcode format.")
}

compute_signature <- function(expr_mat, gene_list, samples) {
  avail <- intersect(gene_list, rownames(expr_mat))
  if (length(avail) < 5) {
    warning("Only ", length(avail), " genes available")
    return(rep(NA_real_, length(samples)))
  }
  sub_mat <- expr_mat[avail, samples, drop = FALSE]
  # Signature = mean of z-scored log2(count+1) across signature genes
  log_mat <- log2(sub_mat + 1)
  # Z-score each gene across patients
  gene_means <- rowMeans(log_mat, na.rm = TRUE)
  gene_sds <- apply(log_mat, 1, sd, na.rm = TRUE)
  gene_sds[gene_sds < 1e-8] <- 1
  z_mat <- (log_mat - gene_means) / gene_sds
  # Per-patient: mean z-score
  colMeans(z_mat, na.rm = TRUE)
}

matched[, dv_score := compute_signature(expr, dv_only_genes, matched$sample)]
matched[, de_score := compute_signature(expr, de_only_genes, matched$sample)]

# Median split
matched[, dv_group := ifelse(dv_score >= median(dv_score, na.rm = TRUE), "High", "Low")]
matched[, de_group := ifelse(de_score >= median(de_score, na.rm = TRUE), "High", "Low")]

# =============================================================================
# 5. Kaplan-Meier + Cox PH analysis
# =============================================================================
run_survival <- function(dt, group_col, label) {
  dt_use <- dt[!is.na(get(group_col))]
  formula_km <- as.formula(paste0("Surv(os_time, os_event) ~ ", group_col))

  # KM
  fit_km <- survfit(formula_km, data = dt_use)

  # Log-rank test
  lr <- survdiff(formula_km, data = dt_use)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  # Cox PH
  fit_cox <- coxph(formula_km, data = dt_use)
  cox_summary <- summary(fit_cox)
  hr <- cox_summary$conf.int[1, 1]
  hr_lower <- cox_summary$conf.int[1, 3]
  hr_upper <- cox_summary$conf.int[1, 4]
  cox_p <- cox_summary$coefficients[1, 5]

  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Log-rank p = %.2e\n", lr_p))
  cat(sprintf("  Cox HR = %.3f (%.3f - %.3f), p = %.2e\n",
              hr, hr_lower, hr_upper, cox_p))

  # Save KM plot
  plot_file <- file.path(OUT_DIR, paste0("tcga_kirc_survival_", tolower(gsub("[^a-zA-Z]", "_", label)), ".pdf"))
  tryCatch({
    p <- ggsurvplot(fit_km, data = dt_use,
                    pval = TRUE, conf.int = TRUE,
                    risk.table = TRUE,
                    title = label,
                    xlab = "Days", ylab = "Overall Survival")
    pdf(plot_file, width = 8, height = 6)
    print(p)
    dev.off()
    cat(sprintf("  Plot saved: %s\n", plot_file))
  }, error = function(e) {
    cat(sprintf("  Plot error: %s\n", conditionMessage(e)))
    # Fallback: base R plot
    pdf(plot_file, width = 8, height = 6)
    plot(fit_km, col = c("blue", "red"), lwd = 2,
         main = label, xlab = "Days", ylab = "Overall Survival")
    legend("topright", legend = names(fit_km$strata),
           col = c("blue", "red"), lwd = 2)
    dev.off()
    cat(sprintf("  Fallback plot saved: %s\n", plot_file))
  })

  data.table(
    signature = label,
    logrank_p = lr_p,
    cox_hr = hr, cox_hr_lower = hr_lower, cox_hr_upper = hr_upper,
    cox_p = cox_p,
    n_high = sum(dt_use[[group_col]] == "High"),
    n_low = sum(dt_use[[group_col]] == "Low"),
    n_events = sum(dt_use$os_event)
  )
}

results <- list()
results[[1]] <- run_survival(matched, "dv_group", "DV-only genes")
results[[2]] <- run_survival(matched, "de_group", "DE-only genes")

summary_dt <- rbindlist(results)

# Save
fwrite(matched[, .(patient, os_time, os_event, dv_score, de_score, dv_group, de_group)],
       file.path(OUT_DIR, "tcga_kirc_dv_survival.csv"))
fwrite(summary_dt, file.path(OUT_DIR, "tcga_kirc_dv_survival_summary.csv"))

cat("\n=== TCGA-KIRC Survival Summary ===\n")
print(summary_dt)
cat("\nDone.\n")
