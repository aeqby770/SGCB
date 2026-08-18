# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)

BASE <- paste0(PROJECT_ROOT, "/benchmark_v2")

cat("===== Bulk SGCB =====\n")
bulk_sgcb_dir <- file.path(BASE, "01_bulk_rnaseq/sgcb/results")
bulk_csvs     <- list.files(bulk_sgcb_dir, pattern = "\\.csv$", full.names = TRUE)
bulk_csvs     <- bulk_csvs[!grepl("_combined", bulk_csvs)]

cat(sprintf(" %d CSV\n", length(bulk_csvs)))

classic_csvs    <- bulk_csvs[grepl("^classic__", basename(bulk_csvs))]
simulation_csvs <- bulk_csvs[grepl("^simulation__", basename(bulk_csvs))]
large_csvs      <- bulk_csvs[grepl("^(gtex_pair|tcga|pseudobulk)__", basename(bulk_csvs))]
null_csvs       <- bulk_csvs[grepl("^null_fdr__", basename(bulk_csvs))]
seqc_csvs       <- bulk_csvs[grepl("^seqc__", basename(bulk_csvs))]
dv_csvs         <- bulk_csvs[grepl("^dv__", basename(bulk_csvs))]

sgcb_classic <- rbindlist(lapply(classic_csvs, fread), fill = TRUE)

sgcb_sim <- rbindlist(lapply(simulation_csvs, fread), fill = TRUE)

sgcb_large <- rbindlist(lapply(large_csvs, fread), fill = TRUE)

sgcb_seqc <- rbindlist(lapply(seqc_csvs, fread), fill = TRUE)

sgcb_dv <- rbindlist(lapply(dv_csvs, fread), fill = TRUE)

sgcb_null_list <- rbindlist(lapply(null_csvs, fread), fill = TRUE)
sgcb_null <- data.table(
  method     = "SGCB",
  mean_FDR01 = mean(sgcb_null_list$n_sig_fdr1 / pmax(sgcb_null_list$n_genes, 1)),
  sd_FDR01   = sd(sgcb_null_list$n_sig_fdr1 / pmax(sgcb_null_list$n_genes, 1)),
  mean_FDR05 = mean(sgcb_null_list$n_sig_fdr5 / pmax(sgcb_null_list$n_genes, 1)),
  sd_FDR05   = sd(sgcb_null_list$n_sig_fdr5 / pmax(sgcb_null_list$n_genes, 1)),
  mean_FDR10 = mean(sgcb_null_list$n_sig_fdr10 / pmax(sgcb_null_list$n_genes, 1)),
  sd_FDR10   = sd(sgcb_null_list$n_sig_fdr10 / pmax(sgcb_null_list$n_genes, 1))
)

cat("===== other_methods =====\n")
om_dir <- file.path(BASE, "01_bulk_rnaseq/other_methods/results")

drop_diphiseq <- function(dt) dt[method != "dv_DiPhiSeq"]

# Classic
om_classic <- fread(file.path(om_dir, "bulk_classic_metrics.csv"))
om_classic <- drop_diphiseq(om_classic)
combined_classic <- rbindlist(list(sgcb_classic, om_classic), fill = TRUE)
fwrite(combined_classic,
       file.path(bulk_sgcb_dir, "combined_classic.csv"))
cat(sprintf("  combined_classic: %d rows (%d methods x %d datasets)\n",
            nrow(combined_classic), uniqueN(combined_classic$method),
            uniqueN(combined_classic$dataset)))

# Simulation
om_sim <- fread(file.path(om_dir, "bulk_simulation_metrics.csv"))
om_sim <- drop_diphiseq(om_sim)
combined_sim <- rbindlist(list(sgcb_sim, om_sim), fill = TRUE)
fwrite(combined_sim,
       file.path(bulk_sgcb_dir, "combined_simulation.csv"))
cat(sprintf("  combined_simulation: %d rows\n", nrow(combined_sim)))

# Large-scale
om_large <- fread(file.path(om_dir, "bulk_large_scale_metrics.csv"))
om_large <- drop_diphiseq(om_large)
combined_large <- rbindlist(list(sgcb_large, om_large), fill = TRUE)
fwrite(combined_large,
       file.path(bulk_sgcb_dir, "combined_large_scale.csv"))
cat(sprintf("  combined_large_scale: %d rows\n", nrow(combined_large)))

# Null FDR
om_null <- fread(file.path(om_dir, "bulk_null_fdr.csv"))
om_null <- drop_diphiseq(om_null)
combined_null <- rbindlist(list(sgcb_null, om_null), fill = TRUE)
fwrite(combined_null,
       file.path(bulk_sgcb_dir, "combined_null_fdr.csv"))
cat(sprintf("  combined_null_fdr: %d rows\n", nrow(combined_null)))

# SEQC
fwrite(sgcb_seqc, file.path(bulk_sgcb_dir, "combined_seqc.csv"))
cat(sprintf("  SEQC: %d rows\n", nrow(sgcb_seqc)))

# DV
fwrite(sgcb_dv, file.path(bulk_sgcb_dir, "combined_dv.csv"))
cat(sprintf("  DV: %d rows\n", nrow(sgcb_dv)))

cat("\n===== scRNA () =====\n")
scrna_out_dir <- file.path(BASE, "02_scrna_pseudobulk/results")
scrna_dirs <- c(
  file.path(BASE, "02_scrna_pseudobulk/results"),
  file.path(BASE, "02_scrna_pseudobulk/sgcb/results"),
  file.path(BASE, "02_scrna_pseudobulk/other_methods/results")
)
scrna_csvs <- unique(unlist(lapply(scrna_dirs, function(d) {
  if (dir.exists(d)) list.files(d, pattern = "\\.csv$", full.names = TRUE) else character(0)
})))
scrna_csvs <- scrna_csvs[!grepl("combined", scrna_csvs)]
scrna_csvs <- scrna_csvs[!grepl("DiPhiSeq", scrna_csvs)]
cat(sprintf("  Found %d scRNA CSVs across all dirs\n", length(scrna_csvs)))

scrna_null_csvs   <- scrna_csvs[grepl("^null_", basename(scrna_csvs))]
scrna_normal_csvs <- scrna_csvs[!grepl("^null_", basename(scrna_csvs))]

combined_scrna <- rbindlist(lapply(scrna_normal_csvs, fread), fill = TRUE)
combined_scrna[, cohort := sub("^PB_", "", dataset)]
combined_scrna[, is_pb := grepl("^PB_", dataset)]
setorder(combined_scrna, method, cohort, is_pb)
combined_scrna <- unique(combined_scrna, by = c("method", "cohort"), fromLast = FALSE)
combined_scrna[, dataset := cohort]
combined_scrna[, c("cohort", "is_pb") := NULL]
fwrite(combined_scrna, file.path(scrna_out_dir, "combined_scrna_de.csv"))
cat(sprintf("  combined_scrna_de: %d rows (%d methods x %d datasets)\n",
            nrow(combined_scrna), uniqueN(combined_scrna$method),
            uniqueN(combined_scrna$dataset)))
cat("  Methods:", paste(sort(unique(combined_scrna$method)), collapse = ", "), "\n")

scrna_null_all <- rbindlist(lapply(scrna_null_csvs, function(f) {
  x <- fread(f)
  x[, null_seed := sub("^null_([0-9]+)__.*", "\\1", basename(f))]
  x
}), fill = TRUE)
scrna_null_all[, cohort := sub("^PB_", "", dataset)]
scrna_null_all[, is_pb := grepl("^PB_", dataset)]
setorder(scrna_null_all, method, cohort, null_seed, is_pb)
scrna_null_all <- unique(scrna_null_all, by = c("method", "cohort", "null_seed"), fromLast = FALSE)
scrna_null_all[, dataset := cohort]
scrna_null_all[, c("cohort", "is_pb", "null_seed") := NULL]
if ("n_sig_fdr1" %in% names(scrna_null_all)) {
  scrna_null_summary <- scrna_null_all[, .(
    mean_FDR01 = mean(n_sig_fdr1 / pmax(n_genes, 1), na.rm = TRUE),
    sd_FDR01   = sd(n_sig_fdr1 / pmax(n_genes, 1), na.rm = TRUE),
    mean_FDR05 = mean(n_sig_fdr5 / pmax(n_genes, 1), na.rm = TRUE),
    sd_FDR05   = sd(n_sig_fdr5 / pmax(n_genes, 1), na.rm = TRUE),
    mean_FDR10 = mean(n_sig_fdr10 / pmax(n_genes, 1), na.rm = TRUE),
    sd_FDR10   = sd(n_sig_fdr10 / pmax(n_genes, 1), na.rm = TRUE),
    n_runs     = .N
  ), by = .(method, dataset)]
} else {
  scrna_null_summary <- data.table()
}
fwrite(scrna_null_summary, file.path(scrna_out_dir, "combined_scrna_null_fdr.csv"))
cat(sprintf("  combined_scrna_null_fdr: %d rows (%d methods)\n",
            nrow(scrna_null_summary), uniqueN(scrna_null_summary$method)))

cat("\n===== Proteomics () =====\n")
prot_out_dir <- file.path(BASE, "03_proteomics/results")
prot_dirs <- c(
  file.path(BASE, "03_proteomics/results"),
  file.path(BASE, "03_proteomics/sgcb/results"),
  file.path(BASE, "03_proteomics/other_methods/results")
)
prot_csvs <- unique(unlist(lapply(prot_dirs, function(d) {
  if (dir.exists(d)) list.files(d, pattern = "\\.csv$", full.names = TRUE) else character(0)
})))
prot_csvs <- prot_csvs[!grepl("combined", prot_csvs)]
prot_csvs <- prot_csvs[!grepl("DiPhiSeq", prot_csvs)]
cat(sprintf("  Found %d proteomics CSVs across all dirs\n", length(prot_csvs)))

prot_null_csvs   <- prot_csvs[grepl("^null_", basename(prot_csvs))]
prot_normal_csvs <- prot_csvs[!grepl("^null_", basename(prot_csvs))]

combined_prot <- rbindlist(lapply(prot_normal_csvs, fread), fill = TRUE)
combined_prot <- unique(combined_prot, by = c("method", "dataset"), fromLast = TRUE)
fwrite(combined_prot, file.path(prot_out_dir, "combined_proteomics_de.csv"))
cat(sprintf("  combined_proteomics_de: %d rows (%d methods x %d datasets)\n",
            nrow(combined_prot), uniqueN(combined_prot$method),
            uniqueN(combined_prot$dataset)))
cat("  Methods:", paste(sort(unique(combined_prot$method)), collapse = ", "), "\n")

prot_null_all <- rbindlist(lapply(prot_null_csvs, fread), fill = TRUE)
n_feat_col <- intersect(c("n_features", "n_genes"), names(prot_null_all))
if (length(n_feat_col) > 0L && "n_sig_fdr1" %in% names(prot_null_all)) {
  if ("n_features" %in% names(prot_null_all) && "n_genes" %in% names(prot_null_all)) {
    prot_null_all[is.na(n_features), n_features := n_genes]
  } else if ("n_genes" %in% names(prot_null_all) && !("n_features" %in% names(prot_null_all))) {
    prot_null_all[, n_features := n_genes]
  }
  prot_null_summary <- prot_null_all[, .(
    mean_FDR01 = mean(n_sig_fdr1 / pmax(n_features, 1), na.rm = TRUE),
    sd_FDR01   = sd(n_sig_fdr1 / pmax(n_features, 1), na.rm = TRUE),
    mean_FDR05 = mean(n_sig_fdr5 / pmax(n_features, 1), na.rm = TRUE),
    sd_FDR05   = sd(n_sig_fdr5 / pmax(n_features, 1), na.rm = TRUE),
    mean_FDR10 = mean(n_sig_fdr10 / pmax(n_features, 1), na.rm = TRUE),
    sd_FDR10   = sd(n_sig_fdr10 / pmax(n_features, 1), na.rm = TRUE),
    n_runs     = .N
  ), by = .(method, dataset)]
} else {
  prot_null_summary <- data.table()
}
fwrite(prot_null_summary, file.path(prot_out_dir, "combined_proteomics_null_fdr.csv"))
cat(sprintf("  combined_proteomics_null_fdr: %d rows (%d methods)\n",
            nrow(prot_null_summary), uniqueN(prot_null_summary$method)))

cat("\n===== =====\n")
