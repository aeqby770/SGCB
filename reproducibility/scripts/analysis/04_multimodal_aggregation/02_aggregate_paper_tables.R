# =============================================================================
# 50_paper_data_aggregate.R
#
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
suppressPackageStartupMessages({
  library(data.table)
})

CFG <- list(
  V2_BULK_SGCB       = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/sgcb/results"),
  V2_BULK_OTHER      = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/other_methods/results"),
  V2_SCRNA_RESULTS   = paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk/results"),
  V2_PROT_RESULTS    = paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics/results"),
  LEGACY_OTHER       = paste0(PROJECT_ROOT, "/benchmark/other_methods_results"),
  LEGACY_PAPER_FIGS  = paste0(PROJECT_ROOT, "/benchmark/output/paper_figures"),
  OUT_DIR            = paste0(PROJECT_ROOT, "/benchmark/output/paper_data_v2")
)
dir.create(CFG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

manifest <- data.table(
  table_name = character(0), n_rows = integer(0), n_cols = integer(0),
  n_methods = integer(0), n_datasets = integer(0),
  methods = character(0), datasets = character(0), note = character(0)
)

add_manifest <- function(tbl, name, note = "") {
  m <- if ("method" %in% names(tbl)) unique(as.character(tbl$method)) else character(0)
  d <- if ("dataset" %in% names(tbl)) unique(as.character(tbl$dataset)) else character(0)
  manifest <<- rbind(manifest, data.table(
    table_name = name,
    n_rows     = nrow(tbl),
    n_cols     = ncol(tbl),
    n_methods  = length(m),
    n_datasets = length(d),
    methods    = paste(sort(m), collapse = "|"),
    datasets   = paste(sort(d), collapse = "|"),
    note       = note
  ))
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
bulk_classic <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_classic.csv"))
fwrite(bulk_classic, file.path(CFG$OUT_DIR, "bulk_classic.csv"))
add_manifest(bulk_classic, "bulk_classic.csv", "v2 ps1: 8 methods x 6 classic datasets")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
bulk_sim <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_simulation.csv"))
fwrite(bulk_sim, file.path(CFG$OUT_DIR, "bulk_simulation.csv"))
add_manifest(bulk_sim, "bulk_simulation.csv", "v2 ps1: 8 methods x 6 NB scenarios")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
bulk_large <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_large_scale.csv"))
tcga_lower_mask <- grepl("^TCGA_[a-z]+$", bulk_large$dataset)
bulk_large[tcga_lower_mask, dataset := sub("^TCGA_", "TCGA_",
                                            sub("^(TCGA_)(.*)$", "\\1\\U\\2", dataset, perl = TRUE))]
bulk_large <- unique(bulk_large, by = c("method", "dataset"), fromLast = TRUE)
fwrite(bulk_large, file.path(CFG$OUT_DIR, "bulk_large_scale.csv"))
add_manifest(bulk_large, "bulk_large_scale.csv",
 "v2 ps1: 8 methods x 16 datasets (5 GTEx + 4 PB + 7 TCGA); TCGA ")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
bulk_null_sum <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_null_fdr.csv"))
fwrite(bulk_null_sum, file.path(CFG$OUT_DIR, "bulk_null_fdr_summary.csv"))
add_manifest(bulk_null_sum, "bulk_null_fdr_summary.csv",
             "v2 ps1: 8 methods mean/sd FDR at 1/5/10%")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
gather_null_per_replicate <- function() {
  out <- list()
  for (i in 1:10) {
    p <- file.path(CFG$V2_BULK_SGCB, sprintf("null_fdr__%d.csv", i))
    if (file.exists(p)) {
      dt <- fread(file = p)
      dt[, replicate := i]
      out[[length(out) + 1L]] <- dt
    }
  }
  for (m in c("DESeq2", "EBSeq", "NOISeq", "edgeR", "glmGamPoi", "limma", "samr")) {
    for (i in 1:10) {
      p <- file.path(CFG$V2_BULK_OTHER, sprintf("%s__null_fdr__%d.csv", m, i))
      if (file.exists(p)) {
        dt <- fread(file = p)
        dt[, replicate := i]
        out[[length(out) + 1L]] <- dt
      }
    }
  }
  if (length(out) == 0L) return(data.table())
  rbindlist(out, fill = TRUE, use.names = TRUE)
}
bulk_null_rep <- gather_null_per_replicate()
fwrite(bulk_null_rep, file.path(CFG$OUT_DIR, "bulk_null_fdr.csv"))
add_manifest(bulk_null_rep, "bulk_null_fdr.csv",
             "8 methods x 10 null replicates (per-replicate)")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
dv_sgcb <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_dv.csv"))
dv_sgcb[, dv_method := "SGCB"]
dv_base_files <- list.files(
  CFG$V2_BULK_OTHER,
  pattern = "^dv_(GAMLSS|MDSeq|clrDV|diffVar)__dv__.*\\.csv$",
  full.names = TRUE
)
dv_baselines <- if (length(dv_base_files) > 0L) {
  rbindlist(lapply(dv_base_files, function(p) {
    dt <- fread(file = p)
    bname <- tools::file_path_sans_ext(basename(p))
    parts <- strsplit(bname, "__", fixed = TRUE)[[1]]
    dt[, dv_method := sub("^dv_", "", parts[1])]
    dt[, dataset := parts[3]]
    dt
  }), fill = TRUE, use.names = TRUE)
} else data.table()
fwrite(dv_sgcb, file.path(CFG$OUT_DIR, "bulk_dv_sgcb.csv"))
add_manifest(dv_sgcb, "bulk_dv_sgcb.csv",
             "v2 ps1: SGCB DV channel on 3 GTEx brain region pairs")
if (nrow(dv_baselines) > 0L) {
  fwrite(dv_baselines, file.path(CFG$OUT_DIR, "bulk_dv_baselines.csv"))
  dv_for_mani <- copy(dv_baselines)[, method := dv_method]
  add_manifest(dv_for_mani, "bulk_dv_baselines.csv",
 "v2 bulk other_methods: GAMLSS/MDSeq/clrDV/diffVar 3 brain pairs; DiPhiSeq excluded")
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
seqc_sgcb <- fread(file = file.path(CFG$V2_BULK_SGCB, "combined_seqc.csv"))
seqc_sgcb_min <- seqc_sgcb[, .(
  method, n_genes, time_sec,
  AUC = AUC_ROC, F1 = F1_fdr5, recall = recall_fdr5, precision = precision_fdr5
)]
seqc_legacy_file <- file.path(CFG$LEGACY_PAPER_FIGS, "test4_seqc_metrics_8methods.csv")
if (file.exists(seqc_legacy_file)) {
  seqc_legacy <- fread(file = seqc_legacy_file)
  seqc_legacy_other <- seqc_legacy[method != "SGCB"]
  bulk_seqc <- rbind(seqc_sgcb_min, seqc_legacy_other, fill = TRUE)
} else {
  bulk_seqc <- seqc_sgcb_min
}
fwrite(bulk_seqc, file.path(CFG$OUT_DIR, "bulk_seqc.csv"))
add_manifest(bulk_seqc, "bulk_seqc.csv",
             "SGCB from v2 ps1; 7 others from legacy test4_seqc_metrics_8methods.csv")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
pb_de   <- fread(file = file.path(CFG$V2_SCRNA_RESULTS, "combined_scrna_de.csv"))
pb_null <- fread(file = file.path(CFG$V2_SCRNA_RESULTS, "combined_scrna_null_fdr.csv"))
fwrite(pb_de,   file.path(CFG$OUT_DIR, "pseudobulk_de.csv"))
fwrite(pb_null, file.path(CFG$OUT_DIR, "pseudobulk_null.csv"))
add_manifest(pb_de,   "pseudobulk_de.csv",
             sprintf("v2: %d methods x %d datasets (all methods, DiPhiSeq excluded)",
                     uniqueN(pb_de$method), uniqueN(pb_de$dataset)))
add_manifest(pb_null, "pseudobulk_null.csv",
             sprintf("v2: %d methods x %d datasets null FDR summary",
                     uniqueN(pb_null$method), uniqueN(pb_null$dataset)))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
prot_de   <- fread(file = file.path(CFG$V2_PROT_RESULTS, "combined_proteomics_de.csv"))
prot_null <- fread(file = file.path(CFG$V2_PROT_RESULTS, "combined_proteomics_null_fdr.csv"))
fwrite(prot_de,   file.path(CFG$OUT_DIR, "proteomics_de.csv"))
fwrite(prot_null, file.path(CFG$OUT_DIR, "proteomics_null.csv"))
add_manifest(prot_de,   "proteomics_de.csv",
             sprintf("v2: %d methods x %d datasets (all methods, DiPhiSeq excluded)",
                     uniqueN(prot_de$method), uniqueN(prot_de$dataset)))
add_manifest(prot_null, "proteomics_null.csv",
             sprintf("v2: %d methods x %d datasets null FDR summary",
                     uniqueN(prot_null$method), uniqueN(prot_null$dataset)))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
fwrite(manifest, file.path(CFG$OUT_DIR, "paper_data_manifest.csv"))
cat("\n=== Paper data aggregation manifest ===\n")
print(manifest[, .(table_name, n_rows, n_cols, n_methods, n_datasets, note)])
cat("\nOutput dir:", CFG$OUT_DIR, "\n")
