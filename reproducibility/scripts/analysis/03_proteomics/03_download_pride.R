# =============================================================================
# =============================================================================
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics/data")
options(timeout = 600)

# =============================================================================
#    Human background + E.coli spiked at 3 ratios
# =============================================================================
cat("=== [1] PXD000279: E.coli spike-in LFQ ===\n")
lfq_dir <- file.path(DATA_DIR, "ecoli_lfq_deqms")
dir.create(lfq_dir, showWarnings = FALSE)

lfq_url <- "https://ftp.ebi.ac.uk/pride-archive/2014/09/PXD000279/proteomebenchmark.zip"
lfq_zip <- file.path(lfq_dir, "proteomebenchmark.zip")
tryCatch({
  download.file(lfq_url, lfq_zip, mode = "wb", quiet = FALSE)
  utils::unzip(lfq_zip, exdir = lfq_dir)
  cat("  Extracted LFQ benchmark data\n")
  cat("  Files:", paste(list.files(lfq_dir), collapse = ", "), "\n")
}, error = function(e) cat("  LFQ download failed:", e$message, "\n"))

writeLines(c(
  "Dataset: E.coli spike-in label-free proteomics (DEqMS core benchmark)",
  "Reference: Zhu et al. (2020) Mol Cell Proteomics 19:1047-1057",
  "PRIDE: PXD000279",
  "Design: Human background + E.coli spiked at multiple ratios",
  "Ground truth: E.coli proteins are DE; human proteins are non-DE",
  "Quantification: MaxQuant LFQ intensities"
), file.path(lfq_dir, "README.txt"))

# =============================================================================
# =============================================================================
cat("\n=== [2] PXD004163: TMT miR proteomics ===\n")
tmt_dir <- file.path(DATA_DIR, "tmt_mir_deqms")
dir.create(tmt_dir, showWarnings = FALSE)

tmt_url <- "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
tmt_file <- file.path(tmt_dir, "Yan_miR_Protein_table.txt")
tryCatch({
  download.file(tmt_url, tmt_file, mode = "wb", quiet = FALSE)
  cat(sprintf("  Downloaded TMT protein table: %.2f MB\n", file.size(tmt_file) / 1e6))
}, error = function(e) cat("  TMT download failed:", e$message, "\n"))

writeLines(c(
  "Dataset: TMT-labeled miR proteomics (DEqMS TMT benchmark)",
  "Reference: Zhu et al. (2020) MCP; Yan et al. (2016)",
  "PRIDE: PXD004163",
  "Design: TMT 10-plex, miR knockdown experiment",
  "Quantification: TMT reporter ion intensities"
), file.path(tmt_dir, "README.txt"))

# =============================================================================
# =============================================================================
cat("\n=== [3] ExperimentHub EH1663: E.coli TMT PSM ===\n")
eh_dir <- file.path(DATA_DIR, "ecoli_tmt_deqms")
dir.create(eh_dir, showWarnings = FALSE)

tryCatch({
  library(ExperimentHub)
  eh <- ExperimentHub()
  dat_psm <- eh[["EH1663"]]
  saveRDS(dat_psm, file.path(eh_dir, "EH1663_ecoli_tmt_psm.rds"))
  cat(sprintf("  Saved EH1663: %d rows x %d cols\n", nrow(dat_psm), ncol(dat_psm)))
}, error = function(e) cat("  ExperimentHub failed:", e$message, "\n"))

writeLines(c(
  "Dataset: E.coli spike-in TMT PSM-level data",
  "Source: ExperimentHub EH1663 (Bioconductor)",
  "Reference: Zhu et al. (2020) MCP",
  "Design: E.coli proteins spiked in human at known TMT ratios",
  "Ground truth: E.coli proteins are DE; human proteins are non-DE"
), file.path(eh_dir, "README.txt"))

# =============================================================================
# =============================================================================
cat("\n=== Cleaning empty directories ===\n")
all_dirs <- list.dirs(DATA_DIR, recursive = FALSE)
invisible(lapply(all_dirs, function(d) {
  data_files <- list.files(d, pattern = "\\.(txt|csv|rds|tsv|zip|xlsx|gz|RData|tab)$")
  readme_only <- length(data_files) == 0 && length(list.files(d)) <= 1
  nm <- basename(d)
  cat(sprintf("  %-30s %2d data files %s\n", nm, length(data_files),
      ifelse(readme_only, "[EMPTY - keep for metadata]", "[OK]")))
}))

# =============================================================================
# =============================================================================
cat("\n========== FINAL INVENTORY ==========\n")
all_dirs <- list.dirs(DATA_DIR, recursive = FALSE)
total_mb <- 0
invisible(lapply(sort(all_dirs), function(d) {
  all_f <- list.files(d, recursive = TRUE, full.names = TRUE)
  sz <- sum(file.size(all_f), na.rm = TRUE) / 1e6
  total_mb <<- total_mb + sz
  data_n <- length(all_f[!grepl("README", all_f)])
  status <- ifelse(data_n > 0, "OK", "--")
  cat(sprintf("  [%2s] %-30s %2d files  (%.1f MB)\n", status, basename(d), data_n, sz))
}))

root_files <- list.files(DATA_DIR, recursive = FALSE, full.names = TRUE)
root_files <- root_files[!file.info(root_files)$isdir]
root_mb <- sum(file.size(root_files), na.rm = TRUE) / 1e6
total_mb <- total_mb + root_mb
cat(sprintf("\n  Root-level files: %d (%.1f MB)\n", length(root_files), root_mb))
cat(sprintf("\n  TOTAL: %.1f MB\n", total_mb))
