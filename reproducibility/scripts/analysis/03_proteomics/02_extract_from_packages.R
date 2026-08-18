# =============================================================================
# =============================================================================
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics/data")

# =============================================================================
# =============================================================================
cat("=== [1] Installing wrProteo for Ramus 2016 data ===\n")
install.packages("wrProteo", repos = "https://cran.r-project.org", quiet = TRUE)
ext <- system.file("extdata", package = "wrProteo")
ramus_dir <- file.path(DATA_DIR, "ramus2016_ups1")
dir.create(ramus_dir, showWarnings = FALSE)
ramus_files <- list.files(ext, pattern = "UPS|Ramus|proteinGroup", full.names = TRUE)
cat("  wrProteo extdata files:\n")
cat(paste("   ", list.files(ext)), sep = "\n")
file.copy(list.files(ext, full.names = TRUE), ramus_dir, overwrite = TRUE)
cat(sprintf("  Copied %d files to %s\n", length(list.files(ramus_dir)), ramus_dir))

# =============================================================================
# =============================================================================
cat("\n=== [2] Installing DEqMS for E.coli spike-in data ===\n")
BiocManager::install("DEqMS", update = FALSE, ask = FALSE, quiet = TRUE)
ext_deqms <- system.file("extdata", package = "DEqMS")
cat("  DEqMS extdata files:\n")
cat(paste("   ", list.files(ext_deqms)), sep = "\n")

lfq_dir <- file.path(DATA_DIR, "ecoli_lfq_deqms")
dir.create(lfq_dir, showWarnings = FALSE)
file.copy(list.files(ext_deqms, full.names = TRUE), lfq_dir, overwrite = TRUE)
cat(sprintf("  Copied %d files to %s\n", length(list.files(lfq_dir)), lfq_dir))

# =============================================================================
# =============================================================================
cat("\n=== [3] Checking msdata for CPTAC ===\n")
tryCatch({
  BiocManager::install("msdata", update = FALSE, ask = FALSE, quiet = TRUE)
  ext_msdata <- system.file("extdata", package = "msdata")
  cat("  msdata extdata:", paste(list.files(ext_msdata), collapse = ", "), "\n")
}, error = function(e) cat("  msdata not available:", e$message, "\n"))

# =============================================================================
# =============================================================================
cat("\n=== [4] Checking proDA for O'Brien 2018 ===\n")
obrien_dir <- file.path(DATA_DIR, "obrien2018_3species")
obrien_file <- file.path(obrien_dir, "proteinGroups.txt")
cat(sprintf("  O'Brien 2018: %s (%.2f MB)\n",
    ifelse(file.exists(obrien_file), "EXISTS", "MISSING"),
    file.size(obrien_file) / 1e6))

# =============================================================================
# =============================================================================
cat("\n=== [5] Additional CPTAC downloads ===\n")
cptac_dir <- file.path(DATA_DIR, "cptac_ups1")

pg_file <- file.path(cptac_dir, "proteinGroups.txt")
cat(sprintf("  CPTAC proteinGroups: %s (%.2f MB)\n",
    ifelse(file.exists(pg_file), "EXISTS", "MISSING"),
    file.size(pg_file) / 1e6))

# =============================================================================
# =============================================================================
cat("\n=== [6] Checking MSstats data ===\n")
tryCatch({
  BiocManager::install("MSstats", update = FALSE, ask = FALSE, quiet = TRUE)
  library(MSstats)
  msstats_dir <- file.path(DATA_DIR, "msstats_benchmark")
  dir.create(msstats_dir, showWarnings = FALSE)

  data("DDARawData", package = "MSstats")
  saveRDS(DDARawData, file.path(msstats_dir, "DDARawData.rds"))
  cat(sprintf("  DDARawData: %d rows\n", nrow(DDARawData)))

  data("DIARawData", package = "MSstats")
  saveRDS(DIARawData, file.path(msstats_dir, "DIARawData.rds"))
  cat(sprintf("  DIARawData: %d rows\n", nrow(DIARawData)))

  data("SRMRawData", package = "MSstats")
  saveRDS(SRMRawData, file.path(msstats_dir, "SRMRawData.rds"))
  cat(sprintf("  SRMRawData: %d rows\n", nrow(SRMRawData)))
}, error = function(e) cat("  MSstats not available:", e$message, "\n"))

# =============================================================================
# =============================================================================
cat("\n========== FINAL DATA INVENTORY ==========\n")
all_dirs <- list.dirs(DATA_DIR, recursive = FALSE)
invisible(lapply(sort(all_dirs), function(d) {
  files <- list.files(d, recursive = FALSE)
  data_files <- files[!grepl("README", files)]
  total_size <- sum(file.size(file.path(d, data_files)), na.rm = TRUE)
  size_mb <- round(total_size / 1e6, 2)
  status <- ifelse(length(data_files) > 0, "OK", "EMPTY")
  cat(sprintf("  [%s] %-30s %2d files  (%s MB)\n",
      status, basename(d), length(data_files), size_mb))
}))
