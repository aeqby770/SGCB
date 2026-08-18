##############################################################################
##############################################################################

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk/data")
set.seed(12345)

library(ExperimentHub)
library(SingleCellExperiment)

eh <- ExperimentHub()

# ====================================================================
# ====================================================================
cat("=== Kang PBMC IFN- ===\n")
library(muscData)
kang <- Kang18_8vs8()
cat(sprintf(" : %d genes %d cells\n", nrow(kang), ncol(kang)))
cat(sprintf(" : %s\n", paste(unique(kang$stim), collapse = ", ")))
cat(sprintf(" : %d\n", length(unique(kang$ind))))
saveRDS(kang, file.path(DATA_DIR, "kang_pbmc_sce.rds"))
cat(" : kang_pbmc_sce.rds\n\n")

# ====================================================================
# ====================================================================
cat("=== Jakel MS snRNA-seq (GSE118257) ===\n")
geo_dir <- file.path(DATA_DIR, "GSE118257_raw")
dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/"
suppl_files <- c(
  "GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz",
  "GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz"
)

options(timeout = 600)
dl_ok <- TRUE
for (f in suppl_files) {
  dest <- file.path(geo_dir, f)
 cat(sprintf(" : %s\n", f))
  tryCatch(
    download.file(paste0(base_url, f), dest, mode = "wb", quiet = TRUE),
 error = function(e) { cat(sprintf(" !! : %s\n", e$message)); dl_ok <<- FALSE }
  )
}

if (dl_ok && file.exists(file.path(geo_dir, suppl_files[1]))) {
 cat(" ...\n")
  expr_mat <- data.table::fread(file.path(geo_dir, suppl_files[1]), header = TRUE)
  genes <- expr_mat[[1]]
  expr_mat <- as.matrix(expr_mat[, -1])
  rownames(expr_mat) <- genes

  anno <- data.table::fread(file.path(geo_dir, suppl_files[2]), header = TRUE)

 cat(sprintf(" : %d genes %d cells\n", nrow(expr_mat), ncol(expr_mat)))
 cat(sprintf(" : %d\n", nrow(anno)))

  sce_ms <- SingleCellExperiment(
    assays = list(counts = expr_mat),
    colData = anno[match(colnames(expr_mat), anno[[1]]), ]
  )
  saveRDS(sce_ms, file.path(DATA_DIR, "jakel_ms_sce.rds"))
 cat(" : jakel_ms_sce.rds\n\n")
} else {
 cat(" !! GSE118257 SCE \n\n")
}

# ====================================================================
# ====================================================================
cat("=== Segerstolpe pancreas ===\n")
seg <- muscData::Segerstolpe18()
cat(sprintf(" : %d genes %d cells\n", nrow(seg), ncol(seg)))
saveRDS(seg, file.path(DATA_DIR, "segerstolpe_pancreas_sce.rds"))
cat(" : segerstolpe_pancreas_sce.rds\n\n")

# ====================================================================
# ====================================================================
cat("===== scRNA =====\n")
rds_files <- list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE)
for (f in rds_files) {
  sz <- file.info(f)$size / 1e6
  cat(sprintf("  %-45s  %.1f MB\n", basename(f), sz))
}
cat("\n")
