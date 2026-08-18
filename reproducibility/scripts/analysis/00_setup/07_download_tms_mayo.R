# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

Sys.setenv(CURLOPT_SSLVERSION = 6L,
           CURLOPT_SSL_VERIFYPEER = 0L,
           CURLOPT_SSL_VERIFYHOST = 0L)
httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L, sslversion = 6L))
options(download.file.method = "libcurl",
        download.file.extra = "--insecure",
        BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

.original_curl_fetch <- curl::curl_fetch_memory
patch_curl <- function(url, handle = curl::new_handle()) {
    curl::handle_setopt(handle, ssl_verifypeer = 0L, ssl_verifyhost = 0L, sslversion = 6L)
    .original_curl_fetch(url, handle)
}
assignInNamespace("curl_fetch_memory", patch_curl, ns = "curl")

library(data.table)
library(SummarizedExperiment)

BULK_DIR <- paste0(PROJECT_ROOT, "/benchmark/data/bulk")
dir.create(BULK_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Tabula Muris Senis Bulk RNA-seq (GSE132040)
# =============================================================================

cat("\n========== 1. Tabula Muris Senis Bulk ==========\n")

library(TabulaMurisSenisData)

tms_bulk <- TabulaMurisSenisBulk()

cat(": ", nrow(tms_bulk), " genes x ", ncol(tms_bulk), " samples\n", sep = "")

counts_tms <- as.matrix(assay(tms_bulk, "counts"))
metadata_tms <- as.data.frame(colData(tms_bulk))

cat(":\n")
print(table(metadata_tms$age))
cat("\n:\n")
print(table(metadata_tms$tissue))

saveRDS(
    list(counts = counts_tms, metadata = metadata_tms),
    file.path(BULK_DIR, "tms_bulk.rds")
)
cat(": tms_bulk.rds (", round(file.size(file.path(BULK_DIR, "tms_bulk.rds")) / 1e6, 1), " MB)\n", sep = "")

# =============================================================================
# =============================================================================

cat("\n========== 2. Tabula Muris Senis FACS ==========\n")

SC_DIR <- paste0(PROJECT_ROOT, "/benchmark/data/single_cell")
dir.create(SC_DIR, showWarnings = FALSE, recursive = TRUE)

tms_facs <- TabulaMurisSenisFACS(tissues = "All", processedCounts = TRUE)

cat(": ", nrow(tms_facs), " genes x ", ncol(tms_facs), " cells\n", sep = "")

saveRDS(tms_facs, file.path(SC_DIR, "tms_facs_all.rds"))
cat(": tms_facs_all.rds\n")

# =============================================================================
#    AMP-AD Synapse syn5550404
#    https://www.synapse.org/#!Synapse:syn5550404
# =============================================================================

cat("\n========== 3. Mayo AD RNA-Seq ==========\n")

library(GEOquery)

# GSE115811: Mayo Clinic RNA-Seq Study (Temporal Cortex)
# AD vs Control, PSP vs Control
mayo_geo <- getGEO("GSE115811", GSEMatrix = TRUE, getGPL = FALSE)

mayo_eset <- mayo_geo[[1]]
cat(": ", nrow(mayo_eset), " features x ", ncol(mayo_eset), " samples\n", sep = "")

mayo_pdata <- pData(mayo_eset)
cat("\n:", paste(colnames(mayo_pdata), collapse = ", "), "\n")

diag_col <- grep("diagnosis|disease|status|group", colnames(mayo_pdata), ignore.case = TRUE, value = TRUE)
cat(":", paste(diag_col, collapse = ", "), "\n")

saveRDS(
    list(expression = exprs(mayo_eset), phenotype = mayo_pdata),
    file.path(BULK_DIR, "mayo_ad_rnaseq.rds")
)
cat(": mayo_ad_rnaseq.rds\n")

cat("\n========== ==========\n")
