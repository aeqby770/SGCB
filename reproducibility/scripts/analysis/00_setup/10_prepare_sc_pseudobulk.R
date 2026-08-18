# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix)
library(data.table)

BULK_DIR <- paste0(PROJECT_ROOT, "/benchmark/data/bulk")
SC_DIR   <- paste0(PROJECT_ROOT, "/benchmark/data/single_cell")

# =============================================================================
# 1. segerstolpe_pancreas: T2D vs Normal (pseudo-bulk by individual)
# =============================================================================

cat("\n=== segerstolpe T2D vs Normal pseudo-bulk ===\n")

seg <- readRDS(file.path(SC_DIR, "segerstolpe_pancreas.rds"))
cd <- as.data.frame(colData(seg))

seg_counts <- counts(seg)
group_keys <- paste(cd$individual, cd$disease, sep = "|")
unique_keys <- sort(unique(group_keys))

cat("", length(unique_keys), " pseudo-bulk...\n")

pb <- sapply(unique_keys, function(k) {
    idx <- which(group_keys == k)
    Matrix::rowSums(seg_counts[, idx])
})

pb_meta <- data.table(
    sample = unique_keys,
    individual = sub("\\|.*", "", unique_keys),
    disease = sub(".*\\|", "", unique_keys)
)

cat("Dim:", nrow(pb), "x", ncol(pb), "\n")
print(table(pb_meta$disease))

saveRDS(
    list(counts = as.matrix(pb), group = factor(pb_meta$disease),
         metadata = pb_meta),
    file.path(BULK_DIR, "segerstolpe_t2d_pseudobulk.rds")
)
cat(": segerstolpe_t2d_pseudobulk.rds\n")
rm(seg, seg_counts, pb); gc()

# =============================================================================
# 2. xin_pancreas: T2D vs Healthy (pseudo-bulk by donor)
# =============================================================================

cat("\n=== xin T2D vs Healthy pseudo-bulk ===\n")

xin <- readRDS(file.path(SC_DIR, "xin_pancreas.rds"))
cd_x <- as.data.frame(colData(xin))

xin_counts <- tryCatch(counts(xin), error = function(e) assay(xin, 1))
group_keys_x <- paste(cd_x$donor.id, cd_x$condition, sep = "|")
unique_keys_x <- sort(unique(group_keys_x))

cat("", length(unique_keys_x), " pseudo-bulk...\n")

pb_x <- sapply(unique_keys_x, function(k) {
    idx <- which(group_keys_x == k)
    Matrix::rowSums(xin_counts[, idx])
})

pb_meta_x <- data.table(
    sample = unique_keys_x,
    donor = sub("\\|.*", "", unique_keys_x),
    condition = sub(".*\\|", "", unique_keys_x)
)

cat("Dim:", nrow(pb_x), "x", ncol(pb_x), "\n")
print(table(pb_meta_x$condition))

saveRDS(
    list(counts = as.matrix(pb_x), group = factor(pb_meta_x$condition),
         metadata = pb_meta_x),
    file.path(BULK_DIR, "xin_t2d_pseudobulk.rds")
)
cat(": xin_t2d_pseudobulk.rds\n")
rm(xin, xin_counts, pb_x); gc()

# =============================================================================
# 3. zilionis_lung: Tumor vs Blood (pseudo-bulk by patient)
# =============================================================================

cat("\n=== zilionis Lung tumor vs blood pseudo-bulk ===\n")

zil <- readRDS(file.path(SC_DIR, "zilionis_lung.rds"))
cd_z <- as.data.frame(colData(zil))

zil_counts <- tryCatch(counts(zil), error = function(e) assay(zil, 1))
group_keys_z <- paste(cd_z$Patient, cd_z$Tissue, sep = "|")
unique_keys_z <- sort(unique(group_keys_z))

cat("", length(unique_keys_z), " pseudo-bulk...\n")

pb_z <- sapply(unique_keys_z, function(k) {
    idx <- which(group_keys_z == k)
    Matrix::rowSums(zil_counts[, idx])
})

pb_meta_z <- data.table(
    sample = unique_keys_z,
    patient = sub("\\|.*", "", unique_keys_z),
    tissue = sub(".*\\|", "", unique_keys_z)
)

cat("Dim:", nrow(pb_z), "x", ncol(pb_z), "\n")
print(table(pb_meta_z$tissue))

saveRDS(
    list(counts = as.matrix(pb_z), group = factor(pb_meta_z$tissue),
         metadata = pb_meta_z),
    file.path(BULK_DIR, "zilionis_lung_pseudobulk.rds")
)
cat(": zilionis_lung_pseudobulk.rds\n")
rm(zil, zil_counts, pb_z); gc()

# =============================================================================
# =============================================================================

cat("\n=== TCGA ===\n")

load(file.path(BULK_DIR, "TCGA_tumor_SE.Rda"))
tumor_se <- se_tumor; rm(se_tumor)
load(file.path(BULK_DIR, "TCGA_normal_SE.Rda"))
normal_se <- se_normal; rm(se_normal)

tumor_barcodes <- colnames(tumor_se)
normal_barcodes <- colnames(normal_se)
tumor_patient <- substr(tumor_barcodes, 1, 12)
normal_patient <- substr(normal_barcodes, 1, 12)
ct_map <- setNames(as.character(colData(tumor_se)$CancerType), tumor_patient)
normal_ct <- ct_map[normal_patient]

tumor_counts <- assay(tumor_se)
normal_counts <- assay(normal_se)

extra_cancers <- c("COAD", "HNSC", "STAD")

sapply(extra_cancers, function(cancer) {
    cat("--- ", cancer, " ---\n", sep = "")

    tumor_idx <- which(colData(tumor_se)$CancerType == cancer)
    normal_idx <- which(normal_ct == cancer)

    n_tumor <- min(length(tumor_idx), 80)
    n_normal <- min(length(normal_idx), 80)

    cat("  Available tumor:", length(tumor_idx), "| normal:", length(normal_idx), "\n")

    sel_tumor <- sample(tumor_idx, n_tumor)
    sel_normal <- sample(normal_idx, n_normal)

    counts_merged <- cbind(tumor_counts[, sel_tumor], normal_counts[, sel_normal])
    group <- factor(c(rep("tumor", n_tumor), rep("normal", n_normal)))

    cat("  Final:", n_tumor, "vs", n_normal, "\n")

    saveRDS(
        list(counts = counts_merged, group = group, cancer = cancer),
        file.path(BULK_DIR, paste0("tcga_", tolower(cancer), "_tumor_vs_normal.rds"))
    )
 cat(" : tcga_", tolower(cancer), "_tumor_vs_normal.rds\n", sep = "")
    NULL
})

rm(tumor_se, normal_se, tumor_counts, normal_counts); gc()

cat("\n", strrep("=", 60), "\n")
cat("!\n")
