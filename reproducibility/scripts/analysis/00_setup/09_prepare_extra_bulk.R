# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(SummarizedExperiment)
library(data.table)

BULK_DIR <- paste0(PROJECT_ROOT, "/benchmark/data/bulk")
SC_DIR   <- paste0(PROJECT_ROOT, "/benchmark/data/single_cell")

# =============================================================================
# =============================================================================

cat("\n=== TCGA tumor vs normal ===\n")

load(file.path(BULK_DIR, "TCGA_tumor_SE.Rda"))
tumor_se <- se_tumor; rm(se_tumor)
load(file.path(BULK_DIR, "TCGA_normal_SE.Rda"))
normal_se <- se_normal; rm(se_normal)

# TCGA barcode: TCGA-XX-XXXX-01 = tumor, -11 = normal
normal_barcodes <- colnames(normal_se)
tumor_barcodes  <- colnames(tumor_se)

tumor_patient <- substr(tumor_barcodes, 1, 12)
normal_patient <- substr(normal_barcodes, 1, 12)

ct_map <- setNames(as.character(colData(tumor_se)$CancerType), tumor_patient)

normal_ct <- ct_map[normal_patient]

target_cancers <- c("BRCA", "LUAD", "KIRC", "LIHC")

cat(" ():\n")
print(table(normal_ct[!is.na(normal_ct)]))

tumor_counts <- assay(tumor_se)
normal_counts <- assay(normal_se)

sapply(target_cancers, function(cancer) {
    cat("\n--- ", cancer, " ---\n", sep = "")
    
    tumor_idx <- which(colData(tumor_se)$CancerType == cancer)
    normal_idx <- which(normal_ct == cancer)
    
    n_tumor <- min(length(tumor_idx), 80)
    n_normal <- min(length(normal_idx), 80)
    
    sel_tumor <- sample(tumor_idx, n_tumor)
    sel_normal <- sample(normal_idx, n_normal)
    
    counts_merged <- cbind(tumor_counts[, sel_tumor], normal_counts[, sel_normal])
    group <- factor(c(rep("tumor", n_tumor), rep("normal", n_normal)))
    
    cat("  Tumor:", n_tumor, "| Normal:", n_normal, "| Genes:", nrow(counts_merged), "\n")
    
    saveRDS(
        list(counts = counts_merged, group = group, cancer = cancer),
        file.path(BULK_DIR, paste0("tcga_", tolower(cancer), "_tumor_vs_normal.rds"))
    )
 cat(" : tcga_", tolower(cancer), "_tumor_vs_normal.rds\n", sep = "")
    NULL
})

rm(tumor_se, normal_se, tumor_counts, normal_counts); gc()

# =============================================================================
# =============================================================================

cat("\n=== muscData Kang18 pseudo-bulk ===\n")

load(file.path(BULK_DIR, "muscData_Kang18.Rda"))

cd <- as.data.frame(colData(sce))
cat(":", length(unique(cd$ind)), "| :", paste(unique(cd$stim), collapse="/"), "\n")
cat(":", length(unique(cd$cell)), "\n")

library(Matrix)
counts_sparse <- counts(sce)
group_keys <- paste(cd$ind, cd$stim, sep = "_")
unique_keys <- sort(unique(group_keys))

cat("", length(unique_keys), " pseudo-bulk ...\n")

pb_counts <- sapply(unique_keys, function(k) {
    idx <- which(group_keys == k)
    Matrix::rowSums(counts_sparse[, idx])
})

pb_meta <- data.table(
    sample = unique_keys,
    ind = sub("_.*", "", unique_keys),
    stim = sub(".*_", "", unique_keys)
)

cat("Pseudo-bulk:", nrow(pb_counts), "genes x", ncol(pb_counts), "samples\n")
print(table(pb_meta$stim))

saveRDS(
    list(counts = as.matrix(pb_counts), group = factor(pb_meta$stim),
         metadata = pb_meta, cell_types = unique(cd$cell)),
    file.path(BULK_DIR, "kang18_pseudobulk.rds")
)
cat(": kang18_pseudobulk.rds\n")
rm(sce, counts_sparse, pb_counts); gc()

# =============================================================================
# =============================================================================

cat("\n=== airway (dexamethasone) ===\n")

library(airway)
data(airway)

airway_counts <- assay(airway)
airway_group <- factor(colData(airway)$dex, levels = c("untrt", "trt"))

cat("Dim:", nrow(airway_counts), "x", ncol(airway_counts), "\n")
print(table(airway_group))

saveRDS(
    list(counts = airway_counts, group = airway_group,
         metadata = as.data.frame(colData(airway))),
    file.path(BULK_DIR, "airway_dex.rds")
)
cat(": airway_dex.rds\n")

# =============================================================================
# =============================================================================

cat("\n=== pasilla (Drosophila RNAi) ===\n")

fn <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla")
pasilla_counts <- as.matrix(read.table(fn, header = TRUE, row.names = 1))
pasilla_group <- factor(c("untreated", "untreated", "untreated", "untreated",
                            "treated", "treated", "treated"))

cat("Dim:", nrow(pasilla_counts), "x", ncol(pasilla_counts), "\n")
print(table(pasilla_group))

saveRDS(
    list(counts = pasilla_counts, group = pasilla_group),
    file.path(BULK_DIR, "pasilla_rnai.rds")
)
cat(": pasilla_rnai.rds\n")

cat("\n", strrep("=", 60), "\n")
cat("!\n")
