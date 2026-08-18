PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = 1L)
data.table::setDTthreads(1L)

PB_DIR  <- paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/05_pseudobulk")
SCE_DIR <- paste0(PROJECT_ROOT, "/benchmark_v2/02_scrna_pseudobulk/data")

library(SingleCellExperiment)
library(Matrix)

# =============================================================================
# 1. jakel MS snRNA-seq (Jakel et al. 2019, GSE118257): donors x MS/control
#    colData: Sample (donor), Condition (MS/control), Celltypes (cell type)
#    Pseudobulk: sum counts per donor within each Condition group
# =============================================================================
cat("=== Creating jakel pseudobulk ===\n")
sce <- readRDS(file.path(SCE_DIR, "jakel_ms_sce.rds"))

cluster_sizes <- sort(table(sce$Celltypes), decreasing = TRUE)
cat("Cluster sizes:\n")
print(cluster_sizes)
top_cluster <- names(cluster_sizes)[1]
cat("Using largest cell type:", top_cluster, "\n")

sce_sub <- sce[, sce$Celltypes == top_cluster]
cat("Subset dim:", nrow(sce_sub), "x", ncol(sce_sub), "\n")

donors <- unique(sce_sub$Sample)
groups <- unique(sce_sub$Condition)
cat("Donors:", paste(donors, collapse=", "), "\n")
cat("Groups:", paste(groups, collapse=", "), "\n")

raw_counts <- counts(sce_sub)
sample_ids <- paste0(sce_sub$Sample, "_", sce_sub$Condition)
unique_samples <- unique(sample_ids)
cat("Unique samples:", length(unique_samples), "\n")

pb_mat <- matrix(0L, nrow = nrow(raw_counts), ncol = length(unique_samples))
rownames(pb_mat) <- rownames(raw_counts)
colnames(pb_mat) <- unique_samples

idx_list <- split(seq_along(sample_ids), sample_ids)
pb_mat <- vapply(unique_samples, function(s) {
  as.integer(rowSums(raw_counts[, idx_list[[s]], drop = FALSE]))
}, integer(nrow(raw_counts)))
rownames(pb_mat) <- rownames(raw_counts)

sample_meta <- data.frame(
  sample = unique_samples,
  donor  = sub("_.*", "", unique_samples),
  group  = sub(".*_", "", unique_samples),
  row.names = unique_samples
)

pb_group <- factor(sample_meta$group)
cat("Pseudobulk dim:", nrow(pb_mat), "x", ncol(pb_mat), "\n")
cat("Group table:\n")
print(table(pb_group))

gene_keep <- rowSums(pb_mat >= 10) >= 2
pb_mat <- pb_mat[gene_keep, ]
cat("After filtering:", nrow(pb_mat), "genes\n")

saveRDS(list(counts = pb_mat, group = pb_group, meta = sample_meta),
        file.path(PB_DIR, "jakel_pseudobulk.rds"))
cat("Saved jakel_pseudobulk.rds\n\n")

# =============================================================================
# 2. Spike-in simulation: DE + DV ground truth
# =============================================================================
cat("=== Creating spike-in simulation pseudobulk ===\n")

kang_pb <- readRDS(file.path(PB_DIR, "kang18_pseudobulk.rds"))
ctrl_idx <- which(kang_pb$group == levels(factor(kang_pb$group))[1])
treat_idx <- which(kang_pb$group == levels(factor(kang_pb$group))[2])
n_ctrl <- length(ctrl_idx)
n_treat <- length(treat_idx)
cat("Base: kang18 pseudobulk,", n_ctrl, "ctrl,", n_treat, "treat\n")

set.seed(12345)
all_ctrl <- c(ctrl_idx, treat_idx)
shuf <- sample(all_ctrl)
grpA_idx <- shuf[seq_len(n_ctrl)]
grpB_idx <- shuf[(n_ctrl + 1):(n_ctrl + n_treat)]

base_counts <- as.matrix(kang_pb$counts)
gene_keep_sim <- rowSums(base_counts >= 5) >= 4
base_counts <- base_counts[gene_keep_sim, ]
n_genes_sim <- nrow(base_counts)
cat("Genes after filtering:", n_genes_sim, "\n")

set.seed(42)
perm <- sample(n_genes_sim)
n_de_only <- round(n_genes_sim * 0.10)
n_dv_only <- round(n_genes_sim * 0.10)
n_de_dv   <- round(n_genes_sim * 0.05)
de_only_idx <- perm[1:n_de_only]
dv_only_idx <- perm[(n_de_only + 1):(n_de_only + n_dv_only)]
de_dv_idx   <- perm[(n_de_only + n_dv_only + 1):(n_de_only + n_dv_only + n_de_dv)]

cat("DE only:", n_de_only, " DV only:", n_dv_only, " DE+DV:", n_de_dv, "\n")

sim_counts <- base_counts

set.seed(123)
de_lfc <- runif(length(c(de_only_idx, de_dv_idx)), min = 1.0, max = 2.0)
de_all_idx <- c(de_only_idx, de_dv_idx)
sim_counts[de_all_idx, grpB_idx] <- round(
  sim_counts[de_all_idx, grpB_idx] * 2^de_lfc
)

set.seed(456)
dv_all_idx <- c(dv_only_idx, de_dv_idx)
dv_scale <- matrix(
  exp(rnorm(length(dv_all_idx) * n_treat, mean = 0, sd = 0.8)),
  nrow = length(dv_all_idx), ncol = n_treat
)
sim_counts[dv_all_idx, grpB_idx] <- round(
  sim_counts[dv_all_idx, grpB_idx] * dv_scale
)

storage.mode(sim_counts) <- "integer"
sim_counts[sim_counts < 0L] <- 0L

sim_group <- factor(c(rep("ctrl", n_ctrl), rep("treat", n_treat)))
colnames(sim_counts) <- paste0("sample_", seq_len(ncol(sim_counts)))

gene_names_sim <- rownames(sim_counts)
category <- rep("ee", n_genes_sim)
category[de_only_idx] <- "de"
category[dv_only_idx] <- "dv"
category[de_dv_idx]   <- "de_dv"

lfc_vec <- rep(0, n_genes_sim)
lfc_vec[de_all_idx] <- de_lfc

truth <- data.frame(
  gene     = gene_names_sim,
  category = category,
  is_de    = category %in% c("de", "de_dv"),
  is_dv    = category %in% c("dv", "de_dv"),
  sim_lfc  = lfc_vec,
  stringsAsFactors = FALSE
)

cat("Truth: DE=", sum(truth$is_de), " DV=", sum(truth$is_dv), " EE=", sum(category == "ee"), "\n")

meanA <- rowMeans(sim_counts[de_all_idx, grpA_idx, drop = FALSE])
meanB <- rowMeans(sim_counts[de_all_idx, grpB_idx, drop = FALSE])
obs_lfc_check <- log2((meanB + 1) / (meanA + 1))
cat("Observed LFC for DE genes - median:", round(median(obs_lfc_check), 2),
    " range:", round(range(obs_lfc_check), 2), "\n")

saveRDS(list(counts = sim_counts, group = sim_group, truth = truth),
        file.path(PB_DIR, "muscat_sim_pseudobulk.rds"))
cat("Saved muscat_sim_pseudobulk.rds\n")
cat("Done.\n")
