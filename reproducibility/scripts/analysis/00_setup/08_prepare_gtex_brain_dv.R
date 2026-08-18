# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

BULK_DIR <- paste0(PROJECT_ROOT, "/benchmark/data/bulk")

d <- readRDS(file.path(BULK_DIR, "gtex_brain.rds"))

cat(" GTEx Brain: ", nrow(d$counts), " genes x ", ncol(d$counts), " samples\n", sep = "")

region <- d$metadata$gtex.smtsd
cat("\n:\n")
print(table(region))

# =============================================================================
# =============================================================================

cat("\n========== : Cerebellum vs Frontal Cortex ==========\n")

cerebellum_idx <- which(region == "Brain - Cerebellum")
cortex_idx <- which(region == "Brain - Frontal Cortex (BA9)")

cat("Cerebellum: ", length(cerebellum_idx), " samples\n", sep = "")
cat("Frontal Cortex: ", length(cortex_idx), " samples\n", sep = "")

sel_idx <- c(cerebellum_idx, cortex_idx)
counts_cb_fc <- d$counts[, sel_idx]
metadata_cb_fc <- d$metadata[sel_idx, ]
group_cb_fc <- factor(
    c(rep("Cerebellum", length(cerebellum_idx)), 
      rep("FrontalCortex", length(cortex_idx)))
)

saveRDS(
    list(counts = counts_cb_fc, metadata = metadata_cb_fc, group = group_cb_fc),
    file.path(BULK_DIR, "gtex_brain_cerebellum_vs_cortex.rds")
)
cat(": gtex_brain_cerebellum_vs_cortex.rds\n")
cat("  ", nrow(counts_cb_fc), " genes x ", ncol(counts_cb_fc), " samples\n", sep = "")

# =============================================================================
# =============================================================================

cat("\n========== : Cerebellar Hemisphere vs Cortex ==========\n")

cerehemi_idx <- which(region == "Brain - Cerebellar Hemisphere")
cortex_broad_idx <- which(region == "Brain - Cortex")

cat("Cerebellar Hemisphere: ", length(cerehemi_idx), " samples\n", sep = "")
cat("Cortex: ", length(cortex_broad_idx), " samples\n", sep = "")

sel_idx2 <- c(cerehemi_idx, cortex_broad_idx)
counts_ch_cx <- d$counts[, sel_idx2]
metadata_ch_cx <- d$metadata[sel_idx2, ]
group_ch_cx <- factor(
    c(rep("CerebellarHemi", length(cerehemi_idx)),
      rep("Cortex", length(cortex_broad_idx)))
)

saveRDS(
    list(counts = counts_ch_cx, metadata = metadata_ch_cx, group = group_ch_cx),
    file.path(BULK_DIR, "gtex_brain_cerehemi_vs_cortex.rds")
)
cat(": gtex_brain_cerehemi_vs_cortex.rds\n")
cat("  ", nrow(counts_ch_cx), " genes x ", ncol(counts_ch_cx), " samples\n", sep = "")

# =============================================================================
# =============================================================================

cat("\n========== : Hippocampus vs Hypothalamus ==========\n")

hippo_idx <- which(region == "Brain - Hippocampus")
hypo_idx <- which(region == "Brain - Hypothalamus")

cat("Hippocampus: ", length(hippo_idx), " samples\n", sep = "")
cat("Hypothalamus: ", length(hypo_idx), " samples\n", sep = "")

sel_idx3 <- c(hippo_idx, hypo_idx)
counts_hi_hy <- d$counts[, sel_idx3]
metadata_hi_hy <- d$metadata[sel_idx3, ]
group_hi_hy <- factor(
    c(rep("Hippocampus", length(hippo_idx)),
      rep("Hypothalamus", length(hypo_idx)))
)

saveRDS(
    list(counts = counts_hi_hy, metadata = metadata_hi_hy, group = group_hi_hy),
    file.path(BULK_DIR, "gtex_brain_hippo_vs_hypo.rds")
)
cat(": gtex_brain_hippo_vs_hypo.rds\n")
cat("  ", nrow(counts_hi_hy), " genes x ", ncol(counts_hi_hy), " samples\n", sep = "")

cat("\n========== GTEx Brain ==========\n")
