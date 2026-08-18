PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = as.character(parallel::detectCores()))
data.table::setDTthreads(parallel::detectCores())

library(data.table)

# =============================================================================
# =============================================================================

DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark_v2/03_proteomics/data")
OUT_FILE <- file.path(DATA_DIR, "sim_proteomics.rds")

pep_dt <- fread(file.path(DATA_DIR, "cptac_lab3_peptides.txt"))
pep_dt <- pep_dt[is.na(Reverse) | Reverse != "+"]
pep_dt <- pep_dt[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]

int_cols <- grep("^Intensity ", names(pep_dt), value = TRUE)
int_mat  <- as.matrix(pep_dt[, ..int_cols])
int_mat[int_mat == 0] <- NA

prot_id <- pep_dt[["Leading razor protein"]]
log2_pep <- log2(int_mat)
pep_tbl  <- data.table(protein = prot_id, log2_pep)
pep_tbl  <- pep_tbl[!is.na(protein) & protein != ""]

pep_counts <- pep_tbl[, .(n_pep = .N), by = protein]

prot_long <- melt(pep_tbl, id.vars = "protein", variable.name = "sample", value.name = "log2int")
prot_wide <- dcast(prot_long, protein ~ sample, value.var = "log2int", fun.aggregate = median, na.rm = TRUE)
prot_mat  <- as.matrix(prot_wide[, -1])
rownames(prot_mat) <- prot_wide$protein
prot_mat[is.nan(prot_mat)] <- NA

n_valid <- rowSums(!is.na(prot_mat))
prot_mat <- prot_mat[n_valid >= 4, ]
n_prot <- nrow(prot_mat)
cat("Proteins after filtering:", n_prot, "\n")

set.seed(42)
n_samp <- ncol(prot_mat)
shuf <- sample(n_samp)
n_A <- n_samp %/% 2
grpA_idx <- shuf[1:n_A]
grpB_idx <- shuf[(n_A + 1):n_samp]
cat("Group A:", n_A, " Group B:", length(grpB_idx), "\n")

set.seed(123)
perm <- sample(n_prot)
n_de_only <- round(n_prot * 0.10)
n_dv_only <- round(n_prot * 0.10)
n_de_dv   <- round(n_prot * 0.05)
de_only_idx <- perm[1:n_de_only]
dv_only_idx <- perm[(n_de_only + 1):(n_de_only + n_dv_only)]
de_dv_idx   <- perm[(n_de_only + n_dv_only + 1):(n_de_only + n_dv_only + n_de_dv)]

cat("DE only:", n_de_only, " DV only:", n_dv_only, " DE+DV:", n_de_dv, "\n")

sim_mat <- prot_mat

set.seed(456)
de_all_idx <- c(de_only_idx, de_dv_idx)
de_lfc <- runif(length(de_all_idx), min = 0.8, max = 1.5)
sim_mat[de_all_idx, grpB_idx] <- sim_mat[de_all_idx, grpB_idx] + de_lfc

set.seed(789)
dv_all_idx <- c(dv_only_idx, de_dv_idx)
dv_noise <- matrix(rnorm(length(dv_all_idx) * length(grpB_idx), 0, 0.6),
                   nrow = length(dv_all_idx), ncol = length(grpB_idx))
sim_mat[dv_all_idx, grpB_idx] <- sim_mat[dv_all_idx, grpB_idx] + dv_noise

sim_group <- factor(c(rep("ctrl", n_A), rep("treat", length(grpB_idx))))
colnames(sim_mat) <- paste0("sample_", seq_len(n_samp))

pep_count_matched <- pep_counts[match(rownames(sim_mat), pep_counts$protein), n_pep]
pep_count_matched[is.na(pep_count_matched)] <- 1L

# Truth table
gene_names <- rownames(sim_mat)
category <- rep("ee", n_prot)
category[de_only_idx] <- "de"
category[dv_only_idx] <- "dv"
category[de_dv_idx]   <- "de_dv"

lfc_vec <- rep(0, n_prot)
lfc_vec[de_all_idx] <- de_lfc

truth <- data.frame(
  gene     = gene_names,
  category = category,
  is_de    = category %in% c("de", "de_dv"),
  is_dv    = category %in% c("dv", "de_dv"),
  sim_lfc  = lfc_vec,
  stringsAsFactors = FALSE
)

cat("Truth: DE=", sum(truth$is_de), " DV=", sum(truth$is_dv), " EE=", sum(category == "ee"), "\n")

meanA <- rowMeans(sim_mat[de_all_idx, grpA_idx, drop = FALSE], na.rm = TRUE)
meanB <- rowMeans(sim_mat[de_all_idx, grpB_idx, drop = FALSE], na.rm = TRUE)
obs_lfc <- meanB - meanA
cat("Observed LFC for DE proteins - median:", round(median(obs_lfc, na.rm = TRUE), 2), "\n")

saveRDS(list(prot_log2 = sim_mat, group = sim_group,
             pep_count = pep_count_matched, truth = truth,
             is_spike = truth$is_de),
        OUT_FILE)
cat("Saved sim_proteomics.rds\n")
