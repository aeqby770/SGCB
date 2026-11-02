# Deterministic input demonstrating the public API.
CONFIG <- list(threads = as.integer(Sys.getenv("SGCB_THREADS", as.character(parallel::detectCores()))),
               n_per_group = 50L, n_genes = 120L,
               output = "SGCB_minimal_results.csv")
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
library(SGCB)

gene <- seq_len(CONFIG$n_genes)
probability <- (seq_len(CONFIG$n_per_group) - 0.5) / CONFIG$n_per_group
mean_a <- exp(seq(log(20), log(2000), length.out = CONFIG$n_genes))
mean_b <- mean_a * ifelse(gene %% 3 == 1, 2, 1)
shape_a <- rep(2, CONFIG$n_genes)
shape_b <- ifelse(gene %% 3 == 2, 0.5, 2)
abundance_a <- matrix(qgamma(rep(probability, each = CONFIG$n_genes),
                             shape = shape_a, scale = mean_a / shape_a), CONFIG$n_genes)
abundance_b <- matrix(qgamma(rep(probability, each = CONFIG$n_genes),
                             shape = shape_b, scale = mean_b / shape_b), CONFIG$n_genes)
abundance <- cbind(abundance_a, abundance_b)
rownames(abundance) <- paste0("example_", gene)
colnames(abundance) <- c(paste0("A", seq_len(CONFIG$n_per_group)),
                        paste0("B", seq_len(CONFIG$n_per_group)))
group <- rep(c("A", "B"), each = CONFIG$n_per_group)

result <- sgcbDE(abundance, group, normalization = "none", offset = 0,
                 min_count = 0, min_samples = 1, alpha = 0.05, threads = CONFIG$threads)
primary <- as.data.frame(result)[, c("gene_id", "log2FoldChange", "pvalue", "padj",
    "ctrl_nonzero_fraction", "treat_nonzero_fraction", "ctrl_Q", "treat_Q",
    "dv_log2_cv2_ratio", "dv_se_log2_cv2", "dv_pvalue", "dv_padj",
    "dv_ci_lower", "dv_ci_upper", "DE_joint_q", "DV_joint_q", "DE", "DV", "channel",
    "fit_converged", "fit_covariance_available", "dv_moments_finite", "dv_inference_available")]
write.csv(primary, CONFIG$output, row.names = FALSE, fileEncoding = "UTF-8")
print(primary[1:6, ])
print(as.data.frame(attr(result, "fit_diagnostics")$ctrl)[1:6, ])
# DE targets unconditional arithmetic means; DV targets abundance CV2.
# dv_ci_lower/upper give the 95% interval. The score regime is recorded in attributes.
# logspread_pvalue is a separate spread test; shape_Q_pvalue is exploratory.
