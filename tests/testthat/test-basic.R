# =============================================================================
# SGCB basic functionality tests
# =============================================================================

test_that("sgcbDE differential analysis works correctly", {
  set.seed(123)
  n_genes <- 500
  n_samples <- 6
  
  # Simulate count data
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("sample_", seq_len(n_samples))
  
  # Add differential expression
  de_idx <- 1:50
  counts[de_idx, 4:6] <- counts[de_idx, 4:6] * 3
  
  group <- rep(c("ctrl", "treat"), each = 3)
  
  # Run analysis
  res <- sgcbDE(counts, group, alpha = 0.1)
  
  # Tests
  expect_s3_class(res, "SGCBResults")
  expect_true(nrow(res) > 0)
  expect_true(nrow(res) <= n_genes)
  expect_true("log2FoldChange" %in% names(res))
  expect_true("padj" %in% names(res))
  expect_true("de_fdr_call" %in% names(res))
  expect_true("dv_fdr_call" %in% names(res))
  expect_true("dg_fdr_call" %in% names(res))
  expect_true("dv_model_fdr_call" %in% names(res))
  expect_true("dv_model_raw_fdr_call" %in% names(res))
  expect_true("dv_bf_fdr_call" %in% names(res))
  expect_true("dv_mix_fdr_call" %in% names(res))
  expect_true("SGCB_Score" %in% names(res))
  expect_true("pvalue_dd" %in% names(res))
  expect_true("p_de_post" %in% names(res))
  expect_true("p_dv_post" %in% names(res))
  expect_true("p_dg_post" %in% names(res))
  expect_true("de_post_call" %in% names(res))
  expect_true("dv_post_call" %in% names(res))
  expect_true("dg_post_call" %in% names(res))
  expect_true("p_de_only_post" %in% names(res))
  expect_true("p_dv_only_post" %in% names(res))
  expect_true("p_dg_only_post" %in% names(res))
  expect_true("p_de_dv_post" %in% names(res))
  expect_true("p_de_dg_post" %in% names(res))
  expect_true("p_dv_dg_post" %in% names(res))
  expect_true("p_de_dv_dg_post" %in% names(res))
  expect_true("model_prob_null" %in% names(res))
  expect_true("model_prob_de_dv_dg" %in% names(res))

  mp <- as.matrix(res[, c(
    "model_prob_null", "model_prob_de", "model_prob_dv", "model_prob_dg",
    "model_prob_de_dv", "model_prob_de_dg", "model_prob_dv_dg", "model_prob_de_dv_dg"
  )])
  expect_true(all(is.finite(mp)))
  expect_true(all(abs(rowSums(mp) - 1) < 1e-8))
  expect_true(all(res$p_de_post >= 0 & res$p_de_post <= 1))
  expect_true(all(res$p_dv_post >= 0 & res$p_dv_post <= 1))
  expect_true(all(res$p_dg_post >= 0 & res$p_dg_post <= 1))
  expect_true(all(res$p_de_only_post >= 0 & res$p_de_only_post <= 1))
  expect_true(all(res$p_dv_only_post >= 0 & res$p_dv_only_post <= 1))
  expect_true(all(res$p_dg_only_post >= 0 & res$p_dg_only_post <= 1))
  expect_true(all(res$p_de_dv_post >= 0 & res$p_de_dv_post <= 1))
  expect_true(all(res$p_de_dg_post >= 0 & res$p_de_dg_post <= 1))
  expect_true(all(res$p_dv_dg_post >= 0 & res$p_dv_dg_post <= 1))
  expect_true(all(res$p_de_dv_dg_post >= 0 & res$p_de_dv_dg_post <= 1))
  expect_true(all(res$de_post_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_post_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dg_post_call %in% c(TRUE, FALSE)))
  expect_true("dv_pvalue_model" %in% names(res))
  expect_true("dv_pvalue_model_raw" %in% names(res))
  expect_true("dv_padj_model" %in% names(res))
  expect_true("dv_padj_model_raw" %in% names(res))
  expect_true("dv_pvalue_bf" %in% names(res))
  expect_true("dv_padj_bf" %in% names(res))
  expect_true("dv_pvalue_mix" %in% names(res))
  expect_true("dv_padj_mix" %in% names(res))
  expect_true(all(res$dv_pvalue_model >= 0 & res$dv_pvalue_model <= 1))
  expect_true(all(res$dv_pvalue_model_raw >= 0 & res$dv_pvalue_model_raw <= 1))
  expect_true(all(res$dv_pvalue_bf >= 0 & res$dv_pvalue_bf <= 1))
  expect_true(all(res$dv_pvalue_mix >= 0 & res$dv_pvalue_mix <= 1))
  expect_true(all(res$de_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dg_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_model_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_model_raw_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_bf_fdr_call %in% c(TRUE, FALSE)))
  expect_true(all(res$dv_mix_fdr_call %in% c(TRUE, FALSE)))
  expect_true(is.numeric(attr(res, "post_dv_realized_fdr")))
  expect_true(attr(res, "post_dv_realized_fdr") >= 0 && attr(res, "post_dv_realized_fdr") <= 0.11)
  expect_true(is.numeric(attr(res, "prior_de_used")))
  expect_true(is.numeric(attr(res, "prior_dv_used")))
  expect_true(is.numeric(attr(res, "prior_dv_upper")))
  expect_true(is.numeric(attr(res, "prior_dg_used")))
  expect_true(attr(res, "prior_de_used") >= 0.01 && attr(res, "prior_de_used") <= 0.99)
  expect_true(attr(res, "prior_dv_used") >= 0.0009 && attr(res, "prior_dv_used") <= 0.99)
  expect_true(attr(res, "prior_dv_upper") >= 0.1 && attr(res, "prior_dv_upper") <= 0.25)
  expect_true(attr(res, "prior_dg_used") >= 0.01 && attr(res, "prior_dg_used") <= 0.99)
  expect_identical(attr(res, "decision_mode"), "frequentist_primary")
  expect_identical(attr(res, "primary_de_channel"), "padj")
  expect_identical(attr(res, "primary_dv_channel"), "dv_padj")
  expect_identical(attr(res, "primary_dg_channel"), "dg_shape_padj")
  
  # Significant genes
  sig <- significantGenes(res, padj_cutoff = 0.1)
  expect_true(nrow(sig) >= 0)
})

test_that("SGCBConfig creation works correctly", {
  cfg <- defaultSGCBConfig()
  expect_s4_class(cfg, "SGCBConfig")
  expect_equal(cfg@bootB, 1000L)
  expect_equal(cfg@seed, 12345L)

  cfg2 <- SGCBConfig(bootB = 500L, seed = 42L)
  expect_equal(cfg2@bootB, 500L)
  expect_equal(cfg2@seed, 42L)
})

test_that("sgcbDE gene filtering works correctly", {
  set.seed(123)
  n_genes <- 100
  n_samples <- 6
  
  counts <- matrix(rpois(n_genes * n_samples, lambda = 100), 
                   nrow = n_genes, ncol = n_samples)
  # Add some low-expression genes
  counts[1:10, ] <- 0
  counts[11:20, ] <- 1
  
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("sample_", seq_len(n_samples))
  
  group <- rep(c("ctrl", "treat"), each = 3)
  
  # Run with default filtering (min_count=10, min_samples=2)
  res <- sgcbDE(counts, group)
  
  # Low-expression genes should be filtered out
  expect_true(nrow(res) < n_genes)
  expect_true(attr(res, "n_filtered") > 0)
})

test_that("GG Fisher information computation works correctly", {
  alpha <- c(2, 3)
  beta <- c(100, 200)
  gamma <- c(1, 1.5)
  info <- fisher_info_diag(alpha, beta, gamma, 10L)
  expect_equal(nrow(info), 2)
  expect_equal(ncol(info), 3)
  expect_true(all(info > 0))
})

test_that("GG log-likelihood computation works correctly", {
  x <- c(1, 2, 3, 4, 5)
  alpha <- rep(2, 5)
  beta <- rep(1, 5)
  gamma <- rep(1, 5)
  
  ll <- gg_loglik_vec(x, alpha, beta, gamma)
  expect_equal(length(ll), 5)
  expect_true(all(is.finite(ll)))
})

test_that("sgcbDReg GLM analysis works correctly", {
  set.seed(123)
  n_genes <- 500
  n_samples <- 6
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  counts[1:50, 4:6] <- counts[1:50, 4:6] * 3
  group <- rep(c("ctrl", "treat"), each = 3)

  res <- sgcbDReg(counts, group = group)
  expect_s3_class(res, "SGCBDRegResults")
  expect_true(nrow(res) > 0)
  expect_true("log2FoldChange" %in% names(res))
  expect_true("padj" %in% names(res))
  expect_true("p_de_post" %in% names(res))
  expect_true("p_de_only_post" %in% names(res))
  expect_true("de_post_call" %in% names(res))
  expect_true("de_fdr_call" %in% names(res))
  expect_true("model_prob_null" %in% names(res))
  expect_true("model_prob_de" %in% names(res))
  expect_true("bf_de" %in% names(res))
  expect_true(all(abs(res$model_prob_null + res$model_prob_de - 1) < 1e-8))
  expect_true(all(res$p_de_post >= 0 & res$p_de_post <= 1))
  expect_true(all(res$p_de_only_post >= 0 & res$p_de_only_post <= 1))
  expect_true(all(res$de_post_call %in% c(TRUE, FALSE)))
  expect_true(is.numeric(attr(res, "prior_de_used")))
  expect_true(attr(res, "prior_de_used") >= 0.01 && attr(res, "prior_de_used") <= 0.99)
  expect_true(sum(res$padj < 0.05, na.rm = TRUE) > 0)
})

test_that("sgcbDReg DV joint posterior works correctly", {
  set.seed(123)
  n_genes <- 300
  n_samples <- 8
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  group <- rep(c("ctrl", "treat"), each = 4)

  res <- sgcbDReg(counts, group = group, design_disp = "auto")
  expect_s3_class(res, "SGCBDRegResults")
  expect_true("p_dv_post" %in% names(res))
  expect_true("dv_post_call" %in% names(res))
  expect_true("dv_fdr_call" %in% names(res))
  expect_true("p_dv_only_post" %in% names(res))
  expect_true("p_de_dv_post" %in% names(res))
  expect_true("model_prob_de_dv" %in% names(res))
  expect_true("ll_null_dv" %in% names(res))
  expect_true("ll_full_dv" %in% names(res))
  expect_true("bf_dv_laplace" %in% names(res))
  expect_true("bf_dv_laplace_prior" %in% names(res))
  mp <- as.matrix(res[, c("model_prob_null", "model_prob_de", "model_prob_dv", "model_prob_de_dv")])
  expect_true(all(is.finite(mp)))
  expect_true(all(abs(rowSums(mp) - 1) < 1e-8))
  expect_true(all(is.finite(res$ll_null_dv)))
  expect_true(all(is.finite(res$ll_full_dv)))
  expect_true(all(res$bf_dv_laplace > 0))
  expect_true(all(res$bf_dv_laplace_prior > 0))
  expect_true(all(res$p_de_post >= 0 & res$p_de_post <= 1))
  expect_true(all(res$p_dv_post >= 0 & res$p_dv_post <= 1))
  expect_true(all(res$p_dv_only_post >= 0 & res$p_dv_only_post <= 1))
  expect_true(all(res$dv_post_call %in% c(TRUE, FALSE)))
  expect_true(all(res$p_de_dv_post >= 0 & res$p_de_dv_post <= 1))
  expect_true(is.numeric(attr(res, "prior_dv_used")))
  expect_true(is.numeric(attr(res, "prior_dv_upper")))
  expect_true(attr(res, "prior_dv_used") >= 0.0009 && attr(res, "prior_dv_used") <= 0.99)
  expect_true(attr(res, "prior_dv_upper") >= 0.1 && attr(res, "prior_dv_upper") <= 0.25)
})

test_that("sgcbContrast wrapper works correctly", {
  set.seed(123)
  n_genes <- 200
  n_samples <- 9
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 120, size = 8), nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("sample_", seq_len(n_samples))
  counts[1:30, 4:6] <- counts[1:30, 4:6] * 2
  sample_data <- data.frame(
    condition = c(rep("A", 3), rep("B", 3), rep("C", 3)),
    batch = rep(c("x", "y", "z"), 3),
    row.names = colnames(counts)
  )
  res <- sgcbContrast(counts, sample_data, group_col = "condition", contrast_levels = c("A", "B"))
  expect_s3_class(res, "SGCBResults")
  expect_identical(attr(res, "deployment_mode"), "contrast_defined_wrapper")
  expect_equal(attr(res, "contrast_levels"), c("A", "B"))
  expect_equal(attr(res, "n_wrapper_samples"), 6)
})

test_that("sgcbPairwise wrapper returns pairwise result list", {
  set.seed(123)
  n_genes <- 150
  n_samples <- 9
  counts <- matrix(rpois(n_genes * n_samples, lambda = 100), nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("sample_", seq_len(n_samples))
  sample_data <- data.frame(
    condition = c(rep("A", 3), rep("B", 3), rep("C", 3)),
    row.names = colnames(counts)
  )
  res_list <- sgcbPairwise(counts, sample_data, group_col = "condition")
  expect_s3_class(res_list, "SGCBContrastList")
  expect_equal(length(res_list), 3)
  expect_true(all(c("A_vs_B", "A_vs_C", "B_vs_C") %in% names(res_list)))
  expect_true(all(vapply(res_list, inherits, logical(1), what = "SGCBResults")))
})
