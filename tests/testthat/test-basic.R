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
  expect_true("SGCB_Score" %in% names(res))
  expect_true("pvalue_dd" %in% names(res))
  
  # Significant genes
  sig <- significantGenes(res, padj_cutoff = 0.1)
  expect_true(nrow(sig) >= 0)
})

test_that("SGCBConfig creation works correctly", {
  # Default configuration
  cfg <- defaultSGCBConfig()
  expect_s4_class(cfg, "SGCBConfig")
  expect_equal(cfg@hiddenWidth, 5000L)
  expect_equal(cfg@learningRate, 0.001)
  
  # Custom configuration
  cfg2 <- SGCBConfig(hiddenWidth = 2000L, maxIter = 100L)
  expect_equal(cfg2@hiddenWidth, 2000L)
  expect_equal(cfg2@maxIter, 100L)
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

test_that("Rcpp functions work correctly", {
  # Test clip_vec
  x <- c(-2, -1, 0, 1, 2, 3)
  clipped <- clip_vec(x, 0, 2)
  expect_equal(clipped, c(0, 0, 0, 1, 2, 2))
  
  # Test softplus_mat
  mat <- matrix(c(0, 1, 2, 20, 30, 40), nrow = 2, ncol = 3)
  sp <- softplus_mat(mat)
  expect_true(all(sp > 0))
  expect_true(all(sp >= mat))
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

test_that("Bootstrap functions work correctly", {
  set.seed(123)
  
  # Test resampling indices
  indices <- resample_indices(10, 5, 3)
  expect_equal(dim(indices), c(5, 3))
  expect_true(all(indices >= 0 & indices < 10))
  
  # Test BH adjustment
  pvals <- c(0.01, 0.02, 0.03, 0.1, 0.5)
  padj <- bh_adjust(pvals)
  expect_equal(length(padj), 5)
  expect_true(all(padj >= pvals))
  expect_true(all(padj <= 1))
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
