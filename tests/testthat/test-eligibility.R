test_that("constant positive values are ineligible for population DV", {
  group <- rep(0:1, each=50L)
  x <- matrix(3, 1L, 100L)
  pure <- SGCB:::joint_population_likelihood_cpp(x, group, 0, 0, threads=1L)
  mixture <- SGCB:::population_mixture_likelihood_cpp(x, group, 0, 0, 0, 0, threads=1L)
  mean <- SGCB:::population_mixture_mean_likelihood_cpp(x, group, 0, 0, 0, 0, threads=1L)
  expect_identical(pure$positive, c(50L, 50L))
  expect_false(any(pure$usable))
  expect_false(any(mixture$usable))
  expect_true(all(mean$usable))
  expect_false(any(mean$cv_usable))
})

test_that("sparse and dense likelihoods use the same exact variation rule", {
  group <- rep(0:1, each=50L)
  x <- rbind(constant=rep(3, 100L), variable=20+seq_len(100L))
  connection <- rawConnection(raw(0), "w+b")
  sparse <- SGCB:::population_mixture_likelihood_cpp(x, group, 0, 0, 0, 0,
    threads=1L, sparse_output=TRUE, connection=connection, sparse_block_rows=2L)
  close(connection)
  dense <- SGCB:::population_mixture_likelihood_cpp(x, group, 0, 0, 0, 0, threads=1L)
  expect_identical(sparse$usable, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(sparse$usable, dense$usable)
})

test_that("unestimable constant inputs return retained unavailable DV rows", {
  group <- rep(0:1, each=50L)
  x <- rbind(constant=rep(20, 100L),
    one_variable_group=c(rep(20, 50L), 20+seq_len(50L)),
    sparse=c(rep(20, 3L), rep(0, 47L), rep(20, 50L)))
  result <- sgcbDE(x, group, endpoint="variability", normalization="none", threads=1L)
  expect_identical(result$gene_id, rownames(x))
  expect_false(any(result$dv_inference_available))
  expect_true(all(is.na(result$dv_pvalue)))
  expect_true(all(is.na(result$dv_padj)))
  expect_false(any(result$DV))
  constant3 <- sgcbDE(matrix(3, 1L, 100L), group, endpoint="variability",
    normalization="none", min_count=1, min_samples=2, threads=1L)
  expect_equal(nrow(constant3), 1L)
  expect_false(constant3$dv_inference_available)
  expect_false(constant3$DV)
})
