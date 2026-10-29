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
  res <- sgcbDE(counts, group, alpha = 0.1, threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  
  expect_s3_class(res, "SGCBResults")
  expect_equal(nrow(res), n_genes)
  expect_equal(res$padj, p.adjust(res$pvalue, "BH", n = nrow(res)))
  posterior <- attr(res,"mean_posterior")
  expect_lte(sum(posterior$posterior_null[res$DE])/max(1,sum(res$DE)),0.1)
  expect_true(all(res$DE == (res$padj < 0.1)))
  expect_true(all(is.finite(res$log2FoldChange)))
  # Three positive observations cannot identify the three-parameter GG fit.
  expect_true(all(!res$dv_inference_available))
  expect_true(all(is.na(res$dv_pvalue)))
  expect_true(all(is.na(res$dv_padj)))
  expect_true(all(!res$DV))
  expect_true(all(is.finite(res$DE_joint_q)))
  expect_true(all(res$DV_joint_q == 1))
  expect_identical(attr(res, "primary_de_channel"), "padj")
  expect_identical(attr(res, "primary_dv_channel"), "dv_padj")
  expect_identical(attr(res, "dv_covariance"), "paired_population_posterior")
  expect_true(all(is.finite(attr(res, "var_shrunk"))))
  # Significant genes
  sig <- significantGenes(res, padj_cutoff = 0.1)
  expect_true(nrow(sig) >= 0)
})

test_that("small feature sets and prepared abundances retain valid mean inference", {
  counts <- outer(c(2,4,7,11),c(1,2,1.5,3,2.5,1.1,2.2,1.7,3.3,2.6),"*")
  group <- factor(rep(c("A","B"),each=5))
  res <- sgcbDE(counts,group,normalization="none",offset=0,min_count=0,min_samples=1,threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  expect_true(all(is.finite(res$pvalue)))
  expect_equal(attr(res,"size_factors"),rep(1,10))
  # The common-level population changes under separate rescaling of each gene.
  # The arithmetic-mean endpoint and posterior expected error rule remain testable.
  expect_equal(log2(res$mean_treat/res$mean_ctrl),
    rep(log2(res$mean_treat[1]/res$mean_ctrl[1]),nrow(res)),tolerance=1e-7)
  expect_true(all(is.finite(res$log2FoldChange) & is.finite(res$lfcSE) & res$lfcSE>0))
  expect_equal(res$log2FC_shrunk,res$log2FoldChange)
  posterior <- attr(res,"mean_posterior")
  expect_true(all(posterior$posterior_null>=0 & posterior$posterior_null<=1))
  expect_lte(sum(posterior$posterior_null[res$DE])/max(1,sum(res$DE)),attr(res,"alpha"))
  expect_error(sgcbDE(counts,replace(group,1,NA),min_count=0))
  expect_error(sgcbDE(counts,group,min_count=1e10))
})

test_that("population CV-squared preserves group direction and the complete BH family", {
  counts <- matrix(exp(sin(seq_len(240))/2+seq_len(240)/1000),nrow=12)
  counts[1,1:7] <- 0
  group <- factor(rep(c("A","B"),each=10))
  # Replay the same existing integration state to test exact group equivariance;
  # this does not generate observations or introduce a seed experiment.
  rng_before <- get(".Random.seed",envir=.GlobalEnv)
  res <- sgcbDE(counts,group,normalization="none",min_count=0,min_samples=1,threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  rng_after <- get(".Random.seed",envir=.GlobalEnv)
  assign(".Random.seed",rng_before,envir=.GlobalEnv)
  reverse <- sgcbDE(counts,factor(group,levels=c("B","A")),normalization="none",min_count=0,min_samples=1,threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  expect_identical(get(".Random.seed",envir=.GlobalEnv),rng_after)
  raw <- attr(res,"distribution_fits")
  expect_equal(res$dv_log2_cv2_ratio,log2(res$dv_cv2_treat/res$dv_cv2_ctrl),tolerance=1e-8)
  expect_equal(res$dv_padj,p.adjust(res$dv_pvalue,"BH",n=nrow(res)))
  expect_false(res$dv_inference_available[1])
  expect_true(is.na(res$dv_pvalue[1]))
  expect_true(sum(res$dv_inference_available)>0)
  expect_equal(res$pvalue,reverse$pvalue,tolerance=1e-7)
  expect_equal(res$log2FoldChange,-reverse$log2FoldChange,tolerance=1e-7)
  pairs <- attr(res,"mean_model")$pairs
  reverse_pairs <- attr(reverse,"mean_model")$pairs
  expect_equal(nrow(unique(pairs[,c("left","right"),with=FALSE])),nrow(pairs))
  expect_equal(pairs$left,reverse_pairs$right)
  expect_equal(pairs$right,reverse_pairs$left)
  expect_equal(pairs$weight,reverse_pairs$weight)
  expect_equal(res$dv_pvalue,reverse$dv_pvalue,tolerance=1e-8)
  expect_equal(res$dv_log2_cv2_ratio,-reverse$dv_log2_cv2_ratio,tolerance=1e-8)
  expect_equal(attr(res,"fit_diagnostics")$ctrl$converged,is.finite(raw$score0) & raw$score0<1e-7)
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
  res <- sgcbDE(counts, group, threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  
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
  res <- sgcbContrast(counts, sample_data, group_col = "condition", contrast_levels = c("A", "B"), threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
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
  res_list <- sgcbPairwise(counts, sample_data, group_col = "condition", threads=as.integer(Sys.getenv("SGCB_THREADS","1")))
  expect_s3_class(res_list, "SGCBContrastList")
  expect_equal(length(res_list), 3)
  expect_true(all(c("A_vs_B", "A_vs_C", "B_vs_C") %in% names(res_list)))
  expect_true(all(vapply(res_list, inherits, logical(1), what = "SGCBResults")))
})

test_that("conditional mean contrasts target arithmetic means and permit an infinite variance", {
  z <- c(-.7,-.5,-.3,-.1,.1,.3,.5,.7)
  x <- matrix(exp(c(z+.3,z-.075)),nrow=1)
  group <- rep(0:1,each=8)
  threads <- as.integer(Sys.getenv("SGCB_THREADS","1"))
  likelihood <- SGCB:::population_mixture_mean_likelihood_cpp(x,group,
    log(c(.5,1)),c(0,0),c(0,0),c(0,0),threads)
  # .3 + .5^2/2 equals -.075 + 1^2/2; the log means differ.
  mean_test <- SGCB:::population_mixture_mean_grid_cpp(x,group,
    log(c(.5,1)),c(0,0),c(0,0),c(0,0),c(.5,.5),
    likelihood$log_likelihood,numeric(),.01,1,1e-11,threads,1L,2L,1)
  expect_equal(unname(mean_test$genes[1,"effect"]),0,tolerance=1e-8)
  expect_equal(unname(mean_test$genes[1,"se"]),sqrt((.5^2+1)/8)/log(2),tolerance=1e-8)
  expect_equal(unname(mean_test$genes[1,"p"]),1,tolerance=1e-8)
  expect_equal(unname(mean_test$genes[1,"available"]),1)

  y <- matrix(rep(exp(z+.3),2),nrow=1)
  heavy <- SGCB:::population_mixture_mean_likelihood_cpp(y,group,
    log(1.05),-.5,0,0,threads)
  expect_true(heavy$mean_finite[1])
  expect_lt(1+2*1.05*(-.5),0)
  heavy_test <- SGCB:::population_mixture_mean_grid_cpp(y,group,
    log(1.05),-.5,0,0,1,heavy$log_likelihood,numeric(),
    .01,1,1e-11,threads,1L,1L,1)
  expect_equal(unname(heavy_test$genes[1,"effect"]),0,tolerance=1e-8)
  expect_equal(unname(heavy_test$genes[1,"p"]),1,tolerance=1e-8)
})


test_that("effect profiling integrates a normalized continuous Gaussian prior", {
  grid <- seq.int(-3000L,3000L)*.01
  effect <- c(0,grid[grid!=0])
  ll <- matrix(dnorm(.125,effect,.2,log=TRUE),nrow=1)
  fit <- SGCB:::mean_effect_profile_cpp(ll,effect,.125,TRUE,1L,.125)
  exact <- dnorm(.125,.125,sqrt(.2^2+.125^2))
  expect_equal(fit$parameters$prior_grid_mass,1,tolerance=1e-10)
  expect_equal(exp(fit$parameters$log_likelihood),exact,tolerance=1e-10)
  wide <- SGCB:::mean_effect_profile_cpp(ll,effect,c(.125,4),TRUE,1L,c(.125,1))
  expect_equal(wide$profile$prior_grid_mass,rep(1,4),tolerance=1e-10)
})


test_that("CV-squared retains the structural all-positive frequency component", {
  threads <- as.integer(Sys.getenv("SGCB_THREADS","1"))
  point <- SGCB:::joint_population_cv_cpp(matrix(0,2L,1L),c(10L,10L),
    c(TRUE,TRUE),log(.5),0,1,c(10L,10L),0,threads=threads,
    moment_increment=0,pair_left=0L,pair_right=0L,pair_weights=1,
    structural_mass=1)
  continuous <- SGCB:::joint_population_cv_cpp(matrix(0,2L,1L),c(10L,10L),
    c(TRUE,TRUE),log(.5),0,1,c(10L,10L),0,threads=threads,
    moment_increment=0,pair_left=0L,pair_right=0L,pair_weights=1,
    structural_mass=0)
  expect_equal(unname(point$genes[1,c("effect","se","p")]),c(0,0,1),tolerance=1e-12)
  expect_gt(continuous$genes[1,"se"],0)
})

test_that("narrow Gaussian mean states retain the analytic null and slab densities", {
  threads <- as.integer(Sys.getenv("SGCB_THREADS","1"))
  sigma <- c(.005,.02)
  standard <- seq(-1,1,length.out=5L)
  x <- rbind(exp(rep(sigma[1]*standard,2L)),
    exp(c(sigma[1]*standard,sigma[1]*standard+sigma[1]/2)))
  group <- rep(0:1,each=5L)
  left <- c(1L,1L,2L,2L)
  right <- c(1L,2L,1L,2L)
  A <- seq(-2,2)*.125
  likelihood <- SGCB:::population_mixture_mean_likelihood_cpp(x,group,log(sigma),
    c(0,0),c(0,0),c(0,0),threads)
  fit <- SGCB:::mean_full_two_group_dictionary_cpp(x,group,log(sigma),
    c(0,0),c(0,0),c(0,0),sigma^2/2,likelihood$log_likelihood,left,right,
    A,.125,.125,4,64L,1e-12,1e-12,1,threads)
  location0 <- rowMeans(log(x[,group==0,drop=FALSE]))/log(2)
  location1 <- rowMeans(log(x[,group==1,drop=FALSE]))/log(2)
  midpoint <- (location0+location1)/2+sigma[1]^2/(2*log(2))
  difference <- location1-location0
  variance <- sigma[1]^2/(5*log(2)^2)
  common <- matrix(dnorm(as.vector(outer(midpoint,A/log(2),"-")),
    sd=sqrt((.125/log(2))^2+variance/2)),nrow=nrow(x))
  loglike <- likelihood$log_likelihood
  peak <- matrixStats::rowMaxs(loglike[1:2,,drop=FALSE])+
    matrixStats::rowMaxs(loglike[3:4,,drop=FALSE])
  probability <- exp(loglike[1:2,1]+loglike[3:4,1]-peak)
  null <- sweep(common,1,probability*dnorm(difference,sd=sqrt(2*variance)),"*")
  slab <- sweep(common,1,probability*dnorm(difference,sd=sqrt(2*variance+16)),"*")
  nA <- length(A)
  nK <- length(left)*nA
  expect_lt(max(abs(fit$likelihood[,seq_len(nA)]-null))/(1+max(null)),1e-12)
  expect_lt(max(abs(fit$likelihood[,nK+seq_len(nA)]-slab))/(1+max(slab)),1e-12)
  expect_lt(max(fit$checks[,"eta_CF_zero_error"]),1e-14)
  expect_lt(max(fit$checks[,"eta_CF_modulus_excess"]),1e-14)
})


test_that("symmetric paired likelihood agrees with its explicit mixture", {
  threads <- as.integer(Sys.getenv("SGCB_THREADS","1"))
  lhs <- matrix(c(.2,.7,.4,.9,.3,.5,.3,.8,.6,.8,.5,.2),nrow=4,byrow=TRUE)
  rhs <- matrix(c(.7,.2,.6,.4,.8,.2,.5,.4,.9,.2,.6,.7),nrow=4,byrow=TRUE)
  pairs <- data.table::CJ(left=1:3,right=1:3)
  ll <- rbind(log(lhs),log(rhs))
  fit <- SGCB:::fit_paired_population_cpp(ll,rep(TRUE,4),rep(1/9,9),
    mixsqp::mixsqp,threads,200L,1e-10,symmetric=TRUE)
  L <- .5*(lhs[,pairs$left]*rhs[,pairs$right]+lhs[,pairs$right]*rhs[,pairs$left])
  reference <- SGCB:::fit_population_dictionary_cpp(L,mixsqp::mixsqp,
    threads,200L,1e-10,reuse_buffer=TRUE,input_is_log=FALSE)
  W <- matrix(fit$weights,3,3,byrow=TRUE)
  expect_true(fit$converged)
  expect_equal(W,t(W),tolerance=1e-12)
  expect_equal(drop(L%*%fit$weights),drop(L%*%reference$weights),tolerance=1e-8)
})
