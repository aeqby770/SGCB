# =============================================================================
# Numerical verification of GG Fisher information matrix
# Strategy: Monte Carlo score covariance must match analytical Fisher
# =============================================================================

test_that("GG score matches finite-difference gradient of log-likelihood", {
  set.seed(12345)
  alpha <- 2.0; beta <- 50.0; gamma <- 1.5
  x <- beta * rgamma(1, shape = alpha)^(1 / gamma)
  x <- max(x, 1e-8)

  ana <- gg_grad_params(x, alpha, beta, gamma)
  eps_fd <- 1e-6

  ll0 <- gg_loglik_vec(x, alpha, beta, gamma)
  fd_a <- (gg_loglik_vec(x, alpha + eps_fd, beta, gamma) - ll0) / eps_fd
  fd_b <- (gg_loglik_vec(x, alpha, beta + eps_fd, gamma) - ll0) / eps_fd
  fd_g <- (gg_loglik_vec(x, alpha, beta, gamma + eps_fd) - ll0) / eps_fd

  expect_equal(ana[1, 1], fd_a, tolerance = 1e-3)
  expect_equal(ana[1, 2], fd_b, tolerance = 1e-3)
  expect_equal(ana[1, 3], fd_g, tolerance = 1e-3)
})

test_that("GG Hessian diagonal matches finite-difference second derivative", {
  set.seed(12345)
  alpha <- 2.0; beta <- 50.0; gamma <- 1.5
  x <- beta * rgamma(1, shape = alpha)^(1 / gamma)
  x <- max(x, 1e-8)

  hess <- gg_hessian_diag(x, alpha, beta, gamma)
  eps_fd <- 1e-5

  g0 <- gg_grad_params(x, alpha, beta, gamma)
  g_a <- gg_grad_params(x, alpha + eps_fd, beta, gamma)
  g_b <- gg_grad_params(x, alpha, beta + eps_fd, gamma)
  g_g <- gg_grad_params(x, alpha, beta, gamma + eps_fd)

  fd_haa <- (g_a[1, 1] - g0[1, 1]) / eps_fd
  fd_hbb <- (g_b[1, 2] - g0[1, 2]) / eps_fd
  fd_hgg <- (g_g[1, 3] - g0[1, 3]) / eps_fd

  expect_equal(hess[1, 1], fd_haa, tolerance = 1e-2)
  expect_equal(hess[1, 2], fd_hbb, tolerance = 1e-2)
  expect_equal(hess[1, 3], fd_hgg, tolerance = 1e-2)
})

test_that("fisher_info_diag matches analytical formula for all 3 diagonal elements", {
  alpha_vals <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  beta_vals  <- c(10, 50, 100, 200, 500)
  gamma_vals <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  n <- 10L

  info <- fisher_info_diag(alpha_vals, beta_vals, gamma_vals, n)

  psi0 <- digamma(alpha_vals)
  psi1 <- trigamma(alpha_vals)

  expected_aa <- n * psi1
  expected_bb <- n * alpha_vals * gamma_vals^2 / beta_vals^2
  expected_gg <- n * (1 + alpha_vals * psi0^2 + 2 * psi0 + alpha_vals * psi1) / gamma_vals^2

  expect_equal(info[, 1], expected_aa, tolerance = 1e-6)
  expect_equal(info[, 2], expected_bb, tolerance = 1e-6)
  expect_equal(info[, 3], expected_gg, tolerance = 1e-6)
})

test_that("Fisher information matches Monte Carlo score covariance", {
  set.seed(12345)
  alpha <- 2.0; beta <- 100.0; gamma <- 1.5
  n_mc <- 50000L

  u <- rgamma(n_mc, shape = alpha, rate = 1)
  x_mc <- beta * u^(1 / gamma)
  x_mc <- pmax(x_mc, 1e-8)

  scores <- gg_grad_params(x_mc,
                           rep(alpha, n_mc),
                           rep(beta, n_mc),
                           rep(gamma, n_mc))

  mc_cov <- cov(scores)

  info_diag <- fisher_info_diag(alpha, beta, gamma, 1L)

  psi0 <- digamma(alpha)
  psi1 <- trigamma(alpha)
  expected_I_aa <- psi1
  expected_I_bb <- alpha * gamma^2 / beta^2
  expected_I_gg <- (1 + alpha * psi0^2 + 2 * psi0 + alpha * psi1) / gamma^2
  expected_I_ab <- gamma / beta
  expected_I_ag <- -psi0 / gamma
  expected_I_bg <- -(1 + alpha * psi0) / beta

  expect_equal(mc_cov[1, 1], expected_I_aa, tolerance = 0.05)
  expect_equal(mc_cov[2, 2], expected_I_bb, tolerance = 0.05)
  expect_equal(mc_cov[3, 3], expected_I_gg, tolerance = 0.05)

  expect_equal(mc_cov[1, 2], expected_I_ab, tolerance = 0.05)
  expect_equal(mc_cov[1, 3], expected_I_ag, tolerance = 0.05)
  expect_equal(mc_cov[2, 3], expected_I_bg, tolerance = 0.05)

  expect_equal(info_diag[1, 1], expected_I_aa, tolerance = 1e-6)
  expect_equal(info_diag[1, 2], expected_I_bb, tolerance = 1e-6)
  expect_equal(info_diag[1, 3], expected_I_gg, tolerance = 1e-6)
})

test_that("Fisher information is correct at Weibull special case (alpha=1)", {
  alpha <- 1.0; beta <- 100.0; gamma <- 2.0
  n <- 1L
  info <- fisher_info_diag(alpha, beta, gamma, n)

  psi0 <- digamma(1)
  psi1 <- trigamma(1)

  expect_equal(info[1, 1], psi1, tolerance = 1e-6)
  expect_equal(info[1, 2], gamma^2 / beta^2, tolerance = 1e-6)
  expect_equal(info[1, 3], (1 + psi0^2 + 2 * psi0 + psi1) / gamma^2, tolerance = 1e-6)
})

test_that("Fisher information is correct at Gamma special case (gamma=1)", {
  alpha <- 3.0; beta <- 50.0; gamma <- 1.0
  n <- 1L
  info <- fisher_info_diag(alpha, beta, gamma, n)

  psi0 <- digamma(3)
  psi1 <- trigamma(3)

  expect_equal(info[1, 1], psi1, tolerance = 1e-6)
  expect_equal(info[1, 2], alpha / beta^2, tolerance = 1e-6)
  expect_equal(info[1, 3], 1 + alpha * psi0^2 + 2 * psi0 + alpha * psi1, tolerance = 1e-6)
})
