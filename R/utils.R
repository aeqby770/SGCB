# =============================================================================
# SGCB: Utility functions module
# Retains only internal DV/DG test functions required by the sgcbDE() pipeline
# =============================================================================

# -----------------------------------------------------------------------------
# GG variance and DV/DG utilities
# -----------------------------------------------------------------------------
.ggVariance <- function(alpha, beta, gamma, eps = 1e-8) {
  a <- pmax(alpha, eps)
  b <- pmax(beta, eps)
  g <- pmax(gamma, eps)
  g0 <- base::gamma(a)
  g1 <- base::gamma(a + 1 / g)
  g2 <- base::gamma(a + 2 / g)
  mean_val <- b * g1 / g0
  var_val <- b^2 * (g2 / g0 - (g1 / g0)^2)
  pmax(var_val, eps)
}

.ggLogVarGrad <- function(alpha, beta, gamma, eps = 1e-6) {
  step_a <- pmax(abs(alpha) * 1e-5, eps)
  step_b <- pmax(abs(beta) * 1e-5, eps)
  step_g <- pmax(abs(gamma) * 1e-5, eps)
  a_lo <- pmax(alpha - step_a, 0.01)
  a_hi <- alpha + step_a
  b_lo <- pmax(beta - step_b, eps)
  b_hi <- beta + step_b
  g_lo <- pmax(gamma - step_g, 0.01)
  g_hi <- gamma + step_g
  log_var <- log(.ggVariance(alpha, beta, gamma, eps))
  grad_a <- (log(.ggVariance(a_hi, beta, gamma, eps)) - log(.ggVariance(a_lo, beta, gamma, eps))) / (a_hi - a_lo)
  grad_b <- (log(.ggVariance(alpha, b_hi, gamma, eps)) - log(.ggVariance(alpha, b_lo, gamma, eps))) / (b_hi - b_lo)
  grad_g <- (log(.ggVariance(alpha, beta, g_hi, eps)) - log(.ggVariance(alpha, beta, g_lo, eps))) / (g_hi - g_lo)
  list(log_var = log_var, grad_alpha = grad_a, grad_beta = grad_b, grad_gamma = grad_g)
}

.ggLogVarSE <- function(alpha, beta, gamma, n_samples, eps = 1e-8) {
  info <- fisher_info_diag(alpha, beta, gamma, as.integer(n_samples))
  var_alpha <- 1 / pmax(info[, 1], eps)
  var_beta <- 1 / pmax(info[, 2], eps)
  var_gamma <- 1 / pmax(info[, 3], eps)
  grad <- .ggLogVarGrad(alpha, beta, gamma, eps)
  var_logvar <- grad$grad_alpha^2 * var_alpha +
    grad$grad_beta^2 * var_beta +
    grad$grad_gamma^2 * var_gamma
  list(log_var = grad$log_var, var_logvar = pmax(var_logvar, eps))
}

.sgcbDVTest <- function(alpha_ctrl, beta_ctrl, gamma_ctrl,
                        alpha_treat, beta_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8) {
  ctrl <- .ggLogVarSE(alpha_ctrl, beta_ctrl, gamma_ctrl, n_ctrl, eps)
  treat <- .ggLogVarSE(alpha_treat, beta_treat, gamma_treat, n_treat, eps)
  diff_log <- treat$log_var - ctrl$log_var
  se_diff <- sqrt(ctrl$var_logvar + treat$var_logvar)
  se_diff <- pmax(se_diff, eps)
  z_stat <- diff_log / se_diff
  p_val <- 2 * stats::pnorm(-abs(z_stat))
  list(
    log2_var_ratio = diff_log / log(2),
    stat = z_stat,
    pvalue = p_val,
    var_ctrl = exp(ctrl$log_var),
    var_treat = exp(treat$log_var)
  )
}

.sgcbDGTest <- function(alpha_ctrl, gamma_ctrl,
                        alpha_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8) {
  info_ctrl <- fisher_info_diag(alpha_ctrl, rep(1, length(alpha_ctrl)), gamma_ctrl, as.integer(n_ctrl))
  info_treat <- fisher_info_diag(alpha_treat, rep(1, length(alpha_treat)), gamma_treat, as.integer(n_treat))
  var_alpha_ctrl <- 1 / pmax(info_ctrl[, 1], eps)
  var_gamma_ctrl <- 1 / pmax(info_ctrl[, 3], eps)
  var_alpha_treat <- 1 / pmax(info_treat[, 1], eps)
  var_gamma_treat <- 1 / pmax(info_treat[, 3], eps)
  log_alpha_ratio <- log(alpha_treat / alpha_ctrl)
  log_gamma_ratio <- log(gamma_treat / gamma_ctrl)
  se_alpha <- sqrt(var_alpha_treat / (alpha_treat^2) + var_alpha_ctrl / (alpha_ctrl^2))
  se_gamma <- sqrt(var_gamma_treat / (gamma_treat^2) + var_gamma_ctrl / (gamma_ctrl^2))
  se_alpha <- pmax(se_alpha, eps)
  se_gamma <- pmax(se_gamma, eps)
  z_alpha <- log_alpha_ratio / se_alpha
  z_gamma <- log_gamma_ratio / se_gamma
  list(
    alpha_log2_ratio = log_alpha_ratio / log(2),
    alpha_stat = z_alpha,
    alpha_pvalue = 2 * stats::pnorm(-abs(z_alpha)),
    gamma_log2_ratio = log_gamma_ratio / log(2),
    gamma_stat = z_gamma,
    gamma_pvalue = 2 * stats::pnorm(-abs(z_gamma))
  )
}
