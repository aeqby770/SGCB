# =============================================================================
# SGCB: Utility functions module
# DV/DG test functions for the sgcbDE() pipeline
#
# v2 (2026-03-26): finite-sample corrections for small-sample DV/DG tests
#   DV: conditional alpha Schur complement after profiling beta and fixing gamma
#   DG-alpha: same conditional alpha Schur complement
#   DG-gamma: conditional diagonal Fisher approximation, empirical-null calibrated
#   Concordance guard: Brown-Forsythe support and MAD-scale null calibration
# =============================================================================

# -----------------------------------------------------------------------------
# GG moment functions (unchanged)
# -----------------------------------------------------------------------------
.ggVariance <- function(alpha, beta, gamma, eps = 1e-8) {
  a <- pmax(alpha, eps)
  b <- pmax(beta, eps)
  g <- pmax(gamma, eps)
  g0 <- base::gamma(a)
  g1 <- base::gamma(a + 1 / g)
  g2 <- base::gamma(a + 2 / g)
  var_val <- b^2 * (g2 / g0 - (g1 / g0)^2)
  pmax(var_val, eps)
}

.ggMean <- function(alpha, beta, gamma, eps = 1e-8) {
  a <- pmax(alpha, eps)
  b <- pmax(beta, eps)
  g <- pmax(gamma, eps)
  g0 <- base::gamma(a)
  g1 <- base::gamma(a + 1 / g)
  pmax(b * g1 / g0, eps)
}

.ggEntropy <- function(alpha, beta, gamma, eps = 1e-8) {
  a <- pmax(alpha, 0.05)
  b <- pmax(beta, eps)
  g <- pmax(gamma, 0.05)
  psi_a <- digamma(a)
  a + log(b) + lgamma(a) - (a - 1 / g) * psi_a - log(g)
}

.shiftedGammaRightTail <- function(stat_sq, eps = 1e-8) {
  x <- pmax(as.numeric(stat_sq), 0)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 25) {
    return(pmin(pmax(stats::pchisq(stat_sq, df = 1, lower.tail = FALSE), eps), 1))
  }

  q_lo <- as.numeric(stats::quantile(x, 0.10, names = FALSE, type = 8, na.rm = TRUE))
  q_hi <- as.numeric(stats::quantile(x, 0.80, names = FALSE, type = 8, na.rm = TRUE))
  bulk <- x[x >= q_lo & x <= q_hi]
  if (length(bulk) < 20) {
    bulk <- x
  }

  m <- mean(bulk)
  v <- stats::var(bulk)
  if (!is.finite(v) || v <= eps || !is.finite(m)) {
    return(pmin(pmax(stats::pchisq(stat_sq, df = 1, lower.tail = FALSE), eps), 1))
  }

  sdv <- sqrt(v)
  skew <- mean((bulk - m)^3) / pmax(sdv^3, eps)

  k_skew <- ifelse(is.finite(skew) & skew > 0.2, (2 / skew)^2, NA_real_)
  k_mom <- pmax((m^2) / pmax(v, eps), 0.2)
  k <- ifelse(is.finite(k_skew) & k_skew > 0.2, k_skew, k_mom)
  k <- pmax(k, 0.2)

  theta <- pmax(sqrt(v / k), eps)
  delta <- m - k * theta
  delta <- pmin(delta, min(bulk) - eps)

  y <- pmax(stat_sq - delta, eps)
  p <- stats::pgamma(y, shape = k, scale = theta, lower.tail = FALSE)
  pmin(pmax(p, eps), 1)
}

.empiricalRightTail <- function(stat_sq, eps = 1e-8) {
  x <- pmax(as.numeric(stat_sq), 0)
  x[!is.finite(x)] <- 0
  n <- length(x)
  rk <- rank(-x, ties.method = "average")
  p <- rk / (n + 1)
  pmin(pmax(p, eps), 1)
}

.scaledChisqRightTail <- function(stat_sq, eps = 1e-8) {
  x_full <- pmax(as.numeric(stat_sq), 0)
  x <- x_full[is.finite(x_full)]
  n <- length(x)
  if (n < 30) {
    return(pmin(pmax(stats::pchisq(stat_sq, df = 1, lower.tail = FALSE), eps), 1))
  }

  qx <- as.numeric(stats::quantile(x, c(0.25, 0.75), names = FALSE, type = 8, na.rm = TRUE))
  ratio_obs <- qx[2] / pmax(qx[1], eps)

  ratio_fun <- function(df) {
    qq <- stats::qchisq(c(0.25, 0.75), df = df)
    qq[2] / pmax(qq[1], eps)
  }

  ratio_lo <- ratio_fun(0.5)
  ratio_hi <- ratio_fun(30)
  ratio_clamped <- min(max(ratio_obs, ratio_hi), ratio_lo)

  obj <- function(df) {
    rr <- ratio_fun(df)
    (rr - ratio_clamped)^2
  }
  fit <- stats::optimize(obj, interval = c(0.5, 30))
  df_hat <- fit$minimum

  qq <- stats::qchisq(c(0.25, 0.75), df = df_hat)
  scale_hat <- (qx[2] - qx[1]) / pmax(qq[2] - qq[1], eps)
  scale_hat <- pmax(scale_hat, eps)

  p <- stats::pchisq(pmax(stat_sq, 0) / scale_hat, df = df_hat, lower.tail = FALSE)
  pmin(pmax(p, eps), 1)
}

# -----------------------------------------------------------------------------
# Analytic gradient of log(CV²) w.r.t. (α, β, γ)
#
# CV² = Var/μ² = Γ(α+2/γ)Γ(α)/Γ(α+1/γ)² − 1
# CV² is β-invariant: its scale derivative is exactly zero.
# -----------------------------------------------------------------------------
.ggLogCV2Grad <- function(alpha, beta, gamma, eps = 1e-6) {
  log_moment_ratio <- lgamma(alpha + 2 / gamma) + lgamma(alpha) -
    2 * lgamma(alpha + 1 / gamma)
  log_cv2 <- log_moment_ratio + log(-expm1(-log_moment_ratio))
  derivative_factor <- 1 / (-expm1(-log_moment_ratio))
  grad_a <- derivative_factor * (digamma(alpha + 2 / gamma) + digamma(alpha) -
                                  2 * digamma(alpha + 1 / gamma))
  grad_g <- derivative_factor * 2 * (digamma(alpha + 1 / gamma) -
                                     digamma(alpha + 2 / gamma)) / gamma^2
  list(log_cv2 = log_cv2, grad_alpha = grad_a,
       grad_beta = rep(0, length(alpha)), grad_gamma = grad_g)
}

# -----------------------------------------------------------------------------
# Model-based delta variance for an interior, regular, three-parameter GG MLE.
# The inverse Schur complement equals the (α,γ) block of the full inverse.
# It propagates estimated γ uncertainty even though the β derivative is zero.
# No generic n/(n-p) multiplier is asserted to be a Bartlett correction.
# This covariance alone does not validate penalized/boundary fits, model
# misspecification, or the downstream empirical calibration in .sgcbDVTest.
# -----------------------------------------------------------------------------
.ggLogCV2SE_v2 <- function(alpha, beta, gamma, n_samples, eps = 1e-8) {
  grad <- .ggLogCV2Grad(alpha, beta, gamma, eps)
  covariance <- fisher_info_schur_ag(alpha, beta, gamma, as.integer(n_samples))
  var_logcv2 <- grad$grad_alpha^2 * covariance[, 1] +
    grad$grad_gamma^2 * covariance[, 2] +
    2 * grad$grad_alpha * grad$grad_gamma * covariance[, 3]
  list(log_cv2 = grad$log_cv2, var_logcv2 = var_logcv2)
}

# -----------------------------------------------------------------------------
# Brown-Forsythe test on log-normalized abundance
# Fully vectorized matrix operations, no loops/apply.
# Compares absolute deviations from group medians, using the approximate
# F reference distribution. This log-scale target is not identical to CV².
# -----------------------------------------------------------------------------
.brownForsytheDV <- function(norm_counts, ctrl_idx, treat_idx,
                             n_ctrl, n_treat, eps = 1e-8) {
  n_genes <- nrow(norm_counts)
  n_total <- n_ctrl + n_treat

  log_counts <- log(norm_counts)
  ctrl_mat <- log_counts[, ctrl_idx, drop = FALSE]
  treat_mat <- log_counts[, treat_idx, drop = FALSE]

  # Brown-Forsythe median centers; means would instead give Levene's test.
  cen_ctrl <- matrixStats::rowMedians(ctrl_mat)
  cen_treat <- matrixStats::rowMedians(treat_mat)

  # Absolute deviations from group center
  dev_ctrl <- abs(ctrl_mat - cen_ctrl)
  dev_treat <- abs(treat_mat - cen_treat)

  # Brown-Forsythe F-statistic (fully vectorized)
  md_ctrl <- rowMeans(dev_ctrl)
  md_treat <- rowMeans(dev_treat)
  grand <- (n_ctrl * md_ctrl + n_treat * md_treat) / n_total

  ss_b <- n_ctrl * (md_ctrl - grand)^2 + n_treat * (md_treat - grand)^2
  ss_w <- rowSums((dev_ctrl - md_ctrl)^2) + rowSums((dev_treat - md_treat)^2)
  ms_w <- ss_w / pmax(n_total - 2, 1)
  f_stat <- ss_b / pmax(ms_w, eps)

  p_bf <- stats::pf(f_stat, 1, n_total - 2, lower.tail = FALSE)
  pmin(pmax(p_bf, 1e-16), 1.0)
}

.dvBootstrapP <- function(norm_counts, ctrl_idx, treat_idx,
                          alpha_ctrl, beta_ctrl, gamma_ctrl,
                          alpha_treat, beta_treat, gamma_treat,
                          n_boot = 40L,
                          eps = 1e-8) {
  n_genes <- nrow(norm_counts)
  n_ctrl <- length(ctrl_idx)
  n_treat <- length(treat_idx)
  n_total <- n_ctrl + n_treat

  log_counts <- log(norm_counts)
  ctrl_mat <- log_counts[, ctrl_idx, drop = FALSE]
  treat_mat <- log_counts[, treat_idx, drop = FALSE]
  cen_ctrl <- matrixStats::rowMedians(ctrl_mat)
  cen_treat <- matrixStats::rowMedians(treat_mat)
  dev_ctrl <- abs(ctrl_mat - cen_ctrl)
  dev_treat <- abs(treat_mat - cen_treat)
  md_ctrl <- rowMeans(dev_ctrl)
  md_treat <- rowMeans(dev_treat)
  grand <- (n_ctrl * md_ctrl + n_treat * md_treat) / n_total
  ss_b <- n_ctrl * (md_ctrl - grand)^2 + n_treat * (md_treat - grand)^2
  ss_w <- rowSums((dev_ctrl - md_ctrl)^2) + rowSums((dev_treat - md_treat)^2)
  ms_w <- ss_w / pmax(n_total - 2, 1)
  f_obs <- ss_b / pmax(ms_w, eps)

  alpha0 <- (alpha_ctrl + alpha_treat) / 2
  beta0 <- (beta_ctrl + beta_treat) / 2
  gamma0 <- (gamma_ctrl + gamma_treat) / 2

  f_null <- matrix(0, nrow = n_genes, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    sim <- gg_sample_mat(n_genes, as.integer(n_total), alpha0, beta0, gamma0, eps)
    sim <- log(sim + eps)
    s_ctrl <- sim[, seq_len(n_ctrl), drop = FALSE]
    s_treat <- sim[, n_ctrl + seq_len(n_treat), drop = FALSE]
    sc <- matrixStats::rowMedians(s_ctrl)
    st <- matrixStats::rowMedians(s_treat)
    dc <- abs(s_ctrl - sc)
    dt <- abs(s_treat - st)
    mc <- rowMeans(dc)
    mt <- rowMeans(dt)
    g <- (n_ctrl * mc + n_treat * mt) / n_total
    ssb <- n_ctrl * (mc - g)^2 + n_treat * (mt - g)^2
    ssw <- rowSums((dc - mc)^2) + rowSums((dt - mt)^2)
    msw <- ssw / pmax(n_total - 2, 1)
    f_null[, b] <- ssb / pmax(msw, eps)
  }

  (rowSums(f_null >= f_obs) + 1) / (n_boot + 1)
}

# -----------------------------------------------------------------------------
# Concordant DV evidence: full GG delta-method uncertainty and median BF.
# The two components target normalized-abundance CV2 and log-scale spread,
# respectively. Both must support a call. For the primary log-spread null,
# max(p_GG, p_BF) >= p_BF preserves the level of the BF anchor regardless of
# GG fit adequacy. It is not a general test of the abundance-CV2 null.
# -----------------------------------------------------------------------------
.sgcbDVTest <- function(alpha_ctrl, beta_ctrl, gamma_ctrl,
                        alpha_treat, beta_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8,
                        norm_counts = NULL, ctrl_idx = NULL, treat_idx = NULL) {
  ctrl <- .ggLogCV2SE_v2(alpha_ctrl, beta_ctrl, gamma_ctrl, n_ctrl, eps)
  treat <- .ggLogCV2SE_v2(alpha_treat, beta_treat, gamma_treat, n_treat, eps)
  diff_log <- treat$log_cv2 - ctrl$log_cv2
  se_model <- sqrt(ctrl$var_logcv2 + treat$var_logcv2)
  z_stat <- diff_log / se_model
  p_model <- pmin(pmax(2 * stats::pnorm(-abs(z_stat)), 1e-16), 1)
  p_bf <- .brownForsytheDV(norm_counts, ctrl_idx, treat_idx, n_ctrl, n_treat, eps)
  p_val <- pmax(p_model, p_bf)
  var_ctrl <- .ggVariance(alpha_ctrl, beta_ctrl, gamma_ctrl, eps)
  var_treat <- .ggVariance(alpha_treat, beta_treat, gamma_treat, eps)
  mean_ctrl <- .ggMean(alpha_ctrl, beta_ctrl, gamma_ctrl, eps)
  mean_treat <- .ggMean(alpha_treat, beta_treat, gamma_treat, eps)
  list(
    log2_var_ratio = (log(var_treat) - log(var_ctrl)) / log(2),
    log2_cv2_ratio = diff_log / log(2),
    se_log2_cv2 = se_model / log(2),
    stat = z_stat, pvalue = p_val,
    p_model_raw = p_model, p_model = p_model, p_bf = p_bf,
    p_mix = p_val,
    var_ctrl = var_ctrl, var_treat = var_treat,
    cv2_ctrl = var_ctrl / pmax(mean_ctrl^2, eps),
    cv2_treat = var_treat / pmax(mean_treat^2, eps)
  )
}

# -----------------------------------------------------------------------------
# Shape diagnostics propagate the full alpha-gamma covariance after
# profiling beta. The gamma coordinate is the remaining GG coordinate when
# mean and CV2 are specified; alpha alone is not a shape-only coordinate.
# Diagnostic p-values use asymptotic working-model inference.
# -----------------------------------------------------------------------------
.sgcbDGTest <- function(alpha_ctrl, gamma_ctrl,
                        alpha_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8) {
  a0 <- pmax(alpha_ctrl, eps)
  a1 <- pmax(alpha_treat, eps)
  u0 <- trigamma(a0) - 1 / a0
  u1 <- trigamma(a1) - 1 / a1
  v0 <- 1 + a0 * trigamma(a0) - 1 / a0
  v1 <- 1 + a1 * trigamma(a1) - 1 / a1
  d0 <- n_ctrl * (u0 * v0 - 1 / a0^2)
  d1 <- n_treat * (u1 * v1 - 1 / a1^2)
  va <- v0 / (a0^2 * d0) + v1 / (a1^2 * d1)
  vg <- u0 / d0 + u1 / d1
  cag <- -1 / (a0^2 * d0) - 1 / (a1^2 * d1)
  da <- log(a1 / a0)
  dg <- log(pmax(gamma_treat, eps) / pmax(gamma_ctrl, eps))
  za <- da / sqrt(va)
  zg <- dg / sqrt(vg)
  joint_stat <- (vg * da^2 + va * dg^2 - 2 * cag * da * dg) / (va * vg - cag^2)
  ent_diff <- .ggEntropy(a1, rep(1, length(a1)), gamma_treat, eps) -
    .ggEntropy(a0, rep(1, length(a0)), gamma_ctrl, eps)
  list(
    alpha_log2_ratio = da / log(2), alpha_se_log2 = sqrt(va) / log(2),
    alpha_stat = za, alpha_pvalue = 2 * stats::pnorm(-abs(za)),
    gamma_log2_ratio = dg / log(2), gamma_se_log2 = sqrt(vg) / log(2),
    gamma_stat = zg, gamma_pvalue = 2 * stats::pnorm(-abs(zg)),
    entropy_diff = ent_diff, entropy_stat = rep(NA_real_, length(da)),
    entropy_pvalue = rep(NA_real_, length(da)),
    shape_stat = zg^2, shape_pvalue = 2 * stats::pnorm(-abs(zg)),
    joint_stat = joint_stat,
    joint_pvalue = stats::pchisq(joint_stat, 2, lower.tail = FALSE)
  )
}

# -----------------------------------------------------------------------------
# SGCB joint posterior helpers
# -----------------------------------------------------------------------------
.sgcbClipProb <- function(x, eps = 1e-12) {
  p <- as.numeric(x)
  p[!is.finite(p)] <- 1
  p[p < eps] <- eps
  p[p > 1 - eps] <- 1 - eps
  p
}

.sgcbClipSE <- function(x, eps = 1e-8) {
  s <- as.numeric(x)
  s[!is.finite(s)] <- NA_real_
  med_s <- stats::median(s[is.finite(s)], na.rm = TRUE)
  med_s <- ifelse(is.finite(med_s), med_s, 1)
  s[!is.finite(s)] <- med_s
  s[s < eps] <- eps
  s
}

.sgcbABF1 <- function(effect, se, prior_sd = 1.0, eps = 1e-12) {
  b <- as.numeric(effect)
  v <- .sgcbClipSE(se, sqrt(eps))^2
  w <- pmax(as.numeric(prior_sd)^2, eps)
  if (length(w) == 1L) {
    w <- rep(w, length(v))
  }
  log_bf <- 0.5 * (log(v) - log(v + w)) + (b^2 * w) / (2 * v * (v + w))
  log_bf <- pmin(log_bf, 700)
  bf <- exp(log_bf)
  bf[!is.finite(bf)] <- exp(700)
  pmax(bf, eps)
}

.sgcbABFDiag <- function(effect_mat, se_mat, prior_sd, eps = 1e-12) {
  b <- as.matrix(effect_mat)
  s <- as.matrix(se_mat)
  stopifnot(nrow(b) == nrow(s), ncol(b) == ncol(s))
  k <- ncol(b)
  psd <- as.numeric(prior_sd)
  if (length(psd) == 1L) {
    psd <- rep(psd, k)
  }
  stopifnot(length(psd) == k)
  v <- pmax(.sgcbClipSE(as.vector(s), sqrt(eps)), sqrt(eps))^2
  v <- matrix(v, nrow = nrow(s), ncol = ncol(s))
  w <- matrix(rep(pmax(psd^2, eps), each = nrow(b)), nrow = nrow(b), ncol = k)
  b2 <- b^2
  b2[!is.finite(b2)] <- 0
  log_bf <- rowSums(0.5 * (log(v) - log(v + w)) + (b2 * w) / (2 * v * (v + w)))
  log_bf <- pmin(log_bf, 700)
  bf <- exp(log_bf)
  bf[!is.finite(bf)] <- exp(700)
  pmax(bf, eps)
}

.sgcbEmpiricalNullInflation <- function(effect, se, eps = 1e-8) {
  b <- as.numeric(effect)
  s <- .sgcbClipSE(se, sqrt(eps))
  z <- b / s
  z[!is.finite(z)] <- 0
  z_med <- stats::median(z, na.rm = TRUE)
  z_mad <- stats::median(abs(z - z_med), na.rm = TRUE) / 0.6745
  infl <- ifelse(is.finite(z_mad), z_mad, 1)
  pmax(infl, 1)
}

.sgcbCalibratedSE <- function(effect, se, eps = 1e-8) {
  s <- .sgcbClipSE(se, sqrt(eps))
  infl <- .sgcbEmpiricalNullInflation(effect, s, eps)
  list(
    se_cal = s * infl,
    inflation = infl
  )
}

.sgcbPosteriorCalls <- function(post_prob, target_fdr = 0.1, eps = 1e-12) {
  p <- .sgcbClipProb(post_prob, eps)
  fdr_target <- pmin(pmax(as.numeric(target_fdr), eps), 1 - eps)
  ord <- order(p, decreasing = TRUE)
  p_ord <- p[ord]
  cum_fdr <- cumsum(1 - p_ord) / seq_along(p_ord)
  idx <- which(cum_fdr <= fdr_target)
  has_sel <- length(idx) > 0
  thr <- ifelse(has_sel, p_ord[max(idx)], 1 + eps)
  call <- p >= thr
  realized <- ifelse(has_sel, cum_fdr[max(idx)], 0)

  list(
    call = call,
    threshold = thr,
    realized_fdr = realized,
    n_call = sum(call)
  )
}

.sgcbEstimatePriorFromP <- function(p, lower = 0.01, upper = 0.99, eps = 1e-12) {
  p_safe <- .sgcbClipProb(p, eps)
  n <- length(p_safe)
  thr <- c(0.001, 0.005, 0.01, 0.02)
  obs <- vapply(thr, function(t) mean(p_safe <= t), numeric(1))
  excess <- pmax(obs - thr, 0)
  pi_excess <- pmin(1, max(excess / pmax(1 - thr, eps)))
  p_bh <- stats::p.adjust(p_safe, method = "BH")
  pi_bh <- mean(p_bh <= 0.1)
  null95 <- stats::qbinom(0.95, size = n, prob = 0.01) / pmax(n, 1)
  pi_binom <- pmax(obs[3] - null95, 0) / 0.99
  pi_binom <- pmin(pmax(pi_binom, 0), 1)
  pi_hat <- max(pi_excess, pi_bh, pi_binom)
  pmin(pmax(pi_hat, lower), upper)
}

.sgcbEstimatePriorFromBF <- function(bf, lower = 0.01, upper = 0.99, eps = 1e-12,
                                     prior_alpha = 1.0, prior_beta = 9.0) {
  bf_safe <- as.numeric(bf)
  bf_safe[!is.finite(bf_safe)] <- 1
  bf_safe <- pmax(bf_safe, eps)
  a0 <- pmax(as.numeric(prior_alpha), eps)
  b0 <- pmax(as.numeric(prior_beta), eps)

  obj <- function(pi0) {
    ll <- sum(log1p(pi0 * (bf_safe - 1)))
    lp <- (a0 - 1) * log(pi0) + (b0 - 1) * log1p(-pi0)
    -(ll + lp)
  }
  fit <- stats::optimize(obj, interval = c(lower, upper))
  pi_hat <- fit$minimum
  pmin(pmax(pi_hat, lower), upper)
}

.sgcbResolvePrior <- function(prior, p = NULL, bf = NULL, lower = 0.01, upper = 0.99, eps = 1e-12,
                              prior_alpha = 1.0, prior_beta = 9.0) {
  is_auto <- is.character(prior) && length(prior) == 1L && identical(prior, "auto")
  prior_raw <- if (is_auto) {
    if (is.null(bf)) {
      .sgcbEstimatePriorFromP(p, lower, upper, eps)
    } else {
      .sgcbEstimatePriorFromBF(
        bf = bf,
        lower = lower,
        upper = upper,
        eps = eps,
        prior_alpha = prior_alpha,
        prior_beta = prior_beta
      )
    }
  } else {
    as.numeric(prior)
  }
  if (!is.finite(prior_raw)) prior_raw <- lower
  pmin(pmax(prior_raw, lower), upper)
}

.sgcbLRBF <- function(lr_stat, n_obs, df, eps = 1e-12) {
  lr <- pmax(as.numeric(lr_stat), 0)
  n <- pmax(as.numeric(n_obs), 2)
  k <- pmax(as.numeric(df), 1)
  bic_delta <- lr - k * log(n)
  log_bf <- 0.5 * bic_delta
  log_bf <- pmin(log_bf, 700)
  bf <- exp(log_bf)
  bf[!is.finite(bf)] <- exp(700)
  pmax(bf, eps)
}

.sgcbLaplaceBF <- function(ll_full, ll_null,
                           logdet_full, logdet_null,
                           df_test,
                           coef_test = NULL,
                           prior_sd = 1.0,
                           eps = 1e-12) {
  ll1 <- as.numeric(ll_full)
  ll0 <- as.numeric(ll_null)
  ld1 <- as.numeric(logdet_full)
  ld0 <- as.numeric(logdet_null)
  k <- as.numeric(df_test)
  n <- length(ll1)

  if (length(ll0) == 1L) ll0 <- rep(ll0, n)
  if (length(ld1) == 1L) ld1 <- rep(ld1, n)
  if (length(ld0) == 1L) ld0 <- rep(ld0, n)
  if (length(k) == 1L) k <- rep(k, n)

  ll0[!is.finite(ll0)] <- ll1[!is.finite(ll0)]
  ld1[!is.finite(ld1)] <- 0
  ld0[!is.finite(ld0)] <- 0
  k[!is.finite(k)] <- 1
  k <- pmax(k, 1)

  # log BF from Laplace approximation:
  #   log BF10 ≈ (ll1-ll0) + 0.5*k*log(2π) - 0.5*(log|H1|-log|H0|) + log π(theta_hat)
  log_bf <- (ll1 - ll0) + 0.5 * k * log(2 * pi) - 0.5 * (ld1 - ld0)
  log_bf[!is.finite(log_bf)] <- log(eps)

  if (!is.null(coef_test)) {
    b <- as.matrix(coef_test)
    stopifnot(nrow(b) == n)
    b[!is.finite(b)] <- 0
    tau <- pmax(as.numeric(prior_sd), sqrt(eps))
    if (length(tau) == 1L) {
      tau <- rep(tau, ncol(b))
    }
    stopifnot(length(tau) == ncol(b))
    tau_mat <- matrix(rep(tau, each = nrow(b)), nrow = nrow(b), ncol = ncol(b))
    quad <- rowSums((b / tau_mat)^2)
    log_prior <- -0.5 * quad - sum(log(tau)) - 0.5 * ncol(b) * log(2 * pi)
    log_bf <- log_bf + log_prior
  }

  log_bf <- pmin(log_bf, 700)
  bf <- exp(log_bf)
  bf[!is.finite(bf)] <- exp(700)
  pmax(bf, eps)
}

.sgcbPseudoBF10 <- function(p, eps = 1e-12) {
  p_safe <- .sgcbClipProb(p, eps)
  thr <- exp(-1)
  bf01_lb <- ifelse(p_safe < thr, -exp(1) * p_safe * log(p_safe), 1)
  1 / pmax(bf01_lb, eps)
}

.sgcbModelPosterior1 <- function(p_de, prior_de = 0.1, bf_de = NULL, eps = 1e-12) {
  pi_de <- pmin(pmax(as.numeric(prior_de), eps), 1 - eps)
  bf_use <- if (is.null(bf_de)) .sgcbPseudoBF10(p_de, eps) else pmax(as.numeric(bf_de), eps)

  lbf_de <- pmin(log(pmax(bf_use, eps)), 700)
  l_null <- log1p(-pi_de)
  l_de <- log(pi_de) + lbf_de
  l_max <- pmax(l_null, l_de)
  w_null <- exp(l_null - l_max)
  w_de <- exp(l_de - l_max)
  den <- pmax(w_null + w_de, eps)

  model_prob <- cbind(
    model_prob_null = w_null / den,
    model_prob_de = w_de / den
  )
  model_prob[!is.finite(model_prob)] <- 0
  model_prob[model_prob < 0] <- 0
  model_prob <- model_prob / pmax(rowSums(model_prob), eps)

  list(
    p_de_post = pmin(pmax(model_prob[, "model_prob_de"], 0), 1),
    p_de_only_post = pmin(pmax(model_prob[, "model_prob_de"], 0), 1),
    bf_de = bf_use,
    model_prob = model_prob
  )
}

.sgcbModelPosterior2 <- function(p_de, p_dv,
                                 prior_de = 0.1, prior_dv = 0.1,
                                 bf_de = NULL, bf_dv = NULL,
                                 eps = 1e-12) {
  pi_de <- pmin(pmax(as.numeric(prior_de), eps), 1 - eps)
  pi_dv <- pmin(pmax(as.numeric(prior_dv), eps), 1 - eps)
  bf_de_use <- if (is.null(bf_de)) .sgcbPseudoBF10(p_de, eps) else pmax(as.numeric(bf_de), eps)
  bf_dv_use <- if (is.null(bf_dv)) .sgcbPseudoBF10(p_dv, eps) else pmax(as.numeric(bf_dv), eps)

  lbf_de <- pmin(log(pmax(bf_de_use, eps)), 700)
  lbf_dv <- pmin(log(pmax(bf_dv_use, eps)), 700)
  l_null <- log1p(-pi_de) + log1p(-pi_dv)
  l_de <- log(pi_de) + lbf_de + log1p(-pi_dv)
  l_dv <- log1p(-pi_de) + log(pi_dv) + lbf_dv
  l_de_dv <- log(pi_de) + lbf_de + log(pi_dv) + lbf_dv
  l_max <- pmax(pmax(l_null, l_de), pmax(l_dv, l_de_dv))
  w_null <- exp(l_null - l_max)
  w_de <- exp(l_de - l_max)
  w_dv <- exp(l_dv - l_max)
  w_de_dv <- exp(l_de_dv - l_max)
  den <- pmax(w_null + w_de + w_dv + w_de_dv, eps)

  model_prob <- cbind(
    model_prob_null = w_null / den,
    model_prob_de = w_de / den,
    model_prob_dv = w_dv / den,
    model_prob_de_dv = w_de_dv / den
  )
  model_prob[!is.finite(model_prob)] <- 0
  model_prob[model_prob < 0] <- 0
  model_prob <- model_prob / pmax(rowSums(model_prob), eps)

  list(
    p_de_post = pmin(pmax(model_prob[, "model_prob_de"] + model_prob[, "model_prob_de_dv"], 0), 1),
    p_dv_post = pmin(pmax(model_prob[, "model_prob_dv"] + model_prob[, "model_prob_de_dv"], 0), 1),
    p_de_only_post = pmin(pmax(model_prob[, "model_prob_de"], 0), 1),
    p_dv_only_post = pmin(pmax(model_prob[, "model_prob_dv"], 0), 1),
    p_de_dv_post = pmin(pmax(model_prob[, "model_prob_de_dv"], 0), 1),
    bf_de = bf_de_use,
    bf_dv = bf_dv_use,
    model_prob = model_prob
  )
}

.sgcbModelPosterior3 <- function(p_de, p_dv, p_dg,
                                 prior_de = 0.1, prior_dv = 0.1, prior_dg = 0.1,
                                 bf_de = NULL, bf_dv = NULL, bf_dg = NULL,
                                 eps = 1e-12) {
  pi_de <- pmin(pmax(as.numeric(prior_de), eps), 1 - eps)
  pi_dv <- pmin(pmax(as.numeric(prior_dv), eps), 1 - eps)
  pi_dg <- pmin(pmax(as.numeric(prior_dg), eps), 1 - eps)
  bf_de_use <- if (is.null(bf_de)) .sgcbPseudoBF10(p_de, eps) else pmax(as.numeric(bf_de), eps)
  bf_dv_use <- if (is.null(bf_dv)) .sgcbPseudoBF10(p_dv, eps) else pmax(as.numeric(bf_dv), eps)
  bf_dg_use <- if (is.null(bf_dg)) .sgcbPseudoBF10(p_dg, eps) else pmax(as.numeric(bf_dg), eps)

  lbf_de <- pmin(log(pmax(bf_de_use, eps)), 700)
  lbf_dv <- pmin(log(pmax(bf_dv_use, eps)), 700)
  lbf_dg <- pmin(log(pmax(bf_dg_use, eps)), 700)

  l_null <- log1p(-pi_de) + log1p(-pi_dv) + log1p(-pi_dg)
  l_de <- log(pi_de) + lbf_de + log1p(-pi_dv) + log1p(-pi_dg)
  l_dv <- log1p(-pi_de) + log(pi_dv) + lbf_dv + log1p(-pi_dg)
  l_dg <- log1p(-pi_de) + log1p(-pi_dv) + log(pi_dg) + lbf_dg
  l_de_dv <- log(pi_de) + lbf_de + log(pi_dv) + lbf_dv + log1p(-pi_dg)
  l_de_dg <- log(pi_de) + lbf_de + log1p(-pi_dv) + log(pi_dg) + lbf_dg
  l_dv_dg <- log1p(-pi_de) + log(pi_dv) + lbf_dv + log(pi_dg) + lbf_dg
  l_de_dv_dg <- log(pi_de) + lbf_de + log(pi_dv) + lbf_dv + log(pi_dg) + lbf_dg

  l_max <- pmax(
    pmax(pmax(l_null, l_de), pmax(l_dv, l_dg)),
    pmax(pmax(l_de_dv, l_de_dg), pmax(l_dv_dg, l_de_dv_dg))
  )
  w_null <- exp(l_null - l_max)
  w_de <- exp(l_de - l_max)
  w_dv <- exp(l_dv - l_max)
  w_dg <- exp(l_dg - l_max)
  w_de_dv <- exp(l_de_dv - l_max)
  w_de_dg <- exp(l_de_dg - l_max)
  w_dv_dg <- exp(l_dv_dg - l_max)
  w_de_dv_dg <- exp(l_de_dv_dg - l_max)
  den <- pmax(w_null + w_de + w_dv + w_dg + w_de_dv + w_de_dg + w_dv_dg + w_de_dv_dg, eps)

  model_prob <- cbind(
    model_prob_null = w_null / den,
    model_prob_de = w_de / den,
    model_prob_dv = w_dv / den,
    model_prob_dg = w_dg / den,
    model_prob_de_dv = w_de_dv / den,
    model_prob_de_dg = w_de_dg / den,
    model_prob_dv_dg = w_dv_dg / den,
    model_prob_de_dv_dg = w_de_dv_dg / den
  )
  model_prob[!is.finite(model_prob)] <- 0
  model_prob[model_prob < 0] <- 0
  model_prob <- model_prob / pmax(rowSums(model_prob), eps)

  list(
    p_de_post = pmin(pmax(model_prob[, "model_prob_de"] + model_prob[, "model_prob_de_dv"] +
      model_prob[, "model_prob_de_dg"] + model_prob[, "model_prob_de_dv_dg"], 0), 1),
    p_dv_post = pmin(pmax(model_prob[, "model_prob_dv"] + model_prob[, "model_prob_de_dv"] +
      model_prob[, "model_prob_dv_dg"] + model_prob[, "model_prob_de_dv_dg"], 0), 1),
    p_dg_post = pmin(pmax(model_prob[, "model_prob_dg"] + model_prob[, "model_prob_de_dg"] +
      model_prob[, "model_prob_dv_dg"] + model_prob[, "model_prob_de_dv_dg"], 0), 1),
    p_de_only_post = pmin(pmax(model_prob[, "model_prob_de"], 0), 1),
    p_dv_only_post = pmin(pmax(model_prob[, "model_prob_dv"], 0), 1),
    p_dg_only_post = pmin(pmax(model_prob[, "model_prob_dg"], 0), 1),
    p_de_dv_post = pmin(pmax(model_prob[, "model_prob_de_dv"], 0), 1),
    p_de_dg_post = pmin(pmax(model_prob[, "model_prob_de_dg"], 0), 1),
    p_dv_dg_post = pmin(pmax(model_prob[, "model_prob_dv_dg"], 0), 1),
    p_de_dv_dg_post = pmin(pmax(model_prob[, "model_prob_de_dv_dg"], 0), 1),
    bf_de = bf_de_use,
    bf_dv = bf_dv_use,
    bf_dg = bf_dg_use,
    model_prob = model_prob
  )
}
