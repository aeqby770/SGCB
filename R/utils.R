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
# Numerical gradient of log(CV²) w.r.t. (α, β, γ)
#
# CV² = Var/μ² = Γ(α+2/γ)Γ(α)/Γ(α+1/γ)² − 1
# Note: CV² is β-invariant (scale-free), so ∂log(CV²)/∂β ≈ 0 analytically.
# We compute all three partial derivatives numerically for safety.
# -----------------------------------------------------------------------------
.ggLogCV2Grad <- function(alpha, beta, gamma, eps = 1e-6) {
  step_a <- pmax(abs(alpha) * 1e-5, eps)
  step_b <- pmax(abs(beta) * 1e-5, eps)
  step_g <- pmax(abs(gamma) * 1e-5, eps)
  a_lo <- pmax(alpha - step_a, 0.01)
  a_hi <- alpha + step_a
  b_lo <- pmax(beta - step_b, eps)
  b_hi <- beta + step_b
  g_lo <- pmax(gamma - step_g, 0.01)
  g_hi <- gamma + step_g

  log_cv2 <- log(.ggVariance(alpha, beta, gamma, eps)) -
    2 * log(.ggMean(alpha, beta, gamma, eps))

  grad_a <- (log(.ggVariance(a_hi, beta, gamma, eps)) - 2 * log(.ggMean(a_hi, beta, gamma, eps)) -
             log(.ggVariance(a_lo, beta, gamma, eps)) + 2 * log(.ggMean(a_lo, beta, gamma, eps))) / (a_hi - a_lo)

  grad_b <- (log(.ggVariance(alpha, b_hi, gamma, eps)) - 2 * log(.ggMean(alpha, b_hi, gamma, eps)) -
             log(.ggVariance(alpha, b_lo, gamma, eps)) + 2 * log(.ggMean(alpha, b_lo, gamma, eps))) / (b_hi - b_lo)

  grad_g <- (log(.ggVariance(alpha, beta, g_hi, eps)) - 2 * log(.ggMean(alpha, beta, g_hi, eps)) -
             log(.ggVariance(alpha, beta, g_lo, eps)) + 2 * log(.ggMean(alpha, beta, g_lo, eps))) / (g_hi - g_lo)

  list(log_cv2 = log_cv2, grad_alpha = grad_a, grad_beta = grad_b, grad_gamma = grad_g)
}

# -----------------------------------------------------------------------------
# [L1] Delta-method variance of log(CV²) using FULL 3×3 Fisher inverse
#
# The correct formula on the GG statistical manifold:
#   Var[log(CV²)] = ∇ᵀ · I(θ)⁻¹ · ∇
#
# where ∇ = (∂log(CV²)/∂α, ∂log(CV²)/∂β, ∂log(CV²)/∂γ)
# and I⁻¹ is the COMPLETE inverse (not diagonal approximation).
#
# The Schur complement shows why diagonal fails:
#   (I⁻¹)_αα = (I_αα − I_αβ²/I_ββ − ...)⁻¹ ≥ 1/I_αα
# with equality only when all off-diagonals vanish.
# For GG, I_αβ = nγ/β ≠ 0, so the diagonal approx underestimates variance.
#
# [L2] Bartlett finite-sample correction:
#   Var_finite ≈ Var_Fisher · n / (n − p)
# where p = number of free parameters (3 for full GG, 2 when γ=1).
# This is the first-order expansion of the exact MLE variance:
#   Var(θ̂_MLE) = I⁻¹(1 + O(1/n)) ≈ I⁻¹ · n/(n−p)
# Reference: Bartlett (1937 JRSS), Cox & Snell (1968 JRSS)
# -----------------------------------------------------------------------------
.ggLogCV2SE_v2 <- function(alpha, beta, gamma, n_samples, eps = 1e-8) {
  # [L1] Conditional profile variance: Var(log CV²|γ̂)
  #
  # Mathematical justification for conditioning on γ:
  #   γ is estimated from POOLED data (both groups, via hierarchical prior)
  #   with effective n ≈ n_ctrl + n_treat, giving O(1/2n) uncertainty.
  #   Per-group α estimates have O(1/n) uncertainty. So γ̂ is √2 more
  #   precise than α̂, and can be treated as approximately known.
  #   This is the "conditional profile likelihood" approach.
  #
  # With γ known, CV² = f(α,γ₀) is a function of α alone.
  # Profiling out β (which CV² doesn't depend on):
  #   Var(α̂|γ) = 1 / (F_αα − F_αβ²/F_ββ)
  #   = 1 / [n(ψ'(α) − 1/α)]   (Stacy-Lawless formula)
  #
  # This is larger than 1/F_αα (the diagonal-only bug) because it
  # accounts for the α-β coupling, but not inflated by α-γ aliasing.

  grad <- .ggLogCV2Grad(alpha, beta, gamma, eps)
  ga <- grad$grad_alpha  # ∂log(CV²)/∂α
  # ∂log(CV²)/∂β ≡ 0 (scale-free), ∂log(CV²)/∂γ absorbed into γ̂ conditioning

  # Fisher info elements (vectorized computation, no C++ call needed)
  ns <- as.double(n_samples)
  psi1_a <- trigamma(pmax(alpha, 0.05))  # ψ'(α)

  # Schur complement of F_ββ in (α,β)|γ submatrix:
  # S_α = F_αα − F_αβ²/F_ββ = n·ψ'(α) − n·(γ/β)² / (n·α·γ²/β²) = n(ψ'(α) − 1/α)
  s_alpha <- ns * (psi1_a - 1 / pmax(alpha, 0.05))
  s_alpha <- pmax(s_alpha, eps)

  # Var(log CV²|γ) = (∂log CV²/∂α)² / S_α
  var_logcv2 <- ga^2 / s_alpha

  # [L2] Bartlett correction: n/(n-1) for the 1-parameter α submodel
  bartlett_factor <- ns / pmax(ns - 1, 1)
  var_logcv2 <- var_logcv2 * bartlett_factor

  list(log_cv2 = grad$log_cv2, var_logcv2 = pmax(var_logcv2, eps))
}

# -----------------------------------------------------------------------------
# Brown-Forsythe DV test: nonparametric safety floor
# Fully vectorized matrix operations, no loops/apply.
# Tests H₀: equal dispersion via absolute deviations from group means.
# Returns F(1, n_total-2) p-values, or rep(1, n_genes) when data unavailable.
# -----------------------------------------------------------------------------
.brownForsytheDV <- function(norm_counts, ctrl_idx, treat_idx,
                             n_ctrl, n_treat, eps = 1e-8) {
  n_genes <- nrow(norm_counts)
  n_total <- n_ctrl + n_treat

  log_counts <- log(norm_counts)
  ctrl_mat <- log_counts[, ctrl_idx, drop = FALSE]
  treat_mat <- log_counts[, treat_idx, drop = FALSE]

  # Group centers (row means — fast, vectorized)
  cen_ctrl <- rowMeans(ctrl_mat)
  cen_treat <- rowMeans(treat_mat)

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
  cen_ctrl <- rowMeans(ctrl_mat)
  cen_treat <- rowMeans(treat_mat)
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
    sc <- rowMeans(s_ctrl)
    st <- rowMeans(s_treat)
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
# DV test v3: MAD-standardized differential CV² test
#
# The GG-model log(CV²) difference is standardized by the genome-wide
# MAD (median absolute deviation), not per-gene SE. This guarantees
# calibration under the empirical null without relying on Fisher
# information SE estimates, which can be pathologically small.
#
# Mathematical background:
#   diff_g = log(CV²_treat) - log(CV²_ctrl) for each gene g
#   Under H₀ (most genes have no DV effect):
#     diff ~ F₀(0, σ²_null) for null genes
#   σ_null = MAD(diff) / 0.6745  (robust estimator of null SD)
#   z_g = diff_g / σ_null
#
# Properties:
#   - Under H₀: z ~ N(0,1) by construction (MAD calibration)
#   - Under H₁: genes with large |diff| → large |z| → small p-value
#   - No per-gene SE needed → no pathological SE outliers
#   - Robust to model misspecification in the SE formula
#
# Trade-off: genes with different true SE are treated equally.
#   Precise genes (large α) lose some power because they share the
#   common σ_null with imprecise genes. But FP control is guaranteed.
#
# The per-gene model-based SE from .ggLogCV2SE_v2() is still computed
# and used as a WEIGHT to combine the MAD null with gene-specific info.
# Specifically: SE_used = max(SE_model, σ_null), ensuring no gene
# has a smaller SE than the genome-wide null SD. This preserves power
# for genes with legitimately large SE (conservative direction) while
# capping precision at the genome-wide level.
# -----------------------------------------------------------------------------
.sgcbDVTest <- function(alpha_ctrl, beta_ctrl, gamma_ctrl,
                        alpha_treat, beta_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8,
                        norm_counts = NULL, ctrl_idx = NULL, treat_idx = NULL) {

  # Compute log(CV²) from GG parameters
  ctrl <- .ggLogCV2SE_v2(alpha_ctrl, beta_ctrl, gamma_ctrl, n_ctrl, eps)
  treat <- .ggLogCV2SE_v2(alpha_treat, beta_treat, gamma_treat, n_treat, eps)
  diff_log <- treat$log_cv2 - ctrl$log_cv2

  # MAD-based null SD: robust genome-wide estimate
  sigma_null <- stats::median(abs(diff_log), na.rm = TRUE) / 0.6745
  sigma_null <- max(sigma_null, eps)

  # Per-gene model SE (for the conservative floor)
  se_model <- sqrt(ctrl$var_logcv2 + treat$var_logcv2)

  # Convex blend of model SE and global null SD (avoid hard max floor)
  n_min <- min(n_ctrl, n_treat)
  w_null_base <- ifelse(n_min <= 4, 0.55,
                        ifelse(n_min <= 6, 0.45,
                               ifelse(n_min <= 10, 0.35,
                                      ifelse(n_min <= 20, 0.25, 0.20))))
  se_ratio <- se_model / sigma_null
  w_null <- w_null_base + 0.15 * pmax(0, 1 - se_ratio)
  w_null <- pmin(pmax(w_null, 0.10), 0.80)
  se_used <- sqrt((1 - w_null) * se_model^2 + w_null * sigma_null^2)

  # =========================================================================
  # Brown-Forsythe concordance weighting
  #
  # BF test on log-counts gives a nonparametric DV signal for each gene.
  # When BF sees no signal (large p_bf), the model-based z is unreliable
  # → inflate the SE to reduce the z-statistic.
  # When BF confirms the signal (small p_bf), keep SE unchanged.
  #
  # Inflation factor: max(1, sqrt(p_bf / 0.05))
  #   p_bf < 0.05: factor = 1 (BF agrees, no penalty)
  #   p_bf = 0.5: factor = √10 ≈ 3.16 (BF skeptical, strong penalty)
  #   p_bf = 0.9: factor = √18 ≈ 4.24 (BF disagrees, very strong penalty)
  #
  # This is an adaptive shrinkage of the z-statistic: genes with concordant
  # parametric and nonparametric evidence retain power, while genes with
  # only parametric evidence (potentially from unstable GG estimates)
  # are penalized proportionally.
  # =========================================================================
  p_bf <- .brownForsytheDV(norm_counts, ctrl_idx, treat_idx, n_ctrl, n_treat, eps)
  bf_inflation <- ifelse(
    p_bf <= 0.05,
    1.00,
    ifelse(
      p_bf <= 0.20,
      1.00 + 0.35 * (p_bf - 0.05) / 0.15,
      ifelse(
        p_bf <= 0.50,
        1.35 + 0.35 * (p_bf - 0.20) / 0.30,
        1.70 + 0.30 * (p_bf - 0.50) / 0.50
      )
    )
  )
  bf_inflation <- pmin(pmax(bf_inflation, 1), 2.0)
  se_used <- se_used * bf_inflation

  z_stat <- diff_log / se_used
  p_model_raw <- stats::pchisq(z_stat^2, df = 1, lower.tail = FALSE)
  p_model_raw <- pmin(pmax(p_model_raw, eps), 1)
  if (n_min <= 6) {
    p_bf <- .dvBootstrapP(
      norm_counts = norm_counts,
      ctrl_idx = ctrl_idx,
      treat_idx = treat_idx,
      alpha_ctrl = alpha_ctrl,
      beta_ctrl = beta_ctrl,
      gamma_ctrl = gamma_ctrl,
      alpha_treat = alpha_treat,
      beta_treat = beta_treat,
      gamma_treat = gamma_treat,
      n_boot = 40L,
      eps = eps
    )
  }
  p_bf <- pmin(pmax(p_bf, eps), 1)

  # Concordance-guarded model channel:
  # require some support from nonparametric BF evidence to prevent
  # anti-conservative parametric tails in finite samples.
  guard_exp <- ifelse(n_min <= 4, 2.8,
                      ifelse(n_min <= 6, 2.6,
                             ifelse(n_min <= 10, 2.2, 2.1)))
  p_model <- pmax(p_model_raw, p_bf^guard_exp)
  p_model <- pmin(pmax(p_model, eps), 1)

  p_mix_mat <- cbind(p_model, p_bf)
  p_mix_mat[p_mix_mat < eps] <- eps
  p_mix_mat[p_mix_mat > 1 - eps] <- 1 - eps
  w_model <- ifelse(n_min <= 4, 1.0,
                    ifelse(n_min <= 6, 1.0,
                           ifelse(n_min <= 10, 1.0, 1.0)))
  w_bf <- ifelse(n_min <= 4, 2.0,
                 ifelse(n_min <= 6, 2.0,
                        ifelse(n_min <= 10, 1.8, 1.5)))
  t_mix <- w_model * tan(pi * (0.5 - p_mix_mat[, 1])) +
    w_bf * tan(pi * (0.5 - p_mix_mat[, 2]))
  p_mix <- 0.5 - atan(t_mix / (w_model + w_bf)) / pi
  p_val <- pmax(p_mix, p_bf)
  p_val <- pmin(pmax(p_val, eps), 1)

  var_ctrl <- .ggVariance(alpha_ctrl, beta_ctrl, gamma_ctrl, eps)
  var_treat <- .ggVariance(alpha_treat, beta_treat, gamma_treat, eps)
  mean_ctrl <- .ggMean(alpha_ctrl, beta_ctrl, gamma_ctrl, eps)
  mean_treat <- .ggMean(alpha_treat, beta_treat, gamma_treat, eps)
  cv2_ctrl <- var_ctrl / pmax(mean_ctrl^2, eps)
  cv2_treat <- var_treat / pmax(mean_treat^2, eps)

  list(
    log2_var_ratio = (log(var_treat) - log(var_ctrl)) / log(2),
    log2_cv2_ratio = diff_log / log(2),
    se_log2_cv2 = se_used / log(2),
    stat = z_stat,
    pvalue = p_val,
    p_model_raw = p_model_raw,
    p_model = p_model,
    p_bf = p_bf,
    p_mix = p_mix,
    var_ctrl = var_ctrl,
    var_treat = var_treat,
    cv2_ctrl = cv2_ctrl,
    cv2_treat = cv2_treat
  )
}

# -----------------------------------------------------------------------------
# DG test v2: Differential Dynamics (α, γ shape parameters)
# with conditional profile variance + Bartlett + EB shrinkage + moderated t
#
# Tests H₀: log(α_treat) = log(α_ctrl) and H₀: log(γ_treat) = log(γ_ctrl)
# on the Lie algebra of (ℝ⁺)³ (log-parameterization).
#
# For log(α): Var[log(α̂)|β,γ] = 1/(α² · S_α) where S_α = n(ψ'(α)−1/α)
# For log(γ): Var[log(γ̂)|α,β] = 1/(γ² · F_γγ) as a conditional diagonal
#   approximation. Full GG Fisher has nonzero α-γ and β-γ coupling, so this
#   channel is diagnostic and empirically calibrated rather than confirmatory.
# -----------------------------------------------------------------------------
.sgcbDGTest <- function(alpha_ctrl, gamma_ctrl,
                        alpha_treat, gamma_treat,
                        n_ctrl, n_treat, eps = 1e-8) {

  # Conditional profile variance for log(α̂): Var = 1/(α²·n(ψ'(α)−1/α))
  psi1_ctrl <- trigamma(pmax(alpha_ctrl, 0.05))
  psi1_treat <- trigamma(pmax(alpha_treat, 0.05))
  s_alpha_ctrl <- as.double(n_ctrl) * (psi1_ctrl - 1 / pmax(alpha_ctrl, 0.05))
  s_alpha_treat <- as.double(n_treat) * (psi1_treat - 1 / pmax(alpha_treat, 0.05))
  var_loga_ctrl <- 1 / (pmax(alpha_ctrl^2, eps) * pmax(s_alpha_ctrl, eps))
  var_loga_treat <- 1 / (pmax(alpha_treat^2, eps) * pmax(s_alpha_treat, eps))

  # Bartlett correction
  bart <- function(n) n / pmax(n - 1, 1)
  var_loga_ctrl <- var_loga_ctrl * bart(n_ctrl)
  var_loga_treat <- var_loga_treat * bart(n_treat)

  # For log(γ̂): use a conditional diagonal Fisher approximation after fixing α,β
  psi_ctrl <- digamma(pmax(alpha_ctrl, 0.05))
  psi_treat <- digamma(pmax(alpha_treat, 0.05))
  fgg_ctrl <- as.double(n_ctrl) * (1 + alpha_ctrl * psi_ctrl^2 + 2 * psi_ctrl + alpha_ctrl * psi1_ctrl) / pmax(gamma_ctrl^2, eps)
  fgg_treat <- as.double(n_treat) * (1 + alpha_treat * psi_treat^2 + 2 * psi_treat + alpha_treat * psi1_treat) / pmax(gamma_treat^2, eps)
  var_logg_ctrl <- 1 / (pmax(gamma_ctrl^2, eps) * pmax(fgg_ctrl, eps))
  var_logg_treat <- 1 / (pmax(gamma_treat^2, eps) * pmax(fgg_treat, eps))
  var_logg_ctrl <- var_logg_ctrl * bart(n_ctrl)
  var_logg_treat <- var_logg_treat * bart(n_treat)

  log_alpha_ratio <- log(pmax(alpha_treat, eps) / pmax(alpha_ctrl, eps))
  log_gamma_ratio <- log(pmax(gamma_treat, eps) / pmax(gamma_ctrl, eps))

  se_alpha <- sqrt(pmax(var_loga_ctrl + var_loga_treat, eps))
  se_gamma <- sqrt(pmax(var_logg_ctrl + var_logg_treat, eps))

  # -----------------------------------------------------------------------
  # Efron empirical-null calibration (fixes boundary-collapse bug)
  #
  # Problem: when many genes have α or γ at boundary, MAD(log_ratio) = 0,
  # destroying the MAD-based SE floor. Model SE alone underestimates noise
  # because Fisher info assumes interior MLE.
  #
  # Fix: compute raw z = log_ratio / model_se, then estimate the empirical
  # null SD from non-boundary genes (where log_ratio ≠ 0). Divide raw z
  # by this empirical SD to get calibrated z ~ N(0,1) under null.
  # -----------------------------------------------------------------------
  z_alpha_raw <- log_alpha_ratio / se_alpha
  z_gamma_raw <- log_gamma_ratio / se_gamma

  # Boundary detection: genes where BOTH groups hit parameter constraints
  at_bnd_alpha <- abs(log_alpha_ratio) < eps
  at_bnd_gamma <- abs(log_gamma_ratio) < eps

  # Empirical null SD from non-boundary genes (Efron 2004 style)
  # Uses Winsorized IQR to be robust against both heavy tails and
  # the boundary-collapse point mass.
  .emp_null_sd <- function(z_raw, at_boundary) {
    free <- !at_boundary & is.finite(z_raw)
    if (sum(free) < 50) return(1.0)
    z_free <- z_raw[free]
    # IQR-based sigma (robust to both tails and boundary mass)
    iqr_sigma <- stats::IQR(z_free, na.rm = TRUE) / 1.349
    # Winsorized SD at 5th/95th percentile (handles heavy tails)
    qq <- stats::quantile(z_free, c(0.05, 0.95), na.rm = TRUE)
    z_wins <- pmin(pmax(z_free, qq[1]), qq[2])
    wins_sd <- stats::sd(z_wins, na.rm = TRUE)
    max(iqr_sigma, wins_sd, 1.0)
  }
  sigma_z_alpha <- .emp_null_sd(z_alpha_raw, at_bnd_alpha)
  sigma_z_gamma <- .emp_null_sd(z_gamma_raw, at_bnd_gamma)

  z_alpha <- z_alpha_raw / sigma_z_alpha
  z_gamma <- z_gamma_raw / sigma_z_gamma

  # Boundary genes: force z = 0 (no evidence of difference)
  z_alpha[at_bnd_alpha] <- 0
  z_gamma[at_bnd_gamma] <- 0

  # entropy channel (objective-Bayes inspired robust shape functional)
  # H(GG) = alpha + log(beta) + lgamma(alpha) - (alpha - 1/gamma)psi(alpha) - log(gamma)
  # use beta=1 to isolate pure shape contribution from (alpha, gamma)
  ent_ctrl_shape <- .ggEntropy(alpha_ctrl, beta = rep(1, length(alpha_ctrl)), gamma_ctrl, eps)
  ent_treat_shape <- .ggEntropy(alpha_treat, beta = rep(1, length(alpha_treat)), gamma_treat, eps)
  ent_diff <- ent_treat_shape - ent_ctrl_shape
  at_bnd_ent <- abs(ent_diff) < eps
  # Entropy has no analytical SE → use IQR-based sigma of non-boundary values
  # then recalibrate using trimmed SD of central z-scores (|z|<3)
  free_ent <- !at_bnd_ent & is.finite(ent_diff)
  iqr_e <- if (sum(free_ent) >= 50) stats::IQR(ent_diff[free_ent], na.rm = TRUE) / 1.349 else 1.0
  iqr_e <- max(iqr_e, eps)
  z_entropy_raw <- ent_diff / iqr_e
  # Trimmed SD recalibration: use central genes (|z|<3) to estimate true null SD
  central <- free_ent & abs(z_entropy_raw) < 3
  sigma_z_ent <- if (sum(central) >= 50) max(stats::sd(z_entropy_raw[central], na.rm = TRUE), 1.0) else 1.0
  z_entropy <- z_entropy_raw / sigma_z_ent
  z_entropy[at_bnd_ent] <- 0

  p_alpha <- 2 * stats::pnorm(-abs(z_alpha))
  p_gamma <- 2 * stats::pnorm(-abs(z_gamma))
  p_entropy <- 2 * stats::pnorm(-abs(z_entropy))

  # shape omnibus: pure (alpha, gamma) chi-sq(2) test
  # Entropy excluded from omnibus — no analytical SE, unreliable at all n.
  shape_stat <- z_alpha^2 + z_gamma^2
  p_shape_core <- stats::pchisq(shape_stat, df = 2, lower.tail = FALSE)
  p_shape_core <- pmin(pmax(p_shape_core, eps), 1)
  # Cauchy combination of alpha and gamma only
  p_cauchy_raw <- 0.5 - atan(
    (tan(pi * (0.5 - p_alpha)) + tan(pi * (0.5 - p_gamma))) / 2
  ) / pi
  p_cauchy <- pmin(pmax(p_cauchy_raw, eps), 1)
  # Bonferroni-Sidak on min(p_alpha, p_gamma)
  p_min <- pmin(p_alpha, p_gamma)
  p_shape_max <- 1 - (1 - pmin(pmax(p_min, eps), 1 - eps))^2
  p_shape_max <- pmin(pmax(p_shape_max, eps), 1)
  p_shape <- pmax(p_shape_core, p_cauchy, p_shape_max)
  p_shape <- pmin(pmax(p_shape, eps), 1)

  list(
    alpha_log2_ratio = log_alpha_ratio / log(2),
    alpha_se_log2 = se_alpha / log(2),
    alpha_stat = z_alpha,
    alpha_pvalue = p_alpha,
    gamma_log2_ratio = log_gamma_ratio / log(2),
    gamma_se_log2 = se_gamma / log(2),
    gamma_stat = z_gamma,
    gamma_pvalue = p_gamma,
    entropy_diff = ent_diff,
    entropy_stat = z_entropy,
    entropy_pvalue = p_entropy,
    shape_stat = shape_stat,
    shape_pvalue = p_shape
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
