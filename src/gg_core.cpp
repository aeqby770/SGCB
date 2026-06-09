// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "sgcb_optimized.h"
#include "fast_special.h"
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace sgcb;

// =============================================================================
// GG log-likelihood (vectorized, OpenMP parallel)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericVector gg_loglik_vec(const Rcpp::NumericVector& x,
                                   const Rcpp::NumericVector& alpha,
                                   const Rcpp::NumericVector& beta,
                                   const Rcpp::NumericVector& gamma,
                                   double eps = 1e-8) {
  const int n = x.size();
  Rcpp::NumericVector out(n);
  
  const double* SGCB_RESTRICT px = &x[0];
  const double* SGCB_RESTRICT pa = &alpha[0];
  const double* SGCB_RESTRICT pb = &beta[0];
  const double* SGCB_RESTRICT pg = &gamma[0];
  double* SGCB_RESTRICT po = &out[0];
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    const double xi = px[i] > eps ? px[i] : eps;
    const double a = pa[i];
    const double b = pb[i];
    const double g = pg[i];
    const double log_xi = std::log(xi);
    const double log_b = std::log(b);
    po[i] = std::log(g) - g * a * log_b - fast_special::lgamma(a) +
            (g * a - 1.0) * log_xi - std::pow(xi / b, g);
  }
  return out;
}

// =============================================================================
// GG parameter gradients (w.r.t. alpha, beta, gamma)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix gg_grad_params(const Rcpp::NumericVector& x,
                                    const Rcpp::NumericVector& alpha,
                                    const Rcpp::NumericVector& beta,
                                    const Rcpp::NumericVector& gamma,
                                    double eps = 1e-8) {
  const int n = x.size();
  Rcpp::NumericMatrix grad(n, 3);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    const double xi = x[i] > eps ? x[i] : eps;
    const double a = alpha[i];
    const double b = beta[i];
    const double g = gamma[i];
    const double log_xi = std::log(xi);
    const double ratio = xi / b;
    const double ratio_g = std::pow(ratio, g);
    const double log_ratio = std::log(ratio);
    
    grad(i, 0) = -g * std::log(b) - fast_special::digamma(a) + g * log_xi;
    grad(i, 1) = -g * a / b + g * ratio_g / b;
    grad(i, 2) = 1.0 / g + a * log_ratio - ratio_g * log_ratio;
  }
  return grad;
}

// =============================================================================
// GG Hessian diagonal elements (for second-order information / confidence intervals)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix gg_hessian_diag(Rcpp::NumericVector x,
                                     Rcpp::NumericVector alpha,
                                     Rcpp::NumericVector beta,
                                     Rcpp::NumericVector gamma,
                                     double eps = 1e-8) {
  int n = x.size();
  Rcpp::NumericMatrix hess(n, 3);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double xi = x[i] + eps;
    double a = alpha[i];
    double b = beta[i];
    double g = gamma[i];
    double log_xi = std::log(xi);
    double ratio = xi / b;
    double ratio_g = std::pow(ratio, g);
    double log_ratio = std::log(ratio);
    double trigamma_a = fast_special::trigamma(a);
    hess(i, 0) = -trigamma_a;
    hess(i, 1) = g * a / (b * b) - g * (g + 1.0) * ratio_g / (b * b);
    hess(i, 2) = -1.0 / (g * g) - ratio_g * log_ratio * log_ratio;
  }
  return hess;
}

// =============================================================================
// GG sampling matrix (OpenMP parallel)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix gg_sample_mat(int n_genes, int n_samples,
                                   Rcpp::NumericVector alpha,
                                   Rcpp::NumericVector beta,
                                   Rcpp::NumericVector gamma,
                                   double eps = 1e-8) {
  Rcpp::NumericMatrix out(n_genes, n_samples);
  for (int g = 0; g < n_genes; g++) {
    double a = alpha[g];
    double b = beta[g];
    double gam = gamma[g];
    for (int s = 0; s < n_samples; s++) {
      double u = R::rgamma(a, 1.0);
      double val = b * std::pow(u + eps, 1.0 / gam);
      out(g, s) = val;
    }
  }
  return out;
}

// =============================================================================
// Compute Fisher information matrix diagonal (for Wald confidence intervals)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix fisher_info_diag(Rcpp::NumericVector alpha,
                                      Rcpp::NumericVector beta,
                                      Rcpp::NumericVector gamma,
                                      int n_samples) {
  int n = alpha.size();
  Rcpp::NumericMatrix info(n, 3);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double a = alpha[i];
    double b = beta[i];
    double g = gamma[i];
    double trigamma_a = R::trigamma(a);
    double digamma_a = R::digamma(a);
    info(i, 0) = n_samples * trigamma_a;
    info(i, 1) = n_samples * a * g * g / (b * b);
    info(i, 2) = n_samples * (1.0 + a * digamma_a * digamma_a + 2.0 * digamma_a + a * trigamma_a) / (g * g);
  }
  return info;
}

// =============================================================================
// Schur complement of I_ββ: conditional (α,γ) precision given β
//
// For functionals that depend only on (α,γ) — such as CV², which is
// scale-free — the correct variance is the CONDITIONAL variance of
// (α̂,γ̂) given β̂, not the marginal variance from the full inverse.
//
// The conditional precision of (α,γ)|β is the Schur complement:
//   S = [I_αα  I_αγ] - (1/I_ββ) [I_αβ] [I_αβ  I_βγ]
//       [I_αγ  I_γγ]             [I_βγ]
//
// S is 2×2 SPD. Its inverse S⁻¹ gives Cov(α̂,γ̂|β̂).
//
// Why this is correct for CV²:
//   CV² = Γ(α+2/γ)Γ(α)/Γ(α+1/γ)² − 1  (no β dependency)
//   The profile likelihood for CV² integrates out β, yielding the
//   conditional Fisher information S as the relevant metric.
//
// Why full I⁻¹ is wrong for CV²:
//   (I⁻¹)_{αα,αγ,γγ} includes the marginal uncertainty from β estimation,
//   which inflates the variance of β-independent functionals.
//   Specifically: Var_marginal(α̂) > Var_conditional(α̂|β̂) when I_αβ ≠ 0.
//
// Returns: n × 3 matrix with columns (S⁻¹_αα, S⁻¹_γγ, S⁻¹_αγ)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix fisher_info_schur_ag(Rcpp::NumericVector alpha,
                                          Rcpp::NumericVector beta,
                                          Rcpp::NumericVector gamma,
                                          int n_samples) {
  const int n = alpha.size();
  const double reg = 1e-10;
  Rcpp::NumericMatrix out(n, 3);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double a = std::max(alpha[i], 0.05);
    double b = std::max(beta[i], reg);
    double g = std::max(gamma[i], 0.05);
    double psi_a  = R::digamma(a);
    double psi1_a = R::trigamma(a);
    int    ns = n_samples;

    // Fisher information matrix elements
    double Faa = ns * psi1_a;
    double Fbb = ns * a * g * g / (b * b);
    double Fgg = ns * (1.0 + a * psi_a * psi_a + 2.0 * psi_a + a * psi1_a) / (g * g);
    double Fab = ns * g / b;
    double Fag = -ns * psi_a / g;
    double Fbg = -ns * (1.0 + a * psi_a) / b;

    // Schur complement: S = [Faa Fag; Fag Fgg] - (1/Fbb)*[Fab; Fbg]*[Fab Fbg]
    double inv_Fbb = 1.0 / std::max(Fbb, reg);
    double Saa = Faa - Fab * Fab * inv_Fbb;
    double Sgg = Fgg - Fbg * Fbg * inv_Fbb;
    double Sag = Fag - Fab * Fbg * inv_Fbb;

    // Invert 2×2 symmetric S → S⁻¹
    double det_S = Saa * Sgg - Sag * Sag;
    det_S = std::max(det_S, reg);
    double inv_det = 1.0 / det_S;

    out(i, 0) = Sgg * inv_det;    // S⁻¹_αα
    out(i, 1) = Saa * inv_det;    // S⁻¹_γγ
    out(i, 2) = -Sag * inv_det;   // S⁻¹_αγ
  }
  return out;
}

// =============================================================================
// Full 3×3 Fisher information inverse for GG distribution
//
// Mathematical background:
//   The GG(α,β,γ) Fisher information matrix is a 3×3 symmetric positive
//   definite matrix with significant off-diagonal coupling:
//     I_αβ = nγ/β,  I_αγ = -nψ(α)/γ,  I_βγ = -n(1+αψ(α))/β
//
//   The DV test computes Var[log(CV²)] = ∇ᵀ I⁻¹ ∇ via the delta method.
//   Using only the diagonal of I⁻¹ (i.e. 1/I_ii) is algebraically wrong:
//     (I⁻¹)_ii ≠ 1/I_ii  when I_ij ≠ 0
//
//   The correct marginal variance requires Cramer's rule on the full matrix:
//     I⁻¹ = adj(I) / det(I)
//
//   For a 3×3 SPD matrix, the Schur complement perspective shows that
//   ignoring the α-β correlation (I_αβ > 0) systematically underestimates
//   the marginal variances, producing inflated z-statistics → false positives.
//
// Returns: n × 6 matrix with columns
//   (inv_aa, inv_bb, inv_gg, inv_ab, inv_ag, inv_bg)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix fisher_info_full_inverse(Rcpp::NumericVector alpha,
                                              Rcpp::NumericVector beta,
                                              Rcpp::NumericVector gamma,
                                              int n_samples) {
  const int n = alpha.size();
  const double reg = 1e-10;
  Rcpp::NumericMatrix out(n, 6);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double a = std::max(alpha[i], 0.05);
    double b = std::max(beta[i], reg);
    double g = std::max(gamma[i], 0.05);
    double psi_a  = R::digamma(a);
    double psi1_a = R::trigamma(a);
    int    ns = n_samples;

    // Fisher information matrix elements (Stacy 1962, Lawless 1980)
    double Faa = ns * psi1_a;
    double Fbb = ns * a * g * g / (b * b);
    double Fgg = ns * (1.0 + a * psi_a * psi_a + 2.0 * psi_a + a * psi1_a) / (g * g);
    double Fab = ns * g / b;
    double Fag = -ns * psi_a / g;
    double Fbg = -ns * (1.0 + a * psi_a) / b;

    // Cramer's rule: cofactors of 3×3 symmetric matrix
    double cof_aa = Fbb * Fgg - Fbg * Fbg;             // minor(0,0)
    double cof_ab = -(Fab * Fgg - Fbg * Fag);           // -minor(0,1)
    double cof_ag = Fab * Fbg - Fbb * Fag;              // minor(0,2)
    double cof_bb = Faa * Fgg - Fag * Fag;              // minor(1,1)
    double cof_bg = -(Faa * Fbg - Fab * Fag);           // -minor(1,2)
    double cof_gg = Faa * Fbb - Fab * Fab;              // minor(2,2)

    double det = Faa * cof_aa + Fab * cof_ab + Fag * cof_ag;

    // Tikhonov regularization for near-singular Fisher matrix (small n)
    double ridge = 1e-6 * std::max(Faa + Fbb + Fgg, 1.0) + reg;
    det = std::max(det, ridge);

    double inv_det = 1.0 / det;

    out(i, 0) = cof_aa * inv_det;   // (I⁻¹)_αα
    out(i, 1) = cof_bb * inv_det;   // (I⁻¹)_ββ
    out(i, 2) = cof_gg * inv_det;   // (I⁻¹)_γγ
    out(i, 3) = cof_ab * inv_det;   // (I⁻¹)_αβ
    out(i, 4) = cof_ag * inv_det;   // (I⁻¹)_αγ
    out(i, 5) = cof_bg * inv_det;   // (I⁻¹)_βγ
  }
  return out;
}

