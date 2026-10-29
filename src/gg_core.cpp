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
// Efficient (profile) information after estimating nuisance scale β.
// S = I_(α,γ) - I_(α,γ),β I_ββ^-1 I_β,(α,γ).
// S^-1 is exactly the (α,γ) block of the full inverse Fisher matrix;
// it is not covariance conditional on a known β. Profiling maximizes
// over β; it does not integrate β out.
// Algebraically simplifying S avoids subtraction of β-dependent terms:
// S/n = [[trigamma(α)-1/α, 1/(αγ)],
//        [1/(αγ), (1+α trigamma(α)-1/α)/γ²]].
//
// Returns: n × 3 matrix with columns (S⁻¹_αα, S⁻¹_γγ, S⁻¹_αγ)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix fisher_info_schur_ag(Rcpp::NumericVector alpha,
                                          Rcpp::NumericVector beta,
                                          Rcpp::NumericVector gamma,
                                          int n_samples) {
  const int n = alpha.size();
  Rcpp::NumericMatrix out(n, 3);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double a = alpha[i];
    double g = gamma[i];
    double tri = R::trigamma(a);
    double u = tri - 1.0 / a;
    double v = 1.0 + a * tri - 1.0 / a;
    double denominator = n_samples * (u * v - 1.0 / (a * a));
    out(i, 0) = v / denominator;
    out(i, 1) = g * g * u / denominator;
    out(i, 2) = -g / (a * denominator);
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
//   A block inverse recovers all six covariance entries without a determinant
//   floor depending on the units of β. A determinant floor is not a ridge.
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
  Rcpp::NumericMatrix out(n, 6);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double a = alpha[i];
    double b = beta[i];
    double g = gamma[i];
    double psi = R::digamma(a);
    double tri = R::trigamma(a);
    double u = tri - 1.0 / a;
    double v = 1.0 + a * tri - 1.0 / a;
    double denominator = n_samples * (u * v - 1.0 / (a * a));
    double cov_aa = v / denominator;
    double cov_gg = g * g * u / denominator;
    double cov_ag = -g / (a * denominator);
    double scale_a = b / (a * g);
    double scale_g = -b * (1.0 + a * psi) / (a * g * g);
    double cov_ab = -(cov_aa * scale_a + cov_ag * scale_g);
    double cov_bg = -(cov_ag * scale_a + cov_gg * scale_g);
    double cov_bb = b * b / (n_samples * a * g * g) -
      scale_a * cov_ab - scale_g * cov_bg;
    out(i, 0) = cov_aa;
    out(i, 1) = cov_bb;
    out(i, 2) = cov_gg;
    out(i, 3) = cov_ab;
    out(i, 4) = cov_ag;
    out(i, 5) = cov_bg;
  }
  return out;
}

