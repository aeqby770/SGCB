// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "sgcb_optimized.hpp"
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
// Parameter clipping
// =============================================================================
// [[Rcpp::export]]
double clip_value(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

// [[Rcpp::export]]
Rcpp::NumericVector clip_vec(Rcpp::NumericVector x, double lo, double hi) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double val = x[i];
    out[i] = (val < lo) ? lo : ((val > hi) ? hi : val);
  }
  return out;
}

// =============================================================================
// Adam optimizer state update (vectorized)
// =============================================================================
// [[Rcpp::export]]
Rcpp::List adam_update(Rcpp::NumericVector param,
                        Rcpp::NumericVector grad,
                        Rcpp::NumericVector m,
                        Rcpp::NumericVector v,
                        double lr,
                        double beta1,
                        double beta2,
                        double eps,
                        int t) {
  int n = param.size();
  Rcpp::NumericVector m_new(n);
  Rcpp::NumericVector v_new(n);
  Rcpp::NumericVector param_new(n);
  double bc1 = 1.0 - std::pow(beta1, t);
  double bc2 = 1.0 - std::pow(beta2, t);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    m_new[i] = beta1 * m[i] + (1.0 - beta1) * grad[i];
    v_new[i] = beta2 * v[i] + (1.0 - beta2) * grad[i] * grad[i];
    double m_hat = m_new[i] / bc1;
    double v_hat = v_new[i] / bc2;
    param_new[i] = param[i] - lr * m_hat / (std::sqrt(v_hat) + eps);
  }
  return Rcpp::List::create(
    Rcpp::Named("param") = param_new,
    Rcpp::Named("m") = m_new,
    Rcpp::Named("v") = v_new
  );
}

// =============================================================================
// Matrix softplus activation and gradient (OpenMP parallel)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix softplus_mat(Rcpp::NumericMatrix x, double threshold = 20.0) {
  int nr = x.nrow();
  int nc = x.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      double val = x(i, j);
      out(i, j) = (val > threshold) ? val : std::log1p(std::exp(val));
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix softplus_grad_mat(Rcpp::NumericMatrix x, double threshold = 20.0) {
  int nr = x.nrow();
  int nc = x.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      double val = x(i, j);
      out(i, j) = (val > threshold) ? 1.0 : 1.0 / (1.0 + std::exp(-val));
    }
  }
  return out;
}

// =============================================================================
// Element-wise matrix clipping
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix clip_mat(Rcpp::NumericMatrix x, double lo, double hi) {
  int nr = x.nrow();
  int nc = x.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      double val = x(i, j);
      out(i, j) = (val < lo) ? lo : ((val > hi) ? hi : val);
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
    info(i, 0) = n_samples * trigamma_a;
    info(i, 1) = n_samples * a * g * g / (b * b);
    info(i, 2) = n_samples / (g * g);
  }
  return info;
}

// =============================================================================
// Compute Log Fold Change and its standard error
// =============================================================================
// [[Rcpp::export]]
Rcpp::List compute_lfc(Rcpp::NumericVector treat_mean,
                        Rcpp::NumericVector control_mean,
                        Rcpp::NumericVector treat_var,
                        Rcpp::NumericVector control_var,
                        int n_treat,
                        int n_control,
                        double eps = 1e-8) {
  int n = treat_mean.size();
  Rcpp::NumericVector lfc(n);
  Rcpp::NumericVector lfc_se(n);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double tm = treat_mean[i] + eps;
    double cm = control_mean[i] + eps;
    lfc[i] = std::log2(tm) - std::log2(cm);
    double var_term = treat_var[i] / (n_treat * tm * tm) + 
                      control_var[i] / (n_control * cm * cm);
    lfc_se[i] = std::sqrt(var_term) / std::log(2.0);
  }
  return Rcpp::List::create(
    Rcpp::Named("lfc") = lfc,
    Rcpp::Named("lfc_se") = lfc_se
  );
}

// =============================================================================
// Compute statistic matrix (batch processing)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_T_stat_batch(const Rcpp::NumericMatrix& X,
                                          const Rcpp::NumericVector& alpha,
                                          const Rcpp::NumericVector& beta,
                                          const Rcpp::NumericVector& gamma,
                                          int n_genes,
                                          int n_samples,
                                          int n_batches,
                                          int batch_size,
                                          double eps = 1e-8) {
  Rcpp::NumericMatrix T_stats(n_genes, n_batches);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) collapse(2)
  #endif
  for (int b = 0; b < n_batches; b++) {
    for (int g = 0; g < n_genes; g++) {
      const int start_col = b * batch_size;
      const int end_col = std::min(start_col + batch_size, n_samples);
      const int actual_size = end_col - start_col;
      
      // Pre-compute GG parameters
      GGParams gg;
      gg.precompute(alpha[g], beta[g], gamma[g]);
      
      double sum_nll = 0.0;
      for (int s = start_col; s < end_col; s++) {
        sum_nll -= gg.loglik(X(g, s));
      }
      T_stats(g, b) = sum_nll / actual_size;
    }
  }
  return T_stats;
}

// =============================================================================
// Calibrated Bootstrap: coverage computation for m search
// =============================================================================
// [[Rcpp::export]]
double compute_coverage_error(Rcpp::NumericMatrix T_null_mat,
                               Rcpp::NumericVector T_obs,
                               double target_coverage = 0.5) {
  int n_genes = T_null_mat.nrow();
  int n_reps = T_null_mat.ncol();
  double total_error = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:total_error) schedule(static)
  #endif
  for (int g = 0; g < n_genes; g++) {
    int count = 0;
    double t_obs = T_obs[g];
    for (int r = 0; r < n_reps; r++) {
      count += (T_null_mat(g, r) <= t_obs) ? 1 : 0;
    }
    double emp_cov = (double)count / n_reps;
    total_error += std::abs(emp_cov - target_coverage);
  }
  return total_error / n_genes;
}

// =============================================================================
// Compute empirical p-values (vectorized)
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericVector compute_empirical_pvalues(Rcpp::NumericMatrix T_null_mat,
                                               Rcpp::NumericVector T_obs) {
  int n_genes = T_null_mat.nrow();
  int B = T_null_mat.ncol();
  Rcpp::NumericVector pvals(n_genes);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int g = 0; g < n_genes; g++) {
    int count = 0;
    double t_obs = T_obs[g];
    for (int b = 0; b < B; b++) {
      count += (T_null_mat(g, b) >= t_obs) ? 1 : 0;
    }
    pvals[g] = (1.0 + count) / (1.0 + B);
  }
  return pvals;
}
