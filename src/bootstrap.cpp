// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "sgcb_optimized.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace sgcb;

// =============================================================================
// Bootstrap-related functions
// =============================================================================

// [[Rcpp::export]]
Rcpp::IntegerMatrix resample_indices(int n, int m, int B, int seed = 12345) {
  // Generate B sets of resampling indices of size m
  Rcpp::IntegerMatrix indices(m, B);
  
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    int thread_id = omp_get_thread_num();
    #else
    int thread_id = 0;
    #endif
    
    std::mt19937 gen(seed + thread_id);
    std::uniform_int_distribution<int> dist(0, n - 1);
    
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int b = 0; b < B; b++) {
      for (int i = 0; i < m; i++) {
        indices(i, b) = dist(gen);
      }
    }
  }
  
  return indices;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrap_resample(Rcpp::NumericMatrix X, 
                                         Rcpp::IntegerVector indices) {
  int nr = X.nrow();
  int nc = indices.size();
  Rcpp::NumericMatrix out(nr, nc);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      out(i, j) = X(i, indices[j]);
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector bootstrap_mean(Rcpp::NumericMatrix X, int B, int seed = 12345) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(B);
  
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(0, nc - 1);
  
  for (int b = 0; b < B; b++) {
    double sum = 0.0;
    for (int j = 0; j < nc; j++) {
      int idx = dist(gen);
      for (int i = 0; i < nr; i++) {
        sum += X(i, idx);
      }
    }
    out[b] = sum / (nr * nc);
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrap_gene_stats(const Rcpp::NumericMatrix& X, 
                                           const Rcpp::NumericVector& alpha,
                                           const Rcpp::NumericVector& beta,
                                           const Rcpp::NumericVector& gamma,
                                           int m, int B,
                                           double eps = 1e-8,
                                           int seed = 12345) {
  const int n_genes = X.nrow();
  const int n_samples = X.ncol();
  Rcpp::NumericMatrix T_boot(n_genes, B);
  
  // Pre-allocate index storage
  std::vector<std::vector<int>> all_indices(B);
  for (int b = 0; b < B; b++) {
    all_indices[b].reserve(m);
  }
  
  // Pre-generate all resampling indices
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(0, n_samples - 1);
  for (int b = 0; b < B; b++) {
    all_indices[b].resize(m);
    for (int i = 0; i < m; i++) {
      all_indices[b][i] = dist(gen);
    }
  }
  
  // Compute all bootstrap statistics in parallel
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) collapse(2)
  #endif
  for (int b = 0; b < B; b++) {
    for (int g = 0; g < n_genes; g++) {
      // Pre-compute GG parameters
      GGParams gg;
      gg.precompute(alpha[g], beta[g], gamma[g]);
      
      double sum_nll = 0.0;
      const auto& indices = all_indices[b];
      for (int i = 0; i < m; i++) {
        sum_nll -= gg.loglik(X(g, indices[i]));
      }
      T_boot(g, b) = sum_nll / m;
    }
  }
  
  return T_boot;
}

// [[Rcpp::export]]
Rcpp::List calibrated_bootstrap_search(Rcpp::NumericMatrix X,
                                         Rcpp::NumericVector alpha,
                                         Rcpp::NumericVector beta,
                                         Rcpp::NumericVector gamma,
                                         Rcpp::NumericVector T_obs,
                                         Rcpp::NumericVector m_fracs,
                                         int cb_reps,
                                         double target_coverage = 0.5,
                                         double eps = 1e-8,
                                         int seed = 12345) {
  
  int n_genes = X.nrow();
  int n_samples = X.ncol();
  int n_m = m_fracs.size();
  
  Rcpp::NumericVector coverage_errors(n_m);
  
  for (int mi = 0; mi < n_m; mi++) {
    double m_frac = m_fracs[mi];
    int m = std::max(1, (int)(n_samples * m_frac));
    
    // Generate bootstrap statistics
    Rcpp::NumericMatrix T_null = bootstrap_gene_stats(
      X, alpha, beta, gamma, m, cb_reps, eps, seed + mi);
    
    // Compute coverage error
    double total_error = 0.0;
    
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:total_error) schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
      int count = 0;
      double t_obs = T_obs[g];
      for (int b = 0; b < cb_reps; b++) {
        if (T_null(g, b) <= t_obs) count++;
      }
      double emp_cov = (double)count / cb_reps;
      total_error += std::abs(emp_cov - target_coverage);
    }
    
    coverage_errors[mi] = total_error / n_genes;
  }
  
  // Find optimal m
  int best_idx = 0;
  double min_error = coverage_errors[0];
  for (int i = 1; i < n_m; i++) {
    if (coverage_errors[i] < min_error) {
      min_error = coverage_errors[i];
      best_idx = i;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("coverage_errors") = coverage_errors,
    Rcpp::Named("best_m_frac") = m_fracs[best_idx],
    Rcpp::Named("best_m") = std::max(1, (int)(n_samples * m_fracs[best_idx])),
    Rcpp::Named("min_error") = min_error
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector permutation_pvalues(Rcpp::NumericVector T_obs,
                                          Rcpp::NumericMatrix T_null) {
  int n = T_obs.size();
  int B = T_null.ncol();
  Rcpp::NumericVector pvals(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    int count = 0;
    double t_obs = T_obs[i];
    for (int b = 0; b < B; b++) {
      if (T_null(i, b) >= t_obs) count++;
    }
    pvals[i] = (1.0 + count) / (1.0 + B);
  }
  
  return pvals;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrap_confidence_intervals(const Rcpp::NumericMatrix& T_null,
                                                     double alpha = 0.05) {
  const int n = T_null.nrow();
  const int B = T_null.ncol();
  Rcpp::NumericMatrix ci(n, 3);
  
  const int lo_idx = static_cast<int>((alpha / 2.0) * (B - 1));
  const int med_idx = B / 2;
  const int hi_idx = static_cast<int>((1.0 - alpha / 2.0) * (B - 1));
  
  #ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<double> vals(B);  // thread-local buffer
    #pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
      for (int b = 0; b < B; b++) {
        vals[b] = T_null(i, b);
      }
      // Use nth_element instead of full sort (O(n) vs O(n log n))
      std::nth_element(vals.begin(), vals.begin() + hi_idx, vals.end());
      std::nth_element(vals.begin(), vals.begin() + med_idx, vals.begin() + hi_idx);
      std::nth_element(vals.begin(), vals.begin() + lo_idx, vals.begin() + med_idx);
      
      ci(i, 0) = vals[lo_idx];
      ci(i, 1) = vals[med_idx];
      ci(i, 2) = vals[hi_idx];
    }
  }
  #else
  std::vector<double> vals(B);
  for (int i = 0; i < n; i++) {
    for (int b = 0; b < B; b++) {
      vals[b] = T_null(i, b);
    }
    std::nth_element(vals.begin(), vals.begin() + hi_idx, vals.end());
    std::nth_element(vals.begin(), vals.begin() + med_idx, vals.begin() + hi_idx);
    std::nth_element(vals.begin(), vals.begin() + lo_idx, vals.begin() + med_idx);
    ci(i, 0) = vals[lo_idx];
    ci(i, 1) = vals[med_idx];
    ci(i, 2) = vals[hi_idx];
  }
  #endif
  
  return ci;
}

// =============================================================================
// Multiple comparison correction
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector bh_adjust(Rcpp::NumericVector pvals) {
  int n = pvals.size();
  Rcpp::NumericVector padj(n);
  
  // Create index-value pairs
  std::vector<std::pair<double, int>> pairs(n);
  for (int i = 0; i < n; i++) {
    pairs[i] = std::make_pair(pvals[i], i);
  }
  
  // Sort by p-value
  std::sort(pairs.begin(), pairs.end());
  
  // BH adjustment
  double cummin = 1.0;
  for (int i = n - 1; i >= 0; i--) {
    double p = pairs[i].first;
    double adj = p * n / (i + 1);
    adj = std::min(adj, 1.0);
    cummin = std::min(cummin, adj);
    padj[pairs[i].second] = cummin;
  }
  
  return padj;
}

// [[Rcpp::export]]
Rcpp::NumericVector bonferroni_adjust(Rcpp::NumericVector pvals) {
  int n = pvals.size();
  Rcpp::NumericVector padj(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    padj[i] = std::min(pvals[i] * n, 1.0);
  }
  
  return padj;
}

// [[Rcpp::export]]
Rcpp::NumericVector holm_adjust(Rcpp::NumericVector pvals) {
  int n = pvals.size();
  Rcpp::NumericVector padj(n);
  
  // Create index-value pairs
  std::vector<std::pair<double, int>> pairs(n);
  for (int i = 0; i < n; i++) {
    pairs[i] = std::make_pair(pvals[i], i);
  }
  
  // Sort by p-value
  std::sort(pairs.begin(), pairs.end());
  
  // Holm adjustment
  double cummax = 0.0;
  for (int i = 0; i < n; i++) {
    double p = pairs[i].first;
    double adj = p * (n - i);
    adj = std::min(adj, 1.0);
    cummax = std::max(cummax, adj);
    padj[pairs[i].second] = cummax;
  }
  
  return padj;
}

// =============================================================================
// Effect size computation
// =============================================================================

// [[Rcpp::export]]
Rcpp::List effect_sizes(Rcpp::NumericMatrix control, Rcpp::NumericMatrix treat,
                         double eps = 1e-8) {
  int n_genes = control.nrow();
  int n_ctrl = control.ncol();
  int n_trt = treat.ncol();
  
  Rcpp::NumericVector lfc(n_genes);
  Rcpp::NumericVector lfc_se(n_genes);
  Rcpp::NumericVector cohens_d(n_genes);
  Rcpp::NumericVector hedges_g(n_genes);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int g = 0; g < n_genes; g++) {
    // Compute means
    double sum_ctrl = 0.0, sum_trt = 0.0;
    double sum_sq_ctrl = 0.0, sum_sq_trt = 0.0;
    
    for (int j = 0; j < n_ctrl; j++) {
      double val = control(g, j);
      sum_ctrl += val;
      sum_sq_ctrl += val * val;
    }
    for (int j = 0; j < n_trt; j++) {
      double val = treat(g, j);
      sum_trt += val;
      sum_sq_trt += val * val;
    }
    
    double mean_ctrl = sum_ctrl / n_ctrl;
    double mean_trt = sum_trt / n_trt;
    double var_ctrl = sum_sq_ctrl / n_ctrl - mean_ctrl * mean_ctrl;
    double var_trt = sum_sq_trt / n_trt - mean_trt * mean_trt;
    
    // Log fold change
    double log2e = 1.0 / std::log(2.0);
    lfc[g] = std::log(mean_trt + eps) * log2e - std::log(mean_ctrl + eps) * log2e;
    
    // LFC standard error
    double var_term = var_trt / (n_trt * (mean_trt + eps) * (mean_trt + eps)) + 
                      var_ctrl / (n_ctrl * (mean_ctrl + eps) * (mean_ctrl + eps));
    lfc_se[g] = std::sqrt(var_term) * log2e;
    
    // Cohen's d
    double pooled_var = ((n_ctrl - 1) * var_ctrl + (n_trt - 1) * var_trt) / 
                         (n_ctrl + n_trt - 2);
    double pooled_sd = std::sqrt(pooled_var + eps);
    cohens_d[g] = (mean_trt - mean_ctrl) / pooled_sd;
    
    // Hedges' g (small-sample bias correction)
    double correction = 1.0 - 3.0 / (4.0 * (n_ctrl + n_trt - 2) - 1);
    hedges_g[g] = cohens_d[g] * correction;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("log2FoldChange") = lfc,
    Rcpp::Named("lfcSE") = lfc_se,
    Rcpp::Named("cohensD") = cohens_d,
    Rcpp::Named("hedgesG") = hedges_g
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector shrink_lfc_normal(Rcpp::NumericVector lfc,
                                        Rcpp::NumericVector lfc_se,
                                        double prior_scale = 1.0) {
  int n = lfc.size();
  Rcpp::NumericVector shrunk(n);
  double prior_var = prior_scale * prior_scale;
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double se2 = lfc_se[i] * lfc_se[i];
    double shrinkage = prior_var / (prior_var + se2);
    shrunk[i] = lfc[i] * shrinkage;
  }
  
  return shrunk;
}

// [[Rcpp::export]]
Rcpp::NumericVector shrink_lfc_cauchy(Rcpp::NumericVector lfc,
                                        Rcpp::NumericVector lfc_se,
                                        double prior_scale = 1.0) {
  int n = lfc.size();
  Rcpp::NumericVector shrunk(n);
  double prior_var = prior_scale * prior_scale;
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double se2 = lfc_se[i] * lfc_se[i];
    double u = lfc[i] / prior_scale;
    double weight = 1.0 / (1.0 + u * u);
    double shrinkage = se2 / (se2 + prior_var);
    shrunk[i] = lfc[i] * (1.0 - weight * shrinkage);
  }
  
  return shrunk;
}
