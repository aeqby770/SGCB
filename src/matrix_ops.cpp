// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// =============================================================================
// High-performance matrix operations
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix matmul_parallel(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int m = A.nrow();
  int k = A.ncol();
  int n = B.ncol();
  
  Rcpp::NumericMatrix C(m, n);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int l = 0; l < k; l++) {
        sum += A(i, l) * B(l, j);
      }
      C(i, j) = sum;
    }
  }
  
  return C;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix tcrossprod_parallel(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int m = A.nrow();
  int k = A.ncol();
  int n = B.nrow();
  
  Rcpp::NumericMatrix C(m, n);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int l = 0; l < k; l++) {
        sum += A(i, l) * B(j, l);
      }
      C(i, j) = sum;
    }
  }
  
  return C;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix crossprod_parallel(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int m = A.ncol();
  int k = A.nrow();
  int n = B.ncol();
  
  Rcpp::NumericMatrix C(m, n);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int l = 0; l < k; l++) {
        sum += A(l, i) * B(l, j);
      }
      C(i, j) = sum;
    }
  }
  
  return C;
}

// =============================================================================
// Row/column statistics
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector rowMeans_parallel(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(nr);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    double sum = 0.0;
    for (int j = 0; j < nc; j++) {
      sum += X(i, j);
    }
    out[i] = sum / nc;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowVars_parallel(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(nr);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int j = 0; j < nc; j++) {
      double val = X(i, j);
      sum += val;
      sum_sq += val * val;
    }
    double mean = sum / nc;
    out[i] = sum_sq / nc - mean * mean;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowSums_parallel(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(nr);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    double sum = 0.0;
    for (int j = 0; j < nc; j++) {
      sum += X(i, j);
    }
    out[i] = sum;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector colSums_parallel(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(nc);
  
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:sum) schedule(static)
    #endif
    for (int i = 0; i < nr; i++) {
      sum += X(i, j);
    }
    out[j] = sum;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowMedians_parallel(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericVector out(nr);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    std::vector<double> row_vals(nc);
    for (int j = 0; j < nc; j++) {
      row_vals[j] = X(i, j);
    }
    std::sort(row_vals.begin(), row_vals.end());
    
    if (nc % 2 == 0) {
      out[i] = (row_vals[nc / 2 - 1] + row_vals[nc / 2]) / 2.0;
    } else {
      out[i] = row_vals[nc / 2];
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rowQuantiles_parallel(Rcpp::NumericMatrix X, 
                                           Rcpp::NumericVector probs) {
  int nr = X.nrow();
  int nc = X.ncol();
  int np = probs.size();
  Rcpp::NumericMatrix out(nr, np);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    std::vector<double> row_vals(nc);
    for (int j = 0; j < nc; j++) {
      row_vals[j] = X(i, j);
    }
    std::sort(row_vals.begin(), row_vals.end());
    
    for (int p = 0; p < np; p++) {
      double prob = probs[p];
      double idx_f = prob * (nc - 1);
      int idx_lo = (int)std::floor(idx_f);
      int idx_hi = (int)std::ceil(idx_f);
      double frac = idx_f - idx_lo;
      
      if (idx_lo == idx_hi) {
        out(i, p) = row_vals[idx_lo];
      } else {
        out(i, p) = row_vals[idx_lo] * (1 - frac) + row_vals[idx_hi] * frac;
      }
    }
  }
  
  return out;
}

// =============================================================================
// Matrix transformations
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix log2_transform(Rcpp::NumericMatrix X, double pseudo = 1.0) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  double log2_const = std::log(2.0);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      out(i, j) = std::log(X(i, j) + pseudo) / log2_const;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix zscore_rows(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  double eps = 1e-8;
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    // Compute mean and standard deviation
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int j = 0; j < nc; j++) {
      double val = X(i, j);
      sum += val;
      sum_sq += val * val;
    }
    double mean = sum / nc;
    double var = sum_sq / nc - mean * mean;
    double sd = std::sqrt(var + eps);
    
    // Standardize
    for (int j = 0; j < nc; j++) {
      out(i, j) = (X(i, j) - mean) / sd;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_rows(Rcpp::NumericMatrix X, Rcpp::NumericVector stats,
                                 int op = 1) {
  // op: 1=subtract, 2=divide, 3=add, 4=multiply
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  double eps = 1e-8;
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      double s = stats[i];
      switch (op) {
        case 1: out(i, j) = X(i, j) - s; break;
        case 2: out(i, j) = X(i, j) / (s + eps); break;
        case 3: out(i, j) = X(i, j) + s; break;
        case 4: out(i, j) = X(i, j) * s; break;
        default: out(i, j) = X(i, j) - s;
      }
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_cols(Rcpp::NumericMatrix X, Rcpp::NumericVector stats,
                                 int op = 1) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  double eps = 1e-8;
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      double s = stats[j];
      switch (op) {
        case 1: out(i, j) = X(i, j) - s; break;
        case 2: out(i, j) = X(i, j) / (s + eps); break;
        case 3: out(i, j) = X(i, j) + s; break;
        case 4: out(i, j) = X(i, j) * s; break;
        default: out(i, j) = X(i, j) - s;
      }
    }
  }
  
  return out;
}

// =============================================================================
// Distance computation
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix euclidean_dist_rows(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix D(nr, nr);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    for (int j = i; j < nr; j++) {
      double sum_sq = 0.0;
      for (int k = 0; k < nc; k++) {
        double diff = X(i, k) - X(j, k);
        sum_sq += diff * diff;
      }
      double dist = std::sqrt(sum_sq);
      D(i, j) = dist;
      D(j, i) = dist;
    }
  }
  
  return D;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix correlation_matrix(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::NumericMatrix R(nc, nc);
  double eps = 1e-8;
  
  // Compute column means and standard deviations
  std::vector<double> col_mean(nc);
  std::vector<double> col_sd(nc);
  
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += X(i, j);
      sum_sq += X(i, j) * X(i, j);
    }
    col_mean[j] = sum / nr;
    col_sd[j] = std::sqrt(sum_sq / nr - col_mean[j] * col_mean[j] + eps);
  }
  
  // Compute correlation coefficients
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int j1 = 0; j1 < nc; j1++) {
    for (int j2 = j1; j2 < nc; j2++) {
      double sum_prod = 0.0;
      for (int i = 0; i < nr; i++) {
        sum_prod += (X(i, j1) - col_mean[j1]) * (X(i, j2) - col_mean[j2]);
      }
      double cor = sum_prod / (nr * col_sd[j1] * col_sd[j2]);
      R(j1, j2) = cor;
      R(j2, j1) = cor;
    }
  }
  
  return R;
}

// =============================================================================
// Sorting and ranking
// =============================================================================

// [[Rcpp::export]]
Rcpp::IntegerMatrix rank_rows(Rcpp::NumericMatrix X) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::IntegerMatrix out(nr, nc);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    // Create index-value pairs
    std::vector<std::pair<double, int>> pairs(nc);
    for (int j = 0; j < nc; j++) {
      pairs[j] = std::make_pair(X(i, j), j);
    }
    
    // Sort
    std::sort(pairs.begin(), pairs.end());
    
    // Assign ranks
    for (int j = 0; j < nc; j++) {
      out(i, pairs[j].second) = j + 1;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix order_rows(Rcpp::NumericMatrix X, bool decreasing = false) {
  int nr = X.nrow();
  int nc = X.ncol();
  Rcpp::IntegerMatrix out(nr, nc);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nr; i++) {
    std::vector<std::pair<double, int>> pairs(nc);
    for (int j = 0; j < nc; j++) {
      pairs[j] = std::make_pair(X(i, j), j);
    }
    
    if (decreasing) {
      std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double, int>>());
    } else {
      std::sort(pairs.begin(), pairs.end());
    }
    
    for (int j = 0; j < nc; j++) {
      out(i, j) = pairs[j].second + 1; // 1-based index
    }
  }
  
  return out;
}

// =============================================================================
// Special functions
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector digamma_vec(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    out[i] = R::digamma(x[i]);
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector trigamma_vec(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    out[i] = R::trigamma(x[i]);
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector lgamma_vec(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    out[i] = std::lgamma(x[i]);
  }
  
  return out;
}

// =============================================================================
// Negative binomial distribution functions (compatibility with other methods)
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector dnbinom_vec(Rcpp::NumericVector x, Rcpp::NumericVector size,
                                 Rcpp::NumericVector mu, bool log_p = false) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    out[i] = R::dnbinom_mu(x[i], size[i], mu[i], log_p ? 1 : 0);
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rnbinom_mat(int n_genes, int n_samples,
                                 Rcpp::NumericVector size, Rcpp::NumericVector mu) {
  Rcpp::NumericMatrix out(n_genes, n_samples);
  
  for (int g = 0; g < n_genes; g++) {
    for (int s = 0; s < n_samples; s++) {
      out(g, s) = R::rnbinom(size[g], size[g] / (size[g] + mu[g]));
    }
  }
  
  return out;
}

// =============================================================================
// Window/sliding functions
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector rolling_mean(Rcpp::NumericVector x, int window) {
  int n = x.size();
  int half_win = window / 2;
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    int start = std::max(0, i - half_win);
    int end = std::min(n - 1, i + half_win);
    double sum = 0.0;
    int count = 0;
    for (int j = start; j <= end; j++) {
      sum += x[j];
      count++;
    }
    out[i] = sum / count;
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector loess_smooth(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                   double span = 0.3) {
  int n = x.size();
  int n_neighbors = std::max(3, (int)(n * span));
  Rcpp::NumericVector out(n);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    // Compute distances
    std::vector<std::pair<double, int>> dists(n);
    for (int j = 0; j < n; j++) {
      dists[j] = std::make_pair(std::abs(x[j] - x[i]), j);
    }
    
    // Partial sort to get nearest neighbors
    std::partial_sort(dists.begin(), dists.begin() + n_neighbors, dists.end());
    
    // Tricube kernel weights
    double max_dist = dists[n_neighbors - 1].first + 1e-8;
    double sum_w = 0.0;
    double sum_wy = 0.0;
    
    for (int k = 0; k < n_neighbors; k++) {
      double u = dists[k].first / max_dist;
      double w = std::pow(1 - u * u * u, 3);
      int idx = dists[k].second;
      sum_w += w;
      sum_wy += w * y[idx];
    }
    
    out[i] = sum_wy / sum_w;
  }
  
  return out;
}
