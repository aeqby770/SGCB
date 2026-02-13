// =============================================================================
// SGCB parameter estimation acceleration module
// Rcpp-accelerated versions of hotspot functions: local regression, TMM normalization, etc.
// =============================================================================

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "sgcb_optimized.hpp"
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace sgcb;

// -----------------------------------------------------------------------------
// Local weighted regression (simplified LOESS) - parallel-accelerated
// Used for dispersion trend fitting and tagwise estimation
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector local_weighted_regression_cpp(const NumericVector& x, const NumericVector& y,
                                             double span = 0.3, int n_threads = 0) {
    const int n = x.size();
    const int n_neighbors = std::max(3, static_cast<int>(n * span));
    NumericVector fitted(n);
    
    // Raw data pointers
    const double* SGCB_RESTRICT px = &x[0];
    const double* SGCB_RESTRICT py = &y[0];
    double* SGCB_RESTRICT pf = &fitted[0];
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        // Thread-local buffer
        std::vector<std::pair<double, int>> dists(n);
        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; i++) {
            const double xi = px[i];
            
            // Compute distances
            for (int j = 0; j < n; j++) {
                dists[j] = {std::abs(xi - px[j]), j};
            }
            
            // nth_element is faster than partial_sort (O(n) vs O(n log k))
            std::nth_element(dists.begin(), dists.begin() + n_neighbors, dists.end());
            
            // Tricube weights
            const double max_dist = dists[n_neighbors - 1].first + EPS;
            double sum_wy = 0.0, sum_w = 0.0;
            
            for (int k = 0; k < n_neighbors; k++) {
                const int idx = dists[k].second;
                const double u = dists[k].first / max_dist;
                const double w = cube(1.0 - cube(u));
                sum_wy += w * py[idx];
                sum_w += w;
            }
            
            pf[i] = sum_wy / (sum_w + EPS);
        }
    }
    #else
    std::vector<std::pair<double, int>> dists(n);
    for (int i = 0; i < n; i++) {
        const double xi = px[i];
        for (int j = 0; j < n; j++) {
            dists[j] = {std::abs(xi - px[j]), j};
        }
        std::nth_element(dists.begin(), dists.begin() + n_neighbors, dists.end());
        const double max_dist = dists[n_neighbors - 1].first + EPS;
        double sum_wy = 0.0, sum_w = 0.0;
        for (int k = 0; k < n_neighbors; k++) {
            const int idx = dists[k].second;
            const double u = dists[k].first / max_dist;
            const double w = cube(1.0 - cube(u));
            sum_wy += w * py[idx];
            sum_w += w;
        }
        pf[i] = sum_wy / (sum_w + EPS);
    }
    #endif
    
    return fitted;
}

// -----------------------------------------------------------------------------
// Batch column median computation - parallel-accelerated
// Used for apply(x, 2, median) in normalization
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector column_medians_cpp(const NumericMatrix& X, int n_threads = 0) {
    const int nr = X.nrow();
    const int nc = X.ncol();
    const int mid = nr / 2;
    NumericVector medians(nc);
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        std::vector<double> col(nr);  // thread-local buffer
        #pragma omp for schedule(static)
        for (int j = 0; j < nc; j++) {
            for (int i = 0; i < nr; i++) {
                col[i] = X(i, j);
            }
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            medians[j] = col[mid];
        }
    }
    #else
    std::vector<double> col(nr);
    for (int j = 0; j < nc; j++) {
        for (int i = 0; i < nr; i++) {
            col[i] = X(i, j);
        }
        std::nth_element(col.begin(), col.begin() + mid, col.end());
        medians[j] = col[mid];
    }
    #endif
    
    return medians;
}

// [[Rcpp::export]]
NumericVector row_medians_cpp(const NumericMatrix& X, int n_threads = 0) {
    const int nr = X.nrow();
    const int nc = X.ncol();
    const int mid = nc / 2;
    NumericVector medians(nr);
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        std::vector<double> row(nc);  // thread-local buffer
        #pragma omp for schedule(static)
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                row[j] = X(i, j);
            }
            std::nth_element(row.begin(), row.begin() + mid, row.end());
            medians[i] = row[mid];
        }
    }
    #else
    std::vector<double> row(nc);
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            row[j] = X(i, j);
        }
        std::nth_element(row.begin(), row.begin() + mid, row.end());
        medians[i] = row[mid];
    }
    #endif
    
    return medians;
}

// -----------------------------------------------------------------------------
// TMM normalization factor computation - batch parallel version
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector tmm_factors_cpp(NumericMatrix counts, int ref_idx,
                               double trim_m = 0.3, double trim_a = 0.05,
                               double eps = 1e-8, int n_threads = 0) {
    int n_genes = counts.nrow();
    int n_samples = counts.ncol();
    
    // Compute library sizes
    NumericVector lib_sizes(n_samples);
    for (int j = 0; j < n_samples; j++) {
        double sum = 0.0;
        for (int i = 0; i < n_genes; i++) {
            sum += counts(i, j);
        }
        lib_sizes[j] = sum;
    }
    
    double lib_ref = lib_sizes[ref_idx];
    NumericVector tmm_factors(n_samples);
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int s = 0; s < n_samples; s++) {
        if (s == ref_idx) {
            tmm_factors[s] = 1.0;
            continue;
        }
        
        double lib_obs = lib_sizes[s];
        
        // Collect M and A values for valid genes
        std::vector<double> M_vals, A_vals, W_vals;
        M_vals.reserve(n_genes);
        A_vals.reserve(n_genes);
        W_vals.reserve(n_genes);
        
        for (int g = 0; g < n_genes; g++) {
            double obs = counts(g, s);
            double ref = counts(g, ref_idx);
            
            if (obs <= 0 || ref <= 0) continue;
            
            double obs_norm = obs / lib_obs;
            double ref_norm = ref / lib_ref;
            
            double M = std::log2(obs_norm / ref_norm + eps);
            double A = 0.5 * std::log2(obs_norm * ref_norm + eps);
            double W = (lib_obs - obs) / (lib_obs * obs + eps) +
                       (lib_ref - ref) / (lib_ref * ref + eps);
            
            M_vals.push_back(M);
            A_vals.push_back(A);
            W_vals.push_back(W);
        }
        
        int n_valid = M_vals.size();
        if (n_valid < 10) {
            tmm_factors[s] = 1.0;
            continue;
        }
        
        // Trimming
        int n_trim_m = (int)(n_valid * trim_m);
        int n_trim_a = (int)(n_valid * trim_a);
        
        // Get sorted indices for M and A
        std::vector<int> order_m(n_valid), order_a(n_valid);
        for (int i = 0; i < n_valid; i++) {
            order_m[i] = i;
            order_a[i] = i;
        }
        
        std::sort(order_m.begin(), order_m.end(),
                  [&M_vals](int a, int b) { return M_vals[a] < M_vals[b]; });
        std::sort(order_a.begin(), order_a.end(),
                  [&A_vals](int a, int b) { return A_vals[a] < A_vals[b]; });
        
        // Find retained indices
        std::vector<bool> keep_m(n_valid, false), keep_a(n_valid, false);
        for (int i = n_trim_m; i < n_valid - n_trim_m; i++) {
            keep_m[order_m[i]] = true;
        }
        for (int i = n_trim_a; i < n_valid - n_trim_a; i++) {
            keep_a[order_a[i]] = true;
        }
        
        // Weighted mean
        double sum_wm = 0.0, sum_w = 0.0;
        for (int i = 0; i < n_valid; i++) {
            if (keep_m[i] && keep_a[i]) {
                double w = 1.0 / (W_vals[i] + eps);
                sum_wm += w * M_vals[i];
                sum_w += w;
            }
        }
        
        double tmm = (sum_w > eps) ? (sum_wm / sum_w) : 0.0;
        tmm_factors[s] = std::pow(2.0, tmm);
    }
    
    return tmm_factors;
}

// -----------------------------------------------------------------------------
// Dispersion trend fitting - parametric version (vectorized)
// log(disp) = a + b/mean
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List fit_dispersion_trend_parametric_cpp(NumericVector means, NumericVector dispersions,
                                          double eps = 1e-8) {
    int n = means.size();
    
    // Filter valid values
    std::vector<int> valid_idx;
    valid_idx.reserve(n);
    for (int i = 0; i < n; i++) {
        if (std::isfinite(means[i]) && std::isfinite(dispersions[i]) &&
            means[i] > eps && dispersions[i] > eps) {
            valid_idx.push_back(i);
        }
    }
    int n_valid = valid_idx.size();
    
    // Design matrix: [1, 1/mean]
    // Normal equations: (X'X)^{-1} X'y
    double sum_1 = n_valid;
    double sum_inv = 0.0, sum_inv2 = 0.0;
    double sum_y = 0.0, sum_inv_y = 0.0;
    
    for (int i = 0; i < n_valid; i++) {
        int idx = valid_idx[i];
        double inv_mean = 1.0 / (means[idx] + eps);
        double log_disp = std::log(dispersions[idx] + eps);
        
        sum_inv += inv_mean;
        sum_inv2 += inv_mean * inv_mean;
        sum_y += log_disp;
        sum_inv_y += inv_mean * log_disp;
    }
    
    // 2x2 matrix inversion
    double det = sum_1 * sum_inv2 - sum_inv * sum_inv;
    double a = (sum_inv2 * sum_y - sum_inv * sum_inv_y) / (det + eps);
    double b = (sum_1 * sum_inv_y - sum_inv * sum_y) / (det + eps);
    
    // Compute fitted values
    NumericVector fitted(n);
    for (int i = 0; i < n; i++) {
        double log_fitted = a + b / (means[i] + eps);
        fitted[i] = std::exp(log_fitted);
    }
    
    return List::create(
        Named("fitted") = fitted,
        Named("coef_a") = a,
        Named("coef_b") = b
    );
}

// -----------------------------------------------------------------------------
// Tagwise dispersion estimation - with empirical Bayes shrinkage
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List estimate_disp_tagwise_cpp(NumericMatrix counts, double span = 0.3,
                                double prior_df = 10.0, double eps = 1e-8,
                                int n_threads = 0) {
    int n_genes = counts.nrow();
    int n_samples = counts.ncol();
    
    // Compute gene-level statistics
    NumericVector means(n_genes), vars(n_genes), disp_gene(n_genes);
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum = 0.0, sum_sq = 0.0;
        for (int j = 0; j < n_samples; j++) {
            double val = counts(g, j);
            sum += val;
            sum_sq += val * val;
        }
        double mean = sum / n_samples;
        double var = sum_sq / n_samples - mean * mean;
        means[g] = mean;
        vars[g] = var;
        
        // NB dispersion estimate
        double excess_var = std::max(var - mean, eps);
        disp_gene[g] = std::max(excess_var / (mean * mean + eps), eps);
    }
    
    // Log-transform
    NumericVector log_means(n_genes), log_disp(n_genes);
    for (int g = 0; g < n_genes; g++) {
        log_means[g] = std::log(means[g] + eps);
        log_disp[g] = std::log(disp_gene[g] + eps);
    }
    
    // Local weighted regression to fit prior
    NumericVector fitted_log_disp = local_weighted_regression_cpp(log_means, log_disp, span, n_threads);
    
    // Empirical Bayes shrinkage
    double shrinkage = prior_df / (prior_df + n_samples - 1);
    NumericVector disp_shrunk(n_genes);
    
    for (int g = 0; g < n_genes; g++) {
        double log_shrunk = shrinkage * fitted_log_disp[g] + (1.0 - shrinkage) * log_disp[g];
        disp_shrunk[g] = std::exp(log_shrunk);
    }
    
    return List::create(
        Named("raw") = disp_gene,
        Named("shrunk") = disp_shrunk,
        Named("prior") = NumericVector::create(std::exp(fitted_log_disp[0])),
        Named("log_means") = log_means,
        Named("log_disp") = log_disp,
        Named("fitted_log_disp") = fitted_log_disp
    );
}

// -----------------------------------------------------------------------------
// Fast quantile computation (batch)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix row_quantiles_cpp(NumericMatrix X, NumericVector probs, int n_threads = 0) {
    int nr = X.nrow();
    int nc = X.ncol();
    int n_probs = probs.size();
    NumericMatrix out(nr, n_probs);
    
    #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < nr; i++) {
        std::vector<double> row(nc);
        for (int j = 0; j < nc; j++) {
            row[j] = X(i, j);
        }
        std::sort(row.begin(), row.end());
        
        for (int p = 0; p < n_probs; p++) {
            double idx_d = probs[p] * (nc - 1);
            int idx_lo = (int)idx_d;
            int idx_hi = std::min(idx_lo + 1, nc - 1);
            double frac = idx_d - idx_lo;
            out(i, p) = row[idx_lo] * (1.0 - frac) + row[idx_hi] * frac;
        }
    }
    
    return out;
}
