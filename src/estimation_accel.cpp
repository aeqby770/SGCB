// =============================================================================
// SGCB estimation utilities: TMM normalization + row quantiles
// =============================================================================

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// TMM normalization factor computation (Robinson & Oshlack 2010)
// Used by both sgcbDE and sgcbDReg pipelines
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
// Fast row quantile computation (used by sgcbDE bootstrap path)
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
