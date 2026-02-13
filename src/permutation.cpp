// =============================================================================
// SGCB permutation test core - Rcpp-accelerated
// =============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

inline double safe_log(double x, double eps = 1e-10) {
    return std::log(std::max(x, eps));
}

inline double safe_div(double a, double b, double eps = 1e-10) {
    return a / std::max(std::abs(b), eps);
}

// Compute row means
// [[Rcpp::export]]
NumericVector row_means_cpp(NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    NumericVector means(n_rows);
    
    #pragma omp parallel for
    for (int i = 0; i < n_rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < n_cols; j++) {
            sum += X(i, j);
        }
        means[i] = sum / n_cols;
    }
    return means;
}

// Compute row variances
// [[Rcpp::export]]
NumericVector row_vars_cpp(NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    NumericVector vars(n_rows);
    
    #pragma omp parallel for
    for (int i = 0; i < n_rows; i++) {
        double mean = 0.0;
        for (int j = 0; j < n_cols; j++) {
            mean += X(i, j);
        }
        mean /= n_cols;
        
        double ss = 0.0;
        for (int j = 0; j < n_cols; j++) {
            double diff = X(i, j) - mean;
            ss += diff * diff;
        }
        vars[i] = ss / (n_cols - 1);
    }
    return vars;
}

// -----------------------------------------------------------------------------
// Permutation test core
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix permutation_t_stats_cpp(
    NumericMatrix X,           // n_genes x n_samples
    IntegerVector group,       // 0/1 group labels
    int n_perm,                // number of permutations
    int seed = 12345
) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    
    // Count samples per group
    int n0 = 0, n1 = 0;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) n0++;
        else n1++;
    }
    
    // Output matrix: n_genes x (1 + n_perm)
    // First column is observed statistic, remaining columns are permutation statistics
    NumericMatrix T_all(n_genes, 1 + n_perm);
    
    // Create sample indices
    std::vector<int> indices(n_samples);
    for (int j = 0; j < n_samples; j++) indices[j] = j;
    
    // Set random seed
    std::mt19937 rng(seed);
    
    // Compute observed and permutation statistics
    #pragma omp parallel
    {
        std::vector<double> mean0(n_genes), mean1(n_genes);
        std::vector<double> var0(n_genes), var1(n_genes);
        std::vector<int> perm_indices(n_samples);
        
        #ifdef _OPENMP
        std::mt19937 local_rng(seed + omp_get_thread_num());
        #else
        std::mt19937 local_rng(seed);
        #endif
        
        #pragma omp for
        for (int p = 0; p <= n_perm; p++) {
            // Get current permutation indices
            if (p == 0) {
                // Observed data, use original group labels
                for (int j = 0; j < n_samples; j++) perm_indices[j] = j;
            } else {
                // Random permutation
                for (int j = 0; j < n_samples; j++) perm_indices[j] = j;
                std::shuffle(perm_indices.begin(), perm_indices.end(), local_rng);
            }
            
            // Compute statistic for each gene
            for (int i = 0; i < n_genes; i++) {
                double sum0 = 0.0, sum1 = 0.0;
                int cnt0 = 0, cnt1 = 0;
                
                // Compute means using permuted indices
                for (int j = 0; j < n_samples; j++) {
                    int idx = perm_indices[j];
                    double val = X(i, idx);
                    if (group[j] == 0) {
                        sum0 += val;
                        cnt0++;
                    } else {
                        sum1 += val;
                        cnt1++;
                    }
                }
                
                double m0 = sum0 / cnt0;
                double m1 = sum1 / cnt1;
                
                // Compute variance
                double ss0 = 0.0, ss1 = 0.0;
                for (int j = 0; j < n_samples; j++) {
                    int idx = perm_indices[j];
                    double val = X(i, idx);
                    if (group[j] == 0) {
                        ss0 += (val - m0) * (val - m0);
                    } else {
                        ss1 += (val - m1) * (val - m1);
                    }
                }
                
                // Pooled variance
                double var_pooled = (ss0 + ss1) / (n_samples - 2);
                double se = std::sqrt(var_pooled * (1.0/cnt0 + 1.0/cnt1) + 1e-10);
                
                // t-statistic
                T_all(i, p) = (m1 - m0) / se;
            }
        }
    }
    
    return T_all;
}

// -----------------------------------------------------------------------------
// Compute permutation p-values
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector permutation_pvalues_cpp(NumericMatrix T_all) {
    int n_genes = T_all.nrow();
    int n_total = T_all.ncol();  // 1 + n_perm
    int n_perm = n_total - 1;
    
    NumericVector pvals(n_genes);
    
    #pragma omp parallel for
    for (int i = 0; i < n_genes; i++) {
        double t_obs = std::abs(T_all(i, 0));
        int count = 0;
        
        for (int p = 1; p < n_total; p++) {
            if (std::abs(T_all(i, p)) >= t_obs) count++;
        }
        
        pvals[i] = (1.0 + count) / (1.0 + n_perm);
    }
    
    return pvals;
}

// -----------------------------------------------------------------------------
// Dropout permutation test (small-sample enhancement)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix dropout_permutation_cpp(
    NumericMatrix X,           // n_genes x n_samples
    IntegerVector group,       // 0/1 group labels
    int n_perm,                // permutations per round
    int n_dropout,             // number of dropout rounds
    double dropout_rate,       // dropout fraction (0.1-0.3)
    int seed = 12345
) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    int n_keep = (int)(n_genes * (1.0 - dropout_rate));
    
    // Output: n_genes x (1 + n_dropout * n_perm)
    // First column is observed statistic from complete data
    NumericMatrix T_all(n_genes, 1 + n_dropout * n_perm);
    
    std::mt19937 rng(seed);
    
    // First column: observed statistic from complete data
    NumericMatrix T_obs = permutation_t_stats_cpp(X, group, 0, seed);
    for (int i = 0; i < n_genes; i++) {
        T_all(i, 0) = T_obs(i, 0);
    }
    
    // Dropout rounds
    std::vector<int> gene_indices(n_genes);
    for (int i = 0; i < n_genes; i++) gene_indices[i] = i;
    
    for (int d = 0; d < n_dropout; d++) {
        // Randomly select genes to keep
        std::shuffle(gene_indices.begin(), gene_indices.end(), rng);
        
        // Create post-dropout matrix
        NumericMatrix X_drop(n_keep, n_samples);
        std::vector<int> kept_genes(n_keep);
        
        for (int i = 0; i < n_keep; i++) {
            kept_genes[i] = gene_indices[i];
            for (int j = 0; j < n_samples; j++) {
                X_drop(i, j) = X(gene_indices[i], j);
            }
        }
        
        // Perform permutation test on dropout data
        NumericMatrix T_drop = permutation_t_stats_cpp(X_drop, group, n_perm, seed + d);
        
        // Map results back to full gene set
        int col_offset = 1 + d * n_perm;
        for (int p = 1; p <= n_perm; p++) {
            // Genes not selected are set to 0 (neutral value)
            for (int i = 0; i < n_genes; i++) {
                T_all(i, col_offset + p - 1) = 0.0;
            }
            // Set statistics for selected genes
            for (int i = 0; i < n_keep; i++) {
                T_all(kept_genes[i], col_offset + p - 1) = T_drop(i, p);
            }
        }
    }
    
    return T_all;
}

// -----------------------------------------------------------------------------
// Dropout p-value computation (weighted)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dropout_pvalues_cpp(
    NumericMatrix T_all,
    int n_dropout,
    double dropout_rate
) {
    int n_genes = T_all.nrow();
    int n_total = T_all.ncol();
    int n_perm_per_dropout = (n_total - 1) / n_dropout;
    
    NumericVector pvals(n_genes);
    double keep_rate = 1.0 - dropout_rate;
    
    #pragma omp parallel for
    for (int i = 0; i < n_genes; i++) {
        double t_obs = std::abs(T_all(i, 0));
        double weighted_count = 0.0;
        double total_weight = 0.0;
        
        for (int d = 0; d < n_dropout; d++) {
            int col_offset = 1 + d * n_perm_per_dropout;
            
            for (int p = 0; p < n_perm_per_dropout; p++) {
                double t_perm = std::abs(T_all(i, col_offset + p));
                
                // Only non-zero statistics count (i.e., gene was retained in this dropout round)
                if (t_perm != 0.0 || t_obs == 0.0) {
                    total_weight += 1.0;
                    if (t_perm >= t_obs) weighted_count += 1.0;
                }
            }
        }
        
        // If gene was never retained, use conservative estimate
        if (total_weight < 1.0) {
            pvals[i] = 1.0;
        } else {
            pvals[i] = (1.0 + weighted_count) / (1.0 + total_weight);
        }
    }
    
    return pvals;
}

// -----------------------------------------------------------------------------
// Dropout stability computation (fast version)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix dropout_lfc_cpp(
    NumericMatrix X,           // n_genes x n_samples (normalized)
    IntegerVector group,       // 0/1 group labels
    int n_dropout,             // number of dropout rounds
    double dropout_rate,       // dropout fraction
    int seed = 12345
) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    double eps = 1e-8;
    
    // Compute sample counts and indices per group
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    int n_ctrl = ctrl_idx.size();
    int n_treat = treat_idx.size();
    
    // Output: n_genes x n_dropout
    NumericMatrix lfc_all(n_genes, n_dropout);
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    for (int d = 0; d < n_dropout; d++) {
        // Random dropout
        std::vector<int> keep_ctrl, keep_treat;
        
        for (int i = 0; i < n_ctrl; i++) {
            if (unif(rng) > dropout_rate) keep_ctrl.push_back(ctrl_idx[i]);
        }
        for (int i = 0; i < n_treat; i++) {
            if (unif(rng) > dropout_rate) keep_treat.push_back(treat_idx[i]);
        }
        
        // Ensure at least 2 samples per group
        while (keep_ctrl.size() < 2 && keep_ctrl.size() < ctrl_idx.size()) {
            int idx = rng() % n_ctrl;
            if (std::find(keep_ctrl.begin(), keep_ctrl.end(), ctrl_idx[idx]) == keep_ctrl.end()) {
                keep_ctrl.push_back(ctrl_idx[idx]);
            }
        }
        while (keep_treat.size() < 2 && keep_treat.size() < treat_idx.size()) {
            int idx = rng() % n_treat;
            if (std::find(keep_treat.begin(), keep_treat.end(), treat_idx[idx]) == keep_treat.end()) {
                keep_treat.push_back(treat_idx[idx]);
            }
        }
        
        int n_keep_ctrl = keep_ctrl.size();
        int n_keep_treat = keep_treat.size();
        
        // Compute LFC
        #pragma omp parallel for
        for (int g = 0; g < n_genes; g++) {
            double sum_ctrl = 0, sum_treat = 0;
            
            for (int i = 0; i < n_keep_ctrl; i++) {
                sum_ctrl += X(g, keep_ctrl[i]);
            }
            for (int i = 0; i < n_keep_treat; i++) {
                sum_treat += X(g, keep_treat[i]);
            }
            
            double mean_ctrl = sum_ctrl / n_keep_ctrl + eps;
            double mean_treat = sum_treat / n_keep_treat + eps;
            
            lfc_all(g, d) = std::log2(mean_treat / mean_ctrl);
        }
    }
    
    return lfc_all;
}

// Compute LFC variance (stability metric)
// [[Rcpp::export]]
NumericVector lfc_variance_cpp(NumericMatrix lfc_all) {
    int n_genes = lfc_all.nrow();
    int n_dropout = lfc_all.ncol();
    
    NumericVector vars(n_genes);
    
    #pragma omp parallel for
    for (int g = 0; g < n_genes; g++) {
        double mean = 0;
        for (int d = 0; d < n_dropout; d++) {
            mean += lfc_all(g, d);
        }
        mean /= n_dropout;
        
        double ss = 0;
        for (int d = 0; d < n_dropout; d++) {
            double diff = lfc_all(g, d) - mean;
            ss += diff * diff;
        }
        vars[g] = ss / (n_dropout - 1);
    }
    
    return vars;
}

// -----------------------------------------------------------------------------
// Log2 Fold Change
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector log2fc_cpp(
    NumericMatrix X,
    IntegerVector group,
    double eps = 1e-8
) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    
    NumericVector lfc(n_genes);
    
    #pragma omp parallel for
    for (int i = 0; i < n_genes; i++) {
        double sum0 = 0.0, sum1 = 0.0;
        int cnt0 = 0, cnt1 = 0;
        
        for (int j = 0; j < n_samples; j++) {
            if (group[j] == 0) {
                sum0 += X(i, j);
                cnt0++;
            } else {
                sum1 += X(i, j);
                cnt1++;
            }
        }
        
        double mean0 = sum0 / cnt0 + eps;
        double mean1 = sum1 / cnt1 + eps;
        lfc[i] = std::log2(mean1 / mean0);
    }
    
    return lfc;
}

// -----------------------------------------------------------------------------
// Variance shrinkage (limma eBayes-style)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector shrink_variance_cpp(
    NumericVector vars,
    double prior_df = 3.0,
    double shrinkage = 0.5
) {
    int n = vars.size();
    NumericVector shrunk(n);
    
    // Compute median as prior variance
    std::vector<double> sorted_vars(n);
    for (int i = 0; i < n; i++) sorted_vars[i] = vars[i];
    std::sort(sorted_vars.begin(), sorted_vars.end());
    double prior_var = sorted_vars[n / 2];
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        shrunk[i] = shrinkage * prior_var + (1.0 - shrinkage) * vars[i];
    }
    
    return shrunk;
}
