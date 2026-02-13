// =============================================================================
// SGCB parametric/smoothed bootstrap
// Addresses p-value discreteness for small samples (n<5)
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// -----------------------------------------------------------------------------
// GG distribution sampling (inverse transform method)
// -----------------------------------------------------------------------------

inline double sample_gg(double alpha, double beta, double gamma, std::mt19937& rng) {
    // GG(alpha, beta, gamma) = (Gamma(alpha, 1))^(1/gamma) * beta
    std::gamma_distribution<double> gamma_dist(alpha, 1.0);
    double g = gamma_dist(rng);
    return beta * std::pow(g, 1.0 / gamma);
}

// -----------------------------------------------------------------------------
// Parametric bootstrap
// Generate new samples from the estimated GG distribution
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix parametric_bootstrap_cpp(NumericVector alpha,   // n_genes
                                        NumericVector beta,    // n_genes
                                        NumericVector gamma,   // n_genes
                                        int n_samples,
                                        int n_boot,
                                        int seed = 12345) {
    const int n_genes = alpha.size();
    
    // Output: per-gene mean for each bootstrap replicate
    NumericMatrix boot_means(n_genes, n_boot);
    
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        #else
        int thread_id = 0;
        #endif
        
        std::mt19937 rng(seed + thread_id);
        
        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int b = 0; b < n_boot; b++) {
            for (int g = 0; g < n_genes; g++) {
                double sum = 0.0;
                for (int j = 0; j < n_samples; j++) {
                    sum += sample_gg(alpha[g], beta[g], gamma[g], rng);
                }
                boot_means(g, b) = sum / n_samples;
            }
        }
    }
    
    return boot_means;
}


// -----------------------------------------------------------------------------
// Smoothed bootstrap
// Add small noise to resampled data
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix smoothed_bootstrap_cpp(NumericMatrix X,         // n_genes x n_samples
                                      int n_boot,
                                      double jitter_scale = 0.1,  // noise scale
                                      int seed = 12345) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    // Compute per-gene standard deviation
    std::vector<double> gene_sd(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum = 0, sum_sq = 0;
        for (int j = 0; j < n_samples; j++) {
            double val = X(g, j);
            sum += val;
            sum_sq += val * val;
        }
        double mean = sum / n_samples;
        double var = sum_sq / n_samples - mean * mean;
        gene_sd[g] = std::sqrt(std::max(var, 1e-8));
    }
    
    // Output: per-gene mean for each bootstrap replicate
    NumericMatrix boot_means(n_genes, n_boot);
    
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        #else
        int thread_id = 0;
        #endif
        
        std::mt19937 rng(seed + thread_id);
        std::uniform_int_distribution<int> sample_dist(0, n_samples - 1);
        
        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int b = 0; b < n_boot; b++) {
            for (int g = 0; g < n_genes; g++) {
                double sum = 0.0;
                double jitter_sd = jitter_scale * gene_sd[g];
                std::normal_distribution<double> noise_dist(0, jitter_sd);
                
                for (int j = 0; j < n_samples; j++) {
                    int idx = sample_dist(rng);
                    double val = X(g, idx) + noise_dist(rng);
                    sum += std::max(val, 0.0);  // keep non-negative
                }
                boot_means(g, b) = sum / n_samples;
            }
        }
    }
    
    return boot_means;
}


// -----------------------------------------------------------------------------
// Adaptive bootstrap selection
// n < 5: use smoothed/parametric; otherwise use standard bootstrap
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List adaptive_bootstrap_lfc_cpp(NumericMatrix X,
                                 IntegerVector group,  // 0 = ctrl, 1 = treat
                                 NumericVector alpha,  // GG parameters (ctrl estimate)
                                 NumericVector beta,
                                 NumericVector gamma,
                                 int n_boot = 200,
                                 double jitter_scale = 0.1,
                                 int seed = 12345) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    // Separate two groups
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    
    const int n_ctrl = ctrl_idx.size();
    const int n_treat = treat_idx.size();
    const int n_min = std::min(n_ctrl, n_treat);
    
    // Decide which bootstrap to use
    // n < 5: use smoothed bootstrap
    // n >= 5: use standard bootstrap
    const bool use_smoothed = (n_min < 5);
    
    // Compute per-group gene standard deviation (for smoothed bootstrap)
    std::vector<double> ctrl_sd(n_genes), treat_sd(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum_c = 0, sum_sq_c = 0;
        for (int idx : ctrl_idx) {
            double val = X(g, idx);
            sum_c += val;
            sum_sq_c += val * val;
        }
        double mean_c = sum_c / n_ctrl;
        double var_c = sum_sq_c / n_ctrl - mean_c * mean_c;
        ctrl_sd[g] = std::sqrt(std::max(var_c, 1e-8));
        
        double sum_t = 0, sum_sq_t = 0;
        for (int idx : treat_idx) {
            double val = X(g, idx);
            sum_t += val;
            sum_sq_t += val * val;
        }
        double mean_t = sum_t / n_treat;
        double var_t = sum_sq_t / n_treat - mean_t * mean_t;
        treat_sd[g] = std::sqrt(std::max(var_t, 1e-8));
    }
    
    // Bootstrap LFC
    NumericMatrix boot_lfc(n_genes, n_boot);
    const double eps = 1e-8;
    
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        #else
        int thread_id = 0;
        #endif
        
        std::mt19937 rng(seed + thread_id);
        std::uniform_int_distribution<int> ctrl_dist(0, n_ctrl - 1);
        std::uniform_int_distribution<int> treat_dist(0, n_treat - 1);
        
        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int b = 0; b < n_boot; b++) {
            for (int g = 0; g < n_genes; g++) {
                double sum_ctrl = 0, sum_treat = 0;
                
                if (use_smoothed) {
                    // Smoothed bootstrap: resampling + noise
                    std::normal_distribution<double> ctrl_noise(0, jitter_scale * ctrl_sd[g]);
                    std::normal_distribution<double> treat_noise(0, jitter_scale * treat_sd[g]);
                    
                    for (int j = 0; j < n_ctrl; j++) {
                        int idx = ctrl_idx[ctrl_dist(rng)];
                        double val = X(g, idx) + ctrl_noise(rng);
                        sum_ctrl += std::max(val, eps);
                    }
                    for (int j = 0; j < n_treat; j++) {
                        int idx = treat_idx[treat_dist(rng)];
                        double val = X(g, idx) + treat_noise(rng);
                        sum_treat += std::max(val, eps);
                    }
                } else {
                    // Standard bootstrap: pure resampling
                    for (int j = 0; j < n_ctrl; j++) {
                        int idx = ctrl_idx[ctrl_dist(rng)];
                        sum_ctrl += X(g, idx);
                    }
                    for (int j = 0; j < n_treat; j++) {
                        int idx = treat_idx[treat_dist(rng)];
                        sum_treat += X(g, idx);
                    }
                }
                
                double mean_ctrl = sum_ctrl / n_ctrl + eps;
                double mean_treat = sum_treat / n_treat + eps;
                boot_lfc(g, b) = std::log2(mean_treat / mean_ctrl);
            }
        }
    }
    
    // Compute statistics
    NumericVector lfc_mean(n_genes);
    NumericVector lfc_se(n_genes);
    NumericVector lfc_ci_lo(n_genes);
    NumericVector lfc_ci_hi(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        std::vector<double> vals(n_boot);
        double sum = 0;
        
        for (int b = 0; b < n_boot; b++) {
            vals[b] = boot_lfc(g, b);
            sum += vals[b];
        }
        
        lfc_mean[g] = sum / n_boot;
        
        double sum_sq = 0;
        for (int b = 0; b < n_boot; b++) {
            double diff = vals[b] - lfc_mean[g];
            sum_sq += diff * diff;
        }
        lfc_se[g] = std::sqrt(sum_sq / (n_boot - 1));
        
        // Quantiles
        std::sort(vals.begin(), vals.end());
        lfc_ci_lo[g] = vals[(int)(0.025 * n_boot)];
        lfc_ci_hi[g] = vals[(int)(0.975 * n_boot)];
    }
    
    return List::create(
        Named("boot_lfc") = boot_lfc,
        Named("lfc_mean") = lfc_mean,
        Named("lfc_se") = lfc_se,
        Named("lfc_ci_lo") = lfc_ci_lo,
        Named("lfc_ci_hi") = lfc_ci_hi,
        Named("use_smoothed") = use_smoothed,
        Named("n_ctrl") = n_ctrl,
        Named("n_treat") = n_treat
    );
}


// -----------------------------------------------------------------------------
// Parametric bootstrap differential test
// Generate virtual samples from GG distribution for testing
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List parametric_de_test_cpp(NumericVector alpha_ctrl,
                             NumericVector beta_ctrl,
                             NumericVector gamma_ctrl,
                             NumericVector alpha_treat,
                             NumericVector beta_treat,
                             NumericVector gamma_treat,
                             int n_ctrl,
                             int n_treat,
                             int n_boot = 500,
                             int seed = 12345) {
    const int n_genes = alpha_ctrl.size();
    const double eps = 1e-8;
    
    // Observed LFC
    NumericVector obs_lfc(n_genes);
    for (int g = 0; g < n_genes; g++) {
        obs_lfc[g] = std::log2(beta_treat[g] / beta_ctrl[g]);
    }
    
    // Null distribution: assume identical parameters in both groups (use ctrl parameters)
    NumericMatrix null_lfc(n_genes, n_boot);
    
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        #else
        int thread_id = 0;
        #endif
        
        std::mt19937 rng(seed + thread_id);
        
        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int b = 0; b < n_boot; b++) {
            for (int g = 0; g < n_genes; g++) {
                // Generate two groups of samples from the same distribution
                double sum_ctrl = 0, sum_treat = 0;
                
                for (int j = 0; j < n_ctrl; j++) {
                    sum_ctrl += sample_gg(alpha_ctrl[g], beta_ctrl[g], gamma_ctrl[g], rng);
                }
                for (int j = 0; j < n_treat; j++) {
                    sum_treat += sample_gg(alpha_ctrl[g], beta_ctrl[g], gamma_ctrl[g], rng);
                }
                
                double mean_ctrl = sum_ctrl / n_ctrl + eps;
                double mean_treat = sum_treat / n_treat + eps;
                null_lfc(g, b) = std::log2(mean_treat / mean_ctrl);
            }
        }
    }
    
    // Compute p-values
    NumericVector pvalue(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double obs = std::abs(obs_lfc[g]);
        int count = 0;
        for (int b = 0; b < n_boot; b++) {
            if (std::abs(null_lfc(g, b)) >= obs) {
                count++;
            }
        }
        pvalue[g] = (double)(count + 1) / (n_boot + 1);  // +1 to prevent p=0
    }
    
    return List::create(
        Named("log2FC") = obs_lfc,
        Named("pvalue") = pvalue,
        Named("null_lfc") = null_lfc
    );
}
