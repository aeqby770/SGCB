// =============================================================================
// SGCB wide-shallow autoencoder - Rcpp-accelerated
// Architecture Section 3: z = softplus(W1' * x + b1), theta = W2' * z + b2
// =============================================================================

#include <Rcpp.h>
#include "fast_special.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// -----------------------------------------------------------------------------
// Softplus activation function
// -----------------------------------------------------------------------------

inline double softplus(double x) {
    return (x > 20.0) ? x : std::log1p(std::exp(x));
}

inline double softplus_grad(double x) {
    return (x > 20.0) ? 1.0 : 1.0 / (1.0 + std::exp(-x));
}

// -----------------------------------------------------------------------------
// Clipping function
// -----------------------------------------------------------------------------

inline double clip(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

// -----------------------------------------------------------------------------
// GG log-likelihood and gradient (single observation)
// -----------------------------------------------------------------------------

inline double gg_loglik(double x, double alpha, double beta, double gamma, double eps) {
    x = std::max(x, eps);
    return std::log(gamma) - alpha * std::log(beta) - fast_special::lgamma(alpha) +
           (gamma * alpha - 1.0) * std::log(x) - std::pow(x / beta, gamma);
}

inline void gg_grad(double x, double alpha, double beta, double gamma, double eps,
                    double& grad_alpha, double& grad_beta, double& grad_gamma) {
    x = std::max(x, eps);
    double log_x = std::log(x);
    double ratio = x / beta;
    double ratio_g = std::pow(ratio, gamma);
    double log_ratio = std::log(ratio);
    
    grad_alpha = -std::log(beta) - fast_special::digamma(alpha) + gamma * log_x;
    grad_beta = -alpha / beta + gamma * ratio_g / beta;
    grad_gamma = 1.0 / gamma + alpha * log_x - ratio_g * log_ratio;
}

// -----------------------------------------------------------------------------
// Wide-shallow autoencoder training (Rcpp-accelerated)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List fit_wide_shallow_ae_cpp(NumericMatrix X, int n_iter = 100, int h_ratio = 2, 
                              double lr = 0.01, double eps = 1e-8) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    int h = std::max(n_genes * h_ratio, 100);
    
    // Parameter initialization
    NumericVector b2(3 * n_genes);
    
    // Data-driven initialization
    NumericVector ctrl_mean(n_genes), ctrl_var(n_genes);
    
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
        ctrl_mean[g] = sum / n_samples + eps;
        ctrl_var[g] = sum_sq / n_samples - ctrl_mean[g] * ctrl_mean[g] + eps;
    }
    
    // Initialize GG parameters (log space)
    for (int g = 0; g < n_genes; g++) {
        double init_alpha = clip(ctrl_mean[g] * ctrl_mean[g] / ctrl_var[g], 0.1, 100);
        double init_beta = clip(ctrl_mean[g], eps, 1e6);
        double init_gamma = 1.0;
        
        b2[3 * g] = std::log(init_alpha);
        b2[3 * g + 1] = std::log(init_beta);
        b2[3 * g + 2] = std::log(init_gamma);
    }
    
    double log_alpha_lo = std::log(0.1), log_alpha_hi = std::log(100);
    double log_beta_lo = std::log(eps), log_beta_hi = std::log(1e6);
    double log_gamma_lo = std::log(0.5), log_gamma_hi = std::log(5);
    
    // Training loop
    NumericVector loss_history(n_iter);
    
    for (int iter = 0; iter < n_iter; iter++) {
        double total_loss = 0;
        
        // Gradient accumulator
        NumericVector grad_b2(3 * n_genes);
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:total_loss) schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            // Current parameters
            double log_alpha = clip(b2[3 * g], log_alpha_lo, log_alpha_hi);
            double log_beta = clip(b2[3 * g + 1], log_beta_lo, log_beta_hi);
            double log_gamma = clip(b2[3 * g + 2], log_gamma_lo, log_gamma_hi);
            
            double alpha = std::exp(log_alpha);
            double beta = std::exp(log_beta);
            double gamma = std::exp(log_gamma);
            
            double sum_grad_alpha = 0, sum_grad_beta = 0, sum_grad_gamma = 0;
            double sum_loss = 0;
            
            for (int j = 0; j < n_samples; j++) {
                double x = X(g, j) + eps;
                
                // Log-likelihood
                sum_loss += gg_loglik(x, alpha, beta, gamma, eps);
                
                // Gradient
                double ga, gb, gg;
                gg_grad(x, alpha, beta, gamma, eps, ga, gb, gg);
                
                sum_grad_alpha += ga * alpha;  // chain rule
                sum_grad_beta += gb * beta;
                sum_grad_gamma += gg * gamma;
            }
            
            total_loss -= sum_loss;  // negative log-likelihood
            
            // Average gradient
            grad_b2[3 * g] = sum_grad_alpha / n_samples;
            grad_b2[3 * g + 1] = sum_grad_beta / n_samples;
            grad_b2[3 * g + 2] = sum_grad_gamma / n_samples;
        }
        
        loss_history[iter] = total_loss / n_genes;
        
        // Gradient ascent update
        for (int i = 0; i < 3 * n_genes; i++) {
            b2[i] += lr * grad_b2[i];
        }
    }
    
    // Extract final parameters
    NumericVector alpha_out(n_genes), beta_out(n_genes), gamma_out(n_genes);
    
    for (int g = 0; g < n_genes; g++) {
        alpha_out[g] = std::exp(clip(b2[3 * g], log_alpha_lo, log_alpha_hi));
        beta_out[g] = std::exp(clip(b2[3 * g + 1], log_beta_lo, log_beta_hi));
        gamma_out[g] = std::exp(clip(b2[3 * g + 2], log_gamma_lo, log_gamma_hi));
    }
    
    return List::create(
        Named("alpha") = alpha_out,
        Named("beta") = beta_out,
        Named("gamma") = gamma_out,
        Named("loss_history") = loss_history
    );
}

// -----------------------------------------------------------------------------
// SDE dropout stability computation (Rcpp-accelerated)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List sde_dropout_stability_cpp(NumericMatrix X, IntegerVector group, 
                                int n_dropout = 50, double dropout_rate = 0.2, 
                                int seed = 12345) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    double eps = 1e-8;
    
    // Group indices
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    int n_ctrl = ctrl_idx.size();
    int n_treat = treat_idx.size();
    
    // LFC matrix
    NumericMatrix lfc_mat(n_genes, n_dropout);
    
    // Dropout loop
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    for (int d = 0; d < n_dropout; d++) {
        // Dropout sampling
        std::vector<int> keep_ctrl, keep_treat;
        
        for (int i = 0; i < n_ctrl; i++) {
            if (unif(rng) > dropout_rate) {
                keep_ctrl.push_back(ctrl_idx[i]);
            }
        }
        for (int i = 0; i < n_treat; i++) {
            if (unif(rng) > dropout_rate) {
                keep_treat.push_back(treat_idx[i]);
            }
        }
        
        // Ensure at least 2 samples per group
        while (keep_ctrl.size() < 2 && keep_ctrl.size() < ctrl_idx.size()) {
            int idx = (int)(unif(rng) * ctrl_idx.size());
            if (std::find(keep_ctrl.begin(), keep_ctrl.end(), ctrl_idx[idx]) == keep_ctrl.end()) {
                keep_ctrl.push_back(ctrl_idx[idx]);
            }
        }
        while (keep_treat.size() < 2 && keep_treat.size() < treat_idx.size()) {
            int idx = (int)(unif(rng) * treat_idx.size());
            if (std::find(keep_treat.begin(), keep_treat.end(), treat_idx[idx]) == keep_treat.end()) {
                keep_treat.push_back(treat_idx[idx]);
            }
        }
        
        int kc = keep_ctrl.size();
        int kt = keep_treat.size();
        
        // Compute LFC
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            double sum_ctrl = 0, sum_treat = 0;
            
            for (int i = 0; i < kc; i++) {
                sum_ctrl += X(g, keep_ctrl[i]);
            }
            for (int i = 0; i < kt; i++) {
                sum_treat += X(g, keep_treat[i]);
            }
            
            double mean_ctrl = sum_ctrl / kc + eps;
            double mean_treat = sum_treat / kt + eps;
            
            lfc_mat(g, d) = std::log2(mean_treat / mean_ctrl);
        }
    }
    
    // Compute summary statistics
    NumericVector lfc_mean(n_genes), lfc_var(n_genes), stability(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum = 0, sum_sq = 0;
        for (int d = 0; d < n_dropout; d++) {
            double val = lfc_mat(g, d);
            sum += val;
            sum_sq += val * val;
        }
        lfc_mean[g] = sum / n_dropout;
        lfc_var[g] = (sum_sq - sum * sum / n_dropout) / (n_dropout - 1);
        stability[g] = 1.0 / (1.0 + lfc_var[g]);
    }
    
    return List::create(
        Named("lfc_stable") = lfc_mean,
        Named("lfc_var") = lfc_var,
        Named("stability") = stability,
        Named("lfc_matrix") = lfc_mat
    );
}

// -----------------------------------------------------------------------------
// Variance shrinkage (empirical Bayes) - Rcpp-accelerated
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List squeeze_var_cpp(NumericVector var, int df, double eps = 1e-15) {
    int n = var.size();
    
    // Clip extreme values
    NumericVector var_copy = clone(var);
    for (int i = 0; i < n; i++) {
        var_copy[i] = std::max(var_copy[i], eps);
    }
    
    // Sort to obtain quantiles
    NumericVector sorted = clone(var_copy);
    std::sort(sorted.begin(), sorted.end());
    
    int lo_idx = (int)(0.01 * n);
    int hi_idx = (int)(0.99 * n);
    double q01 = sorted[lo_idx];
    double q99 = sorted[hi_idx];
    
    // Compute median and MAD
    std::vector<double> trimmed;
    for (int i = 0; i < n; i++) {
        if (var_copy[i] > q01 && var_copy[i] < q99) {
            trimmed.push_back(var_copy[i]);
        }
    }
    
    std::sort(trimmed.begin(), trimmed.end());
    double m = trimmed[trimmed.size() / 2];
    
    // MAD
    std::vector<double> abs_dev(trimmed.size());
    for (size_t i = 0; i < trimmed.size(); i++) {
        abs_dev[i] = std::abs(trimmed[i] - m);
    }
    std::sort(abs_dev.begin(), abs_dev.end());
    double mad = abs_dev[abs_dev.size() / 2];
    double v = mad * mad;
    
    // Prior parameters
    double d0 = std::max(1.0, std::min(100.0, 2 * m * m / (v + 1e-10)));
    double s0_sq = m;
    double df_total = d0 + df;
    
    // Posterior variance
    NumericVector var_post(n);
    for (int i = 0; i < n; i++) {
        var_post[i] = (d0 * s0_sq + df * var_copy[i]) / df_total;
    }
    
    return List::create(
        Named("var_post") = var_post,
        Named("var_prior") = s0_sq,
        Named("df_prior") = d0,
        Named("df_total") = df_total
    );
}

// -----------------------------------------------------------------------------
// Complete SGCB inference pipeline (Rcpp-accelerated)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List sgcb_inference_cpp(NumericMatrix X, IntegerVector group,
                         int n_boot = 200, int n_dropout = 50,
                         double dropout_rate = 0.2, int ae_iter = 50,
                         double ae_lr = 0.01, int seed = 12345) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    double eps = 1e-8;
    
    // Group assignment
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    int n_ctrl = ctrl_idx.size();
    int n_treat = treat_idx.size();
    
    // Extract control group matrix
    NumericMatrix ctrl_mat(n_genes, n_ctrl);
    for (int g = 0; g < n_genes; g++) {
        for (int j = 0; j < n_ctrl; j++) {
            ctrl_mat(g, j) = X(g, ctrl_idx[j]);
        }
    }
    
    // Step 1: autoencoder fitting
    List ae_result = fit_wide_shallow_ae_cpp(ctrl_mat, ae_iter, 2, ae_lr, eps);
    NumericVector ae_alpha = ae_result["alpha"];
    NumericVector ae_beta = ae_result["beta"];
    NumericVector ae_gamma = ae_result["gamma"];
    
    // Step 2: SDE dropout stability
    List sde_result = sde_dropout_stability_cpp(X, group, n_dropout, dropout_rate, seed);
    NumericVector lfc_stable = sde_result["lfc_stable"];
    NumericVector lfc_var = sde_result["lfc_var"];
    NumericVector stability = sde_result["stability"];
    
    // Step 3: basic statistics
    NumericVector mean_ctrl(n_genes), mean_treat(n_genes);
    NumericVector log2FC_raw(n_genes), baseMean(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sc = 0, st = 0;
        for (int j = 0; j < n_ctrl; j++) sc += X(g, ctrl_idx[j]);
        for (int j = 0; j < n_treat; j++) st += X(g, treat_idx[j]);
        
        mean_ctrl[g] = sc / n_ctrl + eps;
        mean_treat[g] = st / n_treat + eps;
        log2FC_raw[g] = std::log2(mean_treat[g] / mean_ctrl[g]);
        baseMean[g] = (mean_ctrl[g] + mean_treat[g]) / 2;
    }
    
    return List::create(
        Named("ae_alpha") = ae_alpha,
        Named("ae_beta") = ae_beta,
        Named("ae_gamma") = ae_gamma,
        Named("lfc_stable") = lfc_stable,
        Named("lfc_var") = lfc_var,
        Named("stability") = stability,
        Named("log2FC_raw") = log2FC_raw,
        Named("baseMean") = baseMean,
        Named("mean_ctrl") = mean_ctrl,
        Named("mean_treat") = mean_treat
    );
}
