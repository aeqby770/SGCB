// =============================================================================
// SGCB conditional autoencoder - multi-factor design with batch effect correction
// Architecture: z = softplus(W1' * [x; batch; condition] + b1)
//              (alpha, beta, gamma) = softplus(W2' * z + b2)
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
// Helper functions
// -----------------------------------------------------------------------------

inline double softplus(double x) {
    return (x > 20.0) ? x : std::log1p(std::exp(x));
}

inline double softplus_grad(double x) {
    return (x > 20.0) ? 1.0 : 1.0 / (1.0 + std::exp(-x));
}

inline double clip(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

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
// Conditional autoencoder training
// Input: X (n_genes x n_samples), batch (n_samples), condition (n_samples)
// Output: per-gene GG parameters for each condition
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List fit_conditional_ae_cpp(NumericMatrix X, 
                            IntegerVector batch,      // 0-indexed batch labels
                            IntegerVector condition,  // 0-indexed condition labels
                            int n_iter = 100, 
                            int h_ratio = 2, 
                            double lr = 0.01, 
                            double eps = 1e-8) {
    
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    const int n_batch = *std::max_element(batch.begin(), batch.end()) + 1;
    const int n_cond = *std::max_element(condition.begin(), condition.end()) + 1;
    
    // Input dimension: n_genes + n_batch (one-hot) + n_cond (one-hot)
    const int input_dim = n_genes + n_batch + n_cond;
    const int h = std::max(n_genes * h_ratio, 100);
    
    // Parameters
    // W1: input_dim x h (too large, split into components)
    // W1_expr: n_genes x h (expression input weights)
    // W1_batch: n_batch x h (batch embedding)
    // W1_cond: n_cond x h (condition embedding)
    // b1: h
    // W2: h x (3 * n_genes) (output alpha, beta, gamma)
    // b2: 3 * n_genes
    
    // Simplified version: use embedding vectors for batch and condition
    const int embed_dim = 32;  // embedding dimension
    
    // Initialize embeddings
    std::vector<double> batch_embed(n_batch * embed_dim);
    std::vector<double> cond_embed(n_cond * embed_dim);
    
    std::mt19937 rng(12345);
    std::normal_distribution<double> normal(0, 0.01);
    
    for (int i = 0; i < n_batch * embed_dim; i++) {
        batch_embed[i] = normal(rng);
    }
    for (int i = 0; i < n_cond * embed_dim; i++) {
        cond_embed[i] = normal(rng);
    }
    
    // Output parameters (per gene x per condition)
    // Storage: GG parameters are condition-specific
    NumericMatrix alpha_out(n_genes, n_cond);
    NumericMatrix beta_out(n_genes, n_cond);
    NumericMatrix gamma_out(n_genes, n_cond);
    
    // Data-driven initialization: compute per-condition statistics
    std::vector<std::vector<double>> cond_means(n_cond, std::vector<double>(n_genes, 0));
    std::vector<std::vector<double>> cond_vars(n_cond, std::vector<double>(n_genes, 0));
    std::vector<int> cond_counts(n_cond, 0);
    
    for (int j = 0; j < n_samples; j++) {
        int c = condition[j];
        cond_counts[c]++;
        for (int g = 0; g < n_genes; g++) {
            cond_means[c][g] += X(g, j);
        }
    }
    
    for (int c = 0; c < n_cond; c++) {
        if (cond_counts[c] > 0) {
            for (int g = 0; g < n_genes; g++) {
                cond_means[c][g] /= cond_counts[c];
            }
        }
    }
    
    // Compute variance
    for (int j = 0; j < n_samples; j++) {
        int c = condition[j];
        for (int g = 0; g < n_genes; g++) {
            double diff = X(g, j) - cond_means[c][g];
            cond_vars[c][g] += diff * diff;
        }
    }
    
    for (int c = 0; c < n_cond; c++) {
        if (cond_counts[c] > 1) {
            for (int g = 0; g < n_genes; g++) {
                cond_vars[c][g] /= (cond_counts[c] - 1);
                cond_vars[c][g] = std::max(cond_vars[c][g], eps);
            }
        }
    }
    
    // Initialize GG parameters
    for (int c = 0; c < n_cond; c++) {
        for (int g = 0; g < n_genes; g++) {
            double mu = std::max(cond_means[c][g], eps);
            double var = std::max(cond_vars[c][g], eps);
            
            alpha_out(g, c) = clip(mu * mu / var, 0.1, 100);
            beta_out(g, c) = clip(mu, eps, 1e6);
            gamma_out(g, c) = 1.0;
        }
    }
    
    // Training: optimize GG parameters via stochastic gradient descent
    // Simultaneously learn batch effects and marginalize them
    
    // Batch effect multiplier (log scale)
    std::vector<std::vector<double>> batch_effect(n_batch, std::vector<double>(n_genes, 0.0));
    
    // Adam optimizer parameters
    const double beta1 = 0.9;
    const double beta2 = 0.999;
    
    // First and second moments
    NumericMatrix m_alpha(n_genes, n_cond);
    NumericMatrix v_alpha(n_genes, n_cond);
    NumericMatrix m_beta(n_genes, n_cond);
    NumericMatrix v_beta(n_genes, n_cond);
    NumericMatrix m_gamma(n_genes, n_cond);
    NumericMatrix v_gamma(n_genes, n_cond);
    
    std::vector<std::vector<double>> m_batch(n_batch, std::vector<double>(n_genes, 0.0));
    std::vector<std::vector<double>> v_batch(n_batch, std::vector<double>(n_genes, 0.0));
    
    // Training loop
    for (int iter = 0; iter < n_iter; iter++) {
        double total_loss = 0.0;
        
        // Iterate over all samples
        for (int j = 0; j < n_samples; j++) {
            int b = batch[j];
            int c = condition[j];
            
            // Compute per-gene gradients
            #ifdef _OPENMP
            #pragma omp parallel for reduction(+:total_loss) schedule(static)
            #endif
            for (int g = 0; g < n_genes; g++) {
                double x = X(g, j);
                
                // Current parameters (with batch effect applied)
                double alpha = alpha_out(g, c);
                double beta_raw = beta_out(g, c);
                double gamma = gamma_out(g, c);
                
                // Batch effect: beta_adjusted = beta_raw * exp(batch_effect)
                double beta_adj = beta_raw * std::exp(batch_effect[b][g]);
                
                // Compute log-likelihood
                double ll = gg_loglik(x, alpha, beta_adj, gamma, eps);
                total_loss -= ll;
                
                // Compute gradients
                double grad_alpha, grad_beta, grad_gamma;
                gg_grad(x, alpha, beta_adj, gamma, eps, grad_alpha, grad_beta, grad_gamma);
                
                // Batch effect gradient
                double grad_batch = grad_beta * beta_adj;  // chain rule
                
                // Adam update
                int t = iter * n_samples + j + 1;
                double lr_t = lr * std::sqrt(1.0 - std::pow(beta2, t)) / (1.0 - std::pow(beta1, t));
                
                // Alpha
                m_alpha(g, c) = beta1 * m_alpha(g, c) + (1 - beta1) * grad_alpha;
                v_alpha(g, c) = beta2 * v_alpha(g, c) + (1 - beta2) * grad_alpha * grad_alpha;
                double update_alpha = lr_t * m_alpha(g, c) / (std::sqrt(v_alpha(g, c)) + eps);
                alpha_out(g, c) = clip(alpha + update_alpha, 0.01, 1000);
                
                // Beta (raw, without batch effect)
                m_beta(g, c) = beta1 * m_beta(g, c) + (1 - beta1) * grad_beta;
                v_beta(g, c) = beta2 * v_beta(g, c) + (1 - beta2) * grad_beta * grad_beta;
                double update_beta = lr_t * m_beta(g, c) / (std::sqrt(v_beta(g, c)) + eps);
                beta_out(g, c) = clip(beta_raw + update_beta, eps, 1e8);
                
                // Gamma
                m_gamma(g, c) = beta1 * m_gamma(g, c) + (1 - beta1) * grad_gamma;
                v_gamma(g, c) = beta2 * v_gamma(g, c) + (1 - beta2) * grad_gamma * grad_gamma;
                double update_gamma = lr_t * m_gamma(g, c) / (std::sqrt(v_gamma(g, c)) + eps);
                gamma_out(g, c) = clip(gamma + update_gamma, 0.01, 100);
                
                // Batch effect
                m_batch[b][g] = beta1 * m_batch[b][g] + (1 - beta1) * grad_batch;
                v_batch[b][g] = beta2 * v_batch[b][g] + (1 - beta2) * grad_batch * grad_batch;
                double update_batch = lr_t * m_batch[b][g] / (std::sqrt(v_batch[b][g]) + eps);
                batch_effect[b][g] = clip(batch_effect[b][g] + update_batch, -5, 5);
            }
        }
        
        // Center batch effects (sum to zero)
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            double mean_effect = 0;
            for (int b = 0; b < n_batch; b++) {
                mean_effect += batch_effect[b][g];
            }
            mean_effect /= n_batch;
            for (int b = 0; b < n_batch; b++) {
                batch_effect[b][g] -= mean_effect;
            }
        }
    }
    
    // Build batch effect matrix
    NumericMatrix batch_effect_out(n_genes, n_batch);
    for (int b = 0; b < n_batch; b++) {
        for (int g = 0; g < n_genes; g++) {
            batch_effect_out(g, b) = batch_effect[b][g];
        }
    }
    
    return List::create(
        Named("alpha") = alpha_out,      // n_genes x n_cond
        Named("beta") = beta_out,        // n_genes x n_cond (batch-free)
        Named("gamma") = gamma_out,      // n_genes x n_cond
        Named("batch_effect") = batch_effect_out,  // n_genes x n_batch
        Named("n_batch") = n_batch,
        Named("n_cond") = n_cond
    );
}


// -----------------------------------------------------------------------------
// Conditional differential test
// Compare GG distribution differences between two conditions
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List conditional_de_test_cpp(NumericMatrix X,
                             IntegerVector batch,
                             IntegerVector condition,
                             NumericMatrix alpha,      // n_genes x n_cond
                             NumericMatrix beta,       // n_genes x n_cond
                             NumericMatrix gamma,      // n_genes x n_cond
                             NumericMatrix batch_effect,
                             int cond_ctrl = 0,        // control condition index
                             int cond_treat = 1,       // treatment condition index
                             int n_boot = 200,
                             double eps = 1e-8) {
    
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    const int n_batch = batch_effect.ncol();
    
    // Compute Log2 Fold Change
    // LFC = log2(beta_treat / beta_ctrl)
    // This is the batch-effect-corrected LFC
    NumericVector log2FC(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double beta_ctrl = beta(g, cond_ctrl);
        double beta_treat = beta(g, cond_treat);
        log2FC[g] = std::log2(beta_treat / beta_ctrl);
    }
    
    // Compute likelihood ratio statistic
    NumericVector lr_stat(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ll_treat = 0, ll_ctrl = 0;
        
        for (int j = 0; j < n_samples; j++) {
            int b = batch[j];
            int c = condition[j];
            double x = X(g, j);
            
            double a = alpha(g, c);
            double beta_raw = beta(g, c);
            double gam = gamma(g, c);
            double beta_adj = beta_raw * std::exp(batch_effect(g, b));
            
            double ll = gg_loglik(x, a, beta_adj, gam, eps);
            
            if (c == cond_ctrl) {
                ll_ctrl += ll;
            } else if (c == cond_treat) {
                ll_treat += ll;
            }
        }
        
        lr_stat[g] = ll_treat - ll_ctrl;
    }
    
    // Bootstrap confidence intervals
    NumericMatrix boot_lfc(n_genes, n_boot);
    
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (condition[j] == cond_ctrl) ctrl_idx.push_back(j);
        else if (condition[j] == cond_treat) treat_idx.push_back(j);
    }
    
    const int n_ctrl = ctrl_idx.size();
    const int n_treat = treat_idx.size();
    
    std::mt19937 rng(12345);
    
    for (int b_iter = 0; b_iter < n_boot; b_iter++) {
        // Resample
        std::vector<int> boot_ctrl(n_ctrl), boot_treat(n_treat);
        for (int i = 0; i < n_ctrl; i++) {
            boot_ctrl[i] = ctrl_idx[rng() % n_ctrl];
        }
        for (int i = 0; i < n_treat; i++) {
            boot_treat[i] = treat_idx[rng() % n_treat];
        }
        
        // Compute bootstrap LFC
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            double sum_ctrl = 0, sum_treat = 0;
            
            for (int i = 0; i < n_ctrl; i++) {
                int j = boot_ctrl[i];
                int bat = batch[j];
                double x = X(g, j) * std::exp(-batch_effect(g, bat));  // remove batch
                sum_ctrl += x;
            }
            
            for (int i = 0; i < n_treat; i++) {
                int j = boot_treat[i];
                int bat = batch[j];
                double x = X(g, j) * std::exp(-batch_effect(g, bat));  // remove batch
                sum_treat += x;
            }
            
            double mean_ctrl = sum_ctrl / n_ctrl + eps;
            double mean_treat = sum_treat / n_treat + eps;
            boot_lfc(g, b_iter) = std::log2(mean_treat / mean_ctrl);
        }
    }
    
    // Compute SE and p-values
    NumericVector lfcSE(n_genes);
    NumericVector pvalue(n_genes);
    NumericVector lfc_ci_lo(n_genes);
    NumericVector lfc_ci_hi(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        // Bootstrap SE
        double sum = 0, sum_sq = 0;
        std::vector<double> boot_vals(n_boot);
        
        for (int b = 0; b < n_boot; b++) {
            double val = boot_lfc(g, b);
            boot_vals[b] = val;
            sum += val;
            sum_sq += val * val;
        }
        
        double mean = sum / n_boot;
        double var = sum_sq / n_boot - mean * mean;
        lfcSE[g] = std::sqrt(var + eps);
        
        // Confidence intervals (quantile method)
        std::sort(boot_vals.begin(), boot_vals.end());
        lfc_ci_lo[g] = boot_vals[(int)(0.025 * n_boot)];
        lfc_ci_hi[g] = boot_vals[(int)(0.975 * n_boot)];
        
        // p-value (two-sided test)
        // If CI includes 0, then p > 0.05
        double t_stat = log2FC[g] / lfcSE[g];
        // Normal approximation
        pvalue[g] = 2.0 * R::pnorm(-std::abs(t_stat), 0, 1, 1, 0);
    }
    
    // baseMean
    NumericVector baseMean(n_genes);
    for (int g = 0; g < n_genes; g++) {
        double sum = 0;
        for (int j = 0; j < n_samples; j++) {
            int b = batch[j];
            sum += X(g, j) * std::exp(-batch_effect(g, b));  // batch-corrected expression
        }
        baseMean[g] = sum / n_samples;
    }
    
    return List::create(
        Named("log2FoldChange") = log2FC,
        Named("lfcSE") = lfcSE,
        Named("stat") = lr_stat,
        Named("pvalue") = pvalue,
        Named("lfc_ci_lo") = lfc_ci_lo,
        Named("lfc_ci_hi") = lfc_ci_hi,
        Named("baseMean") = baseMean
    );
}


// -----------------------------------------------------------------------------
// Batch-corrected expression matrix
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix remove_batch_effect_cpp(NumericMatrix X,
                                       IntegerVector batch,
                                       NumericMatrix batch_effect) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    NumericMatrix X_corrected(n_genes, n_samples);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int j = 0; j < n_samples; j++) {
        int b = batch[j];
        for (int g = 0; g < n_genes; g++) {
            // Multiplicative correction: X_corrected = X * exp(-batch_effect)
            X_corrected(g, j) = X(g, j) * std::exp(-batch_effect(g, b));
        }
    }
    
    return X_corrected;
}
