// =============================================================================
// SGCB Calibrated Bootstrap - Rcpp-accelerated
// Implementation of Architecture Section 4: m-out-of-n bootstrap + GG likelihood ratio
// =============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// -----------------------------------------------------------------------------
// GG log-likelihood (single gene)
// -----------------------------------------------------------------------------

inline double gg_loglik_single(double x, double alpha, double beta, double gamma, double eps = 1e-8) {
    x = std::max(x, eps);
    return std::log(gamma) - gamma * alpha * std::log(beta) - std::lgamma(alpha) +
           (gamma * alpha - 1.0) * std::log(x) - std::pow(x / beta, gamma);
}

// -----------------------------------------------------------------------------
// GG parameter method-of-moments estimation (single gene)
// -----------------------------------------------------------------------------

inline void estimate_gg_params(const std::vector<double>& x, double& alpha, double& beta, double& gamma, double eps = 1e-8) {
    int n = x.size();
    double m1 = 0, m2 = 0, m3 = 0;
    
    for (int i = 0; i < n; i++) {
        double xi = std::max(x[i], eps);
        m1 += xi;
        m2 += xi * xi;
        m3 += xi * xi * xi;
    }
    m1 /= n;
    m2 /= n;
    m3 /= n;
    
    // Variance and standard deviation
    double var = m2 - m1 * m1;
    double sd = std::sqrt(std::max(var, eps));
    double cv = sd / std::max(m1, eps);
    
    // Skewness
    double skew = (m3 - 3 * m1 * m2 + 2 * m1 * m1 * m1) / std::max(sd * sd * sd, eps);
    
    // gamma estimate (based on skewness)
    gamma = std::max(0.5, std::min(5.0, 2.0 / (std::abs(skew) + 0.1)));
    
    // alpha estimate (based on CV)
    alpha = std::max(0.1, std::min(100.0, 1.0 / (cv * cv + eps)));
    
    // beta estimate
    double gamma_ratio = std::tgamma(alpha + 1.0 / gamma) / std::max(std::tgamma(alpha), eps);
    beta = std::max(eps, m1 / std::max(gamma_ratio, eps));
}

// -----------------------------------------------------------------------------
// GG likelihood ratio statistic (single gene)
// -----------------------------------------------------------------------------

inline double gg_likelihood_ratio(const std::vector<double>& x_ctrl, 
                                   const std::vector<double>& x_treat, double eps = 1e-8) {
    // Estimate per-group parameters
    double a_ctrl, b_ctrl, g_ctrl;
    double a_treat, b_treat, g_treat;
    double a_null, b_null, g_null;
    
    estimate_gg_params(x_ctrl, a_ctrl, b_ctrl, g_ctrl, eps);
    estimate_gg_params(x_treat, a_treat, b_treat, g_treat, eps);
    
    // Estimate null parameters from pooled data
    std::vector<double> x_all;
    x_all.reserve(x_ctrl.size() + x_treat.size());
    x_all.insert(x_all.end(), x_ctrl.begin(), x_ctrl.end());
    x_all.insert(x_all.end(), x_treat.begin(), x_treat.end());
    estimate_gg_params(x_all, a_null, b_null, g_null, eps);
    
    // Compute likelihoods
    double ll_alt = 0, ll_null = 0;
    
    for (size_t i = 0; i < x_ctrl.size(); i++) {
        ll_alt += gg_loglik_single(x_ctrl[i], a_ctrl, b_ctrl, g_ctrl, eps);
        ll_null += gg_loglik_single(x_ctrl[i], a_null, b_null, g_null, eps);
    }
    for (size_t i = 0; i < x_treat.size(); i++) {
        ll_alt += gg_loglik_single(x_treat[i], a_treat, b_treat, g_treat, eps);
        ll_null += gg_loglik_single(x_treat[i], a_null, b_null, g_null, eps);
    }
    
    return std::max(0.0, 2.0 * (ll_alt - ll_null));
}

// -----------------------------------------------------------------------------
// Calibrated bootstrap p-value (single gene) - Architecture Section 4
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List calibrated_bootstrap_gene_cpp(NumericVector x_ctrl, NumericVector x_treat,
                                    int n_boot = 200, double m_ratio = 0.7, int seed = 12345) {
    int n_ctrl = x_ctrl.size();
    int n_treat = x_treat.size();
    double eps = 1e-8;
    
    // m-out-of-n
    int m_ctrl = std::max(2, (int)(n_ctrl * m_ratio));
    int m_treat = std::max(2, (int)(n_treat * m_ratio));
    
    // Convert to vector
    std::vector<double> ctrl(x_ctrl.begin(), x_ctrl.end());
    std::vector<double> treat(x_treat.begin(), x_treat.end());
    
    // Observed LR statistic
    double lr_obs = gg_likelihood_ratio(ctrl, treat, eps);
    
    // Pool data for permutation
    std::vector<double> x_all;
    x_all.reserve(n_ctrl + n_treat);
    x_all.insert(x_all.end(), ctrl.begin(), ctrl.end());
    x_all.insert(x_all.end(), treat.begin(), treat.end());
    int n_all = x_all.size();
    
    // Bootstrap
    std::mt19937 rng(seed);
    NumericVector lr_boot(n_boot);
    
    for (int b = 0; b < n_boot; b++) {
        // Permute
        std::shuffle(x_all.begin(), x_all.end(), rng);
        
        // m-out-of-n resampling
        std::vector<double> boot_ctrl(m_ctrl), boot_treat(m_treat);
        std::uniform_int_distribution<int> dist_ctrl(0, n_ctrl - 1);
        std::uniform_int_distribution<int> dist_treat(n_ctrl, n_all - 1);
        
        for (int i = 0; i < m_ctrl; i++) {
            boot_ctrl[i] = x_all[dist_ctrl(rng)];
        }
        for (int i = 0; i < m_treat; i++) {
            boot_treat[i] = x_all[dist_treat(rng)];
        }
        
        lr_boot[b] = gg_likelihood_ratio(boot_ctrl, boot_treat, eps);
    }
    
    // Empirical p-value (continuity correction)
    int count = 0;
    for (int b = 0; b < n_boot; b++) {
        if (lr_boot[b] >= lr_obs) count++;
    }
    double pval = (1.0 + count) / (1.0 + n_boot);
    
    return List::create(
        Named("pvalue") = pval,
        Named("lr_obs") = lr_obs,
        Named("lr_boot") = lr_boot
    );
}

// -----------------------------------------------------------------------------
// Batch calibrated bootstrap (all genes) - OpenMP parallel
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List calibrated_bootstrap_all_cpp(NumericMatrix X, IntegerVector group,
                                   int n_boot = 200, double m_ratio = 0.7, int seed = 12345) {
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
    
    // m-out-of-n
    int m_ctrl = std::max(2, (int)(n_ctrl * m_ratio));
    int m_treat = std::max(2, (int)(n_treat * m_ratio));
    
    // Output
    NumericVector pvalues(n_genes);
    NumericVector lr_stats(n_genes);
    NumericVector gg_alpha(n_genes), gg_beta(n_genes), gg_gamma(n_genes);
    NumericVector negll_stats(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::mt19937 rng(seed + thread_id);
        
        #pragma omp for schedule(dynamic)
        for (int g = 0; g < n_genes; g++) {
    #else
    {
        std::mt19937 rng(seed);
        for (int g = 0; g < n_genes; g++) {
    #endif
            // Extract data
            std::vector<double> x_ctrl(n_ctrl), x_treat(n_treat);
            for (int i = 0; i < n_ctrl; i++) {
                x_ctrl[i] = X(g, ctrl_idx[i]);
            }
            for (int i = 0; i < n_treat; i++) {
                x_treat[i] = X(g, treat_idx[i]);
            }
            
            // Estimate control group GG parameters
            double a, b, gam;
            estimate_gg_params(x_ctrl, a, b, gam, eps);
            gg_alpha[g] = a;
            gg_beta[g] = b;
            gg_gamma[g] = gam;
            
            // Architecture Section 4.3: T_g = negative log-likelihood of treatment under control model
            double negll = 0;
            for (int i = 0; i < n_treat; i++) {
                negll -= gg_loglik_single(x_treat[i], a, b, gam, eps);
            }
            negll_stats[g] = negll / n_treat;
            
            // Observed LR statistic
            double lr_obs = gg_likelihood_ratio(x_ctrl, x_treat, eps);
            lr_stats[g] = lr_obs;
            
            // Bootstrap null distribution
            std::vector<double> x_all;
            x_all.reserve(n_ctrl + n_treat);
            x_all.insert(x_all.end(), x_ctrl.begin(), x_ctrl.end());
            x_all.insert(x_all.end(), x_treat.begin(), x_treat.end());
            int n_all = x_all.size();
            
            int count = 0;
            for (int bb = 0; bb < n_boot; bb++) {
                // Permute
                std::shuffle(x_all.begin(), x_all.end(), rng);
                
                // m-out-of-n resampling
                std::vector<double> boot_ctrl(m_ctrl), boot_treat(m_treat);
                std::uniform_int_distribution<int> dist_ctrl(0, n_ctrl - 1);
                std::uniform_int_distribution<int> dist_treat(n_ctrl, n_all - 1);
                
                for (int i = 0; i < m_ctrl; i++) {
                    boot_ctrl[i] = x_all[dist_ctrl(rng)];
                }
                for (int i = 0; i < m_treat; i++) {
                    boot_treat[i] = x_all[dist_treat(rng)];
                }
                
                double lr_boot = gg_likelihood_ratio(boot_ctrl, boot_treat, eps);
                if (lr_boot >= lr_obs) count++;
            }
            
            // Empirical p-value
            pvalues[g] = (1.0 + count) / (1.0 + n_boot);
        }
    }
    
    return List::create(
        Named("pvalue") = pvalues,
        Named("lr_stat") = lr_stats,
        Named("negll_stat") = negll_stats,
        Named("gg_alpha") = gg_alpha,
        Named("gg_beta") = gg_beta,
        Named("gg_gamma") = gg_gamma
    );
}

// -----------------------------------------------------------------------------
// GG parameter MLE estimation (Newton-Raphson) - single gene
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List gg_mle_single_cpp(NumericVector x, int max_iter = 50, double tol = 1e-6) {
    int n = x.size();
    double eps = 1e-8;
    
    std::vector<double> data(x.begin(), x.end());
    
    // Method-of-moments initial values
    double alpha, beta, gamma;
    estimate_gg_params(data, alpha, beta, gamma, eps);
    
    // Newton-Raphson iterations
    for (int iter = 0; iter < max_iter; iter++) {
        // Compute gradient and Hessian diagonal
        double grad_a = 0, grad_b = 0, grad_g = 0;
        double hess_a = 0, hess_b = 0, hess_g = 0;
        
        for (int i = 0; i < n; i++) {
            double xi = std::max(data[i], eps);
            double log_xi = std::log(xi);
            double ratio = xi / beta;
            double ratio_g = std::pow(ratio, gamma);
            double log_ratio = std::log(ratio);
            
            // Gradient
            grad_a += -gamma * std::log(beta) - R::digamma(alpha) + gamma * log_xi;
            grad_b += -gamma * alpha / beta + gamma * ratio_g / beta;
            grad_g += 1.0 / gamma + alpha * log_ratio - ratio_g * log_ratio;
            
            // Hessian diagonal
            hess_a += -R::trigamma(alpha);
            hess_b += gamma * alpha / (beta * beta) - gamma * (gamma + 1.0) * ratio_g / (beta * beta);
            hess_g += -1.0 / (gamma * gamma) - ratio_g * log_ratio * log_ratio;
        }
        
        // Update (negative Hessian)
        double delta_a = grad_a / std::max(std::abs(hess_a), eps);
        double delta_b = grad_b / std::max(std::abs(hess_b), eps);
        double delta_g = grad_g / std::max(std::abs(hess_g), eps);
        
        alpha = std::max(0.1, std::min(100.0, alpha - delta_a));
        beta = std::max(eps, beta - delta_b);
        gamma = std::max(0.5, std::min(5.0, gamma - delta_g));
        
        // Convergence check
        if (std::abs(delta_a) < tol && std::abs(delta_b) < tol && std::abs(delta_g) < tol) {
            break;
        }
    }
    
    // Compute log-likelihood
    double loglik = 0;
    for (int i = 0; i < n; i++) {
        loglik += gg_loglik_single(data[i], alpha, beta, gamma, eps);
    }
    
    return List::create(
        Named("alpha") = alpha,
        Named("beta") = beta,
        Named("gamma") = gamma,
        Named("loglik") = loglik
    );
}

// -----------------------------------------------------------------------------
// Batch GG MLE estimation
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List gg_mle_all_cpp(NumericMatrix X, int max_iter = 50, double tol = 1e-6) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    
    NumericVector alpha_out(n_genes), beta_out(n_genes), gamma_out(n_genes);
    NumericVector loglik_out(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int g = 0; g < n_genes; g++) {
        std::vector<double> data(n_samples);
        for (int j = 0; j < n_samples; j++) {
            data[j] = X(g, j);
        }
        
        double eps = 1e-8;
        double alpha, beta, gamma;
        estimate_gg_params(data, alpha, beta, gamma, eps);
        
        // Newton-Raphson
        for (int iter = 0; iter < max_iter; iter++) {
            double grad_a = 0, grad_b = 0, grad_g = 0;
            double hess_a = 0, hess_b = 0, hess_g = 0;
            
            for (int i = 0; i < n_samples; i++) {
                double xi = std::max(data[i], eps);
                double log_xi = std::log(xi);
                double ratio = xi / beta;
                double ratio_g = std::pow(ratio, gamma);
                double log_ratio = std::log(ratio);
                
                grad_a += -gamma * std::log(beta) - R::digamma(alpha) + gamma * log_xi;
                grad_b += -gamma * alpha / beta + gamma * ratio_g / beta;
                grad_g += 1.0 / gamma + alpha * log_ratio - ratio_g * log_ratio;
                
                hess_a += -R::trigamma(alpha);
                hess_b += gamma * alpha / (beta * beta) - gamma * (gamma + 1.0) * ratio_g / (beta * beta);
                hess_g += -1.0 / (gamma * gamma) - ratio_g * log_ratio * log_ratio;
            }
            
            double delta_a = grad_a / std::max(std::abs(hess_a), eps);
            double delta_b = grad_b / std::max(std::abs(hess_b), eps);
            double delta_g = grad_g / std::max(std::abs(hess_g), eps);
            
            alpha = std::max(0.1, std::min(100.0, alpha - delta_a));
            beta = std::max(eps, beta - delta_b);
            gamma = std::max(0.5, std::min(5.0, gamma - delta_g));
            
            if (std::abs(delta_a) < tol && std::abs(delta_b) < tol && std::abs(delta_g) < tol) {
                break;
            }
        }
        
        double loglik = 0;
        for (int i = 0; i < n_samples; i++) {
            loglik += gg_loglik_single(data[i], alpha, beta, gamma, eps);
        }
        
        alpha_out[g] = alpha;
        beta_out[g] = beta;
        gamma_out[g] = gamma;
        loglik_out[g] = loglik;
    }
    
    return List::create(
        Named("alpha") = alpha_out,
        Named("beta") = beta_out,
        Named("gamma") = gamma_out,
        Named("loglik") = loglik_out
    );
}
