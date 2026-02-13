// =============================================================================
// SGCB dynamic sample-size-adaptive inference
// Core idea: borrow more cross-gene information for small samples,
// increase prior weight
// =============================================================================

#include <Rcpp.h>
#include "fast_special.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <numeric>

using namespace Rcpp;

// Exact theoretical variance of the GG distribution
// Var[X] = beta^2 * (Gamma(alpha + 2/gamma)/Gamma(alpha) - (Gamma(alpha + 1/gamma)/Gamma(alpha))^2)
inline double gg_theoretical_var(double alpha, double beta, double gamma, double eps = 1e-8) {
    double log_gamma_a = fast_special::lgamma(alpha);
    double log_gamma_a_1g = fast_special::lgamma(alpha + 1.0/gamma);
    double log_gamma_a_2g = fast_special::lgamma(alpha + 2.0/gamma);
    double E_X = beta * std::exp(log_gamma_a_1g - log_gamma_a);
    double E_X2 = beta * beta * std::exp(log_gamma_a_2g - log_gamma_a);
    return std::max(E_X2 - E_X * E_X, eps);
}

// -----------------------------------------------------------------------------
// Adaptive variance shrinkage (sample-size-aware)
// Small samples: stronger prior weight (larger d0)
// Large samples: more data-dependent (smaller d0)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List adaptive_squeeze_var_cpp(NumericVector var, int df, int n_samples, 
                               double eps = 1e-15) {
    int n = var.size();
    
    // =========================================================================
    // limma-style fitFDist: estimate scaled F-distribution hyperparameters (d0, s0^2)
    // Model: s^2_g ~ s0^2 * F(df, d0)
    // Reference: Smyth (2004) Stat. Appl. Genet. Mol. Biol.
    // =========================================================================
    
    // 1. Log-transform and winsorize (1%–99%)
    std::vector<double> log_var(n);
    for (int i = 0; i < n; i++) {
        log_var[i] = std::log(std::max((double)var[i], eps));
    }
    
    std::vector<double> sorted_lv(log_var);
    std::sort(sorted_lv.begin(), sorted_lv.end());
    int lo_idx = std::max(0, (int)(0.01 * n));
    int hi_idx = std::min(n - 1, (int)(0.99 * n));
    double lo_val = sorted_lv[lo_idx];
    double hi_val = sorted_lv[hi_idx];
    
    for (int i = 0; i < n; i++) {
        log_var[i] = std::max(lo_val, std::min(hi_val, log_var[i]));
    }
    
    // 2. Remove known df contribution: e = log(s²) - ψ(df/2) + log(df/2)
    double half_df1 = std::max(df / 2.0, 0.5);
    double psi_hd1 = fast_special::digamma(half_df1);
    double log_hd1 = std::log(half_df1);
    
    double sum_e = 0;
    for (int i = 0; i < n; i++) {
        sum_e += log_var[i] - psi_hd1 + log_hd1;
    }
    double emean = sum_e / n;
    
    // 3. Variance equation: Var[e] - ψ'(df/2) = ψ'(d0/2)
    double ss_e = 0;
    for (int i = 0; i < n; i++) {
        double e_i = log_var[i] - psi_hd1 + log_hd1;
        double diff = e_i - emean;
        ss_e += diff * diff;
    }
    double evar = ss_e / std::max(n - 1, 1) - fast_special::trigamma(half_df1);
    
    double d0, s0_sq;
    if (evar <= 0) {
        d0 = 1e6;
        s0_sq = std::exp(emean);
    } else {
        double half_df2 = fast_special::trigamma_inverse(evar);
        d0 = std::max(0.5, std::min(1e6, 2.0 * half_df2));
        double hd2 = d0 / 2.0;
        s0_sq = std::exp(emean + fast_special::digamma(hd2) - std::log(hd2));
    }
    s0_sq = std::max(s0_sq, eps);
    
    double df_total = d0 + df;
    
    // 4. Posterior variance: (d0·s0² + df·s²_g) / (d0 + df)
    NumericVector var_post(n);
    for (int i = 0; i < n; i++) {
        var_post[i] = (d0 * s0_sq + df * std::max((double)var[i], eps)) / df_total;
    }
    
    return List::create(
        Named("var_post") = var_post,
        Named("var_prior") = s0_sq,
        Named("df_prior") = d0,
        Named("df_total") = df_total,
        Named("adaptive_factor") = d0 / df_total
    );
}

// -----------------------------------------------------------------------------
// GG-prior-based adaptive t-test
// Uses autoencoder-learned GG parameters as prior information
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List adaptive_moderated_t_cpp(NumericMatrix X, IntegerVector group,
                               NumericVector ae_alpha, NumericVector ae_beta,
                               NumericVector ae_gamma, NumericVector stability,
                               double eps = 1e-8) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    
    // Group assignment
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    int n_ctrl = ctrl_idx.size();
    int n_treat = treat_idx.size();
    int n_min = std::min(n_ctrl, n_treat);
    
    // Sample-size adaptive coefficient (moderate version)
    double adaptive_factor = std::max(1.0, 2.0 / std::sqrt((double)n_min));
    
    // Compute basic statistics
    NumericVector mean_ctrl(n_genes), mean_treat(n_genes);
    NumericVector var_ctrl(n_genes), var_treat(n_genes);
    NumericVector log2FC(n_genes);
    
    for (int g = 0; g < n_genes; g++) {
        double sc = 0, st = 0, ssc = 0, sst = 0;
        
        for (int j = 0; j < n_ctrl; j++) {
            double val = std::log2(X(g, ctrl_idx[j]) + 0.5);
            sc += val;
            ssc += val * val;
        }
        for (int j = 0; j < n_treat; j++) {
            double val = std::log2(X(g, treat_idx[j]) + 0.5);
            st += val;
            sst += val * val;
        }
        
        mean_ctrl[g] = sc / n_ctrl;
        mean_treat[g] = st / n_treat;
        var_ctrl[g] = (ssc - sc * sc / n_ctrl) / std::max(n_ctrl - 1, 1);
        var_treat[g] = (sst - st * st / n_treat) / std::max(n_treat - 1, 1);
        log2FC[g] = mean_treat[g] - mean_ctrl[g];
    }
    
    // Pooled variance
    int df_residual = std::max(n_ctrl + n_treat - 2, 1);
    NumericVector var_pooled(n_genes);
    for (int g = 0; g < n_genes; g++) {
        var_pooled[g] = ((n_ctrl - 1) * std::max(var_ctrl[g], eps) + 
                         (n_treat - 1) * std::max(var_treat[g], eps)) / df_residual;
        var_pooled[g] = std::max(var_pooled[g], eps);
    }
    
    // Adaptive variance shrinkage
    List sq = adaptive_squeeze_var_cpp(var_pooled, df_residual, n_min, eps);
    NumericVector var_post = sq["var_post"];
    double d0 = sq["df_prior"];
    double df_total = sq["df_total"];
    
    // ==========================================================================
    // Core innovation: GG prior weighting
    // Uses autoencoder-learned GG parameter distribution to adjust variance estimation
    // ==========================================================================
    
    // Compute global distribution characteristics of GG parameters
    double alpha_mean = 0, beta_mean = 0;
    for (int g = 0; g < n_genes; g++) {
        alpha_mean += ae_alpha[g];
        beta_mean += ae_beta[g];
    }
    alpha_mean /= n_genes;
    beta_mean /= n_genes;
    
    // GG exact theoretical variance: Var(X) = beta^2*(Gamma(alpha+2/gamma)/Gamma(alpha) - (Gamma(alpha+1/gamma)/Gamma(alpha))^2)
    // No longer using 1/alpha approximation; computing exactly
    NumericVector gg_theoretical_variance(n_genes);
    for (int g = 0; g < n_genes; g++) {
        gg_theoretical_variance[g] = gg_theoretical_var(
            ae_alpha[g], ae_beta[g], ae_gamma[g], eps);
    }
    
    // Fused variance estimate: more reliant on GG prior for small samples
    // weight_gg = adaptive_factor / (adaptive_factor + 1)
    double weight_gg = adaptive_factor / (adaptive_factor + 1.0);
    double weight_data = 1.0 - weight_gg;
    
    // Gene-specific GG-aware variance shrinkage
    // The closer sample variance is to GG theoretical variance, the more we trust the GG prior
    NumericVector var_fused(n_genes);
    for (int g = 0; g < n_genes; g++) {
        // Compute agreement: log-ratio of sample variance to GG theoretical variance
        double log_ratio = std::abs(std::log(var_post[g] + eps) - 
                                    std::log(gg_theoretical_variance[g] + eps));
        double agreement = std::exp(-log_ratio);  // 0-1, closer to 1 means better GG fit
        
        // Dynamic weight: base weight * agreement adjustment
        double gene_weight_gg = weight_gg * (0.5 + 0.5 * agreement);
        
        // Shrinkage: sample variance -> GG theoretical variance
        var_fused[g] = (1.0 - gene_weight_gg) * var_post[g] + 
                       gene_weight_gg * gg_theoretical_variance[g];
        var_fused[g] = std::max(var_fused[g], eps);
    }
    
    // t-statistics and p-values
    NumericVector se(n_genes), t_stat(n_genes), pvalue(n_genes);
    for (int g = 0; g < n_genes; g++) {
        se[g] = std::sqrt(var_fused[g] * (1.0 / n_ctrl + 1.0 / n_treat));
        se[g] = std::max(se[g], eps);
        t_stat[g] = log2FC[g] / se[g];
        pvalue[g] = 2.0 * R::pt(-std::abs(t_stat[g]), df_total, 1, 0);
    }
    
    // ==========================================================================
    // Stability-adaptive penalty (stricter version, FDR control)
    // Small samples: more stringent penalty for unstable genes
    // ==========================================================================
    
    // Stability threshold adapts: moderately raised for small samples
    double stability_threshold = 0.5 + 0.1 / std::sqrt((double)n_min);
    stability_threshold = std::min(stability_threshold, 0.70);
    
    // n<=4: dropout delete-1 unreliable at n=2, disable stability penalty (Politis-Romano-Wolf)
    double penalty_factor = (n_min <= 4) ? 1.0 : (1.0 + 0.3 * adaptive_factor);
    
    for (int g = 0; g < n_genes; g++) {
        if (n_min > 4 && stability[g] < stability_threshold) {
            pvalue[g] = std::min(pvalue[g] * penalty_factor, 1.0);
        }
    }
    
    return List::create(
        Named("log2FC") = log2FC,
        Named("se") = se,
        Named("t_stat") = t_stat,
        Named("pvalue") = pvalue,
        Named("df_total") = df_total,
        Named("adaptive_factor") = adaptive_factor,
        Named("weight_gg") = weight_gg,
        Named("stability_threshold") = stability_threshold
    );
}

// -----------------------------------------------------------------------------
// Complete adaptive SGCB inference
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List sgcb_adaptive_inference_cpp(NumericMatrix X, IntegerVector group,
                                  NumericVector ae_alpha, NumericVector ae_beta,
                                  NumericVector ae_gamma,
                                  NumericVector lfc_stable, NumericVector lfc_var,
                                  NumericVector stability,
                                  double eps = 1e-8) {
    int n_genes = X.nrow();
    int n_samples = X.ncol();
    
    // Group assignment
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    int n_ctrl = ctrl_idx.size();
    int n_treat = treat_idx.size();
    int n_min = std::min(n_ctrl, n_treat);
    
    // Adaptive moderated t-test
    List t_result = adaptive_moderated_t_cpp(X, group, ae_alpha, ae_beta, 
                                              ae_gamma, stability, eps);
    
    NumericVector log2FC = t_result["log2FC"];
    NumericVector se = t_result["se"];
    NumericVector t_stat = t_result["t_stat"];
    NumericVector pvalue = t_result["pvalue"];
    double df_total = t_result["df_total"];
    double adaptive_factor = t_result["adaptive_factor"];
    double weight_gg = t_result["weight_gg"];
    double stability_threshold = t_result["stability_threshold"];
    
    // ==========================================================================
    // Additional innovation: LFC stability-based variance correction
    // Small samples: use dropout variance to supplement information
    // ==========================================================================
    
    // LFC variance adaptive fusion
    // n<=4: dropout variance unreliable, use only moderated t SE
    NumericVector lfc_se_fused(n_genes);
    double weight_dropout = (n_min <= 4) ? 0.0 : adaptive_factor / (adaptive_factor + 2.0);
    
    for (int g = 0; g < n_genes; g++) {
        double se_t = se[g];
        double se_dropout = std::sqrt(lfc_var[g] + eps);
        
        lfc_se_fused[g] = (1.0 - weight_dropout) * se_t + weight_dropout * se_dropout;
        lfc_se_fused[g] = std::max(lfc_se_fused[g], eps);
    }
    
    // Recompute t-statistics and p-values using fused SE (optional)
    NumericVector t_stat_fused(n_genes), pvalue_fused(n_genes);
    for (int g = 0; g < n_genes; g++) {
        t_stat_fused[g] = lfc_stable[g] / lfc_se_fused[g];
        pvalue_fused[g] = 2.0 * R::pt(-std::abs(t_stat_fused[g]), df_total, 1, 0);
    }
    
    // Final p-value: Cauchy combination test (avoids FDR inflation from double-dipping)
    // Cauchy combination: T = sum w_i tan(pi(0.5 - p_i)), p_combined = 0.5 - atan(T/sum w_i)/pi
    // For two equal-weight p-values: T = tan(pi(0.5-p1)) + tan(pi(0.5-p2))
    // More conservative than geometric mean, correctly controls FDR
    NumericVector pvalue_final(n_genes);
    const double pi = 3.14159265358979323846;
    for (int g = 0; g < n_genes; g++) {
        // Cauchy combination test
        double p1 = std::max(eps, std::min(1.0 - eps, pvalue[g]));
        double p2 = std::max(eps, std::min(1.0 - eps, pvalue_fused[g]));
        double t1 = std::tan(pi * (0.5 - p1));
        double t2 = std::tan(pi * (0.5 - p2));
        double T_cauchy = (t1 + t2) / 2.0;  // equal-weight average
        pvalue_final[g] = 0.5 - std::atan(T_cauchy) / pi;
        pvalue_final[g] = std::max(pvalue_final[g], eps);
        pvalue_final[g] = std::min(pvalue_final[g], 1.0);
    }
    
    return List::create(
        Named("log2FC") = lfc_stable,
        Named("log2FC_raw") = log2FC,
        Named("lfcSE") = lfc_se_fused,
        Named("t_stat") = t_stat,
        Named("t_stat_fused") = t_stat_fused,
        Named("pvalue") = pvalue_final,
        Named("pvalue_t") = pvalue,
        Named("pvalue_fused") = pvalue_fused,
        Named("df_total") = df_total,
        Named("adaptive_factor") = adaptive_factor,
        Named("weight_gg") = weight_gg,
        Named("weight_dropout") = weight_dropout,
        Named("stability_threshold") = stability_threshold,
        Named("n_ctrl") = n_ctrl,
        Named("n_treat") = n_treat
    );
}
