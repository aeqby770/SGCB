// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include "fast_special.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// =============================================================================
// Information Geometry Optimizer
// Natural gradient descent with log-space reparameterization
//
// Mathematical formulation:
//   Original space: θ = (α, β, γ)
//   Log space: φ = (log α, log β, log γ)
//   Jacobian: J = ∂θ/∂φ = diag(α, β, γ)
//   Log-space Fisher: F_φ = J^T F_θ J = diag(α², β², γ²) ⊙ F_θ
//   Log-space natural gradient: Δφ = F_φ^{-1} ∇_φ L
// =============================================================================

// -----------------------------------------------------------------------------
// Fisher information matrix for the GG distribution (full 3x3 matrix)
// Based on derivations by Stacy (1962) and Lawless (1980)
// -----------------------------------------------------------------------------
struct FisherMatrix3x3 {
    double I_aa, I_ab, I_ag;  // first row
    double I_bb, I_bg;        // second row (diagonal and upper)
    double I_gg;              // third row (diagonal)
    
    // Compute diagonal elements of full 3x3 inverse via Cramer's rule
    // Resolves the zigzag optimization path caused by strong α-β coupling
    void compute_full_inverse(double& inv_aa, double& inv_bb, double& inv_gg,
                               double& inv_ab) const {
        double reg = 1e-8;
        
        // Determinant of 3x3 symmetric matrix
        // |F| = I_aa(I_bb*I_gg - I_bg^2) - I_ab(I_ab*I_gg - I_bg*I_ag) + I_ag(I_ab*I_bg - I_bb*I_ag)
        double minor_aa = I_bb * I_gg - I_bg * I_bg;
        double minor_ab = I_ab * I_gg - I_bg * I_ag;
        double minor_ag = I_ab * I_bg - I_bb * I_ag;
        
        double det = I_aa * minor_aa - I_ab * minor_ab + I_ag * minor_ag;
        det = std::max(det, reg);  // prevent singularity
        
        double inv_det = 1.0 / det;
        
        // Diagonal elements of inverse (cofactors / determinant)
        inv_aa = minor_aa * inv_det;
        inv_bb = (I_aa * I_gg - I_ag * I_ag) * inv_det;
        inv_gg = (I_aa * I_bb - I_ab * I_ab) * inv_det;
        
        // α-β cross-term (key: captures parameter coupling)
        inv_ab = -minor_ab * inv_det;
    }
    
    // Compute 2x2 α-β block inverse (simplified when γ is independent)
    // Key to resolving the "banana-shaped" contour problem
    void compute_block_inverse_ab(double& inv_aa, double& inv_bb, double& inv_ab) const {
        double reg = 1e-8;
        double det_ab = I_aa * I_bb - I_ab * I_ab;
        det_ab = std::max(det_ab, reg);
        
        double inv_det = 1.0 / det_ab;
        inv_aa = I_bb * inv_det;
        inv_bb = I_aa * inv_det;
        inv_ab = -I_ab * inv_det;
    }
};

// -----------------------------------------------------------------------------
// Compute the full Fisher information matrix for the GG distribution
// Parameters: alpha (shape), beta (scale), gamma (power parameter)
// n: sample size
// -----------------------------------------------------------------------------
FisherMatrix3x3 compute_fisher_gg(double alpha, double beta, double gamma, int n) {
    FisherMatrix3x3 F;
    
    // Use fast approximations instead of R's special functions
    double psi1_a = fast_special::trigamma(alpha);  // ψ'(α)
    double psi_a = fast_special::digamma(alpha);    // ψ(α)
    
    const double euler_gamma = 0.5772156649015329;
    
    // Fisher information matrix elements (Lawless 1980)
    F.I_aa = n * psi1_a;
    F.I_bb = n * alpha * gamma * gamma / (beta * beta);
    F.I_gg = n * (1.0 + alpha * psi1_a) / (gamma * gamma);
    F.I_ab = n * gamma / beta;
    F.I_ag = -n * (psi_a - euler_gamma) / gamma;
    F.I_bg = -n * alpha / beta;
    
    return F;
}

// =============================================================================
// Jeffreys/Firth penalty: computation of (1/2) log|I(θ)|
// Small-sample MLE regularization to stabilize likelihood surface singularities
// Reference: Firth (1993) Biometrika, Prentice (1974)
// =============================================================================
inline double compute_log_det_fisher_unit(double alpha, double beta, double gamma, double eps = 1e-8) {
    FisherMatrix3x3 F = compute_fisher_gg(alpha, beta, gamma, 1);
    double minor_aa = F.I_bb * F.I_gg - F.I_bg * F.I_bg;
    double minor_ab = F.I_ab * F.I_gg - F.I_bg * F.I_ag;
    double minor_ag = F.I_ab * F.I_bg - F.I_bb * F.I_ag;
    double det = F.I_aa * minor_aa - F.I_ab * minor_ab + F.I_ag * minor_ag;
    return std::log(std::max(det, eps));
}

inline void invert_sym_3x3(double aa, double ab, double ag,
                            double bb, double bg,
                            double gg,
                            double& inv_aa, double& inv_ab, double& inv_ag,
                            double& inv_bb, double& inv_bg,
                            double& inv_gg) {
     double reg = 1e-10;
     double tr = aa + bb + gg;
     double ridge = 1e-6 * std::max(tr, 0.0) + reg;
     aa += ridge;
     bb += ridge;
     gg += ridge;
     double minor_aa = bb * gg - bg * bg;
     double minor_ab = ab * gg - bg * ag;
     double minor_ag = ab * bg - bb * ag;
     double det = aa * minor_aa - ab * minor_ab + ag * minor_ag;
     det = std::max(det, reg);
     double inv_det = 1.0 / det;
     inv_aa = minor_aa * inv_det;
     inv_bb = (aa * gg - ag * ag) * inv_det;
     inv_gg = (aa * bb - ab * ab) * inv_det;
     inv_ab = -minor_ab * inv_det;
     inv_ag = minor_ag * inv_det;
     inv_bg = (ab * ag - aa * bg) * inv_det;
 }

// -----------------------------------------------------------------------------
// Log-space reparameterized GG distribution
// Optimizes θ = (log α, log β, log γ) instead of (α, β, γ)
// Avoids hard clipping and respects the Lie group structure of positive reals
// -----------------------------------------------------------------------------
struct LogSpaceGGParams {
    double log_alpha, log_beta, log_gamma;
    
    double alpha() const { return std::exp(log_alpha); }
    double beta() const { return std::exp(log_beta); }
    double gamma() const { return std::exp(log_gamma); }
    
    // GG log-likelihood in log space (using fast lgamma)
    double loglik(double x, double eps = 1e-8) const {
        double a = alpha();
        double b = beta();
        double g = gamma();
        double xi = std::max(x, eps);
        double log_xi = std::log(xi);
        double log_b = std::log(b);
        return std::log(g) - a * log_b - fast_special::lgamma(a) +
               (g * a - 1.0) * log_xi - std::pow(xi / b, g);
    }
    
    // Gradient in log space (w.r.t. log_alpha, log_beta, log_gamma)
    // Chain rule: ∂L/∂(log θ) = ∂L/∂θ · θ
    void gradient(double x, double& grad_la, double& grad_lb, double& grad_lg, 
                  double eps = 1e-8) const {
        double a = alpha();
        double b = beta();
        double g = gamma();
        double xi = std::max(x, eps);
        double log_xi = std::log(xi);
        double ratio = xi / b;
        double ratio_g = std::pow(ratio, g);
        double log_ratio = std::log(ratio);
        
        // Original-space gradient
        double grad_a = -std::log(b) - fast_special::digamma(a) + g * log_xi;
        double grad_b = -a / b + g * ratio_g / b;
        double grad_g = 1.0 / g + a * log_xi - ratio_g * log_ratio;
        
        // Chain rule transformation to log space
        grad_la = grad_a * a;
        grad_lb = grad_b * b;
        grad_lg = grad_g * g;
    }
};

// -----------------------------------------------------------------------------
// Riemannian momentum state (Polyak heavy ball on GG manifold)
// Natural gradient F^{-1}∇ℓ already incorporates curvature; first-order momentum only
// Reference: Bécigneul & Ganea (2019) ICLR
// -----------------------------------------------------------------------------
struct NaturalAdamState {
    std::vector<double> m_la, m_lb, m_lg;  // first-order momentum (log space)
    int t;                                   // time step
    
    void init(int n_genes) {
        m_la.resize(n_genes, 0.0);
        m_lb.resize(n_genes, 0.0);
        m_lg.resize(n_genes, 0.0);
        t = 0;
    }
};

// -----------------------------------------------------------------------------
// Adaptive configuration (sample-size-dependent information geometry settings)
// -----------------------------------------------------------------------------
struct InfoGeomConfig {
    int n_iter;
    double lr;
    double beta1, beta2;
    double weight_decay;
    double fisher_reg;      // Fisher matrix regularization
    bool use_natural_grad;  // whether to use natural gradient
    
    static InfoGeomConfig adaptive(int n_samples) {
        InfoGeomConfig cfg;
        cfg.beta1 = 0.9;
        cfg.beta2 = 0.999;
        cfg.fisher_reg = 1e-4;
        cfg.use_natural_grad = true;
        
        if (n_samples <= 5) {
            // Very small sample: strong regularization, more iterations, smaller learning rate
            cfg.n_iter = 300;
            cfg.lr = 0.005;
            cfg.weight_decay = 0.05;
        } else if (n_samples <= 10) {
            // Small sample
            cfg.n_iter = 200;
            cfg.lr = 0.01;
            cfg.weight_decay = 0.02;
        } else if (n_samples <= 30) {
            // Moderate sample
            cfg.n_iter = 150;
            cfg.lr = 0.02;
            cfg.weight_decay = 0.01;
        } else {
            // Large sample: weak regularization, fewer iterations
            cfg.n_iter = 100;
            cfg.lr = 0.03;
            cfg.weight_decay = 0.005;
        }
        return cfg;
    }
};

// =============================================================================
// Core function: natural gradient GG parameter fitting
// =============================================================================
// [[Rcpp::export]]
List fit_gg_natural_gradient_cpp(NumericMatrix X, 
                                  int n_iter = -1,
                                  double lr = -1,
                                  double weight_decay = -1,
                                  bool use_natural_grad = true,
                                  double eps = 1e-8) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    // Adaptive configuration
    InfoGeomConfig cfg = InfoGeomConfig::adaptive(n_samples);
    if (n_iter > 0) cfg.n_iter = n_iter;
    if (lr > 0) cfg.lr = lr;
    if (weight_decay > 0) cfg.weight_decay = weight_decay;
    cfg.use_natural_grad = use_natural_grad;
    
    // Small-sample submodel reduction + Firth penalty (Prentice 1974, Firth 1993)
    // n<=5: fix γ=1 (reduce to standard Gamma) to avoid non-identifiability with 3 params vs <=5 data points
    // n=6-10: add Jeffreys penalty +(1/2)log|I(θ)| to stabilize MLE
    bool fix_gamma = (n_samples <= 5);
    bool use_firth = (n_samples <= 10);
    
    // Log-space parameter initialization
    std::vector<double> log_alpha(n_genes), log_beta(n_genes), log_gamma(n_genes);
    
    // Data-driven initialization
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum_x = 0, sum_x2 = 0, sum_logx = 0;
        for (int j = 0; j < n_samples; j++) {
            double xij = X(g, j) + eps;
            sum_x += xij;
            sum_x2 += xij * xij;
            sum_logx += std::log(xij);
        }
        double mean_x = sum_x / n_samples;
        double var_x = sum_x2 / n_samples - mean_x * mean_x;
        double mean_logx = sum_logx / n_samples;
        
        // Method-of-moments initialization
        double cv2 = var_x / (mean_x * mean_x + eps);
        double init_alpha = std::max(0.5, std::min(10.0, 1.0 / cv2));
        double init_gamma = fix_gamma ? 1.0 : std::max(0.5, std::min(3.0, 1.0 / std::sqrt(cv2 + eps)));
        double init_beta = std::max(eps, mean_x / std::pow(init_alpha, 1.0 / init_gamma));
        
        log_alpha[g] = std::log(init_alpha);
        log_beta[g] = std::log(init_beta);
        log_gamma[g] = std::log(init_gamma);
    }
    
    // Adam state
    NaturalAdamState state;
    state.init(n_genes);
    
    // Training loop
    NumericVector loss_history(cfg.n_iter);
    
    for (int iter = 0; iter < cfg.n_iter; iter++) {
        state.t++;
        double total_loss = 0;
        
        // Compute gradients
        std::vector<double> grad_la(n_genes), grad_lb(n_genes), grad_lg(n_genes);
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:total_loss) schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            LogSpaceGGParams params;
            params.log_alpha = log_alpha[g];
            params.log_beta = log_beta[g];
            params.log_gamma = log_gamma[g];
            
            double sum_grad_la = 0, sum_grad_lb = 0, sum_grad_lg = 0;
            double sum_loss = 0;
            
            for (int j = 0; j < n_samples; j++) {
                double xgj = X(g, j);
                sum_loss += params.loglik(xgj, eps);
                
                double gla, glb, glg;
                params.gradient(xgj, gla, glb, glg, eps);
                sum_grad_la += gla;
                sum_grad_lb += glb;
                sum_grad_lg += glg;
            }
            
            total_loss -= sum_loss;  // negative log-likelihood
            
            // Average gradient (maximize likelihood = gradient ascent = negate for descent)
            grad_la[g] = -sum_grad_la / n_samples;
            grad_lb[g] = -sum_grad_lb / n_samples;
            grad_lg[g] = -sum_grad_lg / n_samples;
            
            // Jeffreys/Firth penalty: -(1/2) d/dφ log|I(θ)| / n
            // Numerical differentiation of log|det(F_unit)| w.r.t. log-parameters
            if (use_firth) {
                double a = params.alpha(), b = params.beta(), gam = params.gamma();
                double h = 1e-4;
                double firth_scale = 0.5 / n_samples;
                
                double ldf_ap = compute_log_det_fisher_unit(a * std::exp(h), b, gam, eps);
                double ldf_am = compute_log_det_fisher_unit(a * std::exp(-h), b, gam, eps);
                grad_la[g] -= firth_scale * (ldf_ap - ldf_am) / (2.0 * h);
                
                double ldf_bp = compute_log_det_fisher_unit(a, b * std::exp(h), gam, eps);
                double ldf_bm = compute_log_det_fisher_unit(a, b * std::exp(-h), gam, eps);
                grad_lb[g] -= firth_scale * (ldf_bp - ldf_bm) / (2.0 * h);
                
                if (!fix_gamma) {
                    double ldf_gp = compute_log_det_fisher_unit(a, b, gam * std::exp(h), eps);
                    double ldf_gm = compute_log_det_fisher_unit(a, b, gam * std::exp(-h), eps);
                    grad_lg[g] -= firth_scale * (ldf_gp - ldf_gm) / (2.0 * h);
                }
            }
            
            // Submodel reduction: zero out γ gradient when fixed
            if (fix_gamma) {
                grad_lg[g] = 0.0;
            }
            
            // Natural gradient correction (full Fisher matrix, resolves α-β coupling)
            // Derivation:
            //   Log-space Fisher: F_φ = J^T F_θ J, J = diag(α,β,γ)
            //   Natural gradient: Δφ = F_φ^{-1} ∇_φ L
            //   Uses full F^{-1} instead of diagonal approximation to capture α-β banana-shaped contours
            if (cfg.use_natural_grad) {
                FisherMatrix3x3 F = compute_fisher_gg(
                    params.alpha(), params.beta(), params.gamma(), n_samples);
                
                double a = params.alpha();
                double b = params.beta();
                double gam = params.gamma();
                
                // Compute full Fisher inverse (including α-β cross-terms)
                double inv_aa, inv_bb, inv_gg, inv_ab;
                F.compute_full_inverse(inv_aa, inv_bb, inv_gg, inv_ab);
                
                // Euclidean gradient in original space (grad_la = α·∇_α, divide by α to recover)
                double grad_a = grad_la[g] / a;
                double grad_b = grad_lb[g] / b;
                double grad_g = grad_lg[g] / gam;
                
                // Full natural gradient (F^{-1} · ∇) in original parameter space
                // d_theta = F_theta^{-1} · grad_theta
                double nat_grad_a = inv_aa * grad_a + inv_ab * grad_b;
                double nat_grad_b = inv_ab * grad_a + inv_bb * grad_b;
                double nat_grad_g = inv_gg * grad_g;
                
                // Transform back to log space: d_phi = d_theta / theta
                // Since phi = log(theta), d_phi/d_theta = 1/theta
                // So d_phi = d_theta / theta
                grad_la[g] = nat_grad_a / a;  // divide by a, not multiply!
                grad_lb[g] = nat_grad_b / b;
                grad_lg[g] = nat_grad_g / gam;
            }
            
            // L2 regularization
            grad_la[g] += cfg.weight_decay * log_alpha[g];
            grad_lb[g] += cfg.weight_decay * log_beta[g];
            grad_lg[g] += cfg.weight_decay * log_gamma[g];
        }
        
        loss_history[iter] = total_loss / n_genes;
        
        // Riemannian momentum update (Polyak heavy ball, no Adam second moment)
        // Natural gradient already corrects curvature via Fisher inverse; Adam v_t diagonal scaling is redundant
        // Reference: Bécigneul & Ganea (2019) "Riemannian Adaptive Optimization Methods"
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            // First-order momentum (exponential moving average)
            state.m_la[g] = cfg.beta1 * state.m_la[g] + (1.0 - cfg.beta1) * grad_la[g];
            state.m_lb[g] = cfg.beta1 * state.m_lb[g] + (1.0 - cfg.beta1) * grad_lb[g];
            state.m_lg[g] = cfg.beta1 * state.m_lg[g] + (1.0 - cfg.beta1) * grad_lg[g];
            
            // Parameter update (gradient descent + momentum)
            log_alpha[g] -= cfg.lr * state.m_la[g];
            log_beta[g] -= cfg.lr * state.m_lb[g];
            log_gamma[g] -= cfg.lr * state.m_lg[g];
            
            // Soft constraints (reasonable range in log space)
            log_alpha[g] = std::max(-2.3, std::min(4.6, log_alpha[g]));
            log_beta[g] = std::max(-10.0, std::min(15.0, log_beta[g]));
            log_gamma[g] = std::max(-0.7, std::min(1.6, log_gamma[g]));
            
            // Submodel reduction: enforce γ=1
            if (fix_gamma) {
                log_gamma[g] = 0.0;
            }
        }
    }
    
    // Transform back to original space
    NumericVector alpha_out(n_genes), beta_out(n_genes), gamma_out(n_genes);
    for (int g = 0; g < n_genes; g++) {
        alpha_out[g] = std::exp(log_alpha[g]);
        beta_out[g] = std::exp(log_beta[g]);
        gamma_out[g] = std::exp(log_gamma[g]);
    }
    
    return List::create(
        Named("alpha") = alpha_out,
        Named("beta") = beta_out,
        Named("gamma") = gamma_out,
        Named("loss_history") = loss_history,
        Named("n_iter") = cfg.n_iter,
        Named("lr") = cfg.lr,
        Named("use_natural_grad") = cfg.use_natural_grad
    );
}

// =============================================================================
// Spike-and-slab style hierarchical Gaussian mixture prior
// Two components:
//   - Background (spike, narrow variance): captures the majority of non-DE genes
//   - Signal (slab, wide variance): captures outlier DE genes
// Prevents true DE genes from being over-shrunk toward the background mean
// =============================================================================
struct HierarchicalPrior {
    // Background component (spike)
    double mu_bg, sigma_bg;
    // Signal component (slab)
    double mu_sig, sigma_sig;
    // Mixture weight (probability of background component)
    double pi_bg;
    
    // Independently estimated for each parameter
    double mu_log_alpha, sigma_log_alpha, sigma_wide_alpha;
    double mu_log_beta, sigma_log_beta, sigma_wide_beta;
    double mu_log_gamma, sigma_log_gamma, sigma_wide_gamma;
    double mix_weight;  // background component weight
    
    // Estimate mixture prior hyperparameters from data
    void estimate_from_data(const std::vector<double>& log_alpha,
                            const std::vector<double>& log_beta,
                            const std::vector<double>& log_gamma) {
        int n = log_alpha.size();
        
        // Fix B: use median for location estimation to prevent DE gene contamination of hyperprior
        {
            std::vector<double> tmp(n);
            for (int i = 0; i < n; i++) tmp[i] = log_alpha[i];
            std::nth_element(tmp.begin(), tmp.begin() + n/2, tmp.end());
            mu_log_alpha = tmp[n/2];
            for (int i = 0; i < n; i++) tmp[i] = log_beta[i];
            std::nth_element(tmp.begin(), tmp.begin() + n/2, tmp.end());
            mu_log_beta = tmp[n/2];
            for (int i = 0; i < n; i++) tmp[i] = log_gamma[i];
            std::nth_element(tmp.begin(), tmp.begin() + n/2, tmp.end());
            mu_log_gamma = tmp[n/2];
        }
        
        double ss_la = 0, ss_lb = 0, ss_lg = 0;
        for (int i = 0; i < n; i++) {
            ss_la += (log_alpha[i] - mu_log_alpha) * (log_alpha[i] - mu_log_alpha);
            ss_lb += (log_beta[i] - mu_log_beta) * (log_beta[i] - mu_log_beta);
            ss_lg += (log_gamma[i] - mu_log_gamma) * (log_gamma[i] - mu_log_gamma);
        }
        double var_la = ss_la / n + 1e-6;
        double var_lb = ss_lb / n + 1e-6;
        double var_lg = ss_lg / n + 1e-6;
        
        // Background component: robust standard deviation via median absolute deviation (MAD)
        std::vector<double> abs_dev_la(n), abs_dev_lb(n), abs_dev_lg(n);
        for (int i = 0; i < n; i++) {
            abs_dev_la[i] = std::abs(log_alpha[i] - mu_log_alpha);
            abs_dev_lb[i] = std::abs(log_beta[i] - mu_log_beta);
            abs_dev_lg[i] = std::abs(log_gamma[i] - mu_log_gamma);
        }
        std::sort(abs_dev_la.begin(), abs_dev_la.end());
        std::sort(abs_dev_lb.begin(), abs_dev_lb.end());
        std::sort(abs_dev_lg.begin(), abs_dev_lg.end());
        
        double mad_la = abs_dev_la[n/2] * 1.4826;  // MAD-to-SD conversion factor
        double mad_lb = abs_dev_lb[n/2] * 1.4826;
        double mad_lg = abs_dev_lg[n/2] * 1.4826;
        
        // Background component (narrow, MAD-based)
        sigma_log_alpha = std::max(0.1, mad_la);
        sigma_log_beta = std::max(0.1, mad_lb);
        sigma_log_gamma = std::max(0.05, mad_lg);
        
        // Signal component (wide, 3-5x background)
        sigma_wide_alpha = std::max(sigma_log_alpha * 4.0, std::sqrt(var_la) * 2.0);
        sigma_wide_beta = std::max(sigma_log_beta * 4.0, std::sqrt(var_lb) * 2.0);
        sigma_wide_gamma = std::max(sigma_log_gamma * 3.0, std::sqrt(var_lg) * 1.5);
        
        // Mixture weight π₀: estimated via simple EM from log_alpha distribution (Efron 2004 two-groups)
        // Data-driven rather than hard-coded
        mix_weight = 0.9;  // EM initial value
        for (int em_iter = 0; em_iter < 5; em_iter++) {
            double sum_resp = 0;
            for (int i = 0; i < n; i++) {
                double z_n = (log_alpha[i] - mu_log_alpha) / sigma_log_alpha;
                double z_w = (log_alpha[i] - mu_log_alpha) / sigma_wide_alpha;
                double p_n = mix_weight * std::exp(-0.5 * z_n * z_n) / sigma_log_alpha;
                double p_w = (1.0 - mix_weight) * std::exp(-0.5 * z_w * z_w) / sigma_wide_alpha;
                double resp = p_n / (p_n + p_w + 1e-10);
                sum_resp += resp;
            }
            mix_weight = std::max(0.5, std::min(0.99, sum_resp / n));
        }
    }
    
    // Negative log-density of the mixture prior
    // -log(π·N(μ,σ²) + (1-π)·N(μ,σ_wide²))
    double neg_log_prior(double la, double lb, double lg) const {
        // Compute mixture density for each parameter
        auto mix_nll = [this](double x, double mu, double sig_narrow, double sig_wide) {
            double z_narrow = (x - mu) / sig_narrow;
            double z_wide = (x - mu) / sig_wide;
            
            // log-sum-exp for numerical stability
            double log_p_narrow = -0.5 * z_narrow * z_narrow - std::log(sig_narrow);
            double log_p_wide = -0.5 * z_wide * z_wide - std::log(sig_wide);
            
            double log_mix = std::log(mix_weight) + log_p_narrow;
            double log_slab = std::log(1.0 - mix_weight) + log_p_wide;
            
            // log(exp(a) + exp(b)) = a + log(1 + exp(b-a))
            double max_log = std::max(log_mix, log_slab);
            double log_sum = max_log + std::log(std::exp(log_mix - max_log) + 
                                                 std::exp(log_slab - max_log));
            return -log_sum;
        };
        
        return mix_nll(la, mu_log_alpha, sigma_log_alpha, sigma_wide_alpha) +
               mix_nll(lb, mu_log_beta, sigma_log_beta, sigma_wide_beta) +
               mix_nll(lg, mu_log_gamma, sigma_log_gamma, sigma_wide_gamma);
    }
    
    // Mixture prior gradient (soft assignment)
    void prior_gradient(double la, double lb, double lg,
                        double& grad_la, double& grad_lb, double& grad_lg) const {
        // Compute posterior probability of each parameter belonging to background
        auto compute_grad = [this](double x, double mu, double sig_narrow, double sig_wide) {
            double z_narrow = (x - mu) / sig_narrow;
            double z_wide = (x - mu) / sig_wide;
            
            double p_narrow = mix_weight * std::exp(-0.5 * z_narrow * z_narrow) / sig_narrow;
            double p_wide = (1.0 - mix_weight) * std::exp(-0.5 * z_wide * z_wide) / sig_wide;
            double p_total = p_narrow + p_wide + 1e-10;
            
            // Posterior-probability-weighted gradient
            double w_narrow = p_narrow / p_total;
            double w_wide = p_wide / p_total;
            
            // Gradient = w_narrow * (x-μ)/σ_narrow² + w_wide * (x-μ)/σ_wide²
            return w_narrow * (x - mu) / (sig_narrow * sig_narrow) +
                   w_wide * (x - mu) / (sig_wide * sig_wide);
        };
        
        grad_la = compute_grad(la, mu_log_alpha, sigma_log_alpha, sigma_wide_alpha);
        grad_lb = compute_grad(lb, mu_log_beta, sigma_log_beta, sigma_wide_beta);
        grad_lg = compute_grad(lg, mu_log_gamma, sigma_log_gamma, sigma_wide_gamma);
    }
};

// =============================================================================
// Natural gradient fitting with hierarchical Bayesian prior
// =============================================================================
// [[Rcpp::export]]
List fit_gg_hierarchical_bayes_cpp(NumericMatrix X,
                                    int n_iter = -1,
                                    double lr = -1,
                                    double prior_strength = -1,
                                    bool use_natural_grad = true,
                                    double eps = 1e-8) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    // Adaptive configuration
    InfoGeomConfig cfg = InfoGeomConfig::adaptive(n_samples);
    if (n_iter > 0) cfg.n_iter = n_iter;
    if (lr > 0) cfg.lr = lr;
    cfg.use_natural_grad = use_natural_grad;
    
    // Prior strength adapts to sample size
    double prior_weight = prior_strength > 0 ? prior_strength : 
                          std::max(0.01, 2.0 / std::sqrt(n_samples + 1.0));
    
    // Small-sample submodel reduction + Firth penalty (Prentice 1974, Firth 1993)
    bool fix_gamma = (n_samples <= 5);
    bool use_firth = (n_samples <= 10);
    
    // Stage 1: fast MLE initialization
    std::vector<double> log_alpha(n_genes), log_beta(n_genes), log_gamma(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum_x = 0, sum_x2 = 0;
        for (int j = 0; j < n_samples; j++) {
            double xij = X(g, j) + eps;
            sum_x += xij;
            sum_x2 += xij * xij;
        }
        double mean_x = sum_x / n_samples;
        double var_x = sum_x2 / n_samples - mean_x * mean_x;
        double cv2 = var_x / (mean_x * mean_x + eps);
        
        log_alpha[g] = std::log(std::max(0.5, std::min(10.0, 1.0 / cv2)));
        log_gamma[g] = fix_gamma ? 0.0 : std::log(std::max(0.5, std::min(3.0, 1.0 / std::sqrt(cv2 + eps))));
        log_beta[g] = std::log(std::max(eps, mean_x));
    }
    
    // Estimate hierarchical prior
    HierarchicalPrior prior;
    prior.estimate_from_data(log_alpha, log_beta, log_gamma);
    
    // Adam state
    NaturalAdamState state;
    state.init(n_genes);
    
    // Stage 2: MAP estimation with prior
    NumericVector loss_history(cfg.n_iter);
    
    for (int iter = 0; iter < cfg.n_iter; iter++) {
        state.t++;
        double total_loss = 0;
        
        std::vector<double> grad_la(n_genes), grad_lb(n_genes), grad_lg(n_genes);
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:total_loss) schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            LogSpaceGGParams params;
            params.log_alpha = log_alpha[g];
            params.log_beta = log_beta[g];
            params.log_gamma = log_gamma[g];
            
            double sum_grad_la = 0, sum_grad_lb = 0, sum_grad_lg = 0;
            double sum_loss = 0;
            
            for (int j = 0; j < n_samples; j++) {
                double xgj = X(g, j);
                sum_loss += params.loglik(xgj, eps);
                
                double gla, glb, glg;
                params.gradient(xgj, gla, glb, glg, eps);
                sum_grad_la += gla;
                sum_grad_lb += glb;
                sum_grad_lg += glg;
            }
            
            // Likelihood gradient (maximization)
            grad_la[g] = -sum_grad_la / n_samples;
            grad_lb[g] = -sum_grad_lb / n_samples;
            grad_lg[g] = -sum_grad_lg / n_samples;
            
            // Jeffreys/Firth penalty
            if (use_firth) {
                double a = params.alpha(), b = params.beta(), gam = params.gamma();
                double h = 1e-4;
                double firth_scale = 0.5 / n_samples;
                
                double ldf_ap = compute_log_det_fisher_unit(a * std::exp(h), b, gam, eps);
                double ldf_am = compute_log_det_fisher_unit(a * std::exp(-h), b, gam, eps);
                grad_la[g] -= firth_scale * (ldf_ap - ldf_am) / (2.0 * h);
                
                double ldf_bp = compute_log_det_fisher_unit(a, b * std::exp(h), gam, eps);
                double ldf_bm = compute_log_det_fisher_unit(a, b * std::exp(-h), gam, eps);
                grad_lb[g] -= firth_scale * (ldf_bp - ldf_bm) / (2.0 * h);
                
                if (!fix_gamma) {
                    double ldf_gp = compute_log_det_fisher_unit(a, b, gam * std::exp(h), eps);
                    double ldf_gm = compute_log_det_fisher_unit(a, b, gam * std::exp(-h), eps);
                    grad_lg[g] -= firth_scale * (ldf_gp - ldf_gm) / (2.0 * h);
                }
            }
            
            if (fix_gamma) {
                grad_lg[g] = 0.0;
            }
            
            // Add hierarchical Bayesian prior gradient
            double prior_gla, prior_glb, prior_glg;
            prior.prior_gradient(log_alpha[g], log_beta[g], log_gamma[g],
                                prior_gla, prior_glb, prior_glg);
            
            grad_la[g] += prior_weight * prior_gla;
            grad_lb[g] += prior_weight * prior_glb;
            grad_lg[g] += prior_weight * prior_glg;
            
            // Natural gradient correction (full Fisher matrix, including α-β coupling)
            if (cfg.use_natural_grad) {
                FisherMatrix3x3 F = compute_fisher_gg(
                    params.alpha(), params.beta(), params.gamma(), n_samples);
                
                double a = params.alpha();
                double b = params.beta();
                double gam = params.gamma();
                
                double inv_aa, inv_bb, inv_gg, inv_ab;
                F.compute_full_inverse(inv_aa, inv_bb, inv_gg, inv_ab);
                
                double grad_a = grad_la[g] / a;
                double grad_b = grad_lb[g] / b;
                double grad_g = grad_lg[g] / gam;
                
                double nat_grad_a = inv_aa * grad_a + inv_ab * grad_b;
                double nat_grad_b = inv_ab * grad_a + inv_bb * grad_b;
                double nat_grad_g = inv_gg * grad_g;
                
                // Transform back to log space: d_phi = d_theta / theta
                grad_la[g] = nat_grad_a / a;
                grad_lb[g] = nat_grad_b / b;
                grad_lg[g] = nat_grad_g / gam;
            }
            
            total_loss -= sum_loss;
            total_loss += prior_weight * prior.neg_log_prior(log_alpha[g], log_beta[g], log_gamma[g]);
        }
        
        loss_history[iter] = total_loss / n_genes;
        
        // Riemannian momentum update (consistent with fit_gg_natural_gradient_cpp)
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            state.m_la[g] = cfg.beta1 * state.m_la[g] + (1.0 - cfg.beta1) * grad_la[g];
            state.m_lb[g] = cfg.beta1 * state.m_lb[g] + (1.0 - cfg.beta1) * grad_lb[g];
            state.m_lg[g] = cfg.beta1 * state.m_lg[g] + (1.0 - cfg.beta1) * grad_lg[g];
            
            log_alpha[g] -= cfg.lr * state.m_la[g];
            log_beta[g] -= cfg.lr * state.m_lb[g];
            log_gamma[g] -= cfg.lr * state.m_lg[g];
            
            log_alpha[g] = std::max(-2.3, std::min(4.6, log_alpha[g]));
            log_beta[g] = std::max(-10.0, std::min(15.0, log_beta[g]));
            log_gamma[g] = std::max(-0.7, std::min(1.6, log_gamma[g]));
            
            if (fix_gamma) {
                log_gamma[g] = 0.0;
            }
        }
        
        // Update prior (empirical Bayes)
        if ((iter + 1) % 20 == 0) {
            prior.estimate_from_data(log_alpha, log_beta, log_gamma);
        }
    }
    
    // Output
    NumericVector alpha_out(n_genes), beta_out(n_genes), gamma_out(n_genes);
    for (int g = 0; g < n_genes; g++) {
        alpha_out[g] = std::exp(log_alpha[g]);
        beta_out[g] = std::exp(log_beta[g]);
        gamma_out[g] = std::exp(log_gamma[g]);
    }
    
    return List::create(
        Named("alpha") = alpha_out,
        Named("beta") = beta_out,
        Named("gamma") = gamma_out,
        Named("loss_history") = loss_history,
        Named("n_iter") = cfg.n_iter,
        Named("lr") = cfg.lr,
        Named("prior_weight") = prior_weight,
        Named("prior_mu_log_alpha") = prior.mu_log_alpha,
        Named("prior_mu_log_beta") = prior.mu_log_beta,
        Named("prior_mu_log_gamma") = prior.mu_log_gamma,
        Named("prior_sigma_log_alpha") = prior.sigma_log_alpha,
        Named("prior_sigma_log_beta") = prior.sigma_log_beta,
        Named("prior_sigma_log_gamma") = prior.sigma_log_gamma
    );
}

// =============================================================================
// Fix A: GG MAP fitting with shared prior (internal function, not exported)
// Same logic as fit_gg_hierarchical_bayes_cpp but uses externally given prior without re-estimation
// =============================================================================
static List fit_gg_with_shared_prior(NumericMatrix X,
                                      const HierarchicalPrior& shared_prior,
                                      bool use_natural_grad,
                                      double eps) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    InfoGeomConfig cfg = InfoGeomConfig::adaptive(n_samples);
    cfg.use_natural_grad = use_natural_grad;
    
    double prior_weight = std::max(0.01, 2.0 / std::sqrt(n_samples + 1.0));
    bool fix_gamma = (n_samples <= 5);
    bool use_firth = (n_samples <= 10);
    
    // Fast MLE initialization
    std::vector<double> log_alpha(n_genes), log_beta(n_genes), log_gamma(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double sum_x = 0, sum_x2 = 0;
        for (int j = 0; j < n_samples; j++) {
            double xij = X(g, j) + eps;
            sum_x += xij;
            sum_x2 += xij * xij;
        }
        double mean_x = sum_x / n_samples;
        double var_x = sum_x2 / n_samples - mean_x * mean_x;
        double cv2 = var_x / (mean_x * mean_x + eps);
        
        log_alpha[g] = std::log(std::max(0.5, std::min(10.0, 1.0 / cv2)));
        log_gamma[g] = fix_gamma ? 0.0 : std::log(std::max(0.5, std::min(3.0, 1.0 / std::sqrt(cv2 + eps))));
        log_beta[g] = std::log(std::max(eps, mean_x));
    }
    
    // Use shared prior (no re-estimation)
    HierarchicalPrior prior = shared_prior;
    
    NaturalAdamState state;
    state.init(n_genes);
    
    NumericVector loss_history(cfg.n_iter);
    
    for (int iter = 0; iter < cfg.n_iter; iter++) {
        state.t++;
        double total_loss = 0;
        std::vector<double> grad_la(n_genes), grad_lb(n_genes), grad_lg(n_genes);
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:total_loss) schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            LogSpaceGGParams params;
            params.log_alpha = log_alpha[g];
            params.log_beta = log_beta[g];
            params.log_gamma = log_gamma[g];
            
            double sum_grad_la = 0, sum_grad_lb = 0, sum_grad_lg = 0;
            double sum_loss = 0;
            
            for (int j = 0; j < n_samples; j++) {
                double xgj = X(g, j);
                sum_loss += params.loglik(xgj, eps);
                double gla, glb, glg;
                params.gradient(xgj, gla, glb, glg, eps);
                sum_grad_la += gla;
                sum_grad_lb += glb;
                sum_grad_lg += glg;
            }
            
            grad_la[g] = -sum_grad_la / n_samples;
            grad_lb[g] = -sum_grad_lb / n_samples;
            grad_lg[g] = -sum_grad_lg / n_samples;
            
            if (use_firth) {
                double a = params.alpha(), b = params.beta(), gam = params.gamma();
                double h = 1e-4;
                double firth_scale = 0.5 / n_samples;
                
                double ldf_ap = compute_log_det_fisher_unit(a * std::exp(h), b, gam, eps);
                double ldf_am = compute_log_det_fisher_unit(a * std::exp(-h), b, gam, eps);
                grad_la[g] -= firth_scale * (ldf_ap - ldf_am) / (2.0 * h);
                
                double ldf_bp = compute_log_det_fisher_unit(a, b * std::exp(h), gam, eps);
                double ldf_bm = compute_log_det_fisher_unit(a, b * std::exp(-h), gam, eps);
                grad_lb[g] -= firth_scale * (ldf_bp - ldf_bm) / (2.0 * h);
                
                if (!fix_gamma) {
                    double ldf_gp = compute_log_det_fisher_unit(a, b, gam * std::exp(h), eps);
                    double ldf_gm = compute_log_det_fisher_unit(a, b, gam * std::exp(-h), eps);
                    grad_lg[g] -= firth_scale * (ldf_gp - ldf_gm) / (2.0 * h);
                }
            }
            
            if (fix_gamma) grad_lg[g] = 0.0;
            
            double prior_gla, prior_glb, prior_glg;
            prior.prior_gradient(log_alpha[g], log_beta[g], log_gamma[g],
                                prior_gla, prior_glb, prior_glg);
            grad_la[g] += prior_weight * prior_gla;
            grad_lb[g] += prior_weight * prior_glb;
            grad_lg[g] += prior_weight * prior_glg;
            
            if (cfg.use_natural_grad) {
                FisherMatrix3x3 F = compute_fisher_gg(
                    params.alpha(), params.beta(), params.gamma(), n_samples);
                double a = params.alpha();
                double b = params.beta();
                double gam = params.gamma();
                double inv_aa, inv_bb, inv_gg, inv_ab;
                F.compute_full_inverse(inv_aa, inv_bb, inv_gg, inv_ab);
                double grad_a = grad_la[g] / a;
                double grad_b = grad_lb[g] / b;
                double grad_g = grad_lg[g] / gam;
                double nat_grad_a = inv_aa * grad_a + inv_ab * grad_b;
                double nat_grad_b = inv_ab * grad_a + inv_bb * grad_b;
                double nat_grad_g = inv_gg * grad_g;
                grad_la[g] = nat_grad_a / a;
                grad_lb[g] = nat_grad_b / b;
                grad_lg[g] = nat_grad_g / gam;
            }
            
            total_loss -= sum_loss;
            total_loss += prior_weight * prior.neg_log_prior(log_alpha[g], log_beta[g], log_gamma[g]);
        }
        
        loss_history[iter] = total_loss / n_genes;
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            state.m_la[g] = cfg.beta1 * state.m_la[g] + (1.0 - cfg.beta1) * grad_la[g];
            state.m_lb[g] = cfg.beta1 * state.m_lb[g] + (1.0 - cfg.beta1) * grad_lb[g];
            state.m_lg[g] = cfg.beta1 * state.m_lg[g] + (1.0 - cfg.beta1) * grad_lg[g];
            log_alpha[g] -= cfg.lr * state.m_la[g];
            log_beta[g] -= cfg.lr * state.m_lb[g];
            log_gamma[g] -= cfg.lr * state.m_lg[g];
            log_alpha[g] = std::max(-2.3, std::min(4.6, log_alpha[g]));
            log_beta[g] = std::max(-10.0, std::min(15.0, log_beta[g]));
            log_gamma[g] = std::max(-0.7, std::min(1.6, log_gamma[g]));
            if (fix_gamma) log_gamma[g] = 0.0;
        }
        // Do not re-estimate prior — use shared prior
    }
    
    NumericVector alpha_out(n_genes), beta_out(n_genes), gamma_out(n_genes);
    for (int g = 0; g < n_genes; g++) {
        alpha_out[g] = std::exp(log_alpha[g]);
        beta_out[g] = std::exp(log_beta[g]);
        gamma_out[g] = std::exp(log_gamma[g]);
    }
    
    return List::create(
        Named("alpha") = alpha_out,
        Named("beta") = beta_out,
        Named("gamma") = gamma_out,
        Named("loss_history") = loss_history
    );
}

// =============================================================================
// Theoretical variance of the GG distribution
// Var[X] = β² * (Γ(α + 2/γ)/Γ(α) - (Γ(α + 1/γ)/Γ(α))²)
// =============================================================================
inline double gg_theoretical_variance(double alpha, double beta, double gamma, double eps = 1e-8) {
    // Use lgamma to avoid overflow
    double log_gamma_a = fast_special::lgamma(alpha);
    double log_gamma_a_1g = fast_special::lgamma(alpha + 1.0/gamma);
    double log_gamma_a_2g = fast_special::lgamma(alpha + 2.0/gamma);
    
    // E[X] = β * Γ(α + 1/γ) / Γ(α)
    double E_X = beta * std::exp(log_gamma_a_1g - log_gamma_a);
    
    // E[X²] = β² * Γ(α + 2/γ) / Γ(α)
    double E_X2 = beta * beta * std::exp(log_gamma_a_2g - log_gamma_a);
    
    // Var[X] = E[X²] - E[X]²
    double var = E_X2 - E_X * E_X;
    return std::max(var, eps);
}

// =============================================================================
// fitFDist: limma-style scaled F-distribution hyperparameter estimation
// Given gene-wise sample variances s²_g (with df1 degrees of freedom), estimate prior (d0, s0²)
// Model: s²_g ~ s0² · F(df1, d0)
// Reference: Smyth (2004) Stat. Appl. Genet. Mol. Biol.
// =============================================================================
struct FDistFit {
    double scale;   // s0² (prior variance)
    double df2;     // d0  (prior degrees of freedom)
};

FDistFit fit_f_dist(const NumericVector& vars, int df1, double eps = 1e-8) {
    FDistFit result;
    const int n = vars.size();
    
    // 1. Take log and winsorize (1%-99%) to remove extremes
    std::vector<double> log_var(n);
    for (int i = 0; i < n; i++) {
        log_var[i] = std::log(std::max((double)vars[i], eps));
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
    
    // 2. Remove known df1 contribution: e_g = log(s²_g) - ψ(df1/2) + log(df1/2)
    double half_df1 = std::max(df1 / 2.0, 0.5);
    double psi_hd1 = fast_special::digamma(half_df1);
    double log_hd1 = std::log(half_df1);
    
    double sum_e = 0;
    for (int i = 0; i < n; i++) {
        sum_e += log_var[i] - psi_hd1 + log_hd1;
    }
    double emean = sum_e / n;
    
    // 3. Variance equation: Var[e] = ψ'(df1/2) + ψ'(d0/2)
    double ss_e = 0;
    for (int i = 0; i < n; i++) {
        double e_i = log_var[i] - psi_hd1 + log_hd1;
        double diff = e_i - emean;
        ss_e += diff * diff;
    }
    double evar = ss_e / std::max(n - 1, 1) - fast_special::trigamma(half_df1);
    
    if (evar <= 0) {
        // No prior variance structure → d0 → ∞ (maximum shrinkage)
        result.df2 = 1e6;
        result.scale = std::exp(emean);
    } else {
        // 4. Estimate d0 via trigamma inverse
        double half_df2 = fast_special::trigamma_inverse(evar);
        result.df2 = 2.0 * half_df2;
        result.df2 = std::max(result.df2, 0.5);
        result.df2 = std::min(result.df2, 1e6);
        
        // 5. Estimate s0²
        half_df2 = result.df2 / 2.0;
        double psi_hd2 = fast_special::digamma(half_df2);
        double log_hd2 = std::log(half_df2);
        result.scale = std::exp(emean + psi_hd2 - log_hd2);
    }
    
    result.scale = std::max(result.scale, eps);
    return result;
}

// =============================================================================
// GG-informed variance shrinkage (using GG theoretical variance as prior target)
// Resolves the "split-brain" problem: shares information between GG parameter estimates and t-test variances
// =============================================================================
// [[Rcpp::export]]
List adaptive_squeeze_var_gg_informed_cpp(NumericVector vars,
                                           NumericVector gg_alpha,
                                           NumericVector gg_beta,
                                           NumericVector gg_gamma,
                                           int n_samples,
                                           double eps = 1e-8) {
    const int n = vars.size();
    
    // Compute per-gene GG theoretical variance as prior target
    NumericVector prior_var_gg(n);
    double sum_log_prior = 0;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:sum_log_prior)
    #endif
    for (int i = 0; i < n; i++) {
        prior_var_gg[i] = gg_theoretical_variance(gg_alpha[i], gg_beta[i], gg_gamma[i], eps);
        sum_log_prior += std::log(prior_var_gg[i] + eps);
    }
    
    // Estimate dispersion of GG theoretical variances (determines shrinkage strength)
    double mu_log_prior = sum_log_prior / n;
    double ss_log_prior = 0;
    for (int i = 0; i < n; i++) {
        double diff = std::log(prior_var_gg[i] + eps) - mu_log_prior;
        ss_log_prior += diff * diff;
    }
    double sigma_log_prior = std::sqrt(ss_log_prior / n + eps);
    
    // Adaptive shrinkage: sample variance -> GG theoretical variance
    // Shrinkage strength depends on: (1) sample size (2) GG fit reliability
    double prior_precision = 1.0 / (sigma_log_prior * sigma_log_prior + eps);
    double data_precision = (n_samples - 1.0) / 2.0;
    double base_shrink = prior_precision / (prior_precision + data_precision + eps);
    
    NumericVector var_shrunk(n);
    NumericVector shrink_factor(n);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < n; i++) {
        // Gene-specific shrinkage: the closer sample variance is to GG theoretical variance, the more we trust GG
        double log_ratio = std::abs(std::log(vars[i] + eps) - std::log(prior_var_gg[i] + eps));
        double agreement = std::exp(-log_ratio);  // 0-1, closer to 1 indicates better GG fit
        
        // Shrinkage factor = base shrinkage * agreement adjustment
        double gene_shrink = base_shrink * (0.5 + 0.5 * agreement);
        
        // Shrunk variance = (1-w) * sample variance + w * GG theoretical variance
        var_shrunk[i] = (1.0 - gene_shrink) * vars[i] + gene_shrink * prior_var_gg[i];
        shrink_factor[i] = gene_shrink;
    }
    
    // Estimate prior degrees of freedom
    double df_prior = 2.0 * prior_precision;
    double df_total = df_prior + n_samples - 1.0;
    
    return List::create(
        Named("var_shrunk") = var_shrunk,
        Named("shrink_factor") = shrink_factor,
        Named("prior_var_gg") = prior_var_gg,
        Named("df_prior") = df_prior,
        Named("df_total") = df_total
    );
}

// =============================================================================
// limma-style variance shrinkage (squeezeVar): fitFDist replaces ad-hoc hierarchical formula
// Model: s²_g ~ s0² · F(df_residual, d0), posterior = (d0·s0² + df·s²_g)/(d0+df)
// Reference: Smyth (2004) Stat. Appl. Genet. Mol. Biol., limma::squeezeVar()
// =============================================================================
// [[Rcpp::export]]
List adaptive_squeeze_var_hierarchical_cpp(NumericVector vars,
                                            NumericVector gene_means,
                                            int n_samples,
                                            double eps = 1e-8) {
    const int n = vars.size();
    int df_residual = std::max(n_samples - 2, 1);
    
    // =====================================================================
    // Step 1: compute variance trend (correct order for limma-trend)
    // Must be done before fitFDist; otherwise trend inflates Var[log(s²)] → d0 is suppressed
    // Reference: internal logic of limma::eBayes(trend=TRUE)
    // =====================================================================
    const int N_BINS = std::min(50, std::max(5, n / 50));
    
    std::vector<double> log_bm(n);
    for (int i = 0; i < n; i++) {
        log_bm[i] = std::log(std::max((double)gene_means[i], eps));
    }
    
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return log_bm[a] < log_bm[b]; });
    
    std::vector<double> bin_center(N_BINS), bin_value(N_BINS);
    int genes_per_bin = n / N_BINS;
    for (int b = 0; b < N_BINS; b++) {
        int bstart = b * genes_per_bin;
        int bend = (b == N_BINS - 1) ? n : (b + 1) * genes_per_bin;
        int bsize = bend - bstart;
        
        double sum_lbm = 0;
        std::vector<double> bvars(bsize);
        for (int i = 0; i < bsize; i++) {
            int g = order[bstart + i];
            sum_lbm += log_bm[g];
            bvars[i] = std::max((double)vars[g], eps);
        }
        bin_center[b] = sum_lbm / bsize;
        
        int bmid = bsize / 2;
        std::nth_element(bvars.begin(), bvars.begin() + bmid, bvars.end());
        bin_value[b] = bvars[bmid];
    }
    
    // Linear interpolation to obtain gene-specific variance trend
    std::vector<double> var_trend(n);
    for (int i = 0; i < n; i++) {
        double x = log_bm[i];
        if (x <= bin_center[0]) {
            var_trend[i] = bin_value[0];
        } else if (x >= bin_center[N_BINS - 1]) {
            var_trend[i] = bin_value[N_BINS - 1];
        } else {
            int lo = 0, hi = N_BINS - 1;
            while (hi - lo > 1) {
                int m = (lo + hi) / 2;
                if (bin_center[m] <= x) lo = m;
                else hi = m;
            }
            double frac = (x - bin_center[lo]) / (bin_center[hi] - bin_center[lo] + eps);
            var_trend[i] = (1.0 - frac) * bin_value[lo] + frac * bin_value[hi];
        }
        var_trend[i] = std::max(var_trend[i], eps);
    }
    
    // =====================================================================
    // Step 2: detrend → fitFDist on detrended variances
    // detrended_var[g] = vars[g] / trend[g], should approximate s0² · F(df, d0)
    // After detrending, Var[log(detrended)] is no longer inflated by mean-var relationship → d0 is correct
    // =====================================================================
    NumericVector detrended_vars(n);
    for (int i = 0; i < n; i++) {
        detrended_vars[i] = std::max((double)vars[i], eps) / var_trend[i];
    }
    
    FDistFit fdist = fit_f_dist(detrended_vars, df_residual, eps);
    
    double d0 = fdist.df2;
    double scale = fdist.scale;   // s0² on the detrended scale (should be close to 1.0)
    double df_total = d0 + df_residual;
    double sw = d0 / df_total;
    
    // =====================================================================
    // Step 3: gene-specific prior = trend[g] * scale
    // s0²_global for output (median of trend * scale)
    // =====================================================================
    NumericVector var_prior_trend(n);
    std::vector<double> trend_vals(n);
    for (int i = 0; i < n; i++) {
        var_prior_trend[i] = var_trend[i] * scale;
        trend_vals[i] = var_prior_trend[i];
    }
    
    int med_idx = n / 2;
    std::nth_element(trend_vals.begin(), trend_vals.begin() + med_idx, trend_vals.end());
    double s0_sq = trend_vals[med_idx];
    
    // Posterior variance: (d0·s0_trend_g + df·s²_g) / (d0 + df)
    NumericVector var_shrunk(n);
    NumericVector shrink_factor(n);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < n; i++) {
        var_shrunk[i] = (d0 * var_prior_trend[i] + df_residual * std::max((double)vars[i], eps)) / df_total;
        shrink_factor[i] = sw;
    }
    
    double mu_log_var = std::log(std::max(s0_sq, eps));
    double sigma_log_var = (d0 > eps) ? std::sqrt(fast_special::trigamma(d0 / 2.0)) : 0.0;
    
    return List::create(
        Named("var_shrunk") = var_shrunk,
        Named("shrink_factor") = shrink_factor,
        Named("prior_var") = s0_sq,
        Named("var_prior_trend") = var_prior_trend,
        Named("df_prior") = d0,
        Named("df_total") = df_total,
        Named("mu_log_var") = mu_log_var,
        Named("sigma_log_var") = sigma_log_var
    );
}

// Rao's geodesic distance on the GG manifold
// Captures the overall displacement of (α,β,γ) triplets, avoiding over-normalization
// =============================================================================
inline double gg_geodesic_distance(double a1, double b1, double g1,
                                    double a2, double b2, double g2,
                                    double eps = 1e-8) {
    // Approximation of Fisher-Rao distance: sqrt(Δθ^T F Δθ)
    // For the GG distribution, uses a weighted combination of log-space distances
    double a1p = std::max(a1, eps), a2p = std::max(a2, eps);
    double b1p = std::max(b1, eps), b2p = std::max(b2, eps);
    double g1p = std::max(g1, eps), g2p = std::max(g2, eps);

    double log_a1 = std::log(a1p), log_a2 = std::log(a2p);
    double log_b1 = std::log(b1p), log_b2 = std::log(b2p);
    double log_g1 = std::log(g1p), log_g2 = std::log(g2p);

    double delta_a = log_a2 - log_a1;
    double delta_b = log_b2 - log_b1;
    double delta_g = log_g2 - log_g1;

    double a_mid = std::sqrt(a1p * a2p);
    double b_mid = std::sqrt(b1p * b2p);
    double g_mid = std::sqrt(g1p * g2p);

    FisherMatrix3x3 F = compute_fisher_gg(a_mid, b_mid, g_mid, 1);

    double L_aa = F.I_aa * a_mid * a_mid;
    double L_bb = F.I_bb * b_mid * b_mid;
    double L_gg = F.I_gg * g_mid * g_mid;
    double L_ab = F.I_ab * a_mid * b_mid;
    double L_ag = F.I_ag * a_mid * g_mid;
    double L_bg = F.I_bg * b_mid * g_mid;

    double dist_sq = L_aa * delta_a * delta_a +
                     L_bb * delta_b * delta_b +
                     L_gg * delta_g * delta_g +
                     2.0 * L_ab * delta_a * delta_b +
                     2.0 * L_ag * delta_a * delta_g +
                     2.0 * L_bg * delta_b * delta_g;

    return std::sqrt(std::max(dist_sq, eps));
}

 inline void estimate_scaled_chisq_df_scale_from_iqr(double chi_q25, double chi_q75,
                                                     double& df_hat, double& scale_hat,
                                                     double eps = 1e-8) {
     double ratio_obs = (chi_q75 + eps) / std::max(chi_q25, eps);
     double lo = 0.5;
     double hi = 30.0;

     for (int it = 0; it < 60; it++) {
         double mid = 0.5 * (lo + hi);
         double ref_q25 = R::qchisq(0.25, mid, 1, 0);
         double ref_q75 = R::qchisq(0.75, mid, 1, 0);
         double ratio_mid = ref_q75 / std::max(ref_q25, eps);

         if (ratio_mid > ratio_obs) {
             lo = mid;
         } else {
             hi = mid;
         }
     }

     df_hat = 0.5 * (lo + hi);

     double chi_ref_q25 = R::qchisq(0.25, df_hat, 1, 0);
     double chi_ref_q75 = R::qchisq(0.75, df_hat, 1, 0);
     double denom = std::max(chi_ref_q75 - chi_ref_q25, eps);
     scale_hat = (chi_q75 - chi_q25) / denom;
     scale_hat = std::max(scale_hat, eps);
 }

// Manifold distance test
// Converts GG parameter differences into p-values, capturing "steady-state dysregulation"
// =============================================================================
// [[Rcpp::export]]
List manifold_distance_test_cpp(NumericVector alpha_ctrl, NumericVector beta_ctrl, 
                                 NumericVector gamma_ctrl,
                                 NumericVector alpha_treat, NumericVector beta_treat,
                                 NumericVector gamma_treat,
                                 int n_ctrl, int n_treat,
                                 double eps = 1e-8) {
    const int n_genes = alpha_ctrl.size();
    
    NumericVector geodesic_dist(n_genes);
    NumericVector pvalue_manifold(n_genes);
    std::vector<double> dist_vec(n_genes);
    std::vector<double> chi_ab_vec(n_genes);
    std::vector<double> chi_g_vec(n_genes);
    std::vector<double> chi_sq_vec(n_genes);
    
    // Estimate distance distribution under the null hypothesis (chi-squared approximation)
    // Under H0, distance² ~ χ²(3) / n (3 parameters)
    double n_eff = 2.0 * n_ctrl * n_treat / (n_ctrl + n_treat);  // effective sample size
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ac = std::max(alpha_ctrl[g], eps);
        double bc = std::max(beta_ctrl[g], eps);
        double gc = std::max(gamma_ctrl[g], eps);
        double at = std::max(alpha_treat[g], eps);
        double bt = std::max(beta_treat[g], eps);
        double gt = std::max(gamma_treat[g], eps);

        dist_vec[g] = gg_geodesic_distance(
            ac, bc, gc,
            at, bt, gt, eps);

        double d_a = std::log(at) - std::log(ac);
        double d_b = std::log(bt) - std::log(bc);
        double d_g = std::log(gt) - std::log(gc);

        FisherMatrix3x3 Fc = compute_fisher_gg(ac, bc, gc, n_ctrl);
        FisherMatrix3x3 Ft = compute_fisher_gg(at, bt, gt, n_treat);

        const double reg = 1e-12;

        double Fca = Fc.I_aa * ac * ac;
        double Fcb = Fc.I_bb * bc * bc;
        double Fcab = Fc.I_ab * ac * bc;
        double trFc = Fca + Fcb;
        double ridgeFc = 1e-12 * std::max(trFc, 0.0) + reg;
        double Fca_r = Fca + ridgeFc;
        double Fcb_r = Fcb + ridgeFc;
        double detFc = Fca_r * Fcb_r - Fcab * Fcab;
        detFc = std::max(detFc, reg);
        double invc_aa = Fcb_r / detFc;
        double invc_bb = Fca_r / detFc;
        double invc_ab = -Fcab / detFc;

        double Fta = Ft.I_aa * at * at;
        double Ftb = Ft.I_bb * bt * bt;
        double Ftab = Ft.I_ab * at * bt;
        double trFt = Fta + Ftb;
        double ridgeFt = 1e-12 * std::max(trFt, 0.0) + reg;
        double Fta_r = Fta + ridgeFt;
        double Ftb_r = Ftb + ridgeFt;
        double detFt = Fta_r * Ftb_r - Ftab * Ftab;
        detFt = std::max(detFt, reg);
        double invt_aa = Ftb_r / detFt;
        double invt_bb = Fta_r / detFt;
        double invt_ab = -Ftab / detFt;

        double Vaa = invc_aa + invt_aa;
        double Vab = invc_ab + invt_ab;
        double Vbb = invc_bb + invt_bb;
        double trV = Vaa + Vbb;
        double ridgeV = 1e-12 * std::max(trV, 0.0) + reg;
        double Vaa_r = Vaa + ridgeV;
        double Vbb_r = Vbb + ridgeV;
        double detV = Vaa_r * Vbb_r - Vab * Vab;
        detV = std::max(detV, reg);
        double Paa = Vbb_r / detV;
        double Pbb = Vaa_r / detV;
        double Pab = -Vab / detV;
        double chi_ab = d_a * (Paa * d_a + Pab * d_b) + d_b * (Pab * d_a + Pbb * d_b);

        double varc_g = 1.0 / std::max(Fc.I_gg * gc * gc, reg);
        double vart_g = 1.0 / std::max(Ft.I_gg * gt * gt, reg);
        double var_g = varc_g + vart_g;
        double chi_g = (d_g * d_g) / std::max(var_g, eps);

        chi_ab_vec[g] = std::max(chi_ab, 0.0);
        chi_g_vec[g] = std::max(chi_g, 0.0);
        chi_sq_vec[g] = std::max(chi_ab + chi_g, 0.0);
    }

    for (int g = 0; g < n_genes; g++) geodesic_dist[g] = dist_vec[g];

    std::vector<double> chi_tmp(n_genes);
    for (int g = 0; g < n_genes; g++) chi_tmp[g] = chi_sq_vec[g];

    int idx25 = n_genes / 4;
    int idx75 = (3 * n_genes) / 4;

    std::nth_element(chi_tmp.begin(), chi_tmp.begin() + idx25, chi_tmp.end());
    double chi_q25 = chi_tmp[idx25];

    std::nth_element(chi_tmp.begin(), chi_tmp.begin() + idx75, chi_tmp.end());
    double chi_q75 = chi_tmp[idx75];

    // IQR-based df and scale estimation (Genomic Control style)
    double df_hat, scale_hat;
    estimate_scaled_chisq_df_scale_from_iqr(chi_q25, chi_q75, df_hat, scale_hat, eps);
    // Fallback: if estimation gives degenerate values, use df=3 with median calibration
    if (df_hat < 0.5 || df_hat > 30.0 || scale_hat < eps) {
        df_hat = 3.0;
        std::vector<double> chi_tmp2(chi_tmp);
        int idx50 = n_genes / 2;
        std::nth_element(chi_tmp2.begin(), chi_tmp2.begin() + idx50, chi_tmp2.end());
        double chi_q50 = chi_tmp2[idx50];
        double chi_ref_q50 = R::qchisq(0.5, df_hat, 1, 0);
        scale_hat = std::max(chi_q50 / std::max(chi_ref_q50, eps), eps);
    }

    for (int g = 0; g < n_genes; g++) {
        double chi_adj = chi_sq_vec[g] / scale_hat;
        pvalue_manifold[g] = R::pchisq(chi_adj, df_hat, 0, 0);
        pvalue_manifold[g] = std::max(eps, std::min(1.0, pvalue_manifold[g]));
    }
    
    return List::create(
        Named("geodesic_dist") = geodesic_dist,
        Named("pvalue_manifold") = pvalue_manifold,
        Named("df") = df_hat,
        Named("scale") = scale_hat,
        Named("n_eff") = n_eff
    );
}

// =============================================================================
// Small dense linear system solver (Gaussian elimination with partial pivoting)
// A: row-major k×k, b: k×1; A and b are modified in-place, b stores the solution
// =============================================================================
inline bool solve_dense_inplace(std::vector<double>& A, std::vector<double>& b, int k) {
    for (int j = 0; j < k; j++) {
        int piv = j;
        double maxval = std::abs(A[j*k+j]);
        for (int i = j+1; i < k; i++) {
            if (std::abs(A[i*k+j]) > maxval) { maxval = std::abs(A[i*k+j]); piv = i; }
        }
        if (piv != j) {
            for (int c = 0; c < k; c++) std::swap(A[j*k+c], A[piv*k+c]);
            std::swap(b[j], b[piv]);
        }
        double diag = A[j*k+j];
        if (std::abs(diag) < 1e-15) return false;
        for (int i = j+1; i < k; i++) {
            double f = A[i*k+j] / diag;
            for (int c = j; c < k; c++) A[i*k+c] -= f * A[j*k+c];
            b[i] -= f * b[j];
        }
    }
    for (int j = k-1; j >= 0; j--) {
        for (int c = j+1; c < k; c++) b[j] -= A[j*k+c] * b[c];
        b[j] /= A[j*k+j];
    }
    return true;
}

// =============================================================================
// Route A: WLS + Empirical Bayes moderated t-test (general design matrix)
//
// y:        genes × samples (log-expression, e.g. logCPM)
// design:   samples × n_coef
// contrast: n_coef (contrast vector specifying the test)
// weights:  genes × samples (precision weights, optional; 0-length = uniform weights)
// gene_means: genes (for trend-based EB shrinkage)
//
// Math: β̂_g = (X'W_gX)⁻¹ X'W_g y_g
//       s²_g = RSS_g / df_resid
//       s̃²_g = (d₀s₀² + df·s²_g) / (d₀+df)   [Smyth 2004]
//       t_g = c'β̂_g / sqrt(s̃²_g · c'(X'W_gX)⁻¹c)
// =============================================================================
// [[Rcpp::export]]
List sgcb_wls_ebayes_cpp(NumericMatrix y,
                          NumericMatrix design,
                          NumericVector contrast,
                          NumericVector gene_means,
                          NumericVector weights,
                          double eps = 1e-8) {
    const int n_genes = y.nrow();
    const int n_samples = y.ncol();
    const int k = design.ncol();
    const int df_residual = std::max(n_samples - k, 1);
    const bool has_weights = (weights.size() == n_genes * n_samples);
    
    NumericVector log2FC_out(n_genes), sigma2_out(n_genes), var_unscaled(n_genes);
    
    // If no weights, precompute (X'X)⁻¹ and (X'X)⁻¹·c (shared across all genes)
    // XtX_inv: k×k row-major
    std::vector<double> XtX_inv_shared(k*k, 0.0);
    double vc_shared = 0.0;  // c'(X'X)⁻¹c
    std::vector<double> XtX_inv_c_shared(k, 0.0);
    
    if (!has_weights) {
        std::vector<double> XtX(k*k, 0.0);
        for (int j = 0; j < n_samples; j++) {
            for (int a = 0; a < k; a++) {
                for (int b = a; b < k; b++) {
                    XtX[a*k+b] += design(j,a) * design(j,b);
                }
            }
        }
        for (int a = 0; a < k; a++)
            for (int b = 0; b < a; b++)
                XtX[a*k+b] = XtX[b*k+a];
        
        // Inversion: solve XtX · col_i = e_i for each i
        for (int col = 0; col < k; col++) {
            std::vector<double> A_copy(XtX);
            std::vector<double> e(k, 0.0);
            e[col] = 1.0;
            solve_dense_inplace(A_copy, e, k);
            for (int r = 0; r < k; r++) XtX_inv_shared[r*k+col] = e[r];
        }
        
        // c'(X'X)⁻¹c
        for (int a = 0; a < k; a++) {
            double tmp = 0;
            for (int b = 0; b < k; b++) tmp += XtX_inv_shared[a*k+b] * contrast[b];
            XtX_inv_c_shared[a] = tmp;
        }
        for (int a = 0; a < k; a++) vc_shared += contrast[a] * XtX_inv_c_shared[a];
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        std::vector<double> beta(k, 0.0);
        double vc_g = 0.0;
        
        if (!has_weights) {
            // β̂ = (X'X)⁻¹ X'y
            std::vector<double> Xty(k, 0.0);
            for (int j = 0; j < n_samples; j++) {
                double yj = y(g, j);
                for (int a = 0; a < k; a++) Xty[a] += design(j, a) * yj;
            }
            for (int a = 0; a < k; a++) {
                double s = 0;
                for (int b = 0; b < k; b++) s += XtX_inv_shared[a*k+b] * Xty[b];
                beta[a] = s;
            }
            vc_g = vc_shared;
        } else {
            // Per-gene weighted: β̂ = (X'WX)⁻¹ X'Wy
            std::vector<double> XtWX(k*k, 0.0);
            std::vector<double> XtWy(k, 0.0);
            for (int j = 0; j < n_samples; j++) {
                double w = weights[g + j * n_genes];  // column-major
                double yj = y(g, j);
                for (int a = 0; a < k; a++) {
                    double xw = design(j, a) * w;
                    XtWy[a] += xw * yj;
                    for (int b = a; b < k; b++) {
                        XtWX[a*k+b] += xw * design(j, b);
                    }
                }
            }
            for (int a = 0; a < k; a++)
                for (int b = 0; b < a; b++)
                    XtWX[a*k+b] = XtWX[b*k+a];
            
            // Invert and compute β̂ and c'(X'WX)⁻¹c
            // First solve (X'WX) β̂ = X'Wy
            std::vector<double> A_work(XtWX);
            beta = XtWy;
            solve_dense_inplace(A_work, beta, k);
            
            // Then solve (X'WX) z = c → vc_g = c'z
            std::vector<double> A_work2(XtWX);
            std::vector<double> c_copy(k);
            for (int a = 0; a < k; a++) c_copy[a] = contrast[a];
            solve_dense_inplace(A_work2, c_copy, k);
            for (int a = 0; a < k; a++) vc_g += contrast[a] * c_copy[a];
        }
        
        // log2FC = c'β̂
        double lfc = 0;
        for (int a = 0; a < k; a++) lfc += contrast[a] * beta[a];
        log2FC_out[g] = lfc;
        
        // Residual variance s²
        double rss = 0;
        for (int j = 0; j < n_samples; j++) {
            double fitted = 0;
            for (int a = 0; a < k; a++) fitted += design(j, a) * beta[a];
            double resid = y(g, j) - fitted;
            double w = has_weights ? weights[g + j * n_genes] : 1.0;
            rss += w * resid * resid;
        }
        sigma2_out[g] = std::max(rss / df_residual, eps);
        var_unscaled[g] = vc_g;
    }
    
    // EB shrinkage (trend-based)
    List shrink = adaptive_squeeze_var_hierarchical_cpp(sigma2_out, gene_means, n_samples, eps);
    double d0 = as<double>(shrink["df_prior"]);
    double df_total = d0 + df_residual;
    NumericVector var_shrunk_out = shrink["var_shrunk"];
    NumericVector var_prior_trend = shrink["var_prior_trend"];
    
    // Moderated t-test
    NumericVector t_stat_out(n_genes), pvalue_out(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double se = std::sqrt(std::max(var_shrunk_out[g] * var_unscaled[g], eps));
        t_stat_out[g] = log2FC_out[g] / (se + eps);
        pvalue_out[g] = 2.0 * R::pt(-std::abs(t_stat_out[g]), df_total, 1, 0);
        pvalue_out[g] = std::max(eps, std::min(1.0, pvalue_out[g]));
    }
    
    return List::create(
        Named("log2FC") = log2FC_out,
        Named("sigma2") = sigma2_out,
        Named("var_shrunk") = var_shrunk_out,
        Named("var_prior_trend") = var_prior_trend,
        Named("var_unscaled") = var_unscaled,
        Named("t_stat") = t_stat_out,
        Named("pvalue") = pvalue_out,
        Named("df_prior") = d0,
        Named("df_total") = df_total,
        Named("df_residual") = df_residual
    );
}

// =============================================================================
// Complete information-geometric SGCB inference
// =============================================================================
// [[Rcpp::export]]
List sgcb_info_geom_inference_cpp(NumericMatrix X,
                                   IntegerVector group,
                                   int n_dropout = 50,
                                   double dropout_rate = 0.2,
                                   bool use_natural_grad = true,
                                   bool use_hierarchical = true,
                                   double eps = 1e-8) {
    const int n_genes = X.nrow();
    const int n_samples = X.ncol();
    
    // Group assignment
    std::vector<int> ctrl_idx, treat_idx;
    for (int j = 0; j < n_samples; j++) {
        if (group[j] == 0) ctrl_idx.push_back(j);
        else treat_idx.push_back(j);
    }
    const int n_ctrl = ctrl_idx.size();
    const int n_treat = treat_idx.size();
    const int n_min = std::min(n_ctrl, n_treat);
    
    // Extract control group data
    NumericMatrix ctrl_mat(n_genes, n_ctrl);
    for (int g = 0; g < n_genes; g++) {
        for (int j = 0; j < n_ctrl; j++) {
            ctrl_mat(g, j) = X(g, ctrl_idx[j]);
        }
    }
    
    // Extract treatment group data
    NumericMatrix treat_mat(n_genes, n_treat);
    for (int g = 0; g < n_genes; g++) {
        for (int j = 0; j < n_treat; j++) {
            treat_mat(g, j) = X(g, treat_idx[j]);
        }
    }
    
    // Fix A: estimate shared hyperprior from pooled data, used by both groups
    // Avoids ctrl/treat estimating separate priors, which causes DE contamination → implicit batch effect
    List gg_ctrl, gg_treat;
    if (use_hierarchical) {
        // 1) Pooled moment-of-methods initialization → shared prior
        std::vector<double> pool_la(n_genes), pool_lb(n_genes), pool_lg(n_genes);
        bool pool_fix_gamma = (n_min <= 5);
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            double sum_x = 0, sum_x2 = 0;
            for (int j = 0; j < n_samples; j++) {
                double xij = X(g, j) + eps;
                sum_x += xij;
                sum_x2 += xij * xij;
            }
            double mean_x = sum_x / n_samples;
            double var_x = sum_x2 / n_samples - mean_x * mean_x;
            double cv2 = var_x / (mean_x * mean_x + eps);
            pool_la[g] = std::log(std::max(0.5, std::min(10.0, 1.0 / cv2)));
            pool_lg[g] = pool_fix_gamma ? 0.0 : std::log(std::max(0.5, std::min(3.0, 1.0 / std::sqrt(cv2 + eps))));
            pool_lb[g] = std::log(std::max(eps, mean_x));
        }
        
        HierarchicalPrior shared_prior;
        shared_prior.estimate_from_data(pool_la, pool_lb, pool_lg);
        
        // 2) Per-group MAP fitting using shared prior
        gg_ctrl = fit_gg_with_shared_prior(ctrl_mat, shared_prior, use_natural_grad, eps);
        gg_treat = fit_gg_with_shared_prior(treat_mat, shared_prior, use_natural_grad, eps);
    } else {
        gg_ctrl = fit_gg_natural_gradient_cpp(ctrl_mat, -1, -1, -1, use_natural_grad, eps);
        gg_treat = fit_gg_natural_gradient_cpp(treat_mat, -1, -1, -1, use_natural_grad, eps);
    }
    
    NumericVector ae_alpha = gg_ctrl["alpha"];
    NumericVector ae_beta = gg_ctrl["beta"];
    NumericVector ae_gamma = gg_ctrl["gamma"];
    
    NumericVector treat_alpha = gg_treat["alpha"];
    NumericVector treat_beta = gg_treat["beta"];
    NumericVector treat_gamma = gg_treat["gamma"];
    
    // Basic statistics (original scale + log2 scale)
    NumericVector mean_ctrl(n_genes), mean_treat(n_genes);
    NumericVector var_ctrl_log2(n_genes), var_treat_log2(n_genes);
    NumericVector log2FC(n_genes);
    NumericVector baseMean(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        // Original-scale means (for baseMean and log2FC)
        double sum_c = 0, sum_t = 0;
        for (int j = 0; j < n_ctrl; j++) sum_c += X(g, ctrl_idx[j]);
        for (int j = 0; j < n_treat; j++) sum_t += X(g, treat_idx[j]);
        mean_ctrl[g] = sum_c / n_ctrl + eps;
        mean_treat[g] = sum_t / n_treat + eps;
        baseMean[g] = (mean_ctrl[g] + mean_treat[g]) / 2;
        
        // log2-scale mean and variance (computed directly in log2 space, consistent with limma-voom)
        double sum_lc = 0, sum_lc2 = 0;
        for (int j = 0; j < n_ctrl; j++) {
            double lv = std::log2(X(g, ctrl_idx[j]) + 0.5);
            sum_lc += lv;
            sum_lc2 += lv * lv;
        }
        double sum_lt = 0, sum_lt2 = 0;
        for (int j = 0; j < n_treat; j++) {
            double lv = std::log2(X(g, treat_idx[j]) + 0.5);
            sum_lt += lv;
            sum_lt2 += lv * lv;
        }
        // Fix D: log2FC computed as mean difference directly in log2 space, same scale as variance
        log2FC[g] = sum_lt / n_treat - sum_lc / n_ctrl;
        var_ctrl_log2[g] = (sum_lc2 - sum_lc * sum_lc / n_ctrl) / std::max(1, n_ctrl - 1) + eps;
        var_treat_log2[g] = (sum_lt2 - sum_lt * sum_lt / n_treat) / std::max(1, n_treat - 1) + eps;
    }
    
    // =========================================================================
    // GG-Informed Variance Shrinkage v2 (Trend + Fisher df)
    //
    // Fix: old per-gene GG variance was an "intra-gene estimator" (correlated with s_g²/δ̂_g),
    //      violating Smyth (2004)'s independent χ² condition → reference distribution mismatch.
    //
    // Approach 1 (Trend): smooth GG variance along baseMean via binned-median
    //   s_{0,g}² = (1-w)·s₀² + w·v_trend(baseMean_g)
    //   v_trend is a function of an external covariate, approximately independent of s_g² (under normality)
    //   → restores limma-trend's theoretical guarantees
    //
    // Approach 2 (Fisher df): GG prior uncertainty quantified via Fisher information
    //   Var(log v̂) = ∇ᵀ I(θ̂)⁻¹ ∇  (delta method)
    //   d_GG = 2/Var(log v̂)  (χ² CV² correspondence)
    //   d_{0,g} = (1-w)·d₀ + w·d_GG,g  → gene-specific reference df
    //
    // Mathematical basis: for X ~ GG(α,β,γ):
    //   Var[log₂(X)] = ψ'(α) / (γ² · ln²(2))  (closed form)
    // =========================================================================
    
    // Step 1: per-gene GG theoretical variance (log2 scale, closed form)
    const double ln2_sq = std::log(2.0) * std::log(2.0);
    NumericVector var_gg_log2(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ac = std::max(ae_alpha[g], 0.05);
        double gc = std::max(ae_gamma[g], 0.05);
        double at = std::max(treat_alpha[g], 0.05);
        double gt = std::max(treat_gamma[g], 0.05);
        
        double v_ctrl = fast_special::trigamma(ac) / (gc * gc * ln2_sq);
        double v_treat = fast_special::trigamma(at) / (gt * gt * ln2_sq);
        var_gg_log2[g] = std::max((v_ctrl + v_treat) / 2.0, eps);
    }
    
    // Step 2: pooled sample variance + limma fitFDist (global prior s0², d0)
    int df_residual = std::max(n_ctrl + n_treat - 2, 1);
    NumericVector pooled_var_log2(n_genes);
    for (int g = 0; g < n_genes; g++) {
        pooled_var_log2[g] = ((n_ctrl - 1) * var_ctrl_log2[g] + 
                              (n_treat - 1) * var_treat_log2[g]) / df_residual;
    }
    
    List shrink_result = adaptive_squeeze_var_hierarchical_cpp(
        pooled_var_log2, baseMean, n_ctrl + n_treat, eps);
    
    double s0_sq = as<double>(shrink_result["prior_var"]);
    double d0 = as<double>(shrink_result["df_prior"]);
    double df_total = d0 + df_residual;
    NumericVector var_prior_trend = shrink_result["var_prior_trend"];
    
    // =====================================================================
    // Step 3 (Approach 1): Trend-based GG variance
    // Bin by baseMean → take binned median of GG variance → linear interpolation
    // → s_{0,g}² becomes a smooth function of baseMean, approximately independent of s_g²
    // =====================================================================
    const int N_BINS = std::min(50, std::max(5, n_genes / 50));
    
    std::vector<int> order(n_genes);
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> log_bm(n_genes);
    for (int g = 0; g < n_genes; g++) log_bm[g] = std::log(baseMean[g] + eps);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return log_bm[a] < log_bm[b]; });
    
    std::vector<double> bin_center(N_BINS), bin_value(N_BINS);
    int genes_per_bin = n_genes / N_BINS;
    for (int b = 0; b < N_BINS; b++) {
        int bstart = b * genes_per_bin;
        int bend = (b == N_BINS - 1) ? n_genes : (b + 1) * genes_per_bin;
        int bsize = bend - bstart;
        
        double sum_lbm = 0;
        std::vector<double> bgg(bsize);
        for (int i = 0; i < bsize; i++) {
            int g = order[bstart + i];
            sum_lbm += log_bm[g];
            bgg[i] = var_gg_log2[g];
        }
        bin_center[b] = sum_lbm / bsize;
        
        int bmid = bsize / 2;
        std::nth_element(bgg.begin(), bgg.begin() + bmid, bgg.end());
        bin_value[b] = bgg[bmid];
    }
    
    NumericVector var_gg_trend(n_genes);
    for (int g = 0; g < n_genes; g++) {
        double x = log_bm[g];
        if (x <= bin_center[0]) {
            var_gg_trend[g] = bin_value[0];
        } else if (x >= bin_center[N_BINS - 1]) {
            var_gg_trend[g] = bin_value[N_BINS - 1];
        } else {
            int lo = 0, hi = N_BINS - 1;
            while (hi - lo > 1) {
                int m = (lo + hi) / 2;
                if (bin_center[m] <= x) lo = m;
                else hi = m;
            }
            double frac = (x - bin_center[lo]) / (bin_center[hi] - bin_center[lo] + eps);
            var_gg_trend[g] = (1.0 - frac) * bin_value[lo] + frac * bin_value[hi];
        }
        var_gg_trend[g] = std::max(var_gg_trend[g], eps);
    }
    
    {
        std::vector<double> tmp_tr(n_genes), tmp_sp(n_genes);
        for (int g = 0; g < n_genes; g++) {
            tmp_tr[g] = var_gg_trend[g];
            tmp_sp[g] = pooled_var_log2[g];
        }
        int idx_med = n_genes / 2;
        std::nth_element(tmp_tr.begin(), tmp_tr.begin() + idx_med, tmp_tr.end());
        double med_tr = tmp_tr[idx_med];
        std::nth_element(tmp_sp.begin(), tmp_sp.begin() + idx_med, tmp_sp.end());
        double med_sp = tmp_sp[idx_med];
        double trend_scale = med_sp / (med_tr + eps);
        
        for (int g = 0; g < n_genes; g++) {
            var_gg_trend[g] *= trend_scale;
            var_gg_trend[g] = std::max(0.25 * s0_sq, std::min(4.0 * s0_sq, var_gg_trend[g]));
        }
    }
    
    // =====================================================================
    // Step 4 (Approach 2): Fisher-based effective df per gene
    // Delta method: Var(log v̂) = ∇ᵀ I(θ̂)⁻¹ ∇
    //   ∇ = (ψ''(α)/ψ'(α), -2/γ)  (gradient of log v w.r.t. α, γ)
    // Effective df: d_GG = 2/Var(log v̂)  (χ² CV² = 2/ν)
    // =====================================================================
    NumericVector d_gg_eff(n_genes);
    const double h_nd = 1e-5;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ac = std::max(ae_alpha[g], 0.05);
        double gc = std::max(ae_gamma[g], 0.05);
        double bc = std::max((double)ae_beta[g], eps);
        double at = std::max(treat_alpha[g], 0.05);
        double gt = std::max(treat_gamma[g], 0.05);
        double bt = std::max((double)treat_beta[g], eps);
        
        // Ctrl: Fisher 3x3 → full invert → extract (α,γ) entries
        FisherMatrix3x3 F_c = compute_fisher_gg(ac, bc, gc, n_ctrl);
        double c_iaa, c_iab, c_iag, c_ibb, c_ibg, c_igg;
        invert_sym_3x3(F_c.I_aa, F_c.I_ab, F_c.I_ag,
                       F_c.I_bb, F_c.I_bg, F_c.I_gg,
                       c_iaa, c_iab, c_iag, c_ibb, c_ibg, c_igg);
        
        // Treat: same
        FisherMatrix3x3 F_t = compute_fisher_gg(at, bt, gt, n_treat);
        double t_iaa, t_iab, t_iag, t_ibb, t_ibg, t_igg;
        invert_sym_3x3(F_t.I_aa, F_t.I_ab, F_t.I_ag,
                       F_t.I_bb, F_t.I_bg, F_t.I_gg,
                       t_iaa, t_iab, t_iag, t_ibb, t_ibg, t_igg);
        
        // ∂log(v)/∂α = ψ''(α)/ψ'(α), numerical central difference for ψ''
        double psi2_c = (fast_special::trigamma(ac + h_nd) - fast_special::trigamma(ac - h_nd)) / (2.0 * h_nd);
        double ga_c = psi2_c / (fast_special::trigamma(ac) + eps);
        double gg_c = -2.0 / gc;
        double var_lv_c = ga_c * ga_c * c_iaa + gg_c * gg_c * c_igg
                        + 2.0 * ga_c * gg_c * c_iag;
        
        double psi2_t = (fast_special::trigamma(at + h_nd) - fast_special::trigamma(at - h_nd)) / (2.0 * h_nd);
        double ga_t = psi2_t / (fast_special::trigamma(at) + eps);
        double gg_t = -2.0 / gt;
        double var_lv_t = ga_t * ga_t * t_iaa + gg_t * gg_t * t_igg
                        + 2.0 * ga_t * gg_t * t_iag;
        
        double var_lv = (std::max(var_lv_c, 0.0) + std::max(var_lv_t, 0.0)) / 4.0;
        var_lv = std::max(var_lv, 0.01);
        
        d_gg_eff[g] = std::min(200.0, std::max(2.0, 2.0 / var_lv));
    }
    
    // Step 5: mixture weight (sample-size adaptive)
    // n<=5: GG estimates unreliable, fully degrades to limma-trend
    // n>5: gradually introduces GG information (bidirectional fusion)
    double w_base;
    if      (n_min <= 5)  w_base = 0.0;
    else if (n_min <= 10) w_base = 0.20;
    else if (n_min <= 20) w_base = 0.35;
    else if (n_min <= 50) w_base = 0.50;
    else                  w_base = 0.65;
    
    // Step 6: Gene-specific EB shrinkage (Trend + bidirectional GG fusion)
    // Unified formula:
    //   s0_trend_g = var_prior_trend[g] (from P0: mean-dependent trend prior)
    //   w_eff_g = w_base * min(1, d_gg_eff_g / d0)  (GG precision controls weight)
    //   s0_gene = exp((1-w_eff)*log(s0_trend_g) + w_eff*log(v_gg_trend_g))
    //   df_total remains fixed (safety: never modify reference distribution thickness)
    //
    // n<=5: w_base=0 → s0_gene = s0_trend_g (pure trend prior, not fixed s0²)
    // n>5: bidirectional fusion, GG can raise or lower variance, bounded by d_gg_eff
    NumericVector var_shrunk(n_genes);
    NumericVector gg_weight(n_genes);
    NumericVector df_gene(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        df_gene[g] = df_total;
        
        double s0_trend_g = std::max(var_prior_trend[g], eps);
        
        double w_eff = (d0 > eps) ? w_base * std::min(1.0, d_gg_eff[g] / d0) : 0.0;
        gg_weight[g] = w_eff;
        
        double s0_gene;
        if (w_eff < 1e-10) {
            s0_gene = s0_trend_g;
        } else {
            double log_s0t = std::log(s0_trend_g);
            double log_vgg = std::log(std::max((double)var_gg_trend[g], eps));
            s0_gene = std::exp((1.0 - w_eff) * log_s0t + w_eff * log_vgg);
        }
        // clamp: do not allow final prior to deviate too far from trend prior (safety valve)
        s0_gene = std::max(0.1 * s0_trend_g, std::min(10.0 * s0_trend_g, s0_gene));
        
        var_shrunk[g] = (d0 * s0_gene + df_residual * std::max((double)pooled_var_log2[g], eps)) / df_total;
    }
    
    // Step 7: t-test with gene-specific reference distribution
    NumericVector t_stat(n_genes), pvalue_t(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double se = std::sqrt(var_shrunk[g] * (1.0 / n_ctrl + 1.0 / n_treat) + eps);
        t_stat[g] = log2FC[g] / (se + eps);
        pvalue_t[g] = 2.0 * R::pt(-std::abs(t_stat[g]), df_gene[g], 1, 0);
    }
    
    // =========================================================================
    // Step 7b: GG Mean (μ) Wald Test
    // μ = β · Γ(α + 1/γ) / Γ(α)
    // H₀: log(μ_treat) = log(μ_ctrl)
    // Var(log μ̂) via delta method: J' · I⁻¹_logspace · J
    //   J = (∂log μ/∂logα, ∂log μ/∂logβ, ∂log μ/∂logγ)
    //     = (α·[ψ(α+1/γ) - ψ(α)],  1,  -ψ(α+1/γ)/γ)
    // =========================================================================
    NumericVector log2FC_gg(n_genes), pvalue_mu_wald(n_genes);
    std::vector<double> z_mu_raw(n_genes);
    
    // Phase 1: compute raw z-statistics (parallel)
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ac = std::max((double)ae_alpha[g], 0.05);
        double bc = std::max((double)ae_beta[g], eps);
        double gc = std::max((double)ae_gamma[g], 0.05);
        double at = std::max((double)treat_alpha[g], 0.05);
        double bt = std::max((double)treat_beta[g], eps);
        double gt = std::max((double)treat_gamma[g], 0.05);
        
        // GG mean: log(μ) = log(β) + lgamma(α + 1/γ) - lgamma(α)
        double log_mu_c = std::log(bc) + fast_special::lgamma(ac + 1.0/gc) - fast_special::lgamma(ac);
        double log_mu_t = std::log(bt) + fast_special::lgamma(at + 1.0/gt) - fast_special::lgamma(at);
        
        log2FC_gg[g] = (log_mu_t - log_mu_c) / std::log(2.0);
        
        // --- Ctrl: Jacobian + Fisher inverse ---
        double psi_apg_c = fast_special::digamma(ac + 1.0/gc);
        double ja_c = ac * (psi_apg_c - fast_special::digamma(ac));
        double jb_c = 1.0;
        double jg_c = -psi_apg_c / gc;
        
        FisherMatrix3x3 Fc_w = compute_fisher_gg(ac, bc, gc, n_ctrl);
        double fc_aa = Fc_w.I_aa * ac * ac;
        double fc_ab = Fc_w.I_ab * ac * bc;
        double fc_ag = Fc_w.I_ag * ac * gc;
        double fc_bb = Fc_w.I_bb * bc * bc;
        double fc_bg = Fc_w.I_bg * bc * gc;
        double fc_gg = Fc_w.I_gg * gc * gc;
        
        double ic_aa, ic_ab, ic_ag, ic_bb, ic_bg, ic_gg;
        invert_sym_3x3(fc_aa, fc_ab, fc_ag, fc_bb, fc_bg, fc_gg,
                       ic_aa, ic_ab, ic_ag, ic_bb, ic_bg, ic_gg);
        
        double var_c = ja_c*ja_c*ic_aa + jb_c*jb_c*ic_bb + jg_c*jg_c*ic_gg
                     + 2.0*ja_c*jb_c*ic_ab + 2.0*ja_c*jg_c*ic_ag + 2.0*jb_c*jg_c*ic_bg;
        
        // --- Treat: Jacobian + Fisher inverse ---
        double psi_apg_t = fast_special::digamma(at + 1.0/gt);
        double ja_t = at * (psi_apg_t - fast_special::digamma(at));
        double jb_t = 1.0;
        double jg_t = -psi_apg_t / gt;
        
        FisherMatrix3x3 Ft_w = compute_fisher_gg(at, bt, gt, n_treat);
        double ft_aa = Ft_w.I_aa * at * at;
        double ft_ab = Ft_w.I_ab * at * bt;
        double ft_ag = Ft_w.I_ag * at * gt;
        double ft_bb = Ft_w.I_bb * bt * bt;
        double ft_bg = Ft_w.I_bg * bt * gt;
        double ft_gg = Ft_w.I_gg * gt * gt;
        
        double it_aa, it_ab, it_ag, it_bb, it_bg, it_gg;
        invert_sym_3x3(ft_aa, ft_ab, ft_ag, ft_bb, ft_bg, ft_gg,
                       it_aa, it_ab, it_ag, it_bb, it_bg, it_gg);
        
        double var_t = ja_t*ja_t*it_aa + jb_t*jb_t*it_bb + jg_t*jg_t*it_gg
                     + 2.0*ja_t*jb_t*it_ab + 2.0*ja_t*jg_t*it_ag + 2.0*jb_t*jg_t*it_bg;
        
        double delta_logmu = log_mu_t - log_mu_c;
        double var_delta = std::max(var_c, 0.0) + std::max(var_t, 0.0);
        z_mu_raw[g] = delta_logmu / std::sqrt(var_delta + eps);
    }
    
    // Phase 2: Genomic Control calibration (Devlin & Roeder 1999)
    // Under H₀, z² ~ χ²(1), median = 0.4549. Under model misspecification, z² ~ λ·χ²(1)
    // λ̂ = median(z²) / 0.4549, only deflate (λ≥1)
    {
        const double chi1_median = 0.4549364;  // qchisq(0.5, 1)
        std::vector<double> z_sq(n_genes);
        for (int g = 0; g < n_genes; g++) z_sq[g] = z_mu_raw[g] * z_mu_raw[g];
        int idx_med = n_genes / 2;
        std::nth_element(z_sq.begin(), z_sq.begin() + idx_med, z_sq.end());
        double lambda_gc = std::max(1.0, z_sq[idx_med] / chi1_median);
        double sqrt_lambda = std::sqrt(lambda_gc);
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < n_genes; g++) {
            double z_corrected = z_mu_raw[g] / sqrt_lambda;
            pvalue_mu_wald[g] = 2.0 * R::pnorm(-std::abs(z_corrected), 0.0, 1.0, 1, 0);
            pvalue_mu_wald[g] = std::max(eps, std::min(1.0, pvalue_mu_wald[g]));
        }
    }
    
    // =========================================================================
    // Manifold distance test (captures overall displacement of α,β,γ)
    // Core innovation of SGCB: tests not only the mean (β) but also transcriptional dynamics (α,γ)
    // =========================================================================
    NumericVector geodesic_dist(n_genes), pvalue_manifold(n_genes);
    std::vector<double> dist_vec(n_genes);
    std::vector<double> chi_ab_vec(n_genes);
    std::vector<double> chi_g_vec(n_genes);
    std::vector<double> chi_sq_vec(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        double ac = std::max(ae_alpha[g], eps);
        double bc = std::max(ae_beta[g], eps);
        double gc = std::max(ae_gamma[g], eps);
        double at = std::max(treat_alpha[g], eps);
        double bt = std::max(treat_beta[g], eps);
        double gt = std::max(treat_gamma[g], eps);

        dist_vec[g] = gg_geodesic_distance(
            ac, bc, gc,
            at, bt, gt, eps);

        double d_a = std::log(at) - std::log(ac);
        double d_b = std::log(bt) - std::log(bc);
        double d_g = std::log(gt) - std::log(gc);

        FisherMatrix3x3 Fc = compute_fisher_gg(ac, bc, gc, n_ctrl);
        FisherMatrix3x3 Ft = compute_fisher_gg(at, bt, gt, n_treat);

        const double reg = 1e-12;

        double Fca = Fc.I_aa * ac * ac;
        double Fcb = Fc.I_bb * bc * bc;
        double Fcab = Fc.I_ab * ac * bc;
        double trFc = Fca + Fcb;
        double ridgeFc = 1e-12 * std::max(trFc, 0.0) + reg;
        double Fca_r = Fca + ridgeFc;
        double Fcb_r = Fcb + ridgeFc;
        double detFc = Fca_r * Fcb_r - Fcab * Fcab;
        detFc = std::max(detFc, reg);
        double invc_aa = Fcb_r / detFc;
        double invc_bb = Fca_r / detFc;
        double invc_ab = -Fcab / detFc;

        double Fta = Ft.I_aa * at * at;
        double Ftb = Ft.I_bb * bt * bt;
        double Ftab = Ft.I_ab * at * bt;
        double trFt = Fta + Ftb;
        double ridgeFt = 1e-12 * std::max(trFt, 0.0) + reg;
        double Fta_r = Fta + ridgeFt;
        double Ftb_r = Ftb + ridgeFt;
        double detFt = Fta_r * Ftb_r - Ftab * Ftab;
        detFt = std::max(detFt, reg);
        double invt_aa = Ftb_r / detFt;
        double invt_bb = Fta_r / detFt;
        double invt_ab = -Ftab / detFt;

        double Vaa = invc_aa + invt_aa;
        double Vab = invc_ab + invt_ab;
        double Vbb = invc_bb + invt_bb;
        double trV = Vaa + Vbb;
        double ridgeV = 1e-12 * std::max(trV, 0.0) + reg;
        double Vaa_r = Vaa + ridgeV;
        double Vbb_r = Vbb + ridgeV;
        double detV = Vaa_r * Vbb_r - Vab * Vab;
        detV = std::max(detV, reg);
        double Paa = Vbb_r / detV;
        double Pbb = Vaa_r / detV;
        double Pab = -Vab / detV;
        double chi_ab = d_a * (Paa * d_a + Pab * d_b) + d_b * (Pab * d_a + Pbb * d_b);

        double varc_g = 1.0 / std::max(Fc.I_gg * gc * gc, reg);
        double vart_g = 1.0 / std::max(Ft.I_gg * gt * gt, reg);
        double var_g = varc_g + vart_g;
        double chi_g = (d_g * d_g) / std::max(var_g, eps);

        chi_ab_vec[g] = std::max(chi_ab, 0.0);
        chi_g_vec[g] = std::max(chi_g, 0.0);
        chi_sq_vec[g] = std::max(chi_ab + chi_g, 0.0);
    }

    for (int g = 0; g < n_genes; g++) geodesic_dist[g] = dist_vec[g];

    std::vector<double> chi_tmp(n_genes);
    for (int g = 0; g < n_genes; g++) chi_tmp[g] = chi_sq_vec[g];

    int idx25 = n_genes / 4;
    int idx75 = (3 * n_genes) / 4;

    std::nth_element(chi_tmp.begin(), chi_tmp.begin() + idx25, chi_tmp.end());
    double chi_q25 = chi_tmp[idx25];

    std::nth_element(chi_tmp.begin(), chi_tmp.begin() + idx75, chi_tmp.end());
    double chi_q75 = chi_tmp[idx75];

    // IQR-based df and scale estimation (Genomic Control style)
    double df_hat, scale_hat;
    estimate_scaled_chisq_df_scale_from_iqr(chi_q25, chi_q75, df_hat, scale_hat, eps);
    // Fallback: if estimation gives degenerate values, use df=3 with median calibration
    if (df_hat < 0.5 || df_hat > 30.0 || scale_hat < eps) {
        df_hat = 3.0;
        std::vector<double> chi_tmp2(chi_tmp);
        int idx50 = n_genes / 2;
        std::nth_element(chi_tmp2.begin(), chi_tmp2.begin() + idx50, chi_tmp2.end());
        double chi_q50 = chi_tmp2[idx50];
        double chi_ref_q50 = R::qchisq(0.5, df_hat, 1, 0);
        scale_hat = std::max(chi_q50 / std::max(chi_ref_q50, eps), eps);
    }

    // Fix E: disable manifold test when n_min < 6 (asymptotic theory unreliable)
    bool manifold_enabled = (n_min >= 6);
    for (int g = 0; g < n_genes; g++) {
        if (manifold_enabled) {
            double chi_adj = chi_sq_vec[g] / scale_hat;
            pvalue_manifold[g] = R::pchisq(chi_adj, df_hat, 0, 0);
            pvalue_manifold[g] = std::max(eps, std::min(1.0, pvalue_manifold[g]));
        } else {
            pvalue_manifold[g] = 1.0;
        }
    }
    
    // =========================================================================
    // Fix C: primary p-value uses moderated t-test only
    // Manifold test is retained as a diagnostic column, not used for FDR
    // =========================================================================
    NumericVector pvalue(n_genes);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < n_genes; g++) {
        pvalue[g] = std::max(eps, std::min(1.0, pvalue_t[g]));
    }
    
    return List::create(
        Named("log2FC") = log2FC,
        Named("log2FC_gg") = log2FC_gg,
        Named("mean_ctrl") = mean_ctrl,
        Named("mean_treat") = mean_treat,
        Named("t_stat") = t_stat,
        Named("pvalue") = pvalue,
        Named("pvalue_t") = pvalue_t,
        Named("pvalue_mu_wald") = pvalue_mu_wald,
        Named("pvalue_manifold") = pvalue_manifold,
        Named("geodesic_dist") = geodesic_dist,
        Named("manifold_df") = df_hat,
        Named("manifold_scale") = scale_hat,
        Named("var_shrunk") = var_shrunk,
        Named("var_prior_trend") = var_prior_trend,
        Named("var_gg_log2") = var_gg_log2,
        Named("var_gg_trend") = var_gg_trend,
        Named("gg_weight") = gg_weight,
        Named("d_gg_eff") = d_gg_eff,
        Named("df_gene") = df_gene,
        Named("baseMean") = baseMean,
        Named("df_total") = df_total,
        Named("ctrl_alpha") = ae_alpha,
        Named("ctrl_beta") = ae_beta,
        Named("ctrl_gamma") = ae_gamma,
        Named("treat_alpha") = treat_alpha,
        Named("treat_beta") = treat_beta,
        Named("treat_gamma") = treat_gamma,
        Named("use_natural_grad") = use_natural_grad,
        Named("use_hierarchical") = use_hierarchical
    );
}
