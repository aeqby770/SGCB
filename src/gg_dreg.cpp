// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
// gg_dreg.cpp - GG Distributional Regression V1
//
// Y_{gi} ~ GG(alpha_g, beta_{gi}, gamma_g), log(beta_{gi}) = x_i' b_g
// b-block: Fisher scoring (weight gamma^2*alpha = const -> OLS)
// (alpha,gamma)-block: Natural gradient on Fisher-Rao 2x2 sub-manifold
// EB shrinkage: fitFDist on Pearson dispersion, moderated t

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

// ========================= Dense solver =========================
static bool dreg_solve(std::vector<double>& A, std::vector<double>& b, int P) {
    for (int j = 0; j < P; j++) {
        int piv = j; double mx = std::abs(A[j*P+j]);
        for (int i = j+1; i < P; i++) {
            double v = std::abs(A[i*P+j]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (piv != j) {
            for (int c = 0; c < P; c++) std::swap(A[j*P+c], A[piv*P+c]);
            std::swap(b[j], b[piv]);
        }
        if (std::abs(A[j*P+j]) < 1e-15) return false;
        for (int i = j+1; i < P; i++) {
            double f = A[i*P+j] / A[j*P+j];
            for (int c = j; c < P; c++) A[i*P+c] -= f * A[j*P+c];
            b[i] -= f * b[j];
        }
    }
    for (int j = P-1; j >= 0; j--) {
        for (int c = j+1; c < P; c++) b[j] -= A[j*P+c] * b[c];
        b[j] /= A[j*P+j];
    }
    return true;
}

static void dreg_invert(const std::vector<double>& M, std::vector<double>& inv, int P) {
    for (int col = 0; col < P; col++) {
        std::vector<double> A(M);
        std::vector<double> e(P, 0.0);
        e[col] = 1.0;
        dreg_solve(A, e, P);
        for (int r = 0; r < P; r++) inv[r*P + col] = e[r];
    }
}

static double dreg_logdet_lu(std::vector<double> A, int P, double eps = 1e-12) {
    double logdet = 0.0;
    for (int j = 0; j < P; j++) {
        int piv = j;
        double mx = std::abs(A[j*P + j]);
        for (int i = j + 1; i < P; i++) {
            double v = std::abs(A[i*P + j]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (mx < eps) return std::log(eps);
        if (piv != j) {
            for (int c = 0; c < P; c++) std::swap(A[j*P + c], A[piv*P + c]);
        }
        double diag = A[j*P + j];
        if (std::abs(diag) < eps) return std::log(eps);
        logdet += std::log(std::abs(diag));
        for (int i = j + 1; i < P; i++) {
            double f = A[i*P + j] / diag;
            for (int c = j + 1; c < P; c++) A[i*P + c] -= f * A[j*P + c];
        }
    }
    return logdet;
}

// ========================= fitFDist =========================
struct DregFDist { double scale; double df2; };

static DregFDist dreg_fitfdist(const std::vector<double>& vars, int n, int df1,
                                double eps = 1e-8) {
    DregFDist res;
    std::vector<double> lv(n);
    for (int i = 0; i < n; i++) lv[i] = std::log(std::max(vars[i], eps));

    std::vector<double> sv(lv); std::sort(sv.begin(), sv.end());
    double lo = sv[std::max(0, (int)(0.01*n))];
    double hi = sv[std::min(n-1, (int)(0.99*n))];
    for (int i = 0; i < n; i++) lv[i] = std::max(lo, std::min(hi, lv[i]));

    double hd1 = std::max(df1 / 2.0, 0.5);
    double psi1 = fast_special::digamma(hd1);
    double log1 = std::log(hd1);

    double sm = 0;
    for (int i = 0; i < n; i++) sm += lv[i] - psi1 + log1;
    double emean = sm / n;

    double ss = 0;
    for (int i = 0; i < n; i++) {
        double d = (lv[i] - psi1 + log1) - emean; ss += d*d;
    }
    double evar = ss / std::max(n-1,1) - fast_special::trigamma(hd1);

    if (evar <= 0) {
        res.df2 = 1e6; res.scale = std::exp(emean);
    } else {
        double hd2 = fast_special::trigamma_inverse(evar);
        res.df2 = std::max(0.5, std::min(1e6, 2.0*hd2));
        hd2 = res.df2 / 2.0;
        res.scale = std::exp(emean + fast_special::digamma(hd2) - std::log(hd2));
    }
    res.scale = std::max(res.scale, eps);
    return res;
}

// ========================= Main engine =========================
// [[Rcpp::export]]
List sgcb_dreg_v1_cpp(
    NumericMatrix y,        // G x N positive values (TMM-normalized)
    NumericMatrix design,   // N x P design matrix
    NumericVector contrast, // P contrast vector
    int max_outer = 20,
    int max_inner_ag = 40,
    double ag_lr = 0.15,
    double tol = 1e-5,
    double eps = 1e-8)
{
    const int G = y.nrow(), N = y.ncol(), P = design.ncol();
    const int df_resid = std::max(N - P, 1);
    const bool hard_fix_gamma = (df_resid <= 8);

    // ---- Precompute (X'X)^{-1} and hat matrix ----
    std::vector<double> XtX(P*P, 0.0);
    for (int j = 0; j < N; j++)
        for (int a = 0; a < P; a++)
            for (int b = a; b < P; b++)
                XtX[a*P+b] += design(j,a) * design(j,b);
    for (int a = 0; a < P; a++)
        for (int b = 0; b < a; b++) XtX[a*P+b] = XtX[b*P+a];

    std::vector<double> XtX_inv(P*P, 0.0);
    dreg_invert(XtX, XtX_inv, P);

    // c'(X'X)^{-1}c
    std::vector<double> Hc(P, 0.0);
    double vc = 0.0;
    for (int a = 0; a < P; a++)
        for (int b = 0; b < P; b++) Hc[a] += XtX_inv[a*P+b] * contrast[b];
    for (int a = 0; a < P; a++) vc += contrast[a] * Hc[a];

    // (X'X)^{-1} X' : P x N
    std::vector<double> HXt(P*N, 0.0);
    for (int a = 0; a < P; a++)
        for (int j = 0; j < N; j++) {
            double s = 0;
            for (int b = 0; b < P; b++) s += XtX_inv[a*P+b] * design(j,b);
            HXt[a*N+j] = s;
        }

    // ---- Initialize b, alpha, gamma ----
    std::vector<std::vector<double>> coef(G, std::vector<double>(P, 0.0));
    std::vector<double> log_a(G), log_g(G);

    // Harrell-Davis quantile weights for robust initialization
    std::vector<double> hd_w25(N), hd_w50(N), hd_w75(N);
    {
        auto fill_hd = [&](double p, std::vector<double>& out) {
            double aa = (N + 1.0) * p;
            double bb = (N + 1.0) * (1.0 - p);
            double prev = 0.0;
            for (int i = 1; i <= N; i++) {
                double u = static_cast<double>(i) / static_cast<double>(N);
                double cur = R::pbeta(u, aa, bb, 1, 0);
                out[i - 1] = std::max(0.0, cur - prev);
                prev = cur;
            }
            double s = std::accumulate(out.begin(), out.end(), 0.0);
            if (s <= 1e-12) {
                for (int i = 0; i < N; i++) out[i] = 1.0 / static_cast<double>(N);
            } else {
                for (int i = 0; i < N; i++) out[i] /= s;
            }
        };
        fill_hd(0.25, hd_w25);
        fill_hd(0.50, hd_w50);
        fill_hd(0.75, hd_w75);
    }

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < G; g++) {
        // b from OLS on log(y)
        for (int a = 0; a < P; a++) {
            double s = 0;
            for (int j = 0; j < N; j++)
                s += HXt[a*N+j] * std::log(std::max((double)y(g,j), eps));
            coef[g][a] = s;
        }
        // Robust initialization via HD quantile standardization
        std::vector<double> yv(N);
        for (int j = 0; j < N; j++) yv[j] = std::max((double)y(g,j), eps);
        std::sort(yv.begin(), yv.end());

        double q25 = 0.0, q50 = 0.0, q75 = 0.0;
        for (int j = 0; j < N; j++) {
            q25 += hd_w25[j] * yv[j];
            q50 += hd_w50[j] * yv[j];
            q75 += hd_w75[j] * yv[j];
        }
        q50 = std::max(q50, eps);
        double iqr = std::max(q75 - q25, eps);
        double lo = std::max(q25 - 4.0 * iqr, eps);
        double hi = q75 + 4.0 * iqr;

        double sy = 0, sy2 = 0;
        for (int j = 0; j < N; j++) {
            double yw = std::max(lo, std::min(hi, std::max((double)y(g,j), eps)));
            double ys = yw / q50;
            sy += ys;
            sy2 += ys * ys;
        }
        double my = sy / N;
        double cv2_raw = std::max((sy2/N - my*my) / (my*my + eps), 0.05);
        log_a[g] = std::log(std::max(0.5, std::min(15.0, 1.0 / cv2_raw)));
        // Start at gamma=1 (stable Gamma anchor), then free with warm-up schedule
        log_g[g] = 0.0;
    }

    // ---- Hierarchical prior on (log_alpha, log_gamma) ----
    auto update_prior = [&](double& mu_la, double& mu_lg,
                            double& sig_la, double& sig_lg) {
        std::vector<double> la(log_a), lg(log_g);
        std::nth_element(la.begin(), la.begin()+G/2, la.end()); mu_la = la[G/2];
        std::nth_element(lg.begin(), lg.begin()+G/2, lg.end()); mu_lg = lg[G/2];
        std::vector<double> da(G), dg(G);
        for (int g = 0; g < G; g++) {
            da[g] = std::abs(log_a[g] - mu_la);
            dg[g] = std::abs(log_g[g] - mu_lg);
        }
        std::sort(da.begin(), da.end());
        std::sort(dg.begin(), dg.end());
        sig_la = std::max(0.1, da[G/2]*1.4826);
        sig_lg = std::max(0.05, dg[G/2]*1.4826);
    };

    double mu_la, mu_lg, sig_la, sig_lg;
    update_prior(mu_la, mu_lg, sig_la, sig_lg);
    double pw = std::max(0.01, 2.0/std::sqrt((double)N+1.0));

    // ---- Blockwise iteration ----
    NumericVector ll_hist(max_outer);
    int n_done = 0;

    for (int outer = 0; outer < max_outer; outer++) {
        n_done = outer + 1;

        // Gamma warm-up: fix gamma=1 for first 5 outer iters when freed
        bool cur_fix_gamma = hard_fix_gamma || (outer < 5);
        // Gradual gamma bounds: tight initially, widen later
        double lg_lo = -0.5, lg_hi = 1.1;
        if (!hard_fix_gamma && outer < 10) {
            lg_lo = -0.25; lg_hi = 0.35;  // gamma in [0.78, 1.42]
        }

        // --- b-block: OLS on working response ---
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int g = 0; g < G; g++) {
            double ag = std::exp(log_a[g]), gg = std::exp(log_g[g]);
            std::vector<double> z(N);
            for (int j = 0; j < N; j++) {
                double eta = 0;
                for (int a = 0; a < P; a++) eta += design(j,a)*coef[g][a];
                double bj = std::exp(eta);
                double ratio = std::max((double)y(g,j), eps) / bj;
                double log_u = gg * std::log(std::max(ratio, eps));
                log_u = std::max(-20.0, std::min(20.0, log_u));
                double u = std::exp(log_u);
                double adj = (u - ag) / (gg * ag + eps);
                adj = std::max(-5.0, std::min(5.0, adj));
                z[j] = eta + adj;
            }
            for (int a = 0; a < P; a++) {
                double s = 0;
                for (int j = 0; j < N; j++) s += HXt[a*N+j]*z[j];
                coef[g][a] = s;
            }
        }

        // --- (alpha,gamma)-block: natural gradient ---
        int ni = (outer == 0) ? max_inner_ag : std::max(5, max_inner_ag/4);
        double lr = (outer == 0) ? ag_lr : ag_lr * 0.7;
        double total_ll = 0.0;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:total_ll)
        #endif
        for (int g = 0; g < G; g++) {
            std::vector<double> rv(N), lrv(N);
            for (int j = 0; j < N; j++) {
                double eta = 0;
                for (int a = 0; a < P; a++) eta += design(j,a)*coef[g][a];
                rv[j] = std::max((double)y(g,j), eps) / std::exp(eta);
                lrv[j] = std::log(rv[j]);
            }

            double la = log_a[g], lg = log_g[g];
            for (int it = 0; it < ni; it++) {
                double a = std::exp(la), gm = std::exp(lg);
                double slr = 0, srg = 0, srglr = 0;
                for (int j = 0; j < N; j++) {
                    double rg = std::exp(std::min(gm*lrv[j], 500.0));
                    slr += lrv[j]; srg += rg; srglr += rg*lrv[j];
                }
                double sc_a = -N*fast_special::digamma(a) + gm*slr;
                double sc_g = N/gm + a*slr - srglr;
                double gla = a*sc_a - pw*(la-mu_la)/(sig_la*sig_la+eps);
                double glg = gm*sc_g - pw*(lg-mu_lg)/(sig_lg*sig_lg+eps);
                if (cur_fix_gamma) glg = 0.0;

                double p1a = fast_special::trigamma(a);
                double psa = fast_special::digamma(a);
                double Faa = a*a*N*p1a + 1e-6;
                double Fag = -a*N*psa;
                double Fgg = N*(1.0 + a*psa*psa + 2.0*psa + a*p1a) + 1e-6;
                double det = std::max(Faa*Fgg - Fag*Fag, 1e-10);

                double nla, nlg;
                if (cur_fix_gamma) { nla = gla/Faa; nlg = 0; }
                else {
                    nla = (Fgg*gla - Fag*glg) / det;
                    nlg = (-Fag*gla + Faa*glg) / det;
                    nlg *= 0.3;  // dampen gamma to prevent overshooting
                }
                la += lr*nla;  lg += lr*nlg;
                la = std::max(-1.0, std::min(3.5, la));
                lg = std::max(lg_lo, std::min(lg_hi, lg));
                if (cur_fix_gamma) lg = 0.0;
            }
            log_a[g] = la; log_g[g] = lg;

            double a = std::exp(la), gm = std::exp(lg);
            double ll = 0;
            for (int j = 0; j < N; j++) {
                double eta = 0;
                for (int aa = 0; aa < P; aa++) eta += design(j,aa)*coef[g][aa];
                double xi = std::max((double)y(g,j), eps);
                ll += std::log(gm) - fast_special::lgamma(a)
                    - a*gm*eta + (a*gm-1.0)*std::log(xi)
                    - std::exp(std::min(gm*std::log(xi/std::exp(eta)), 500.0));
            }
            total_ll += ll;
        }
        ll_hist[outer] = total_ll;

        // Early termination when log-likelihood converges (uses `tol` parameter)
        if (outer > 0 &&
            std::abs(ll_hist[outer] - ll_hist[outer - 1]) <
                tol * (std::abs(ll_hist[outer]) + eps)) {
            for (int r = outer + 1; r < max_outer; r++) ll_hist[r] = total_ll;
            break;
        }

        if ((outer+1) % 5 == 0) update_prior(mu_la, mu_lg, sig_la, sig_lg);
    }

    // ---- Pearson dispersion + EB shrinkage ----
    std::vector<double> phi(G), bm(G);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < G; g++) {
        double a = std::exp(log_a[g]), gm = std::exp(log_g[g]);
        double R1 = std::exp(fast_special::lgamma(a+1.0/gm) - fast_special::lgamma(a));
        double R2 = std::exp(fast_special::lgamma(a+2.0/gm) - fast_special::lgamma(a));
        double cv2 = std::max(R2/(R1*R1) - 1.0, eps);

        double ps = 0, ms = 0;
        for (int j = 0; j < N; j++) {
            double eta = 0;
            for (int aa = 0; aa < P; aa++) eta += design(j,aa)*coef[g][aa];
            double mu = std::exp(eta) * R1;
            double vr = mu*mu*cv2;
            double res = y(g,j) - mu;
            ps += res*res / (vr + eps);
            ms += mu;
        }
        phi[g] = std::max(ps / df_resid, eps);
        bm[g] = ms / N;
    }

    DregFDist fd = dreg_fitfdist(phi, G, df_resid, eps);
    double d0 = fd.df2, s0sq = fd.scale, df_tot = d0 + df_resid;

    // ---- Inference ----
    NumericVector lfc(G), se(G), tst(G), pv(G);
    NumericVector ao(G), go(G), dsp(G), dsps(G), bmo(G);
    NumericMatrix cfo(G, P);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < G; g++) {
        double a = std::exp(log_a[g]), gm = std::exp(log_g[g]);
        double phi_s = (d0*s0sq + df_resid*phi[g]) / df_tot;

        double cb = 0;
        for (int aa = 0; aa < P; aa++) {
            cb += contrast[aa]*coef[g][aa];
            cfo(g, aa) = coef[g][aa];
        }
        lfc[g] = cb / std::log(2.0);
        double var_cb = phi_s / (gm*gm*a + eps) * vc;
        se[g] = std::sqrt(std::max(var_cb, eps)) / std::log(2.0);
        tst[g] = lfc[g] / (se[g] + eps);
        pv[g] = 2.0 * R::pt(-std::abs(tst[g]), df_tot, 1, 0);
        pv[g] = std::max(eps, std::min(1.0, pv[g]));

        ao[g] = a; go[g] = gm;
        dsp[g] = phi[g]; dsps[g] = phi_s; bmo[g] = bm[g];
    }

    return List::create(
        Named("log2FoldChange") = lfc, Named("lfcSE") = se,
        Named("stat") = tst, Named("pvalue") = pv,
        Named("baseMean") = bmo,
        Named("alpha") = ao, Named("gamma") = go,
        Named("dispersion") = dsp, Named("dispersion_shrunk") = dsps,
        Named("coefficients") = cfo,
        Named("df_prior") = d0, Named("var_prior") = s0sq,
        Named("df_total") = df_tot, Named("df_residual") = df_resid,
        Named("vc") = vc, Named("n_outer") = n_done, Named("loglik") = ll_hist);
}

// =============================================================================
// V2 DV test: Dispersion GLM + Likelihood Ratio test
//
// Given V1 results (b, alpha_shared, gamma), fit a dispersion GLM:
//   log(alpha_{gi}) = w_i' d_g
// then compute LR = 2*(ll_full - ll_null) ~ chi^2(Q-1)
//
// Algorithm:
//   d-block: Fisher scoring (WLS) on dispersion parameters
//   b-block: WLS with sample-specific alpha weights (optional re-fit)
//   gamma: fixed from V1
// =============================================================================

static double dreg_gene_loglik(const double* yg, int N, const double* eta,
                                double alpha, double gamma, double eps) {
    double ll = 0;
    for (int j = 0; j < N; j++) {
        double xi = std::max(yg[j], eps);
        double log_r = std::log(xi) - eta[j];
        double rg = std::exp(std::min(gamma * log_r, 500.0));
        ll += std::log(gamma) - fast_special::lgamma(alpha)
            + alpha * gamma * log_r - std::log(xi) - rg;
    }
    return ll;
}

static double dreg_gene_loglik_varalpha(const double* yg, int N,
    const double* eta, const double* log_alpha_j, double gamma, double eps) {
    double ll = 0;
    for (int j = 0; j < N; j++) {
        double xi = std::max(yg[j], eps);
        double aj = std::exp(log_alpha_j[j]);
        double log_r = std::log(xi) - eta[j];
        double rg = std::exp(std::min(gamma * log_r, 500.0));
        ll += std::log(gamma) - fast_special::lgamma(aj)
            + aj * gamma * log_r - std::log(xi) - rg;
    }
    return ll;
}

// [[Rcpp::export]]
List sgcb_dreg_dv_cpp(
    NumericMatrix y,            // G x N
    NumericMatrix design_mean,  // N x P
    NumericMatrix design_disp,  // N x Q (dispersion design, usually intercept + group)
    NumericVector contrast_disp,// Q (dispersion contrast, e.g. c(0,1))
    NumericMatrix coef_v1,      // G x P (V1 mean coefficients)
    NumericVector alpha_v1,     // G (V1 alpha, shared across samples)
    NumericVector gamma_v1,     // G (V1 gamma)
    int max_outer_dv = 15,
    int max_inner_d = 8,
    double d_lr = 0.5,
    double eps = 1e-8)
{
    const int G = y.nrow(), N = y.ncol();
    const int P = design_mean.ncol(), Q = design_disp.ncol();

    // ---- Precompute mean hat matrix (X'X)^{-1}X' for b re-fit ----
    std::vector<double> XtX(P*P, 0.0);
    for (int j = 0; j < N; j++)
        for (int a = 0; a < P; a++)
            for (int b = a; b < P; b++)
                XtX[a*P+b] += design_mean(j,a) * design_mean(j,b);
    for (int a = 0; a < P; a++)
        for (int b = 0; b < a; b++) XtX[a*P+b] = XtX[b*P+a];

    // ---- Per-gene: null loglik, then d-block optimization ----
    NumericVector ll_null(G), ll_full(G), lr_stat(G), pv_dv(G);
    NumericVector logdet_full_d(G), logdet_null_d(G);
    NumericVector dv_logratio(G), dv_se_log2(G);
    NumericMatrix d_coef(G, Q);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int g = 0; g < G; g++) {
        double gm = gamma_v1[g];
        double la_null = std::log(std::max((double)alpha_v1[g], eps));

        // --- eta from V1 b ---
        std::vector<double> eta(N), log_r(N);
        for (int j = 0; j < N; j++) {
            eta[j] = 0;
            for (int a = 0; a < P; a++) eta[j] += design_mean(j,a) * coef_v1(g,a);
            log_r[j] = std::log(std::max((double)y(g,j), eps)) - eta[j];
        }

        // --- Refine alpha toward pure-likelihood null MLE ---
        // V1's alpha is MAP (includes prior); for a proper LR test the null
        // loglik must be evaluated at the constrained MLE (no prior).
        // A few Newton-Raphson steps on the 1-D profile loglik suffice.
        // Profile score:  d ell / d(log alpha) = alpha * [-N psi(alpha) + gamma * sum(log_r)]
        // Profile Hessian: d^2 / d(log alpha)^2 = -alpha^2 * N * psi'(alpha)
        {
            double la_nr = la_null;
            double slr = 0;
            for (int j = 0; j < N; j++) slr += log_r[j];
            for (int nr = 0; nr < 5; nr++) {
                double a_cur = std::exp(la_nr);
                double score_la = a_cur * (-N * fast_special::digamma(a_cur) + gm * slr);
                double H_la = a_cur * a_cur * N * fast_special::trigamma(a_cur);
                if (H_la < eps) break;
                double step = score_la / (H_la + eps);
                step = std::max(-1.0, std::min(1.0, step));   // guard large jumps
                la_nr += step;
                la_nr = std::max(-2.0, std::min(5.0, la_nr));
                if (std::abs(step) < 1e-8) break;
            }
            la_null = la_nr;
        }
        double alpha_null_mle = std::exp(la_null);
        double info_null_la = static_cast<double>(N) * alpha_null_mle * alpha_null_mle *
                              std::max(fast_special::trigamma(alpha_null_mle), eps);
        double logdet_null_gene = std::log(std::max(info_null_la, eps));

        // --- Null loglik (shared alpha, pure-likelihood MLE) ---
        std::vector<double> yg(N);
        for (int j = 0; j < N; j++) yg[j] = y(g,j);
        ll_null[g] = dreg_gene_loglik(yg.data(), N, eta.data(),
                                       alpha_null_mle, gm, eps);
        logdet_null_d[g] = logdet_null_gene;

        // --- Initialize d: intercept = log(alpha_v1), rest = 0 ---
        std::vector<double> dg(Q, 0.0);
        dg[0] = la_null;

        // --- d-block Fisher scoring (WLS) ---
        for (int outer = 0; outer < max_outer_dv; outer++) {
            for (int it = 0; it < max_inner_d; it++) {
                // Compute alpha_j, score, Fisher weight
                std::vector<double> la_j(N), score(N), omega(N);
                for (int j = 0; j < N; j++) {
                    la_j[j] = 0;
                    for (int q = 0; q < Q; q++) la_j[j] += design_disp(j,q)*dg[q];
                    double aj = std::exp(la_j[j]);
                    double psi_aj = fast_special::digamma(aj);
                    double psi1_aj = fast_special::trigamma(aj);
                    score[j] = -psi_aj + gm * log_r[j];
                    omega[j] = aj * aj * std::max(psi1_aj, eps);
                }

                // WLS: d_new = (W'OmW)^{-1} W'Om z_d
                // z_dj = w_j'd + score_j / (alpha_j * psi1_aj)
                std::vector<double> WtOW(Q*Q, 0.0), WtOz(Q, 0.0);
                for (int j = 0; j < N; j++) {
                    double aj = std::exp(la_j[j]);
                    double psi1 = fast_special::trigamma(aj);
                    double z_d = la_j[j] + score[j] / (aj * std::max(psi1, eps));
                    for (int a = 0; a < Q; a++) {
                        WtOz[a] += omega[j] * design_disp(j,a) * z_d;
                        for (int b = a; b < Q; b++)
                            WtOW[a*Q+b] += omega[j] * design_disp(j,a) * design_disp(j,b);
                    }
                }
                for (int a = 0; a < Q; a++)
                    for (int b = 0; b < a; b++) WtOW[a*Q+b] = WtOW[b*Q+a];

                // Solve
                std::vector<double> d_new(WtOz);
                std::vector<double> A(WtOW);
                dreg_solve(A, d_new, Q);

                // Damped update
                double max_step = 0;
                for (int q = 0; q < Q; q++) {
                    double step = d_new[q] - dg[q];
                    step = std::max(-2.0, std::min(2.0, step));
                    dg[q] += d_lr * step;
                    if (std::abs(step) > max_step) max_step = std::abs(step);
                }
                if (max_step < 1e-6) break;
            }

            // --- b re-fit (WLS with sample-specific alpha) ---
            std::vector<double> la_j(N);
            for (int j = 0; j < N; j++) {
                la_j[j] = 0;
                for (int q = 0; q < Q; q++) la_j[j] += design_disp(j,q)*dg[q];
            }
            std::vector<double> z_b(N), omega_b(N);
            for (int j = 0; j < N; j++) {
                double aj = std::exp(la_j[j]);
                double bj = std::exp(eta[j]);
                double ratio = std::max((double)y(g,j), eps) / bj;
                double log_u = gm * std::log(std::max(ratio, eps));
                log_u = std::max(-20.0, std::min(20.0, log_u));
                double u = std::exp(log_u);
                double adj = (u - aj) / (gm * aj + eps);
                adj = std::max(-5.0, std::min(5.0, adj));
                z_b[j] = eta[j] + adj;
                omega_b[j] = gm * gm * aj;
            }
            // WLS for b
            std::vector<double> XtWX(P*P, 0.0), XtWz(P, 0.0);
            for (int j = 0; j < N; j++) {
                for (int a = 0; a < P; a++) {
                    XtWz[a] += omega_b[j] * design_mean(j,a) * z_b[j];
                    for (int b = a; b < P; b++)
                        XtWX[a*P+b] += omega_b[j] * design_mean(j,a) * design_mean(j,b);
                }
            }
            for (int a = 0; a < P; a++)
                for (int b = 0; b < a; b++) XtWX[a*P+b] = XtWX[b*P+a];
            std::vector<double> b_new(XtWz);
            std::vector<double> Ab(XtWX);
            dreg_solve(Ab, b_new, P);

            // Update eta and log_r
            for (int j = 0; j < N; j++) {
                eta[j] = 0;
                for (int a = 0; a < P; a++) eta[j] += design_mean(j,a) * b_new[a];
                log_r[j] = std::log(std::max((double)y(g,j), eps)) - eta[j];
            }
        }

        // --- Full loglik (sample-specific alpha) ---
        std::vector<double> la_final(N);
        for (int j = 0; j < N; j++) {
            la_final[j] = 0;
            for (int q = 0; q < Q; q++) la_final[j] += design_disp(j,q)*dg[q];
        }
        ll_full[g] = dreg_gene_loglik_varalpha(yg.data(), N, eta.data(),
                                                 la_final.data(), gm, eps);

        // --- LR stat (omnibus: tests all Q-1 non-intercept dispersion parameters) ---
        // NOTE: pvalue_dv is an omnibus test (df = Q-1), while dv_log2ratio below
        //       is contrast-specific. When Q=2 (typical intercept+group) they coincide;
        //       for Q>2 designs the p-value and effect size address different hypotheses.
        double lr = 2.0 * (ll_full[g] - ll_null[g]);
        lr_stat[g] = std::max(lr, 0.0);
        pv_dv[g] = R::pchisq(lr_stat[g], std::max(Q - 1, 1), 0, 0);
        pv_dv[g] = std::max(eps, std::min(1.0, pv_dv[g]));

        // --- DV log-ratio: c'd in log2 scale ---
        double cd = 0;
        for (int q = 0; q < Q; q++) {
            cd += contrast_disp[q] * dg[q];
            d_coef(g, q) = dg[q];
        }
        dv_logratio[g] = cd / std::log(2.0);

        // --- Wald SE for contrast on d-block ---
        // Cov(d_hat) ≈ (W' Ω W)^(-1), Ω_j = alpha_j^2 * trigamma(alpha_j)
        std::vector<double> WtOW(Q*Q, 0.0);
        for (int j = 0; j < N; j++) {
            double aj = std::exp(la_final[j]);
            double omega = aj * aj * std::max(fast_special::trigamma(aj), eps);
            for (int a = 0; a < Q; a++) {
                for (int b = a; b < Q; b++) {
                    WtOW[a*Q+b] += omega * design_disp(j,a) * design_disp(j,b);
                }
            }
        }
        for (int a = 0; a < Q; a++)
            for (int b = 0; b < a; b++) WtOW[a*Q+b] = WtOW[b*Q+a];

        std::vector<double> WtOW_det(WtOW);
        for (int q = 0; q < Q; q++) WtOW_det[q*Q + q] += eps;
        logdet_full_d[g] = dreg_logdet_lu(WtOW_det, Q, eps);

        std::vector<double> rhs(Q, 0.0);
        for (int q = 0; q < Q; q++) rhs[q] = contrast_disp[q];
        std::vector<double> Ase(WtOW);
        bool ok_cov = dreg_solve(Ase, rhs, Q);
        double var_cd = 1e6;
        if (ok_cov) {
            var_cd = 0.0;
            for (int q = 0; q < Q; q++) var_cd += contrast_disp[q] * rhs[q];
            var_cd = std::max(var_cd, eps);
        }
        dv_se_log2[g] = std::sqrt(var_cd) / std::log(2.0);
    }

    return List::create(
        Named("pvalue_dv") = pv_dv,
        Named("lr_stat_dv") = lr_stat,
        Named("ll_null") = ll_null,
        Named("ll_full") = ll_full,
        Named("logdet_full_d") = logdet_full_d,
        Named("logdet_null_d") = logdet_null_d,
        Named("dv_log2ratio") = dv_logratio,
        Named("dv_se_log2") = dv_se_log2,
        Named("d_coefficients") = d_coef);
}
