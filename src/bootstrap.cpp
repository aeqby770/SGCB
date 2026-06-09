// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// =============================================================================
// LFC shrinkage (used by sgcbDE pipeline)
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector shrink_lfc_normal(Rcpp::NumericVector lfc,
                                        Rcpp::NumericVector lfc_se,
                                        double prior_scale = 1.0) {
  int n = lfc.size();
  Rcpp::NumericVector shrunk(n);
  double prior_var = prior_scale * prior_scale;
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double se2 = lfc_se[i] * lfc_se[i];
    double shrinkage = prior_var / (prior_var + se2);
    shrunk[i] = lfc[i] * shrinkage;
  }
  
  return shrunk;
}

// [[Rcpp::export]]
Rcpp::NumericVector shrink_lfc_cauchy(Rcpp::NumericVector lfc,
                                        Rcpp::NumericVector lfc_se,
                                        double prior_scale = 1.0) {
  int n = lfc.size();
  Rcpp::NumericVector shrunk(n);
  double prior_var = prior_scale * prior_scale;
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < n; i++) {
    double se2 = lfc_se[i] * lfc_se[i];
    double u = lfc[i] / prior_scale;
    double weight = 1.0 / (1.0 + u * u);
    double shrinkage = se2 / (se2 + prior_var);
    shrunk[i] = lfc[i] * (1.0 - weight * shrinkage);
  }
  
  return shrunk;
}
