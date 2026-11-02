// Adapted from mixsqp 0.3-54, distributed under the MIT license.
// Copyright (c) 2017-2020 Youngseok Kim, Peter Carbonetto,
// Matthew Stephens, and Mihai Anitescu.
// The terminal Hessian is deferred; exact Hessian cross-products use Eigen.
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Eigen/Core>
#include <cstdlib>
namespace sgcb_mixsqp {


// FUNCTION DECLARATIONS
// ---------------------
void scalecols          (arma::mat& A, const arma::vec& b);
void normalizerows      (arma::mat& A);
void normalizerowsbymax (arma::mat& A);



// FUNCTION DECLARATIONS
// ---------------------
void mixem_update (const arma::mat& L, const arma::vec& w,
		   arma::vec& x, arma::mat& P);



// FUNCTION DECLARATIONS
// ---------------------
double compute_objective (const arma::mat& L, const arma::vec& w,
			  const arma::vec& x, const arma::vec& z,
			  const arma::vec& e);

double compute_objective (const arma::mat& L, const arma::mat& U,
			  const arma::mat& V, const arma::vec& w,
			  const arma::vec& x, const arma::vec& z,
			  const arma::vec& e, bool usesvd);


using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Scale each column A[,i] by b[i].
void scalecols (mat& A, const vec& b) {
  unsigned int n = A.n_cols;
  for (unsigned int i = 0; i < n; i++)
    A.col(i) *= b[i];
}

// Normalize each row of A so that the entries in each row sum to 1.
void normalizerows (mat& A) {
  vec b = sum(A,1);
  A.each_col() /= b;
}

// Scale each row of A so that the largest entry in each row is 1.
void normalizerowsbymax (mat& A) {
  vec b = max(A,1);
  A.each_col() /= b;
}

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double compute_objective_helper (const vec& u, const vec& w, const vec& z);

// FUNCTION DEFINITIONS
// --------------------
// Compute the value of the (unmodified) objective at x.
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& z, const vec& e) {
  vec u = L*x + e;
  return compute_objective_helper(u,w,z);
}

// Compute the value of the (unmodified) objective at x, using either
// the exact matrix L, or a low-rank SVD approximation of L.
double compute_objective (const mat& L, const mat& U, const mat& V,
			  const vec& w, const vec& x, const vec& z,
			  const vec& e, bool usesvd) {
  vec u;
  if (usesvd)
    u = U*(trans(V)*x);
  else
    u = L*x;
  u += e;
  return compute_objective_helper(u,w,z);
}

// Compute the value of the (unmodified) objective at x, after
// precalculating u = L*x + e.
double compute_objective_helper (const vec& u, const vec& w, const vec& z) {
  if (u.min() <= 0)
    stop("Objective is -Inf");
  return -sum(w % (z + log(u)));
}  

using namespace Rcpp;
using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform several EM updates of the mixture weights.
//
List mixem_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& z,
		 const arma::vec& x0, const arma::vec& eps, int numiter,
		 double zerothresholdsolution, bool verbose) {
  vec obj(numiter);
  vec nnz(numiter);
  vec dmax(numiter);
  int m = L.n_cols;
  mat P = L;
  vec x = x0;
  vec xold(m);
  vec d(m);
  for (unsigned int i = 0; i < numiter; i++) {

    // Store the current estimate of the mixture weights.
    xold = x;
    
    // Run a single EM update.
    mixem_update(L,w,x,P);
    
    // Compute the value of the objective at x.
    obj(i) = compute_objective(L,w,x,z,eps);
    
    // Record the algorithm's progress.
    nnz(i)  = sum(x > zerothresholdsolution);
    d       = abs(x - xold);
    dmax(i) = d.max();
    if (verbose)
      Rprintf("%4d %+0.9e  -- EM -- %4d 1.00e+00 %0.2e  --  --\n",
	      i + 1,obj(i),int(nnz(i)),dmax(i));
  }
  
  return List::create(Named("x")         = x,
		      Named("objective") = obj,
		      Named("nnz")       = nnz,
		      Named("max.diff")  = dmax);
}

// Perform a single EM update of the mixture weights.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;
  
  // Compute the n x m matrix of posterior mixture assignment
  // probabilities (L is an n x m matrix). This is the "E step".
  //
  // The equivalent R code when e = 0 is
  //
  //   P <- t(t(L) * x)
  //   P <- P / rowSums(P)
  // 
  P = L;
  scalecols(P,x + e);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = trans(P) * w;
}
// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.


using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_grad (const mat& L, const mat& U, const mat& V, const vec& w,
		   const vec& x, const vec& e, vec& g, mat& H, mat& Z,
		   bool usesvd);
int  activesetqp  (const mat& H, const vec& g, vec& y, int maxiter,
		   double zerosearchdir, double tol, double ainc);
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B,
				  double ainc);
int  backtracking_line_search (double f, const mat& L, const mat& U,
			       const mat& V, const vec& w, const vec& z,
			       const vec& g, const vec& x, const vec& y,
			       const vec& e, bool usesvd, double suffdecr, 
			       double beta, double amin, double& a, 
			       vec& xnew);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
List solve (const arma::mat& L, const arma::mat& U, const arma::mat& V,
		  const arma::vec& w, const arma::vec& z, const arma::vec& x0, 
		  bool usesvd, bool runem, double convtolsqp, 
		  double convtolactiveset, double zerothresholdsolution, 
		  double zerothresholdsearchdir, double suffdecr, 
		  double stepsizereduce, double minstepsize,
		  double identitycontribincrease, const arma::vec& eps,
		  int maxitersqp, int maxiteractiveset, bool verbose) {
  
  // Get the number columns of the data matrix.
  int m = L.n_cols;

  // PREPARE DATA STRUCTURES
  // -----------------------
  // Initialize storage for the outputs obj, gmin, nnz, nqp and dmax.
  vec obj(maxitersqp,fill::zeros);
  vec gmin(maxitersqp,fill::zeros);
  vec nnz(maxitersqp,fill::zeros);
  vec nqp(maxitersqp);
  vec nls(maxitersqp);
  vec stepsize(maxitersqp);
  vec dmax(maxitersqp);
  nqp.fill(-1);
  nls.fill(-1);
  stepsize.fill(-1);
  dmax.fill(-1);

  // Initialize the solution.
  vec    x = x0;
  double status = 1;
  int    i;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  vec  xold(m);
  vec  g(m);
  vec  ghat(m);
  mat  H(m,m);
  uvec j0(m);
  uvec j1(m);
  vec  y(m);
  vec  d(m);
  vec  xnew(m);
  mat  P = L;
  mat  Z;
  if (usesvd)
    Z = L;
  else
    Z = U;
      
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (i = 0; i < maxitersqp; i++) {

    // Store the current estimate of the mixture weights.
    xold = x;
    
    // Perform a single EM update.
    if (runem)
      mixem_update(L,w,x,P);

    // Find the zero (j0) and non-zero (j1) co-ordinates.
    j0 = find(x <= zerothresholdsolution);
    j1 = find(x > zerothresholdsolution);

    // Zero any co-ordinates that are below the specified threshold.
    x(j0).fill(0);
    
    // Compute the value of the objective at x.
    obj(i) = compute_objective(L,U,V,w,x,z,eps,usesvd);

    // Compute the gradient and scaled likelihood.
    compute_grad(L,U,V,w,x,eps,g,H,Z,usesvd);

    // Report on the algorithm's progress. Here we compute: the
    // smallest gradient value, or, equivalently, dual residual,
    // corresponding to the nonzero co-ordinates (gmin), which is used
    // as a convergence criterion; and the number of nonzeros in the
    // solution (nnz). Note that only the dual residuals (gmin's)
    // corresponding to the nonzero co-ordinates are relevant.
    gmin(i) = 1 + g(j1).min();
    nnz(i)  = sum(x > 0);
    if (verbose) {
      if (i == 0)
        Rprintf("%4d %+0.9e %+0.3e%4d  ------   ------   --  --\n",
		i + 1,obj(i),-gmin(i),int(nnz(i)));
      else
        Rprintf("%4d %+0.9e %+0.3e%4d %0.2e %0.2e %3d %3d\n",i + 1,obj(i),
		-gmin(i),(int) nnz(i),stepsize(i-1),dmax(i-1),
		(int) nqp(i-1),(int) nls(i-1));
    }
    
    // Check convergence. Convergence is reached with the maximum dual
    // residual is small. The negative of "gmin" is also the maximum
    // dual residual (denoted as "rdual" on p. 609 of Boyd &
    // Vandenberghe, "Convex Optimization", 2009).
    if (gmin(i) >= -convtolsqp) {
      status = 0;
      i++;
      break;
    }

    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    checkUserInterrupt();
    
    // Form the Hessian only when a quadratic step is required.
    // The complete cross-product is evaluated without an SVD approximation.
    if (usesvd) {
      H = V * (trans(Z) * Z) * trans(V);
    } else {
      const Eigen::Map<const Eigen::MatrixXd> scaled(Z.memptr(),Z.n_rows,Z.n_cols);
      Eigen::Map<Eigen::MatrixXd> hessian(H.memptr(),H.n_rows,H.n_cols);
      hessian.noalias() = scaled.transpose()*scaled;
    }

    // Solve the quadratic subproblem to obtain a search direction.
    ghat   = g - H*x + 1;
    y      = x;
    nqp(i) = (double) activesetqp(H,ghat,y,maxiteractiveset,
				  zerothresholdsearchdir,convtolactiveset,
				  identitycontribincrease);
    // Extra code used to diagnose line search issues:
    //
    //   if (verbose)
    //     Rprintf("%0.2e\n",norm(x - y,"inf"));
    // 

    // Run backtracking line search.
    nls(i) = (double) backtracking_line_search(obj(i),L,U,V,w,z,g,x,y,eps,
					       usesvd,suffdecr,stepsizereduce,
					       minstepsize,stepsize(i),xnew);
    
    // Update the solution, and store the largest change in the
    // mixture weights.
    d       = abs(xnew - xold);
    dmax(i) = d.max();
    x       = xnew;
  }

  // CONSTRUCT OUTPUT
  // ----------------
  return List::create(Named("x")         = x,
		      Named("status")    = status,
		      Named("objective") = obj.head(i),
		      Named("max.rdual") = -gmin.head(i),
		      Named("nnz")       = nnz.head(i),
		      Named("stepsize")  = stepsize.head(i),
		      Named("max.diff")  = dmax.head(i),
		      Named("nqp")       = nqp.head(i),
		      Named("nls")       = nls.head(i));
}

// Return a or b, which ever is smaller.
inline double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}

// Compute the gradient and scaled likelihood of the objective at x.
void compute_grad (const mat& L, const mat& U, const mat& V, const vec& w,
		   const vec& x, const vec& e, vec& g, mat& H, mat& Z, 
		   bool usesvd) {
  vec u;
  if (usesvd) {
    u = U*(trans(V)*x) + e;
    g = -V * (trans(U) * (w/u));
    Z = U;
    Z.each_col() %= (sqrt(w)/u);

  } else {
    u = L*x + e;
    g = -trans(L) * (w/u);
    Z = L;
    Z.each_col() %= (sqrt(w)/u);

  } 
}

// Return the largest step size maintaining feasibility (x >= 0) for
// the given the search direction (p);
inline void feasible_stepsize (const vec& x, const vec& p, int& j, double& a) {
  uvec i = find(p < 0);
  a = 1;
  j = -1;
  if (!i.is_empty()) {
    vec t = -x(i)/p(i);
    j = t.index_min();
    if (t(j) < 1)
      a = t(j);
    j = i(j);
  }
}

// This implements the active-set method from p. 472 of of Nocedal &
// Wright, Numerical Optimization, 2nd ed, 2006.
int activesetqp (const mat& H, const vec& g, vec& y, int maxiter,
		 double zerosearchdir, double tol, double ainc) {
  int    m = g.n_elem;
  int    k;
  int    iter;
  double a;
  vec    b(m);
  vec    p(m);
  vec    bs(m);
  vec    ps(m);
  mat    Hs(m,m);
  mat    Bs(m,m);
  uvec   i(m);
  uvec   j(m);
  bool   add_to_working_set;

  // This vector is used to keep track of the working set; all zero
  // entries of "t" are co-ordinates belonging to the working set.
  uvec t = (y > 0);
  
  // Run active set method to solve the quadratic subproblem.
  for (iter = 0; iter < maxiter; iter++) {

    // Find the co-ordinates inside (j) and outside (i) the working
    // set.
    i = find(t != 0);
    j = find(t == 0);

    // Make sure the co-ordinates in the working set are set to zero.
    y(j).fill(0);
    
    // Define the smaller quadratic subproblem.
    Hs    = H(i,i);
    b     = g;
    b(i) += Hs*y(i);
    bs    = b(i);

    // Solve the quadratic subproblem to obtain a search direction.
    p.fill(0);
    compute_activeset_searchdir(Hs,bs,ps,Bs,ainc);
    p(i) = ps;

    // Reset the step size.
    a = 1;
    add_to_working_set = false;
    
    // Check that the search direction is close to zero.
    if (norm(p,"inf") <= zerosearchdir) {

      // Calculate b for all co-ordinates.
      b = g + H*y;
      
      // If all the Lagrange multiplers in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (j.is_empty()) {
	iter++;
	break;
      } else if (b(j).min() >= -tol) {
	iter++;
	break;
      } else {

        // Find a co-ordinate with the smallest multiplier, and remove
        // it from the working set.
        k    = j(b(j).index_min());
        t(k) = 1;
      }
      
    // In this next part, we consider adding a co-ordinate to the
    // working set (but only if there are two or more non-zero
    // co-ordinates).
    } else {
      
      // Define the step size.
      feasible_stepsize(y,p,k,a);
      if (k >= 0 && a < 1) {

        // Blocking constraint exists; find and add it to the
        // working set (but only if there are two or more non-zero
        // co-ordinates).
	if (i.n_elem > 1)
	  add_to_working_set = true;
      }
      
      // Move to the new iteration along the search direction. The new
      // iterate must have only positive entries.
      y += a*p;
      j  = find(y < 0);
      y(j).fill(0);

      // Make sure all co-ordinates in the working set are set to
      // zero.
      y(j).fill(0);
      if (add_to_working_set) {
	t(k) = 0;
	y(k) = 0;
      }
    }
  }

  return iter;
}

// Get the initial scalar multiplier for the identity matrix based on
// examining the diagonal entries of the Hessian.
inline double init_hessian_correction (const mat& H, double a0) {
  double d = H.diag().min();
  double a;
  if (d > a0)
    a = 0;
  else
    a = a0 - d;
  return a;
}

// This implements Algorithm 3.3, "Cholesky with added multiple of the
// identity", from Nocedal & Wright, 2nd ed, p. 51.
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p, 
				  mat& B, double ainc) {
  double a0   = 1e-15;
  double amax = 1e15;
  int    n    = y.n_elem;
  mat    I(n,n,fill::eye);
  mat    R(n,n);

  // Get the initial scalar multiplier for the identity matrix.
  double a = init_hessian_correction(H,a0);

  // Repeat until a modified Hessian is found that is symmetric
  // positive definite, or until we cannot modify it any longer.
  while (true) {

    // Compute the modified Hessian.
    B = H + a*I;
    
    // Attempt to compute the Cholesky factorization of the modified
    // Hessian. If this fails, increase the contribution of the
    // identity matrix in the modified Hessian.
    if (a*ainc > amax)
      break;
    else if (chol(R,B))
      break;
    else if (a <= 0)
      a = a0;
    else
      a *= ainc;
  }

  // Compute the search direction using the modified Hessian.
  p = solve(B,-y);
}

// This implements the backtracking line search algorithm from p. 37
// of Nocedal & Wright, Numerical Optimization, 2nd ed, 2006.
// the search direction is given by p = y - x.
int backtracking_line_search (double f, const mat& L, const mat& U,
			      const mat& V, const vec& w, const vec& z,
			      const vec& g, const vec& x, const vec& y,
			      const vec& e, bool usesvd, double suffdecr, 
			      double beta, double amin, double& a, 
			      vec& xnew) {
  int    k;
  double afeas;
  double fnew;
  int    nls = 0;
  
  // Determine the largest step size maintaining feasibility; if it is
  // larger than the minimum step size, return the minimum step size
  // that maintains feasibility of the solution. Otherwise, continue
  // to backtracking line search.
  vec p = y - x;
  feasible_stepsize(x,p,k,afeas);
  if (afeas <= amin) {
    a    = afeas;
    xnew = afeas*y + (1 - afeas)*x;
  } else {

    // Set the initial step size.
    a = min(1,afeas);

    // Iteratively reduce the step size until either (1) we can't reduce
    // any more (because we have hit the minimum step size constraint),
    // or (2) the new candidate solution satisfies the "sufficient
    // decrease" condition.
    while (true) {
      xnew = a*y + (1 - a)*x;
      fnew = compute_objective(L,U,V,w,xnew,z,e,usesvd);
      nls++;

      // Check whether the new candidate solution satisfies the
      // sufficient decrease condition, and remains feasible. If so,
      // accept this candidate solution.
      if ((xnew.min() >= 0) &&
  	  (fnew + sum(xnew) <= f + sum(x) + suffdecr*a*dot(y - x,g + 1)))
        break;

      // If we cannot decrease the step size further, terminate the
      // backtracking line search, and set the step size to be the
      // minimum step size.
      else if (a*beta < amin) {

        // We need to terminate backtracking line search because we have
        // arrived at the smallest allowable step size.
        a    = amin;
        xnew = a*y + (1 - a)*x;
        if (xnew.min() < 0) {
  	  a    = 0;
	  xnew = x;
        }
        break;
      }
    
      // The new candidate does not satisfy the sufficient decrease
      // condition, so we need to try again with a smaller step size.
      a *= beta;
    }
  }

  return nls;
}
} // namespace sgcb_mixsqp

// [[Rcpp::export]]
Rcpp::List mixsqp_eigen_rcpp (const arma::mat& L,const arma::mat& U,const arma::mat& V,
    const arma::vec& w,const arma::vec& z,const arma::vec& x0,bool usesvd,bool runem,
    double convtolsqp,double convtolactiveset,double zerothresholdsolution,
    double zerothresholdsearchdir,double suffdecr,double stepsizereduce,
    double minstepsize,double identitycontribincrease,const arma::vec& eps,
    int maxitersqp,int maxiteractiveset,bool verbose) {
    const char* configured=std::getenv("SGCB_THREADS");
    Eigen::setNbThreads(configured ? std::atoi(configured) : 1);
    return sgcb_mixsqp::solve(L,U,V,w,z,x0,usesvd,runem,convtolsqp,
        convtolactiveset,zerothresholdsolution,zerothresholdsearchdir,suffdecr,
        stepsizereduce,minstepsize,identitycontribincrease,eps,maxitersqp,
        maxiteractiveset,verbose);
}
