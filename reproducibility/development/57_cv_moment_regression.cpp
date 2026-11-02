#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

// Log-linear conditional mean and CV-squared. The second equation uses
// squared relative residuals, so the fitted variance remains positive.
// The joint influence function includes uncertainty from estimating the
// mean before donor clustering; no likelihood for squared data is assumed.
static bool chol_inverse(const std::vector<double>& a,std::vector<double>& inverse,int p) {
    std::vector<double> l(p*p,0),linv(p*p,0);
    for (int i=0;i<p;++i) for (int j=0;j<=i;++j) {
        double value=a[i*p+j];
        for (int k=0;k<j;++k) value-=l[i*p+k]*l[j*p+k];
        if (i==j) {
            if (!(value>0) || !std::isfinite(value)) return false;
            l[i*p+j]=std::sqrt(value);
        } else l[i*p+j]=value/l[j*p+j];
    }
    for (int col=0;col<p;++col) for (int i=0;i<p;++i) {
        double value=i==col ? 1 : 0;
        for (int j=0;j<i;++j) value-=l[i*p+j]*linv[j*p+col];
        linv[i*p+col]=value/l[i*p+i];
    }
    inverse.assign(p*p,0);
    for (int i=0;i<p;++i) for (int j=0;j<p;++j)
        for (int k=0;k<p;++k) inverse[i*p+j]+=linv[k*p+i]*linv[k*p+j];
    return true;
}

struct MomentFit {
    std::vector<double> beta,inverse,fitted;
    bool converged;
    int iterations;
};

static double moment_loss(const std::vector<double>& y,const std::vector<double>& x,
    const std::vector<double>& beta,int n,int p) {
    double loss=0;
    for (int i=0;i<n;++i) {
        double eta=0;
        for (int j=0;j<p;++j) eta+=x[i*p+j]*beta[j];
        loss+=std::exp(eta)-y[i]*eta;
    }
    return loss;
}

static MomentFit fit_moment(const std::vector<double>& y,const std::vector<double>& x,int n,int p) {
    MomentFit out;out.beta.assign(p,0);out.inverse.assign(p*p,NA_REAL);
    out.fitted.assign(n,NA_REAL);out.converged=false;out.iterations=0;
    double sum0=0,sum1=0,n0=0,n1=0;
    for (int i=0;i<n;++i) {
        if (x[i*p+1]==0) {sum0+=y[i];++n0;} else {sum1+=y[i];++n1;}
    }
    if (!(sum0>0 && sum1>0)) return out;
    out.beta[0]=std::log(sum0/n0);
    out.beta[1]=std::log(sum1/n1)-out.beta[0];
    double value=moment_loss(y,x,out.beta,n,p);
    for (int iter=0;iter<80;++iter) {
        out.iterations=iter;
        std::vector<double> gradient(p,0),hessian(p*p,0),inverse,step(p,0);
        for (int i=0;i<n;++i) {
            double eta=0;
            for (int j=0;j<p;++j) eta+=x[i*p+j]*out.beta[j];
            const double m=std::exp(eta);
            out.fitted[i]=m;
            for (int j=0;j<p;++j) {
                gradient[j]+=(m-y[i])*x[i*p+j];
                for (int k=0;k<p;++k) hessian[j*p+k]+=m*x[i*p+j]*x[i*p+k];
            }
        }
        if (!chol_inverse(hessian,inverse,p)) return out;
        out.inverse=inverse;
        double decrement=0;
        for (int j=0;j<p;++j) {
            for (int k=0;k<p;++k) step[j]+=inverse[j*p+k]*gradient[k];
            decrement+=step[j]*gradient[j];
        }
        double maximum_step=0;
        for (double s:step) maximum_step=std::max(maximum_step,std::abs(s));
        if (maximum_step<1e-8) {out.converged=true;return out;}
        bool accepted=false;double length=1;
        std::vector<double> candidate(p);
        for (int search=0;search<35;++search) {
            for (int j=0;j<p;++j) candidate[j]=out.beta[j]-length*step[j];
            const double next=moment_loss(y,x,candidate,n,p);
            const double roundoff=8*std::numeric_limits<double>::epsilon()*(1+std::abs(value));
            if (std::isfinite(next) && next<=value-1e-4*length*decrement+roundoff) {
                out.beta=candidate;value=next;accepted=true;break;
            }
            length*=0.5;
        }
        if (!accepted) return out;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cv_moment_regression_cpp(Rcpp::NumericMatrix abundance,
    Rcpp::NumericMatrix design,Rcpp::IntegerVector cluster,int threads=1) {
    const int ng=abundance.nrow(),n=abundance.ncol(),p=design.ncol();
    const int clusters=*std::max_element(cluster.begin(),cluster.end())+1;
    std::vector<double> x(n*p);
    for (int i=0;i<n;++i) for (int j=0;j<p;++j) x[i*p+j]=design(i,j);
    Rcpp::NumericMatrix output(ng,10);
    std::fill(output.begin(),output.end(),NA_REAL);
    Rcpp::colnames(output)=Rcpp::CharacterVector::create("cv2_group0","cv2_group1",
        "log2_cv2_ratio","se_log2_cv2","pvalue","mean0_at_reference","mean1_at_reference",
        "first_moment_converged","variance_model_converged","moments_consistent");
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,16)
    #endif
    for (int gene=0;gene<ng;++gene) {
        std::vector<double> y(n),relative_squared(n),ratio(n);
        double scale=0;
        for (int i=0;i<n;++i) scale+=abundance(gene,i)/n;
        if (!(scale>0)) continue;
        for (int i=0;i<n;++i) y[i]=abundance(gene,i)/scale;
        const MomentFit first=fit_moment(y,x,n,p);
        output(gene,7)=first.converged;
        if (!first.converged) continue;
        for (int i=0;i<n;++i) {
            ratio[i]=y[i]/first.fitted[i];
            relative_squared[i]=(ratio[i]-1)*(ratio[i]-1);
        }
        const MomentFit second=fit_moment(relative_squared,x,n,p);
        output(gene,8)=second.converged;
        if (!second.converged) continue;
        output(gene,0)=std::exp(second.beta[0]);
        output(gene,1)=std::exp(second.beta[0]+second.beta[1]);
        output(gene,5)=scale*std::exp(first.beta[0]);
        output(gene,6)=scale*std::exp(first.beta[0]+first.beta[1]);
        output(gene,9)=1;
        std::vector<double> cross(p*p,0),mean_gradient(p,0),h1(p,0),h2(p,0),score(clusters,0);
        for (int j=0;j<p;++j) h2[j]=second.inverse[j*p+1];
        for (int i=0;i<n;++i) for (int j=0;j<p;++j) for (int k=0;k<p;++k)
            cross[j*p+k]+=2*ratio[i]*(ratio[i]-1)*x[i*p+j]*x[i*p+k];
        for (int j=0;j<p;++j) for (int k=0;k<p;++k) mean_gradient[j]-=cross[k*p+j]*h2[k];
        for (int j=0;j<p;++j) for (int k=0;k<p;++k)
            h1[j]+=first.inverse[j*p+k]*mean_gradient[k];
        for (int i=0;i<n;++i) {
            double a1=0,a2=0;
            for (int j=0;j<p;++j) {a1+=x[i*p+j]*h1[j];a2+=x[i*p+j]*h2[j];}
            score[cluster[i]]+=(y[i]-first.fitted[i])*a1+(relative_squared[i]-second.fitted[i])*a2;
        }
        double variance=0;
        for (double s:score) variance+=s*s;
        variance*=static_cast<double>(clusters)/(clusters-1)*(n-1)/(n-p);
        const double difference=second.beta[1];
        output(gene,2)=difference/std::log(2.0);
        output(gene,3)=std::sqrt(variance)/std::log(2.0);
        output(gene,4)=2*R::pt(-std::abs(difference)/std::sqrt(variance),clusters-1,1,0);
    }
    return output;
}
