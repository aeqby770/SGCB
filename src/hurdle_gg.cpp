#include <Rcpp.h>
#include "source_identity.h"
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
using V3 = std::array<double, 3>;
using M3 = std::array<double, 9>;

// Prentice log-generalized-gamma likelihood, including the Q=0 normal limit.
// Standardized log observations; a weak normal prior regularizes Q only.
static double lgg_constant(double q) {
    const double q2 = q*q;
    if (std::abs(q) < 0.05)
        return -0.9189385332046727 - q2/12.0 + q2*q2*q2/360.0
            - std::pow(q2, 5)/1260.0;
    const double k = 1.0/q2;
    return std::log(std::abs(q)) + k*std::log(k) - R::lgammafn(k) - k;
}

static double lgg_constant_derivative(double q) {
    if (std::abs(q) < 0.05)
        return -q/6.0 + std::pow(q,5)/60.0 - std::pow(q,9)/126.0;
    const double k = 1.0/(q*q);
    return 1.0/q - 2.0/(q*q*q)*(std::log(k)-R::digamma(k));
}

static double point_loss(double y, const V3& p, V3& gradient,
                         double constant, double constant_derivative) {
    const double sigma = std::exp(p[1]);
    const double q = p[2], z = (y-p[0])/sigma, t = q*z;
    double hq2, ezq, dq;
    if (std::abs(t) < 0.005) {
        hq2 = z*z*(0.5+t*(1.0/6.0+t*(1.0/24.0+t*(1.0/120.0+t/720.0))));
        ezq = z*(1.0+t*(0.5+t*(1.0/6.0+t*(1.0/24.0+t/120.0))));
        dq = z*z*z*(1.0/6.0+t*(1.0/12.0+t*(1.0/40.0+t/180.0)));
    } else {
        const double e = std::expm1(t), h = e-t;
        hq2 = h/(q*q);
        ezq = e/q;
        dq = (t*e-2.0*h)/(q*q*q);
    }
    gradient = {-ezq/sigma, 1.0-z*ezq, dq-constant_derivative};
    return p[1]-constant+hq2;
}

static double observation_loss(double y, const V3& p, V3& gradient) {
    return point_loss(y,p,gradient,lgg_constant(p[2]),lgg_constant_derivative(p[2]));
}

static double objective(const std::vector<double>& y, const V3& p, V3& gradient,
                        double prior_weight=0) {
    gradient = {0.0,0.0,0.0};
    double value = 0.0;
    const double constant=lgg_constant(p[2]), derivative=lgg_constant_derivative(p[2]);
    for (double yi : y) {
        V3 gi;
        value += point_loss(yi,p,gi,constant,derivative)/y.size();
        for (int j=0;j<3;++j) gradient[j] += gi[j]/y.size();
    }
    gradient[2] += prior_weight*p[2];
    return value+0.5*prior_weight*p[2]*p[2];
}

static V3 multiply(const M3& a,const V3& b) {
    V3 out = {0,0,0};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) out[i]+=a[3*i+j]*b[j];
    return out;
}
static double dot(const V3& a,const V3& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static double quadratic(const V3& a,const M3& b) { return dot(a,multiply(b,a)); }

static bool inverse(const M3& a,M3& b) {
    const double aa=a[0], ab=a[1], ac=a[2], bb=a[4], bc=a[5], cc=a[8];
    const double det=aa*(bb*cc-bc*bc)-ab*(ab*cc-ac*bc)+ac*(ab*bc-ac*bb);
    if (!(aa>0 && aa*bb-ab*ab>0 && det>0)) return false;
    b={(bb*cc-bc*bc)/det,(ac*bc-ab*cc)/det,(ab*bc-ac*bb)/det,
       (ac*bc-ab*cc)/det,(aa*cc-ac*ac)/det,(ab*ac-aa*bc)/det,
       (ab*bc-ac*bb)/det,(ab*ac-aa*bc)/det,(aa*bb-ab*ab)/det};
    return true;
}

static M3 expected_information(const V3& p) {
    const double sigma=std::exp(p[1]),q=p[2],u=q*q;
    if (q==0) return {1/(sigma*sigma),0,-0.5/sigma,
        0,2,0,-0.5/sigma,0,5.0/12.0};
    double d,ss_minus2,qq;
    if (std::abs(q)<0.2) {
        d=-u/2-u*u/12+std::pow(u,4)/120-std::pow(u,6)/252;
        ss_minus2=u/2+u*u/6-std::pow(u,4)/30+std::pow(u,6)/42+d*d/u+2*d;
        qq=ss_minus2/u+2.0/3-2*u*u/15+2*std::pow(u,4)/21-2*std::pow(u,6)/15;
    } else {
        const double a=1/u,t=R::trigamma(a);
        d=R::digamma(a)-std::log(a);
        ss_minus2=a*(t+d*d)+2*d-1;
        qq=(ss_minus2+4*a*a*t-4*a-2)/u;
    }
    const double mq=-(1+d/u)/sigma,ms=-q*mq,sq=-ss_minus2/q;
    return {1/(sigma*sigma),ms,mq,ms,ss_minus2+2,sq,mq,sq,qq};
}

struct Fit {
    V3 p;
    M3 model, robust, frequentist, expected;
    double loss, start_loss, residual;
    int iterations;
    bool converged, covariance;
};

static Fit fit(const std::vector<double>& raw,double shape_precision=0) {
    const int n=raw.size();
    const double prior_weight=shape_precision/n;
    double center=0, scale=0;
    for (double x : raw) center+=x/n;
    for (double x : raw) scale+=(x-center)*(x-center)/n;
    scale=std::sqrt(scale);
    Fit out;
    out.p={center,std::log(scale),0};
    out.model.fill(NA_REAL); out.robust.fill(NA_REAL); out.frequentist.fill(NA_REAL);
    out.expected.fill(NA_REAL);
    out.converged=false;out.covariance=false;out.iterations=0;
    out.loss=NA_REAL;out.start_loss=NA_REAL;out.residual=NA_REAL;
    if (!(scale>0) || n<4) return out;
    std::vector<double> y(n);
    for (int i=0;i<n;++i) y[i]=(raw[i]-center)/scale;
    V3 p={0,0,0},g;
    M3 b={1,0,0,0,0.5,0,0,0,6};
    double value=objective(y,p,g,prior_weight);
    out.start_loss=value+std::log(scale);
    for (int iteration=0;iteration<400;++iteration) {
        out.residual=std::max({std::abs(g[0]),std::abs(g[1]),std::abs(g[2])});
        if (out.residual<1e-7) {out.converged=true;break;}
        V3 direction=multiply(b,g);
        for (double& d : direction) d=-d;
        double slope=dot(g,direction);
        if (!(slope<0)) {direction={-g[0],-g[1],-g[2]};slope=-dot(g,g);}
        double step=1;
        bool accepted=false;
        V3 candidate,new_gradient;
        double next_value=value;
        for (int backtrack=0;backtrack<35;++backtrack) {
            for (int j=0;j<3;++j) candidate[j]=p[j]+step*direction[j];
            next_value=objective(y,candidate,new_gradient,prior_weight);
            if (std::isfinite(next_value) && next_value<=value+1e-4*step*slope) {
                accepted=true;break;
            }
            step*=0.5;
        }
        if (!accepted) break;
        V3 s,v;
        for (int j=0;j<3;++j) {s[j]=candidate[j]-p[j];v[j]=new_gradient[j]-g[j];}
        const double sv=dot(s,v);
        if (sv>1e-12*std::sqrt(dot(s,s)*dot(v,v))) {
            const V3 bv=multiply(b,v);
            const double factor=(sv+dot(v,bv))/(sv*sv);
            for (int i=0;i<3;++i) for (int j=0;j<3;++j)
                b[3*i+j]+=factor*s[i]*s[j]-(bv[i]*s[j]+s[i]*bv[j])/sv;
        }
        p=candidate;g=new_gradient;value=next_value;
        ++out.iterations;
    }
    out.residual=std::max({std::abs(g[0]),std::abs(g[1]),std::abs(g[2])});
    out.converged=out.residual<1e-7;
    out.p={center+scale*p[0],std::log(scale)+p[1],p[2]};
    out.loss=value+std::log(scale);
    M3 h;
    for (int j=0;j<3;++j) {
        const double delta=1e-4*(1+std::abs(p[j]));
        V3 pp=p,pm=p,gp,gm;
        pp[j]+=delta;pm[j]-=delta;
        objective(y,pp,gp,prior_weight);objective(y,pm,gm,prior_weight);
        for (int i=0;i<3;++i) h[3*i+j]=(gp[i]-gm[i])/(2*delta);
    }
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j)
        h[3*i+j]=h[3*j+i]=(h[3*i+j]+h[3*j+i])/2;
    M3 hi;
    out.covariance=inverse(h,hi);
    if (!out.covariance) return out;
    M3 meat={0,0,0,0,0,0,0,0,0};
    V3 likelihood_gradient;
    objective(y,p,likelihood_gradient);
    for (double yi : y) {
        V3 score;observation_loss(yi,p,score);
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            meat[3*i+j]+=(score[i]-likelihood_gradient[i])*(score[j]-likelihood_gradient[j])/n;
    }
    const V3 transform={scale,1,1};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        double sandwich=0;
        for (int k=0;k<3;++k) for (int l=0;l<3;++l)
            sandwich+=hi[3*i+k]*meat[3*k+l]*hi[3*l+j];
        out.model[3*i+j]=hi[3*i+j]*transform[i]*transform[j]/n;
        // Sampling covariance of the penalized estimator: Hpen^-1 Hlik Hpen^-1.
        // The inverse penalized Hessian alone is a posterior covariance.
        out.frequentist[3*i+j]=(hi[3*i+j]-prior_weight*hi[3*i+2]*hi[6+j])
            *transform[i]*transform[j]/n;
        out.robust[3*i+j]=sandwich*transform[i]*transform[j]/(n-3);
    }
    M3 information=expected_information(out.p),information_inverse;
    information[8]+=prior_weight;
    if (inverse(information,information_inverse)) {
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            out.expected[3*i+j]=(information_inverse[3*i+j]
                -prior_weight*information_inverse[3*i+2]*information_inverse[6+j])/n;
    }
    return out;
}

static double log_moment(double r,const V3& p) {
    const double sigma=std::exp(p[1]),q=p[2],rs=r*sigma;
    if (q==0) return r*p[0]+0.5*rs*rs;
    const double k=1/(q*q),t=rs/q;
    if (!(k+t>0)) return NA_REAL;
    if (std::abs(q)>=0.05)
        return r*p[0]+R::lgammafn(k+t)-R::lgammafn(k)-t*std::log(k);
    const double nodes[8]={0.0198550717512319,0.101666761293187,0.237233795041836,
        0.408282678752175,0.591717321247825,0.762766204958164,0.898333238706813,0.980144928248768};
    const double weights[8]={0.0506142681451881,0.111190517226687,0.156853322938944,
        0.181341891689181,0.181341891689181,0.156853322938944,0.111190517226687,0.0506142681451881};
    double integral=0;
    for (int j=0;j<8;++j) {
        const double v=t*nodes[j]/k,invx=1/(k+t*nodes[j]);
        const double psi_diff=std::log1p(v)-0.5*invx-invx*invx/12
            +std::pow(invx,4)/120-std::pow(invx,6)/252;
        integral+=weights[j]*t*psi_diff;
    }
    return r*p[0]+integral;
}

static double logcv2(const V3& p) {
    const double v=log_moment(2,p)-2*log_moment(1,p);
    return v+std::log(-std::expm1(-v));
}
static double logvariance(const V3& p) {
    const double q=p[2],q2=q*q;
    if (std::abs(q)<0.01)
        return 2*p[1]+std::log1p(q2/2+q2*q2/6-std::pow(q2,4)/30);
    return 2*p[1]+std::log(R::trigamma(1/q2)/q2);
}
static V3 log_moment_gradient(double r,const V3& p) {
    const double sigma=std::exp(p[1]),q=p[2],rs=r*sigma;
    if (q==0) return {r,rs*rs,-rs/2-rs*rs*rs/6};
    const double k=1/(q*q),t=rs/q;
    if (!(k+t>0)) return {r,NA_REAL,NA_REAL};
    if (std::abs(q)>=0.05) {
        const double a=R::digamma(k+t)-std::log(k);
        return {r,t*a,(-2*k*(R::digamma(k+t)-R::digamma(k))-t*a+2*t)/q};
    }
    const double nodes[8]={0.0198550717512319,0.101666761293187,0.237233795041836,
        0.408282678752175,0.591717321247825,0.762766204958164,0.898333238706813,0.980144928248768};
    const double weights[8]={0.0506142681451881,0.111190517226687,0.156853322938944,
        0.181341891689181,0.181341891689181,0.156853322938944,0.111190517226687,0.0506142681451881};
    V3 out={r,0,0};
    for (int j=0;j<8;++j) {
        const double v=nodes[j]*rs*q,inv=q*q/(1+v),w=v/(1+v);
        const double c=-inv/2-inv*inv/12+std::pow(inv,4)/120-std::pow(inv,6)/252;
        const double cp=-0.5-inv/6+std::pow(inv,3)/30-std::pow(inv,5)/42;
        const double w_log=std::abs(v)<0.005 ? v*v*(-0.5+v*(2.0/3+v*(-0.75+
            v*(0.8+v*(-5.0/6+v*(6.0/7-v*0.875)))))) : w-std::log1p(v);
        out[1]+=weights[j]*t*(std::log1p(v)+c+w*(1-inv*cp));
        out[2]+=weights[j]*(t/q)*(w_log-c+inv*cp*(2-w));
    }
    return out;
}
static V3 log_moment_ratio_gradient(const V3& p) {
    const V3 first=log_moment_gradient(1,p),second=log_moment_gradient(2,p);
    return {0,second[1]-2*first[1],second[2]-2*first[2]};
}
static V3 effect_gradient(const V3& p,int effect) {
    if (effect==0) {
        V3 out=log_moment_ratio_gradient(p);
        const double denominator=-std::expm1(2*log_moment(1,p)-log_moment(2,p));
        out[1]/=denominator;out[2]/=denominator;
        return out;
    }
    const double q=p[2],u=q*q;
    if (std::abs(q)<0.01)
        return {0,2,(q+2*q*u/3-4*q*std::pow(u,3)/15)/(1+u/2+u*u/6-std::pow(u,4)/30)};
    const double k=1/u;
    return {0,2,(-2*k*R::psigamma(k,2)/R::trigamma(k)-2)/q};
}

// A point mass at zero and a regularized log-generalized-gamma positive part.
// The Bernoulli frequency and positive-part parameters have orthogonal scores.
// The CV-squared target refers to the entire non-negative abundance distribution.
// [[Rcpp::export]]
Rcpp::NumericMatrix hurdle_gg_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    double shape_precision=1,int threads=1) {
    const int ng=x.nrow(),ns=x.ncol();
    const std::vector<std::string> columns={"mu0","sigma0","Q0","mu1","sigma1","Q1",
        "converged0","converged1","iterations0","iterations1","score0","score1",
        "covariance0","covariance1","positive0","positive1","n_positive0","n_positive1",
        "logcv2_0","logcv2_1","logvar_positive0","logvar_positive1",
        "cv_var0_model","cv_var1_model","cv_var0_sandwich","cv_var1_sandwich",
        "lv_var0_model","lv_var1_model","lv_var0_sandwich","lv_var1_sandwich",
        "shape_var0_model","shape_var1_model","shape_var0_sandwich","shape_var1_sandwich",
        "fitted_mean0","fitted_mean1","cv_var0_penalized","cv_var1_penalized",
        "cv_var0_expected","cv_var1_expected"};
    Rcpp::NumericMatrix out(ng,columns.size());
    out.attr("native_source_id")=SGCB_SOURCE_ID;
    Rcpp::colnames(out)=Rcpp::wrap(columns);std::fill(out.begin(),out.end(),NA_REAL);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,16)
    #endif
    for (int gene=0;gene<ng;++gene) {
        std::vector<double> y[2];int n[2]={0,0};
        for (int j=0;j<ns;++j) {
            const int k=group[j];++n[k];
            if (x(gene,j)>0) y[k].push_back(std::log(x(gene,j)));
        }
        for (int k=0;k<2;++k) {
            const Fit a=fit(y[k],shape_precision);
            const double prob=static_cast<double>(y[k].size())/n[k];
            out(gene,3*k)=a.p[0];out(gene,3*k+1)=std::exp(a.p[1]);out(gene,3*k+2)=a.p[2];
            out(gene,6+k)=a.converged;out(gene,8+k)=a.iterations;out(gene,10+k)=a.residual;
            out(gene,12+k)=a.covariance;out(gene,14+k)=prob;out(gene,16+k)=y[k].size();
            if (!a.converged || !a.covariance) continue;
            const double logR=log_moment(2,a.p)-2*log_moment(1,a.p)-std::log(prob);
            out(gene,18+k)=logR+std::log(-std::expm1(-logR));
            out(gene,20+k)=logvariance(a.p);
            const double factor=1/(-std::expm1(-logR));
            V3 dc=log_moment_ratio_gradient(a.p),dv=effect_gradient(a.p,1);
            for (double& v: dc) v*=factor;
            const double frequency_variance=factor*factor*(1-prob)/(n[k]*prob);
            out(gene,22+k)=quadratic(dc,a.model)+frequency_variance;
            out(gene,24+k)=quadratic(dc,a.robust)+frequency_variance;
            out(gene,26+k)=quadratic(dv,a.model);out(gene,28+k)=quadratic(dv,a.robust);
            out(gene,30+k)=a.model[8];out(gene,32+k)=a.robust[8];
            out(gene,34+k)=prob*std::exp(log_moment(1,a.p));
            out(gene,36+k)=quadratic(dc,a.frequentist)+frequency_variance;
            out(gene,38+k)=quadratic(dc,a.expected)+frequency_variance;
        }
    }
    return out;
}
