// Population likelihoods and zero-aware relative-variability inference.
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
#include <limits>
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

using V5=std::array<double,5>;
using M5=std::array<double,25>;
static constexpr double profile_tolerance=1e-7;
static constexpr double moment_tolerance=1e-12;
static constexpr double allocation_floor=1e-12;
static constexpr double boundary_roundoff=32*std::numeric_limits<double>::epsilon();

static double softplus(double x) {
    return std::max(x,0.0)+std::log1p(std::exp(-std::abs(x)));
}

// K=log(E[X^2|X>0])-2log(E[X|X>0]) increases with positive-part scale.
static double sigma_for_K(double target,double q) {
    if (!(target>0) || !std::isfinite(target)) return NA_REAL;
    if (q==0) return std::sqrt(target);
    double lower=0,upper=q<0 ? std::nextafter(-0.5/q,0.0) : std::sqrt(target);
    for (int i=0;i<60 && q>0;++i) {
        const V3 p={0,std::log(upper),q};
        const double value=log_moment(2,p)-2*log_moment(1,p);
        if (!std::isfinite(value) || value>=target) break;
        upper*=2;
    }
    double sigma=std::min(std::sqrt(target),upper*0.5);
    for (int i=0;i<80;++i) {
        const V3 p={0,std::log(sigma),q};
        const double value=log_moment(2,p)-2*log_moment(1,p);
        if (std::isfinite(value) && std::abs(value-target)<moment_tolerance*(1+target))
            return sigma;
        if (!std::isfinite(value) || value>target) upper=sigma; else lower=sigma;
        const double derivative=log_moment_ratio_gradient(p)[1]/sigma;
        double next=sigma-(value-target)/derivative;
        if (!(next>lower && next<upper) || !std::isfinite(next)) next=(lower+upper)/2;
        sigma=next;
    }
    return NA_REAL;
}

// Location is a nuisance for CV2 and has an analytic positive-likelihood maximum.
struct PositiveProfile {
    V3 p={NA_REAL,NA_REAL,NA_REAL},gradient={NA_REAL,NA_REAL,NA_REAL};
    double loss=NA_REAL;
};
static PositiveProfile positive_profile(const std::vector<double>& y,double sigma,
    double q,double precision) {
    PositiveProfile out;
    double center=0;
    for (double yi:y) center+=yi/y.size();
    double location=center;
    if (q!=0) {
        double maximum=R_NegInf;
        for (double yi:y) maximum=std::max(maximum,q*(yi-center)/sigma);
        double sum=0;
        if (maximum<0.5) {
            for (double yi:y) sum+=std::expm1(q*(yi-center)/sigma)/y.size();
            location+=sigma/q*std::log1p(sum);
        } else {
            for (double yi:y) sum+=std::exp(q*(yi-center)/sigma-maximum)/y.size();
            location+=sigma/q*(maximum+std::log(sum));
        }
    }
    out.p={location,std::log(sigma),q};
    out.loss=y.size()*objective(y,out.p,out.gradient,precision/y.size());
    for (double& v:out.gradient) v*=y.size();
    return out;
}

struct ProfileData { std::vector<double> y[2]; int n[2]; double delta,precision; };
struct ProfileState {
    PositiveProfile positive[2];
    double pi[2]={NA_REAL,NA_REAL},phi[2]={NA_REAL,NA_REAL};
};

// eta=log(CV0^2), CV1^2=exp(eta+delta).  K_g=u_g*log(1+CV_g^2),
// pi_g=exp(-(1-u_g)*log(1+CV_g^2)), 0<u_g<=1.  This enforces the
// exact common-CV2 constraint without forcing equal positivity or GG shape.
static double profile_objective(const ProfileData& data,const V5& p,V5& gradient,
    ProfileState* saved=nullptr) {
    gradient.fill(0);
    ProfileState state;
    double loss=0;
    for (int g=0;g<2;++g) {
        const int ui=1+2*g,qi=2+2*g,m=data.y[g].size(),zeros=data.n[g]-m;
        const double eta=p[0]+g*data.delta,phi=softplus(eta),u=p[ui],q=p[qi];
        if (!(u>=allocation_floor && u<=1)) return R_PosInf;
        const double sigma=sigma_for_K(phi*u,q),logpi=-phi*(1-u);
        if (!std::isfinite(sigma) || (zeros>0 && logpi>=0)) return R_PosInf;
        const PositiveProfile positive=positive_profile(data.y[g],sigma,q,data.precision);
        if (!std::isfinite(positive.loss)) return R_PosInf;
        const V3 dk=log_moment_ratio_gradient(positive.p);
        const double positive_K=positive.gradient[1]/dk[1];
        const double bernoulli_gradient=zeros==0 ? -m : -m+zeros/std::expm1(-logpi);
        loss+=positive.loss-m*logpi;
        if (zeros>0) loss-=zeros*std::log(-std::expm1(logpi));
        gradient[0]+=(positive_K*u+bernoulli_gradient*(u-1))*std::exp(eta-phi);
        gradient[ui]=phi*(positive_K+bernoulli_gradient);
        gradient[qi]=positive.gradient[2]-positive.gradient[1]*dk[2]/dk[1];
        state.positive[g]=positive;state.pi[g]=std::exp(logpi);state.phi[g]=phi;
    }
    const double total=data.n[0]+data.n[1];
    for (double& v:gradient) v/=total;
    if (saved) *saved=state;
    return loss/total;
}

static V5 projected_gradient(const V5& p,const V5& gradient) {
    V5 out=gradient;
    for (int j:{1,3}) {
        if ((p[j]>=1 && gradient[j]<0) ||
            (p[j]<=allocation_floor && gradient[j]>0)) out[j]=0;
    }
    return out;
}
static double dot5(const V5& a,const V5& b) {
    double out=0;for (int j=0;j<5;++j) out+=a[j]*b[j];return out;
}
struct ProfileFit { V5 p; ProfileState state; double loss,residual; int iterations; bool converged; };
static ProfileFit fit_profile(const ProfileData& data,const V5& start) {
    ProfileFit out;out.p=start;out.iterations=0;out.converged=false;
    M5 inverse_hessian={0};for (int j=0;j<5;++j) inverse_hessian[5*j+j]=1;
    V5 gradient;
    double value=profile_objective(data,out.p,gradient);
    for (int iteration=0;iteration<400;++iteration) {
        const V5 pg=projected_gradient(out.p,gradient);
        out.residual=0;for (double v:pg)
            out.residual=std::isfinite(v) ? std::max(out.residual,std::abs(v)) : R_PosInf;
        if (out.residual<profile_tolerance && std::isfinite(value)) break;
        // An active bound is not a free nuisance direction.  Use the reduced
        // inverse-Hessian block instead of mixing its nonzero KKT multiplier.
        for (int j:{1,3}) if (pg[j]==0 && gradient[j]!=0) {
            for (int k=0;k<5;++k) inverse_hessian[5*j+k]=inverse_hessian[5*k+j]=0;
            inverse_hessian[5*j+j]=1;
        }
        V5 direction={0};
        for (int j=0;j<5;++j) for (int k=0;k<5;++k)
            direction[j]-=inverse_hessian[5*j+k]*pg[k];
        for (int j:{1,3}) {
            if ((out.p[j]>=1 && direction[j]>0) ||
                (out.p[j]<=allocation_floor && direction[j]<0)) direction[j]=0;
        }
        double slope=dot5(gradient,direction);
        if (!(slope<0)) {
            for (int j=0;j<5;++j) direction[j]=-pg[j];
            slope=-dot5(pg,pg);
        }
        double step=1;
        for (int j:{1,3}) {
            if (direction[j]>0) step=std::min(step,(1-out.p[j])/direction[j]);
            if (direction[j]<0) step=std::min(step,(allocation_floor-out.p[j])/direction[j]);
        }
        V5 next,new_gradient;double next_value=value;bool accepted=false;
        for (int backtrack=0;backtrack<35;++backtrack) {
            for (int j=0;j<5;++j) next[j]=out.p[j]+step*direction[j];
            for (int j:{1,3}) next[j]=std::min(1.0,std::max(allocation_floor,next[j]));
            for (int g=0;g<2;++g) if (data.n[g]==int(data.y[g].size()) &&
                next[1+2*g]>1-boundary_roundoff) next[1+2*g]=1;
            next_value=profile_objective(data,next,new_gradient);
            if (std::isfinite(next_value) && next_value<=value+1e-4*step*slope) {
                accepted=true;break;
            }
            step*=0.5;
        }
        if (!accepted) break;
        V5 change,gradient_change,hg={0};
        for (int j=0;j<5;++j) {change[j]=next[j]-out.p[j];gradient_change[j]=new_gradient[j]-gradient[j];}
        const V5 next_projected=projected_gradient(next,new_gradient);
        for (int j:{1,3}) if (pg[j]==0 && gradient[j]!=0 &&
            next_projected[j]==0 && new_gradient[j]!=0) gradient_change[j]=0;
        const double sy=dot5(change,gradient_change);
        if (sy>1e-12*std::sqrt(dot5(change,change)*dot5(gradient_change,gradient_change))) {
            for (int j=0;j<5;++j) for (int k=0;k<5;++k)
                hg[j]+=inverse_hessian[5*j+k]*gradient_change[k];
            const double factor=(sy+dot5(gradient_change,hg))/(sy*sy);
            for (int j=0;j<5;++j) for (int k=0;k<5;++k)
                inverse_hessian[5*j+k]+=factor*change[j]*change[k]-(hg[j]*change[k]+change[j]*hg[k])/sy;
        }
        out.p=next;gradient=new_gradient;value=next_value;++out.iterations;
    }
    out.loss=profile_objective(data,out.p,gradient,&out.state)*(data.n[0]+data.n[1]);
    const V5 pg=projected_gradient(out.p,gradient);
    out.residual=0;for (double v:pg)
        out.residual=std::isfinite(v) ? std::max(out.residual,std::abs(v)) : R_PosInf;
    out.converged=std::isfinite(out.loss) && out.residual<profile_tolerance;
    return out;
}

// A quadratic penalized likelihood ratio has scale V_sampling/V_curvature
// in the tested scalar direction.  Evaluate both at the constrained fit.
static double profile_penalty_scale(const ProfileData& data,const ProfileFit& fit) {
    if (!fit.converged) return NA_REAL;
    double sampling=0,curvature=0;
    for (int g=0;g<2;++g) {
        const double m=data.y[g].size(),prior=data.precision/m,pi=fit.state.pi[g];
        const V3 p=fit.state.positive[g].p;
        M3 info=expected_information(p),inv;info[8]+=prior;
        if (!inverse(info,inv)) return NA_REAL;
        V3 dc=log_moment_ratio_gradient(p);
        const double factor=1/(-std::expm1(-fit.state.phi[g]));
        for (double& v:dc) v*=factor;
        const double bernoulli=factor*factor*(1-pi)/(data.n[g]*pi);
        const double curve=quadratic(dc,inv)/m+bernoulli;
        const double shape=dc[0]*inv[2]+dc[1]*inv[5]+dc[2]*inv[8];
        curvature+=curve;sampling+=curve-prior*shape*shape/m;
    }
    return sampling/curvature;
}


struct ScaleKernel {double log_density,scale_score,shape_score;};

// Integrate location with prior d(mu) d(log sigma). Constants independent of
// sigma cancel here. The Q=0 limit is sigma^2 = SS / chi-square_(m-1).
static ScaleKernel marginal_scale_kernel(const std::vector<double>& centered,double s,double q) {
    const double m=centered.size(),inverse_sigma=std::exp(-s);
    double largest=0;
    for (double value:centered) largest=std::max(largest,std::abs(q*value*inverse_sigma));
    if (largest<0.005) {
        double a2=0,a3=0,a4=0,a5=0,a6=0;
        for (double value:centered) {
            const double z=value*inverse_sigma,z2=z*z;
            a2+=z2/m;a3+=z2*z/m;a4+=z2*z2/m;
            a5+=z2*z2*z/m;a6+=z2*z2*z2/m;
        }
        const double c4=a4-3*a2*a2,c5=a5-10*a2*a3;
        const double c6=a6-15*a4*a2-10*a3*a3+30*a2*a2*a2;
        const double ratio=a2/2+q*(a3/6+q*(c4/24+q*(c5/120+q*c6/720)));
        const double scale=a2+q*(a3/2+q*(c4/6+q*(c5/24+q*c6/120)));
        const double shape=a3/6+q*(c4/12+q*(c5/40+q*c6/180));
        return {-(m-1)*s-m*ratio,-(m-1)+m*scale,-m*shape};
    }
    double maximum=R_NegInf;
    for (double value:centered) maximum=std::max(maximum,q*value*inverse_sigma);
    double sum=0,weighted=0,small=0;
    for (double value:centered) {
        const double z=q*value*inverse_sigma,w=std::exp(z-maximum);
        sum+=w;weighted+=w*z;
        if (largest<0.5) small+=std::expm1(z)/m;
    }
    const double logmean=largest<0.5 ? std::log1p(small) : maximum+std::log(sum/m);
    const double ez=weighted/sum,k=1/(q*q);
    return {-(m-1)*s-m*k*logmean,-(m-1)+m*k*ez,m*k*(2*logmean-ez)/q};
}

static void legendre_rule(int n,std::vector<double>& nodes,std::vector<double>& weights) {
    nodes.resize(n);weights.resize(n);
    for (int j=0;j<(n+1)/2;++j) {
        double z=std::cos(3.14159265358979323846*(j+0.75)/(n+0.5)),derivative=0;
        for (int iteration=0;iteration<30;++iteration) {
            double previous=1,current=z;
            for (int k=2;k<=n;++k) {
                const double next=((2*k-1)*z*current-(k-1)*previous)/k;
                previous=current;current=next;
            }
            derivative=n*(z*current-previous)/(z*z-1);
            const double step=current/derivative;z-=step;
            if (std::abs(step)<2e-15) break;
        }
        nodes[j]=-z;nodes[n-1-j]=z;
        weights[j]=weights[n-1-j]=2/((1-z*z)*derivative*derivative);
    }
}

struct ScaleMoments {double mean,variance,shape_derivative,normal_error,log_integral;};
static ScaleMoments marginal_scale_moments(const std::vector<double>& y,double q,
    const std::vector<double>& nodes,const std::vector<double>& weights) {
    const double m=y.size();double center=0,ss=0;
    for (double value:y) center+=value/m;
    std::vector<double> centered(m);
    for (int j=0;j<m;++j) {centered[j]=y[j]-center;ss+=centered[j]*centered[j];}
    double lower=0.5*std::log(ss/(m-1)),upper=lower;
    for (int j=0;j<80 && marginal_scale_kernel(centered,lower,q).scale_score<0;++j) lower-=0.5;
    for (int j=0;j<80 && marginal_scale_kernel(centered,upper,q).scale_score>0;++j) upper+=0.5;
    for (int j=0;j<60;++j) {
        const double middle=(lower+upper)/2;
        if (marginal_scale_kernel(centered,middle,q).scale_score>0) lower=middle; else upper=middle;
    }
    const double mode=(lower+upper)/2,peak=marginal_scale_kernel(centered,mode,q).log_density;
    for (int j=0;j<120 && marginal_scale_kernel(centered,lower,q).log_density>peak-40;++j) lower-=0.5;
    for (int j=0;j<120 && marginal_scale_kernel(centered,upper,q).log_density>peak-40;++j) upper+=0.5;
    const double half=(upper-lower)/2,midpoint=(upper+lower)/2;
    double mass=0,first=0,second=0,qfirst=0,cross=0;
    for (int j=0;j<int(nodes.size());++j) {
        const double s=midpoint+half*nodes[j];
        const ScaleKernel value=marginal_scale_kernel(centered,s,q);
        const double w=weights[j]*std::exp(value.log_density-peak),deviation=s-mode;
        mass+=w;first+=w*deviation;second+=w*deviation*deviation;
        qfirst+=w*value.shape_score;cross+=w*deviation*value.shape_score;
    }
    const double mean=mode+first/mass,variance=second/mass-std::pow(first/mass,2);
    const double log_integral=peak+std::log(half*mass);
    double error=0;
    if (q==0) error=std::max(std::abs(mean-.5*(std::log(ss/2)-R::digamma((m-1)/2))),
        std::abs(variance-.25*R::trigamma((m-1)/2)));
    if (q==0) error=std::max(error,std::abs(log_integral+std::log(2.0)-
        R::lgammafn((m-1)/2)+(m-1)/2*std::log(ss/2)));
    return {mean,variance,cross/mass-first*qfirst/(mass*mass),error,log_integral};
}

struct FixedShapeFit {
    V3 p={NA_REAL,NA_REAL,NA_REAL};
    double residual=NA_REAL;
    bool converged=false;
};

// One conditional scale fit; the starting scale is the stored original fit.
static FixedShapeFit fit_fixed_shape(const std::vector<double>& y,double q,double start_sigma) {
    FixedShapeFit out;
    double s=std::log(start_sigma);
    for (int iteration=0;iteration<150;++iteration) {
        const PositiveProfile at=positive_profile(y,std::exp(s),q,0);
        const double gradient=at.gradient[1]/y.size(),value=at.loss/y.size();
        out.p=at.p;out.residual=std::abs(gradient);
        if (!std::isfinite(value)) return out;
        if (out.residual<profile_tolerance) {out.converged=true;return out;}
        const M3 information=expected_information(at.p);
        const double conditional=information[4]-information[1]*information[1]/information[0];
        const double direction=-gradient/conditional;
        double step=1;
        bool accepted=false;
        for (int backtrack=0;backtrack<35;++backtrack) {
            const double next=s+step*direction;
            const PositiveProfile proposed=positive_profile(y,std::exp(next),q,0);
            if (std::isfinite(proposed.loss) &&
                proposed.loss/y.size()<=value+1e-4*step*gradient*direction) {
                s=next;accepted=true;break;
            }
            step*=0.5;
        }
        if (!accepted) return out;
    }
    return out;
}

// Approximate empirical-Bayes inference for Q, followed by conditional GG fits.
// Each component retains scale/frequency uncertainty and correlated shape
// uncertainty. Mixture tails include uncertainty between sharing components.

static double location_integral_constant(double m,double q) {
    const double c0=-(m-1)/2*std::log(2*3.14159265358979323846)-.5*std::log(m);
    if (std::abs(q)<0.15) {
        const double q2=q*q;
        return c0+(1/m-m)*q2/12+(m-1/std::pow(m,3))*std::pow(q2,3)/360+
            (1/std::pow(m,5)-m)*std::pow(q2,5)/1260;
    }
    const double k=1/(q*q);
    return (m-1)*std::log(std::abs(q))+R::lgammafn(m*k)-m*R::lgammafn(k)-m*k*std::log(m);
}

// Exact integration of location and deterministic quadrature over log scale.
// Every pair with at least four nonconstant positive observations is eligible.
#include "sgcb_pocket_fft.h"
#include "sgcb_sparse_mixture.h"
#include <chrono>
#include <mutex>
#include <map>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <tuple>
#include <Eigen/Core>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

// Dictionary generation needs the density but not its scale or shape derivatives.
static double population_scale_logdensity(const std::vector<double>& centered,double s,double q) {
    const double m=centered.size(),inverse_sigma=std::exp(-s);
    if(q==0) {
        double a2=0;
        for(double value:centered){const double z=value*inverse_sigma;a2+=z*z/m;}
        return -(m-1)*s-m*(a2/2);
    }
    double largest=0,maximum=R_NegInf;
    for(double value:centered) {
        const double z=q*value*inverse_sigma;
        largest=std::max(largest,std::abs(z));maximum=std::max(maximum,z);
    }
    if(largest<.005) {
        double a2=0,a3=0,a4=0,a5=0,a6=0;
        for(double value:centered) {
            const double z=value*inverse_sigma,z2=z*z;
            a2+=z2/m;a3+=z2*z/m;a4+=z2*z2/m;
            a5+=z2*z2*z/m;a6+=z2*z2*z2/m;
        }
        const double c4=a4-3*a2*a2,c5=a5-10*a2*a3,
            c6=a6-15*a4*a2-10*a3*a3+30*a2*a2*a2;
        const double ratio=a2/2+q*(a3/6+q*(c4/24+q*(c5/120+q*c6/720)));
        return -(m-1)*s-m*ratio;
    }
    double sum=0,small=0;
    if(largest<.5)for(double value:centered)small+=std::expm1(q*value*inverse_sigma)/m;
    else for(double value:centered)sum+=std::exp(q*value*inverse_sigma-maximum);
    const double logmean=largest<.5?std::log1p(small):maximum+std::log(sum/m),k=1/(q*q);
    return -(m-1)*s-m*k*logmean;
}

// Likelihood for the joint population of log scale and shape. Location is
// integrated analytically; zero-frequency likelihood remains gene-specific.
// [[Rcpp::export]]
Rcpp::List joint_population_likelihood_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,int threads=1,
    Rcpp::Nullable<Rcpp::NumericMatrix> previous_loglike=R_NilValue) {
    const int ng=x.nrow(),nk=shape.size();
    int totals[2]={0,0};for(int g:group)++totals[g];
    const int maximum_n=std::max(totals[0],totals[1]);
    std::vector<double> distinct(shape.begin(),shape.end());
    std::sort(distinct.begin(),distinct.end());
    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
    const int nq=distinct.size();std::vector<int> shape_index(nk);
    for(int k=0;k<nk;++k)shape_index[k]=std::lower_bound(distinct.begin(),distinct.end(),shape[k])-distinct.begin();
    std::vector<double> constants(size_t(maximum_n+1)*nq);
    for(int m=4;m<=maximum_n;++m)for(int j=0;j<nq;++j)
        constants[size_t(m)*nq+j]=location_integral_constant(m,distinct[j]);
    const int previous_columns=previous_loglike.isNotNull() ? Rcpp::NumericMatrix(previous_loglike).ncol() : 0;
    Rcpp::NumericMatrix loglike(2*ng,previous_columns+nk);Rcpp::IntegerVector count(2*ng);
    Rcpp::LogicalVector usable(2*ng);std::fill(loglike.begin(),loglike.end(),NA_REAL);
    if (previous_loglike.isNotNull()) {
        const Rcpp::NumericMatrix previous(previous_loglike);
        std::copy(previous.begin(),previous.end(),loglike.begin());
    }
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,8)
    #endif
    for (int row=0;row<2*ng;++row) {
        const int gene=row%ng,g=row/ng;std::vector<double> y;
        double center=0,ss=0;
        for (int j=0;j<x.ncol();++j) if (group[j]==g && x(gene,j)>0) y.push_back(std::log(x(gene,j)));
        std::sort(y.begin(),y.end());
        const int m=y.size();count[row]=m;
        for (double value:y) center+=value/m;
        for (double& value:y) {value-=center;ss+=value*value;}
        usable[row]=m>=4 && ss>0;
        if (!usable[row]) continue;
        for (int k=0;k<nk;++k) loglike(row,previous_columns+k)=
            constants[size_t(m)*nq+shape_index[k]]+population_scale_logdensity(y,log_sigma[k],shape[k]);
    }
    return Rcpp::List::create(Rcpp::Named("log_likelihood")=loglike,
        Rcpp::Named("positive")=count,Rcpp::Named("usable")=usable);
}


// Stream complete row blocks through an ordinary R binary connection.
// Only the exact thresholded likelihood columns remain resident.
// [[Rcpp::export(rng = false)]]
Rcpp::List streamed_refinement_dictionary_cpp(Rcpp::NumericMatrix x,
    Rcpp::IntegerVector group,Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,
    Rcpp::NumericMatrix previous_loglike,Rcpp::LogicalVector eligible,SEXP connection,
    int threads=1,double tolerance=1e-6,int block_rows=256) {
    const int ng=x.nrow(),old=previous_loglike.ncol(),nf=shape.size(),nk=old+nf;
    if(previous_loglike.nrow()!=2*ng || eligible.size()!=2*ng || log_sigma.size()!=nf)
        Rcpp::stop("Incompatible refinement likelihood dimensions");
    std::unique_ptr<SparsePopulationBuffer> buffer(new SparsePopulationBuffer);
    for(int row=0;row<2*ng;++row)if(eligible[row])buffer->rows.push_back(row);
    const int n=buffer->rows.size();
    buffer->n=n;buffer->nk=nk;buffer->omitted_bound=tolerance/4;
    buffer->threshold=buffer->omitted_bound/n;
    buffer->row_scale.resize(n);buffer->umbrella.resize(n);
    int totals[2]={0,0};for(int g:group)++totals[g];
    const int maximum_n=std::max(totals[0],totals[1]);
    std::vector<double> distinct(shape.begin(),shape.end());
    std::sort(distinct.begin(),distinct.end());
    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
    const int nq=distinct.size();
    std::vector<int> shape_index(nf);
    for(int k=0;k<nf;++k)
        shape_index[k]=std::lower_bound(distinct.begin(),distinct.end(),shape[k])-distinct.begin();
    std::vector<double> constants(size_t(maximum_n+1)*nq);
    for(int m=4;m<=maximum_n;++m)for(int q=0;q<nq;++q)
        constants[size_t(m)*nq+q]=location_integral_constant(m,distinct[q]);
    std::vector<size_t> count(nk,0);
    Rcpp::Function write=Rcpp::Environment::base_env()["writeBin"];
    Rcpp::Function read=Rcpp::Environment::base_env()["readBin"];
    Rcpp::Function seek=Rcpp::Environment::base_env()["seek"];
    Rcpp::Function flush=Rcpp::Environment::base_env()["flush"];
    Rcpp::Function gc=Rcpp::Environment::base_env()["gc"];
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #endif
    Rcpp::NumericVector block(R_xlen_t(std::min(block_rows,n))*nk);
    for(int begin=0;begin<n;begin+=block_rows) {
        Rcpp::checkUserInterrupt();
        const int nr=std::min(block_rows,n-begin);
        if(block.size()!=R_xlen_t(nr)*nk)block=Rcpp::NumericVector(R_xlen_t(nr)*nk);
        std::vector<std::vector<double>> centered(nr);
        for(int i=0;i<nr;++i) {
            const int row=buffer->rows[begin+i],gene=row%ng,g=row/ng;
            auto& y=centered[i];double center=0;
            for(int j=0;j<x.ncol();++j)
                if(group[j]==g && x(gene,j)>0)y.push_back(std::log(x(gene,j)));
            std::sort(y.begin(),y.end());
            const int m=y.size();
            for(double value:y)center+=value/m;
            for(double& value:y)value-=center;
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<old;++k)
            for(int i=0;i<nr;++i)block[R_xlen_t(k)*nr+i]=previous_loglike(buffer->rows[begin+i],k);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nf;++k)
            for(int i=0;i<nr;++i) {
                const int m=centered[i].size();
                block[R_xlen_t(old+k)*nr+i]=constants[size_t(m)*nq+shape_index[k]]+
                    population_scale_logdensity(centered[i],log_sigma[k],shape[k]);
            }
        for(int i=0;i<nr;++i) {
            double largest=R_NegInf;
            for(int k=0;k<nk;++k)largest=std::max(largest,block[R_xlen_t(k)*nr+i]);
            buffer->row_scale[begin+i]=largest;
        }
        std::vector<double> partial(size_t(threads)*nr,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            int thread=0;
            #ifdef _OPENMP
            thread=omp_get_thread_num();
            #endif
            double* sum=partial.data()+size_t(thread)*nr;
            for(int i=0;i<nr;++i) {
                const double value=std::exp(block[R_xlen_t(k)*nr+i]-buffer->row_scale[begin+i]);
                if(value>buffer->threshold){++count[k];sum[i]+=value/nk;}
            }
        }
        for(int i=0;i<nr;++i) {
            double sum=0;for(int thread=0;thread<threads;++thread)sum+=partial[size_t(thread)*nr+i];
            buffer->umbrella[begin+i]=sum;
        }
        write(block,connection,Rcpp::Named("size")=8,Rcpp::Named("endian")="little");
    }
    flush(connection);
    block=Rcpp::NumericVector(0);
    gc(false);
    const double temporary_bytes=Rcpp::as<double>(seek(connection,
        Rcpp::Named("where")=NA_REAL,Rcpp::Named("rw")="write"));
    seek(connection,Rcpp::Named("where")=0,Rcpp::Named("origin")="start",Rcpp::Named("rw")="read");
    buffer->column_start.assign(size_t(nk)+1,0);
    for(int k=0;k<nk;++k)buffer->column_start[k+1]=buffer->column_start[k]+count[k];
    buffer->column_row.resize(buffer->column_start.back());
    buffer->column_value.resize(buffer->column_start.back());
    std::vector<size_t> at=buffer->column_start;
    for(int begin=0;begin<n;begin+=block_rows) {
        Rcpp::checkUserInterrupt();
        const int nr=std::min(block_rows,n-begin);
        const R_xlen_t length=R_xlen_t(nr)*nk;
        gc(false);
        const Rcpp::NumericVector block=read(connection,Rcpp::Named("what")="double",
            Rcpp::Named("n")=double(length),Rcpp::Named("size")=8,Rcpp::Named("endian")="little");
        if(block.size()!=length)Rcpp::stop("Incomplete temporary refinement likelihood block");
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            size_t entry=at[k];
            for(int i=0;i<nr;++i) {
                const double value=std::exp(block[R_xlen_t(k)*nr+i]-buffer->row_scale[begin+i]);
                if(value>buffer->threshold) {
                    buffer->column_row[entry]=begin+i;
                    buffer->column_value[entry++]=value;
                }
            }
            at[k]=entry;
        }
    }
    const double entries=double(buffer->column_value.size());
    Rcpp::XPtr<SparsePopulationBuffer> pointer(buffer.release(),true);
    return Rcpp::List::create(Rcpp::Named("dictionary")=pointer,
        Rcpp::Named("sparse_entries")=entries,Rcpp::Named("dictionary_rows")=n,
        Rcpp::Named("dictionary_states")=nk,Rcpp::Named("temporary_bytes")=temporary_bytes,
        Rcpp::Named("maximum_block_bytes")=double(std::min(block_rows,n))*nk*8);
}
// End streamed refinement dictionary.

// Complete checkpoint files are immutable; incomplete writes have another suffix.
static std::string population_checkpoint_path(const std::string& directory,
    const std::string& prefix,int index) {
    const std::string number=std::to_string(index);
    return directory+"/"+prefix+std::string(10-number.size(),'0')+number+".rds";
}
static void population_save_checkpoint(const Rcpp::List& state,const std::string& path) {
    const Rcpp::String target(path,CE_UTF8),temporary(path+".writing",CE_UTF8);
    Rcpp::Function save=Rcpp::Environment::base_env()["saveRDS"];
    Rcpp::Function rename=Rcpp::Environment::base_env()["file.rename"];
    save(state,Rcpp::Named("file")=temporary,Rcpp::Named("compress")=false);
    if (!Rcpp::as<bool>(rename(temporary,target))) Rcpp::stop("Checkpoint rename failed");
}

// Exact zero-weight columns do not contribute to the current mixture.
// Remove them from restricted solves; the full dictionary KKT scan still
// checks every atom and can reintroduce any discarded column.
static void prune_population_support(std::vector<int>& support,std::vector<double>& weights) {
    std::vector<int> retained;std::vector<double> kept(1,weights[0]);
    retained.reserve(support.size());kept.reserve(weights.size());
    for (int j=0;j<int(support.size());++j) if (weights[j+1]!=0) {
        retained.push_back(support[j]);kept.push_back(weights[j+1]);
    }
    support.swap(retained);weights.swap(kept);
}

// Form the endpoint from nonnegative mixture contributions. Subtracting the
// leaving column from the current density can create negative endpoint values
// by cancellation, reversing the line-search derivative at the boundary.
static void exchange_population_weight(const Rcpp::NumericMatrix& active,
    std::vector<double>& weights,std::vector<double>& density,int entering,int leaving) {
    const int nr=active.nrow();
    const double leaving_weight=weights[leaving];
    Eigen::VectorXd endpoint_weights=Eigen::Map<Eigen::VectorXd>(weights.data(),weights.size());
    endpoint_weights[entering]+=leaving_weight;endpoint_weights[leaving]=0;
    const Eigen::Map<const Eigen::MatrixXd> likelihood(active.begin(),nr,active.ncol());
    const Eigen::VectorXd endpoint=likelihood*endpoint_weights;
    // Per-unit transfer scores remain representable when the leaving mass
    // is tiny and endpoint-density rounds to zero.
    const Eigen::VectorXd direction=likelihood.col(entering)-likelihood.col(leaving);
    double low=0,high=1,boundary_derivative=0;
    for (int i=0;i<nr;++i) boundary_derivative+=direction[i]/endpoint[i];
    if (boundary_derivative<=0) for (int search=0;search<50;++search) {
        const double fraction=(low+high)/2;double derivative=0;
        for (int i=0;i<nr;++i) derivative+=direction[i]/
            ((1-fraction)*density[i]+fraction*endpoint[i]);
        if (derivative>0) low=fraction;else high=fraction;
    }
    const double fraction=boundary_derivative>0 ? 1 : (low+high)/2;
    const double change=fraction*leaving_weight;
    weights[entering]+=change;weights[leaving]-=change;
    for (int i=0;i<nr;++i)
        density[i]=(1-fraction)*density[i]+fraction*endpoint[i];
}

// r-log(1+r) is second order: its series avoids subtracting its linear term.
static long double population_log1p_remainder(long double r) {
    if (std::abs(r)>=0.125L) return r-std::log1p(r);
    long double term=r*r,value=term/2;
    for (int k=3;k<128;++k) {
        term*=-r;
        const long double next=value+term/k;
        if (next==value) break;
        value=next;
    }
    return value;
}

// Exact difference for -mean(log(L*x))+sum(x); sum(x) is fixed on the simplex.
static long double population_objective_difference(const Eigen::MatrixXd& likelihood,
    const Eigen::VectorXd& density,const Eigen::VectorXd& gradient,
    const Eigen::VectorXd& difference,const Eigen::VectorXd& next_density) {
    const Eigen::VectorXd projected=likelihood*difference;
    long double value=0;
    for (int k=0;k<difference.size();++k)
        value+=static_cast<long double>(gradient[k])*difference[k];
    for (int i=0;i<density.size();++i) {
        const long double r=static_cast<long double>(projected[i])/density[i];
        const long double remainder=std::abs(r)<0.125L ? population_log1p_remainder(r) :
            r-std::log(static_cast<long double>(next_density[i])/density[i]);
        value+=remainder/density.size();
    }
    return value;
}

// Projected Newton-CG for the same nonnegative mixture likelihood.
// Hessian products use L' diag(1/(n*d^2)) L v; no dense Hessian is formed.
// The outer algorithm still tests every candidate against the original KKT bound.
// [[Rcpp::export]]
Rcpp::List population_weights_cpp(Rcpp::NumericMatrix likelihood,Rcpp::NumericVector start,
    double restricted_tolerance=1e-8,int threads=1) {
    const int nr=likelihood.nrow(),nk=likelihood.ncol();
    Rcpp::Environment ns=Rcpp::Environment::namespace_env("mixsqp");
    Rcpp::Function defaults=ns["mixsqp_control_default"];
    const Rcpp::List control=defaults();
    const int maximum_iterations=Rcpp::as<int>(control["maxiter.sqp"]);
    const double sufficient_decrease=Rcpp::as<double>(control["suffdecr.linesearch"]);
    const double reduction=Rcpp::as<double>(control["stepsizereduce"]);
    const double minimum_step=Rcpp::as<double>(control["minstepsize"]);
    const double zero_threshold=1e-12,machine_epsilon=std::numeric_limits<double>::epsilon();
    Eigen::setNbThreads(threads);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #endif
    Eigen::MatrixXd L=Eigen::Map<const Eigen::MatrixXd>(likelihood.begin(),nr,nk);
    const Eigen::VectorXd row_max=L.rowwise().maxCoeff();
    L.array().colwise()/=row_max.array();
    Eigen::VectorXd x=Eigen::Map<const Eigen::VectorXd>(start.begin(),nk);
    x/=x.sum();
    Eigen::VectorXd density=L*x,gradient(nk),diagonal(nk);
    double gap=R_PosInf;
    int iteration=0,cg_iterations=0;
    for (;iteration<maximum_iterations;++iteration) {
        gradient=Eigen::VectorXd::Ones(nk)-L.transpose()*(density.array().inverse()/nr).matrix();
        gap=std::max(0.0,-gradient.minCoeff());
        if (gap<restricted_tolerance) break;
        Rcpp::checkUserInterrupt();
        Eigen::VectorXd free=((x.array()>0)||(gradient.array()<0)).cast<double>().matrix();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int k=0;k<nk;++k)
            diagonal[k]=(L.col(k).array()/density.array()).square().sum()/nr;
        // Damping only regularizes the search direction, not the objective.
        const double ridge=std::sqrt(machine_epsilon)*diagonal.maxCoeff();
        diagonal.array()+=ridge;
        Eigen::VectorXd residual=(-gradient.array()*free.array()).matrix();
        const double target=std::min(.5,std::sqrt(std::max(gap,machine_epsilon)))*residual.norm();
        Eigen::VectorXd z=(residual.array()/diagonal.array()).matrix();
        Eigen::VectorXd direction=z,step_direction=Eigen::VectorXd::Zero(nk);
        double rz=residual.dot(z);
        for (int cg=0;cg<nk && residual.norm()>target;++cg) {
            const Eigen::VectorXd projected=L*direction;
            Eigen::VectorXd product=L.transpose()*
                (projected.array()/density.array().square()/nr).matrix()+ridge*direction;
            product.array()*=free.array();
            const double curvature=direction.dot(product);
            if (!(curvature>0) || !(rz>0)) break;
            double scale=rz/curvature;
            int boundary=-1;
            // Solve the same quadratic direction inside the nonnegative face.
            for(int j=0;j<nk;++j)if(free[j] && direction[j]<0) {
                const double limit=std::max(0.0,x[j]+step_direction[j])/(-direction[j]);
                if(limit<scale){scale=limit;boundary=j;}
            }
            if(boundary>=0)free[boundary]=0;
            step_direction+=scale*direction;
            residual-=scale*product;
            residual.array()*=free.array();
            if(boundary>=0)step_direction[boundary]=-x[boundary];
            ++cg_iterations;
            if (residual.norm()<=target) break;
            z=(residual.array()/diagonal.array()).matrix();
            const double next_rz=residual.dot(z);
            if(boundary>=0)direction=z;
            else direction=z+(next_rz/rz)*direction;
            rz=next_rz;
        }
        for (int k=0;k<nk;++k)
            if (x[k]==0 && step_direction[k]<0) step_direction[k]=0;
        if (!(gradient.dot(step_direction)<0))
            step_direction=(-gradient.array()*free.array()/diagonal.array()).matrix();
        double step=1;
        bool accepted=false;
        for (;step>=minimum_step;step*=reduction) {
            Eigen::VectorXd candidate=(x+step*step_direction).cwiseMax(0);
            for (int k=0;k<nk;++k) if (candidate[k]<zero_threshold) candidate[k]=0;
            const double total=candidate.sum();
            if (!(total>0)) continue;
            candidate/=total;
            const Eigen::VectorXd next_density=L*candidate;
            if (!(next_density.minCoeff()>0)) continue;
            const Eigen::VectorXd difference=candidate-x;
            const double slope=gradient.dot(difference);
            const long double objective_difference=population_objective_difference(
                L,density,gradient,difference,next_density);
            if (slope<0 && objective_difference<=sufficient_decrease*slope) {
                x=std::move(candidate);density=next_density;accepted=true;break;
            }
        }
        if (!accepted) break;
    }
    gradient=Eigen::VectorXd::Ones(nk)-L.transpose()*(density.array().inverse()/nr).matrix();
    gap=std::max(0.0,-gradient.minCoeff());
    Rcpp::NumericVector fitted(nk);
    std::copy(x.data(),x.data()+nk,fitted.begin());
    const std::string status=Rcpp::as<std::string>(
        ns[gap<restricted_tolerance ? "mixsqp.status.converged" : "mixsqp.status.didnotconverge"]);
    return Rcpp::List::create(Rcpp::Named("x")=fitted,Rcpp::Named("status")=status,
        Rcpp::Named("restricted_kkt_gap")=gap,Rcpp::Named("newton_iterations")=iteration,
        Rcpp::Named("cg_iterations")=cg_iterations);
}

// Convex column generation. Each restricted mixture is solved by the existing
// constrained mixture objective with matrix-free restricted Newton steps.
// The first column is the uniform mixture over the entire dictionary, so the
// initial density is positive for every observation. All columns are tested
// against the KKT inequality at each iteration; no genes are selected.
// [[Rcpp::export]]
Rcpp::List fit_population_dictionary_cpp(Rcpp::NumericMatrix loglike,Rcpp::Function solver,
    int threads=1,int maximum_columns=200,double tolerance=1e-6,
    Rcpp::Nullable<Rcpp::NumericVector> start_weights=R_NilValue,bool dense_start=false,
    Rcpp::Nullable<Rcpp::LogicalVector> eligible=R_NilValue,bool reuse_buffer=false,
    std::string checkpoint_directory="",int column_batch=8,bool input_is_log=true) {
    const int nk=loglike.ncol();std::vector<int> rows;
    if (eligible.isNotNull()) {
        const Rcpp::LogicalVector included(eligible);
        for (int i=0;i<loglike.nrow();++i) if (included[i]) rows.push_back(i);
    } else for (int i=0;i<loglike.nrow();++i) rows.push_back(i);
    const int nr=rows.size();
    std::vector<double> likelihood(reuse_buffer ? 0 : size_t(nr)*nk),umbrella(nr,0),density(nr);
    // Refinement can discard its full log-likelihood dictionary after fitting.
    // Reuse that buffer and reconstruct only retained new columns afterwards.
    const int input_rows=loglike.nrow();double* buffer=loglike.begin();
    auto value_at=[&](int i,int k) {
        return reuse_buffer ? buffer[size_t(k)*input_rows+rows[i]] : likelihood[size_t(i)*nk+k];
    };
    std::vector<double> anchor_weight(nk,1.0/nk);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(static)
    #endif
    for (int i=0;i<nr;++i) {
        double largest=0;
        if (input_is_log) {
            largest=R_NegInf;
            for (int k=0;k<nk;++k) largest=std::max(largest,loglike(rows[i],k));
        }
        for (int k=0;k<nk;++k) {
            const double value=input_is_log ? std::exp(loglike(rows[i],k)-largest) : loglike(rows[i],k);
            if (reuse_buffer) buffer[size_t(k)*input_rows+rows[i]]=value;
            else likelihood[size_t(i)*nk+k]=value;
            umbrella[i]+=value/nk;
        }
        density[i]=umbrella[i];
    }
    std::vector<int> support;std::vector<double> weights(1,1);
    if (start_weights.isNotNull() && dense_start) {
        const Rcpp::NumericVector start(start_weights);
        anchor_weight.assign(start.begin(),start.end());
        // Preserve the entire previous density as one convex component.
        // Column generation then improves this same fit without a restart.
        for (int i=0;i<nr;++i) {
            umbrella[i]=0;
            for (int k=0;k<nk;++k) umbrella[i]+=anchor_weight[k]*value_at(i,k);
            density[i]=umbrella[i];
        }
    } else if (start_weights.isNotNull()) {
        const Rcpp::NumericVector start(start_weights);weights[0]=0;
        for (int k=0;k<nk;++k) if (start[k]>0) {support.push_back(k);weights.push_back(start[k]);}
        for (int i=0;i<nr;++i) {
            density[i]=0;
            for (int j=0;j<int(support.size());++j) density[i]+=weights[j+1]*value_at(i,support[j]);
        }
    }
    const int starting_columns=support.size();
    int completed_additions=0;
    if (!checkpoint_directory.empty()) {
        Rcpp::Function list_files=Rcpp::Environment::base_env()["list.files"];
        Rcpp::Function read=Rcpp::Environment::base_env()["readRDS"];
        Rcpp::CharacterVector paths=list_files(Rcpp::Named("path")=Rcpp::String(checkpoint_directory,CE_UTF8),
            Rcpp::Named("pattern")="^weights_[0-9]+[.]rds$",Rcpp::Named("full.names")=true);
        if (paths.size()>0) {
            const Rcpp::List saved=read(Rcpp::CharacterVector::create(paths[paths.size()-1]));
            if (Rcpp::as<int>(saved["atoms"])!=nk || Rcpp::as<int>(saved["rows"])!=nr ||
                Rcpp::as<bool>(saved["dense_start"])!=dense_start) Rcpp::stop("Checkpoint dictionary mismatch");
            support=Rcpp::as<std::vector<int>>(saved["support"]);
            weights=Rcpp::as<std::vector<double>>(saved["restricted_weights"]);
            anchor_weight=Rcpp::as<std::vector<double>>(saved["anchor_weights"]);
            completed_additions=Rcpp::as<int>(saved["completed_additions"]);
            density=Rcpp::as<std::vector<double>>(saved["density"]);
            umbrella=Rcpp::as<std::vector<double>>(saved["umbrella"]);
            Rcpp::Rcout<<"restored dictionary additions="<<completed_additions<<"\n";
        }
    }
    prune_population_support(support,weights);
    double gap=R_PosInf;int iteration=0,added=0,checkpointed=0;std::string solver_status;
    for (;added<maximum_columns && iteration<maximum_columns;++iteration) {
        std::vector<double> partial(threads*nk,0),gradient(nk,0);
        if (reuse_buffer) {
            std::vector<double> factor(nr);
            for (int i=0;i<nr;++i) factor[i]=1/(nr*density[i]);
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int k=0;k<nk;++k) {
                // Match the original static row partitions and addition order,
                // while reading the reused R matrix by contiguous columns.
                for (int t=0;t<threads;++t) {
                    const int begin=t*(nr/threads)+std::min(t,nr%threads);
                    const int end=begin+nr/threads+(t<nr%threads);
                    double score=0;
                    for (int i=begin;i<end;++i) score+=value_at(i,k)*factor[i];
                    partial[t*nk+k]=score;
                }
            }
        } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int i=0;i<nr;++i) {
            int thread=0;
            #ifdef _OPENMP
            thread=omp_get_thread_num();
            #endif
            const double factor=1/(nr*density[i]);
            for (int k=0;k<nk;++k) partial[thread*nk+k]+=likelihood[size_t(i)*nk+k]*factor;
        }
        }
        int best=0;
        for (int k=0;k<nk;++k) {
            for (int t=0;t<threads;++t) gradient[k]+=partial[t*nk+k];
            if (gradient[k]>gradient[best]) best=k;
        }
        gap=gradient[best]-1;
        if (gap<tolerance) break;
        // Add a batch of violated vertices before the restricted solve.
        // Their mean is a feasible mixture direction; the same exact concave
        // line search gives positive starting mass without changing the optimum.
        std::vector<int> entering;
        for (int k=0;k<nk;++k)
            if (gradient[k]>1+tolerance && std::find(support.begin(),support.end(),k)==support.end())
                entering.push_back(k);
        std::stable_sort(entering.begin(),entering.end(),
            [&](int a,int b) {return gradient[a]>gradient[b];});
        entering.resize(std::min({int(entering.size()),column_batch,maximum_columns-added}));
        // A selected vertex can still violate KKT after a restricted solve.
        // Repeat that solve when no additional vertex enters.
        if (!entering.empty()) {
            std::vector<double> candidate(nr,0);
            for (int k:entering) for (int i=0;i<nr;++i)
                candidate[i]+=value_at(i,k)/entering.size();
            double lower_step=0,upper_step=1;
            for (int search=0;search<50;++search) {
                const double step=(lower_step+upper_step)/2;double derivative=0;
                for (int i=0;i<nr;++i) {
                    const double direction=candidate[i]-density[i];
                    derivative+=direction/(density[i]+step*direction);
                }
                if (derivative>0) lower_step=step;else upper_step=step;
            }
            const double step=(lower_step+upper_step)/2;
            for (double& value:weights) value*=1-step;
            for (int k:entering) {support.push_back(k);weights.push_back(step/entering.size());}
            added+=entering.size();
        }
        const int na=support.size();Rcpp::NumericMatrix active(nr,na+1);
        for (int i=0;i<nr;++i) {
            active(i,0)=umbrella[i];
            for (int j=0;j<na;++j) active(i,j+1)=value_at(i,support[j]);
        }
        // Intermediate supports need only reduce the current global violation.
        // The full dictionary KKT test above retains the requested final tolerance.
        const double restricted_tolerance=std::max(tolerance*.1,std::min(.1,gap*.1));
        Rcpp::List result=population_weights_cpp(active,Rcpp::wrap(weights),restricted_tolerance,threads);
        Rcpp::NumericVector fitted=result["x"];weights.assign(fitted.begin(),fitted.end());
        solver_status=Rcpp::as<std::string>(result["status"]);
        for (int i=0;i<nr;++i) {
            density[i]=weights[0]*umbrella[i];
            for (int j=0;j<na;++j) density[i]+=weights[j+1]*value_at(i,support[j]);
        }
        // Check every restricted coordinate, including zero weights. The
        // library can stop with a violated inactive-coordinate inequality.
        // Exact vertex exchanges finish that same convex optimization.
        for (int exchange=0;exchange<2000;++exchange) {
            std::vector<double> score(na+1,0);
            // R matrices are column-major; preserve each column's summation
            // order while reading its observations contiguously.
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int j=0;j<=na;++j) for (int i=0;i<nr;++i)
                score[j]+=active(i,j)/(nr*density[i]);
            int entering=0,leaving=-1;
            for (int j=0;j<=na;++j) {
                if (score[j]>score[entering]) entering=j;
                if (weights[j]>0 && (leaving<0 || score[j]<score[leaving])) leaving=j;
            }
            if (score[entering]-score[leaving]<restricted_tolerance) break;
            exchange_population_weight(active,weights,density,entering,leaving);
        }
        prune_population_support(support,weights);
        if (added/20>checkpointed/20) {
            const int additions=completed_additions+added;
            checkpointed=added;
            if (!checkpoint_directory.empty()) {
                Rcpp::NumericVector state_weights(nk);
                for (int k=0;k<nk;++k) state_weights[k]=weights[0]*anchor_weight[k];
                for (int j=0;j<int(support.size());++j) state_weights[support[j]]+=weights[j+1];
                population_save_checkpoint(Rcpp::List::create(Rcpp::Named("weights")=state_weights,
                    Rcpp::Named("support")=support,Rcpp::Named("restricted_weights")=weights,
                    Rcpp::Named("anchor_weights")=anchor_weight,Rcpp::Named("dense_start")=dense_start,
                    Rcpp::Named("density")=density,Rcpp::Named("umbrella")=umbrella,
                    Rcpp::Named("atoms")=nk,Rcpp::Named("rows")=nr,
                    Rcpp::Named("completed_additions")=additions),
                    population_checkpoint_path(checkpoint_directory,"weights_",additions));
            }
            Rcpp::Rcout<<"dictionary columns="<<additions<<", KKT gap="<<gap<<"\n";
        }
    }
    Rcpp::NumericVector fitted(nk);
    for (int k=0;k<nk;++k) fitted[k]=weights[0]*anchor_weight[k];
    for (int j=0;j<int(support.size());++j) fitted[support[j]]+=weights[j+1];
    return Rcpp::List::create(Rcpp::Named("weights")=fitted,Rcpp::Named("kkt_gap")=gap,
        Rcpp::Named("columns")=starting_columns+completed_additions+added,
        Rcpp::Named("restricted_solves")=iteration,
        Rcpp::Named("retained_columns")=support.size(),Rcpp::Named("converged")=(gap<tolerance),
        Rcpp::Named("uniform_component_weight")=(dense_start ? NA_REAL : weights[0]),
        Rcpp::Named("initial_component_weight")=weights[0],Rcpp::Named("last_solver_status")=solver_status);
}

// A blocked matrix product reuses cache lines for the complete paired-state
// gradient, without forming the feature-by-Cartesian-state likelihood matrix.
static Rcpp::NumericMatrix population_crossprod(const Rcpp::NumericMatrix& left,
    const Rcpp::NumericMatrix& right,int threads) {
    const Eigen::Map<const Eigen::MatrixXd> lhs(left.begin(),left.nrow(),left.ncol());
    const Eigen::Map<const Eigen::MatrixXd> rhs(right.begin(),right.nrow(),right.ncol());
    Rcpp::NumericMatrix result(left.ncol(),right.ncol());
    Eigen::Map<Eigen::MatrixXd> output(result.begin(),result.nrow(),result.ncol());
    Eigen::setNbThreads(threads);
    output.noalias()=lhs.transpose()*rhs;
    return result;
}

// The paired dictionary factorizes as L0(i,k)*L1(i,l). Its complete gradient
// is L0' diag(1/density) L1 / n, so no feature-by-Cartesian-state matrix is
// needed. The optimization and retained states are otherwise unchanged.
// Block maxima bound every pair score; the current eighth score can tighten
// the search threshold without removing any candidate from the KKT problem.
struct PairedPricing {
    int nr,ns,block_size,blocks;
    Eigen::MatrixXd max_left,max_right;
    double inflation;
    PairedPricing(const Rcpp::NumericMatrix& left,const Rcpp::NumericMatrix& right,
        int size,int threads):nr(left.nrow()),ns(left.ncol()),block_size(size),
        blocks((ns+size-1)/size),max_left(nr,blocks),max_right(nr,blocks) {
        const Eigen::Map<const Eigen::MatrixXd> l(left.begin(),nr,ns),r(right.begin(),nr,ns);
        #ifdef _OPENMP
        omp_set_num_threads(threads);
        #pragma omp parallel for schedule(static)
        #endif
        for (int b=0;b<blocks;++b) {
            const int start=b*block_size,count=std::min(block_size,ns-start);
            max_left.col(b)=l.middleCols(start,count).rowwise().maxCoeff();
            max_right.col(b)=r.middleCols(start,count).rowwise().maxCoeff();
        }
        const double eps=std::numeric_limits<double>::epsilon();
        const double gamma=(2.0*nr+4)*eps/(1-(2.0*nr+4)*eps);
        inflation=1/(1-gamma);
    }
    std::vector<R_xlen_t> select(const Rcpp::NumericMatrix& weighted,
        const Rcpp::NumericMatrix& right,const std::vector<double>& density,
        const std::vector<R_xlen_t>& support,int count,double tolerance,int threads,double& gap) {
        Eigen::setNbThreads(threads);
        const Eigen::Map<const Eigen::MatrixXd> w(weighted.begin(),nr,ns),r(right.begin(),nr,ns);
        const Eigen::Map<const Eigen::VectorXd> d(density.data(),nr);
        const Eigen::MatrixXd wm=max_left.array().colwise()/d.array();
        const Eigen::MatrixXd upper=(wm.transpose()*max_right)*(inflation/nr);
        std::vector<int> order(blocks*blocks);std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
            return upper(a/blocks,a%blocks)>upper(b/blocks,b%blocks);
        });
        const std::unordered_set<R_xlen_t> occupied(support.begin(),support.end());
        std::vector<R_xlen_t> chosen;std::vector<double> scores;
        double best=0,remaining=0;
        const double cut=1+tolerance;
        for (int block:order) {
            const double threshold=chosen.size()<size_t(count) ? cut : std::max(cut,scores.back());
            const int a=block/blocks,b=block%blocks;
            if (upper(a,b)<threshold) {remaining=upper(a,b);break;}
            const int i=a*block_size,j=b*block_size;
            const int ni=std::min(block_size,ns-i),nj=std::min(block_size,ns-j);
            const Eigen::MatrixXd tile=w.middleCols(i,ni).transpose()*r.middleCols(j,nj);
            for (int u=0;u<ni;++u) for (int v=0;v<nj;++v) {
                const double score=tile(u,v)/nr;best=std::max(best,score);
                const R_xlen_t index=R_xlen_t(i+u)*ns+j+v;
                if (!(score>cut) || occupied.count(index)) continue;
                int position=0;
                while (position<int(chosen.size()) && (scores[position]>score ||
                    (scores[position]==score && chosen[position]<index))) ++position;
                if (position>=count) continue;
                chosen.insert(chosen.begin()+position,index);scores.insert(scores.begin()+position,score);
                if (chosen.size()>size_t(count)) {chosen.pop_back();scores.pop_back();}
            }
        }
        // When fewer than count violations exist, unvisited upper bounds
        // still certify the complete dictionary, including occupied pairs.
        gap=std::max(best,remaining)-1;
        return chosen;
    }
};
static void prune_priced_support(std::vector<R_xlen_t>& support,std::vector<double>& weights) {
    std::vector<R_xlen_t> retained;std::vector<double> kept(1,weights[0]);
    for (int j=0;j<int(support.size());++j) if (weights[j+1]!=0) {
        retained.push_back(support[j]);kept.push_back(weights[j+1]);
    }
    support.swap(retained);weights.swap(kept);
}

// [[Rcpp::export]]
Rcpp::List fit_paired_population_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::LogicalVector usable_pair,Rcpp::NumericVector start_weights,Rcpp::Function solver,
    int threads=1,int maximum_columns=200,double tolerance=1e-6,int column_batch=8,
    std::string checkpoint_directory="",Rcpp::Nullable<Rcpp::NumericVector> marginal_start=R_NilValue,
    bool sparse_output=false,int block_size=32,bool symmetric=false,bool initialize_support=false) {
    const int ng=loglike.nrow()/2,ns=loglike.ncol();
    const R_xlen_t nk=R_xlen_t(ns)*ns;
    const bool factorized=marginal_start.isNotNull();
    if(initialize_support && (factorized || start_weights.size()!=nk))
        Rcpp::stop("Support initialization requires explicit joint weights");
    const Rcpp::NumericVector marginal=factorized ? Rcpp::NumericVector(marginal_start) : Rcpp::NumericVector();
    std::vector<int> rows;
    for (int i=0;i<ng;++i) if (usable_pair[i]) rows.push_back(i);
    const int nr=rows.size();
    Rcpp::NumericMatrix left(nr,ns),right(nr,ns),weighted(nr,ns),reciprocal(nr,1);
    Rcpp::Environment base=Rcpp::Environment::base_env();
    Rcpp::Function multiply=base["%*%"],crossprod=base["crossprod"];
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(static)
    #endif
    for (int i=0;i<nr;++i) {
        double peak0=R_NegInf,peak1=R_NegInf;
        for (int k=0;k<ns;++k) {
            peak0=std::max(peak0,loglike(rows[i],k));
            peak1=std::max(peak1,loglike(rows[i]+ng,k));
        }
        for (int k=0;k<ns;++k) {
            left(i,k)=std::exp(loglike(rows[i],k)-peak0);
            right(i,k)=std::exp(loglike(rows[i]+ng,k)-peak1);
        }
    }
    auto pair_value=[&](int i,R_xlen_t k) {
        const int a=k/ns,b=k%ns;
        const double value=left(i,a)*right(i,b);
        return symmetric ? .5*(value+left(i,b)*right(i,a)) : value;
    };
    if(symmetric && !factorized) {
        start_weights=Rcpp::clone(start_weights);
        for(int k=0;k<ns;++k)for(int l=k+1;l<ns;++l) {
            const double mass=.5*(start_weights[R_xlen_t(k)*ns+l]+start_weights[R_xlen_t(l)*ns+k]);
            start_weights[R_xlen_t(k)*ns+l]=mass;
            start_weights[R_xlen_t(l)*ns+k]=mass;
        }
    }
    std::vector<double> umbrella(nr,0),density(nr),entering_column(nr);
    std::vector<R_xlen_t> initial_support;
    double initial_total=0;R_xlen_t anchor=-1;
    if(initialize_support) {
        for(R_xlen_t k=0;k<nk;++k)if(start_weights[k]>0) {
            initial_support.push_back(k);initial_total+=start_weights[k];
            if(anchor<0 || start_weights[k]>start_weights[anchor])anchor=k;
        }
        if(!(initial_total>0))Rcpp::stop("Empty joint starting support");
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int i=0;i<nr;++i)for(R_xlen_t k:initial_support)
            umbrella[i]+=start_weights[k]/initial_total*pair_value(i,k);
    } else if (factorized) {
        const Eigen::Map<const Eigen::MatrixXd> l(left.begin(),nr,ns),r(right.begin(),nr,ns);
        const Eigen::Map<const Eigen::VectorXd> p(marginal.begin(),ns);
        const Eigen::VectorXd product=((l*p).array()*(r*p).array()).matrix();
        std::copy(product.data(),product.data()+nr,umbrella.begin());
    } else {
        Rcpp::NumericMatrix mixing(ns,ns);
        for (int k=0;k<ns;++k) for (int l=0;l<ns;++l)
            mixing(k,l)=start_weights[R_xlen_t(k)*ns+l];
        const Rcpp::NumericMatrix projected=Rcpp::as<Rcpp::NumericMatrix>(multiply(left,mixing));
        for (int i=0;i<nr;++i) for (int l=0;l<ns;++l)
            umbrella[i]+=projected(i,l)*right(i,l);
    }
    density=umbrella;
    std::vector<R_xlen_t> support;std::vector<double> weights(1,1);
    if(initialize_support) {
        // One occupied pair anchors column zero; every other initial pair is
        // free to reweight immediately. No frozen compound start remains.
        weights[0]=start_weights[anchor]/initial_total;
        for(R_xlen_t k:initial_support)if(k!=anchor) {
            support.push_back(k);weights.push_back(start_weights[k]/initial_total);
        }
        for(int i=0;i<nr;++i)umbrella[i]=pair_value(i,anchor);
    }
    int completed_additions=0;
    if (!checkpoint_directory.empty()) {
        Rcpp::Function list_files=base["list.files"],read=base["readRDS"];
        Rcpp::CharacterVector paths=list_files(Rcpp::Named("path")=Rcpp::String(checkpoint_directory,CE_UTF8),
            Rcpp::Named("pattern")="^weights_[0-9]+[.]rds$",Rcpp::Named("full.names")=true);
        if (paths.size()>0) {
            const Rcpp::List saved=read(Rcpp::CharacterVector::create(paths[paths.size()-1]));
            const Rcpp::NumericVector previous_start=saved["initial_weights"];
            const bool previous_symmetric=saved.containsElementNamed("symmetric") &&
                Rcpp::as<bool>(saved["symmetric"]);
            if(previous_symmetric!=symmetric)Rcpp::stop("Checkpoint paired symmetry mismatch");
            const bool previous_support=saved.containsElementNamed("initialize_support") &&
                Rcpp::as<bool>(saved["initialize_support"]);
            if(previous_support!=initialize_support)Rcpp::stop("Checkpoint paired initialization mismatch");
            if (Rcpp::as<int>(saved["atoms"])!=ns || Rcpp::as<int>(saved["rows"])!=nr ||
                previous_start.size()!=start_weights.size() ||
                !std::equal(previous_start.begin(),previous_start.end(),start_weights.begin()))
                Rcpp::stop("Checkpoint paired dictionary mismatch");
            if (factorized) {
                const Rcpp::NumericVector previous_marginal=saved["initial_marginal"];
                if (previous_marginal.size()!=marginal.size() ||
                    !std::equal(previous_marginal.begin(),previous_marginal.end(),marginal.begin()))
                    Rcpp::stop("Checkpoint marginal start mismatch");
            }
            support=Rcpp::as<std::vector<R_xlen_t>>(saved["support"]);
            weights=Rcpp::as<std::vector<double>>(saved["restricted_weights"]);
            density=Rcpp::as<std::vector<double>>(saved["density"]);
            umbrella=Rcpp::as<std::vector<double>>(saved["umbrella"]);
            completed_additions=Rcpp::as<int>(saved["completed_additions"]);
            Rcpp::Rcout<<"restored paired additions="<<completed_additions<<"\n";
        }
    }
    std::unique_ptr<PairedPricing> pricing;
    if (sparse_output && !symmetric) pricing=std::make_unique<PairedPricing>(left,right,block_size,threads);
    double gap=R_PosInf;std::string solver_status;int iteration=0,added=0,checkpointed=0;
    for (;added<maximum_columns && iteration<maximum_columns;++iteration) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int k=0;k<ns;++k) for (int i=0;i<nr;++i) weighted(i,k)=left(i,k)/density[i];
        std::vector<R_xlen_t> entering;
        std::vector<R_xlen_t> seeded_support;
        if(anchor>=0){seeded_support=support;seeded_support.push_back(anchor);}
        const auto& occupied=anchor>=0?seeded_support:support;
        if (sparse_output && !symmetric) {
            entering=pricing->select(weighted,right,density,occupied,
                std::min(column_batch,maximum_columns-added),tolerance,threads,gap);
        } else {
            Rcpp::NumericMatrix gradient=population_crossprod(weighted,right,threads);
            if(symmetric)for(int k=0;k<ns;++k)for(int l=k+1;l<ns;++l) {
                const double score=.5*(gradient(k,l)+gradient(l,k));
                gradient(k,l)=score;gradient(l,k)=score;
            }
            R_xlen_t best=0;double best_score=gradient(0,0)/nr;
            for (int k=0;k<ns;++k) for (int l=0;l<ns;++l) {
                const double score=gradient(k,l)/nr;
                if (score>best_score) {best=R_xlen_t(k)*ns+l;best_score=score;}
            }
            gap=best_score-1;
            for (R_xlen_t k=0;k<nk;++k)
                if ((!symmetric || k/ns<=k%ns) && gradient(k/ns,k%ns)/nr>1+tolerance &&
                    std::find(occupied.begin(),occupied.end(),k)==occupied.end()) entering.push_back(k);
            std::stable_sort(entering.begin(),entering.end(),
                [&](R_xlen_t a,R_xlen_t b) {return gradient(a/ns,a%ns)>gradient(b/ns,b%ns);});
            entering.resize(std::min({int(entering.size()),column_batch,maximum_columns-added}));
        }
        if (gap<tolerance) break;
        // Continue the restricted solve when a selected pair still violates KKT.
        if (!entering.empty()) {
            std::fill(entering_column.begin(),entering_column.end(),0);
            for (R_xlen_t k:entering) for (int i=0;i<nr;++i)
                entering_column[i]+=pair_value(i,k)/entering.size();
            double lower_step=0,upper_step=1;
            for (int search=0;search<50;++search) {
                const double step=(lower_step+upper_step)/2;double derivative=0;
                for (int i=0;i<nr;++i) {
                    const double direction=entering_column[i]-density[i];
                    derivative+=direction/(density[i]+step*direction);
                }
                if (derivative>0) lower_step=step;else upper_step=step;
            }
            const double step=(lower_step+upper_step)/2;
            for (double& value:weights) value*=1-step;
            for (R_xlen_t k:entering) {support.push_back(k);weights.push_back(step/entering.size());}
            added+=entering.size();
        }
        const int na=support.size();Rcpp::NumericMatrix active(nr,na+1);
        for (int i=0;i<nr;++i) {
            active(i,0)=umbrella[i];
            for (int j=0;j<na;++j)
                active(i,j+1)=pair_value(i,support[j]);
        }
        // Intermediate supports need only reduce the current global violation.
        // The full dictionary KKT test above retains the requested final tolerance.
        const double restricted_tolerance=std::max(tolerance*.1,std::min(.1,gap*.1));
        Rcpp::List result=population_weights_cpp(active,Rcpp::wrap(weights),restricted_tolerance,threads);
        Rcpp::NumericVector fitted=result["x"];weights.assign(fitted.begin(),fitted.end());
        solver_status=Rcpp::as<std::string>(result["status"]);
        for (int i=0;i<nr;++i) {
            density[i]=0;
            for (int j=0;j<=na;++j) density[i]+=weights[j]*active(i,j);
        }
        for (int exchange=0;exchange<2000;++exchange) {
            for (int i=0;i<nr;++i) reciprocal(i,0)=1/(nr*density[i]);
            Rcpp::NumericMatrix score=population_crossprod(active,reciprocal,threads);
            int entering=0,leaving=-1;
            for (int j=0;j<=na;++j) {
                if (score(j,0)>score(entering,0)) entering=j;
                if (weights[j]>0 && (leaving<0 || score(j,0)<score(leaving,0))) leaving=j;
            }
            if (score(entering,0)-1<restricted_tolerance) break;
            exchange_population_weight(active,weights,density,entering,leaving);
        }
        prune_priced_support(support,weights);
        if (added/20>checkpointed/20) {
            const int additions=completed_additions+added;
            checkpointed=added;
            if (!checkpoint_directory.empty()) {
                population_save_checkpoint(Rcpp::List::create(Rcpp::Named("support")=support,
                    Rcpp::Named("restricted_weights")=weights,Rcpp::Named("density")=density,
                    Rcpp::Named("umbrella")=umbrella,Rcpp::Named("initial_weights")=start_weights,
                    Rcpp::Named("initial_marginal")=marginal,
                    Rcpp::Named("initialize_support")=initialize_support,
                    Rcpp::Named("atoms")=ns,Rcpp::Named("rows")=nr,
                    Rcpp::Named("completed_additions")=additions,Rcpp::Named("symmetric")=symmetric),
                    population_checkpoint_path(checkpoint_directory,"weights_",additions));
            }
            Rcpp::Rcout<<"paired columns="<<additions<<", KKT gap="<<gap<<"\n";
        }
    }
    Rcpp::NumericVector fitted;
    Rcpp::IntegerVector pair_left,pair_right;
    if (sparse_output) {
        std::map<R_xlen_t,double> combined;
        if (weights[0]>0) {
            if(anchor>=0) {
                const R_xlen_t reverse=(anchor%ns)*ns+anchor/ns;
                if(symmetric && anchor!=reverse) {
                    combined[anchor]=.5*weights[0];combined[reverse]=.5*weights[0];
                } else combined[anchor]=weights[0];
            }
            else if (factorized) {
                std::vector<int> occupied;
                for (int k=0;k<ns;++k) if (marginal[k]>0) occupied.push_back(k);
                for (int k:occupied) for (int l:occupied)
                    combined[R_xlen_t(k)*ns+l]=weights[0]*(marginal[k]*marginal[l]);
            } else for (R_xlen_t k=0;k<nk;++k) if (start_weights[k]>0)
                combined[k]=weights[0]*start_weights[k];
        }
        for (int j=0;j<int(support.size());++j) {
            const R_xlen_t k=support[j],reverse=(k%ns)*ns+k/ns;
            if(symmetric && k!=reverse) {
                combined[k]+=.5*weights[j+1];combined[reverse]+=.5*weights[j+1];
            } else combined[k]+=weights[j+1];
        }
        fitted=Rcpp::NumericVector(combined.size());
        pair_left=Rcpp::IntegerVector(combined.size());pair_right=Rcpp::IntegerVector(combined.size());
        int position=0;
        for (const auto& item:combined) {
            pair_left[position]=item.first/ns+1;pair_right[position]=item.first%ns+1;
            fitted[position++]=item.second;
        }
    } else {
        fitted=Rcpp::NumericVector(nk);
        if(anchor>=0) {
            const R_xlen_t reverse=(anchor%ns)*ns+anchor/ns;
            if(symmetric && anchor!=reverse) {
                fitted[anchor]=.5*weights[0];fitted[reverse]=.5*weights[0];
            } else fitted[anchor]=weights[0];
        }
        else for (R_xlen_t k=0;k<nk;++k)
            fitted[k]=weights[0]*(factorized ? marginal[k/ns]*marginal[k%ns] : start_weights[k]);
        for (int j=0;j<int(support.size());++j) {
            const R_xlen_t k=support[j],reverse=(k%ns)*ns+k/ns;
            if(symmetric && k!=reverse) {
                fitted[k]+=.5*weights[j+1];fitted[reverse]+=.5*weights[j+1];
            } else fitted[k]+=weights[j+1];
        }
    }
    return Rcpp::List::create(Rcpp::Named("weights")=fitted,Rcpp::Named("left")=pair_left,Rcpp::Named("right")=pair_right,
        Rcpp::Named("symmetric")=symmetric,Rcpp::Named("kkt_gap")=gap,
        Rcpp::Named("columns")=completed_additions+added,Rcpp::Named("converged")=(gap<tolerance),
        Rcpp::Named("restricted_solves")=iteration,Rcpp::Named("retained_columns")=support.size(),
        Rcpp::Named("initial_support_columns")=initial_support.size(),
        Rcpp::Named("initial_component_is_pair")=initialize_support,
        Rcpp::Named("uniform_component_weight")=NA_REAL,
        Rcpp::Named("initial_component_weight")=weights[0],Rcpp::Named("last_solver_status")=solver_status);
}

static double grid_cdf(const std::vector<double>& cdf,int index) {
    if (index<0) return 0;
    if (index>=int(cdf.size())) return cdf.back();
    return cdf[index];
}

// Full nonlinear CV2 posterior on a finite dictionary, integrating positivity
// by its beta CDF. Undefined second-moment mass is retained. Binned convolution
// gives lower/upper CDF bounds, not a normal approximation or random draws.
// Empirical-Bayes mass of the structural pi=1 component.
// The continuous component retains its complete beta-binomial marginal.
// [[Rcpp::export(rng = false)]]
Rcpp::List population_frequency_structure_cpp(Rcpp::IntegerVector positive,
    Rcpp::IntegerVector totals) {
    const int ng=positive.size()/2,n=positive.size();
    int complete[2]={0,0};
    double beta[2],score0=0;
    for (int g=0;g<2;++g) {
        beta[g]=std::exp(R::lbeta(totals[g]+.5,.5)-R::lbeta(.5,.5));
        for (int i=0;i<ng;++i) complete[g]+=positive[g*ng+i]==totals[g];
        score0+=complete[g]*(1-beta[g])/beta[g];
    }
    const int partial=n-complete[0]-complete[1];
    score0-=partial;
    double mass=0;
    if (partial==0) mass=1;
    else if (score0>0) {
        double lower=0,upper=1;
        for (int iteration=0;iteration<80;++iteration) {
            const double candidate=(lower+upper)/2;
            double score=-partial/(1-candidate);
            for (int g=0;g<2;++g)
                score+=complete[g]*(1-beta[g])/(candidate+(1-candidate)*beta[g]);
            if (score>0) lower=candidate;else upper=candidate;
        }
        mass=(lower+upper)/2;
    }
    return Rcpp::List::create(Rcpp::Named("mass")=mass,
        Rcpp::Named("complete_groups")=Rcpp::IntegerVector::create(complete[0],complete[1]),
        Rcpp::Named("partial_groups")=partial,Rcpp::Named("score_at_zero")=score0);
}


// [[Rcpp::export]]
Rcpp::List joint_population_cv_cpp(Rcpp::NumericMatrix loglike,Rcpp::IntegerVector positive,
    Rcpp::LogicalVector usable,Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,
    Rcpp::NumericVector weight,Rcpp::IntegerVector totals,Rcpp::NumericVector truth,
    double grid_step=0.025,double probability_tolerance=1e-12,double alpha=.05,int threads=1,
    bool shared_shape=false,Rcpp::Nullable<Rcpp::NumericVector> moment_increment=R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> sharing_key=R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> pair_left=R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> pair_right=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> pair_weights=R_NilValue,double structural_mass=.5,
    bool shared_shape_only=false) {
    const int ng=loglike.nrow()/2,nk=weight.size(),maximum_n=std::max(totals[0],totals[1]);
    const bool joint_pair=pair_weights.isNotNull();
    std::vector<int> pair0,pair1;std::vector<double> joint_weight;
    if (joint_pair) {
        pair0=Rcpp::as<std::vector<int>>(pair_left);
        pair1=Rcpp::as<std::vector<int>>(pair_right);
        joint_weight=Rcpp::as<std::vector<double>>(pair_weights);
    }
    std::vector<double> key(shape.begin(),shape.end()),increment(nk,0);
    if (sharing_key.isNotNull()) {
        const Rcpp::IntegerVector supplied(sharing_key);
        for (int k=0;k<nk;++k) key[k]=supplied[k];
    }
    if (moment_increment.isNotNull()) increment=Rcpp::as<std::vector<double>>(moment_increment);
    std::vector<double> distinct_shape=key;
    std::sort(distinct_shape.begin(),distinct_shape.end());
    distinct_shape.erase(std::unique(distinct_shape.begin(),distinct_shape.end()),distinct_shape.end());
    const int nq=distinct_shape.size();
    std::vector<int> shape_index(nk);
    std::vector<double> prior_shape(nq,0),posterior(joint_pair ? 0 : 2*ng*nk,0),
        posterior_shape(joint_pair ? 0 : 2*ng*nq,0),bayes_factor(ng,0);
    for (int k=0;k<nk;++k) {
        shape_index[k]=std::lower_bound(distinct_shape.begin(),distinct_shape.end(),key[k])-distinct_shape.begin();
        prior_shape[shape_index[k]]+=weight[k];
    }
    for (int row=0;row<2*ng;++row) {
        if (joint_pair || !usable[row]) continue;
        double maximum=R_NegInf,sum=0;
        for (int k=0;k<nk;++k) maximum=std::max(maximum,std::log(weight[k])+loglike(row,k));
        for (int k=0;k<nk;++k) {
            posterior[row*nk+k]=std::exp(std::log(weight[k])+loglike(row,k)-maximum);
            sum+=posterior[row*nk+k];
        }
        for (int k=0;k<nk;++k) {
            posterior[row*nk+k]/=sum;
            posterior_shape[row*nq+shape_index[k]]+=posterior[row*nk+k];
        }
    }
    // Both branches have the same marginal population prior. The shared
    // branch couples Q only; scale remains conditional on that Q. Its Bayes
    // factor is sum_q P(q|y0) P(q|y1) / P(q), with locations separate.
    for (int i=0;i<ng;++i) if (!joint_pair && usable[i] && usable[i+ng])
        for (int q=0;q<nq;++q)
            bayes_factor[i]+=posterior_shape[i*nq+q]*posterior_shape[(i+ng)*nq+q]/prior_shape[q];
    double rho=0;
    if (shared_shape && !joint_pair) {
        double lower_rho=0,upper_rho=1;
        for (int iteration=0;iteration<80;++iteration) {
            const double candidate=(lower_rho+upper_rho)/2;double derivative=0;
            for (int i=0;i<ng;++i) if (usable[i] && usable[i+ng])
                derivative+=(bayes_factor[i]-1)/(1-candidate+candidate*bayes_factor[i]);
            if (derivative>0) lower_rho=candidate;else upper_rho=candidate;
        }
        rho=(lower_rho+upper_rho)/2;
    }
    if (shared_shape_only)
        return Rcpp::List::create(Rcpp::Named("shared_shape_probability")=rho);
    std::vector<double> K(nk);double lower=R_PosInf,maximum_K=R_NegInf;
    for (int k=0;k<nk;++k) {
        const V3 p={0,log_sigma[k],shape[k]};K[k]=log_moment(2,p)-2*log_moment(1,p)+increment[k];
        if (std::isfinite(K[k]) && K[k]>0) {
            lower=std::min(lower,K[k]+std::log(-std::expm1(-K[k])));
            maximum_K=std::max(maximum_K,K[k]);
        }
    }
    const double minimum_pi=R::qbeta(probability_tolerance,4.5,maximum_n-3.5,true,false);
    const double maximum_R=maximum_K-std::log(minimum_pi);
    lower=std::floor(lower/grid_step)*grid_step-2*grid_step;
    const double upper=maximum_R+std::log(-std::expm1(-maximum_R))+2*grid_step;
    const int nd=int(std::ceil((upper-lower)/grid_step))+1;
    // Observed positive counts identify identical conditional beta laws.
    // Reuse their CDF grids and moments without grouping gene likelihoods.
    std::vector<int> count_slot(2*(maximum_n+1),-1),count_group,count_value;
    for (int row=0;row<2*ng;++row) if (usable[row])
        count_slot[(row/ng)*(maximum_n+1)+positive[row]]=0;
    for (int g=0;g<2;++g) for (int m=4;m<=totals[g];++m)
        if (count_slot[g*(maximum_n+1)+m]==0) {
            count_slot[g*(maximum_n+1)+m]=count_group.size();
            count_group.push_back(g);count_value.push_back(m);
        }
    const int cache_rows=count_group.size()*nk;
    std::vector<double> cache(size_t(cache_rows)*nd,0),cache_mass(cache_rows,0),
        cache_mean(cache_rows,0),cache_variance(cache_rows,0);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,8)
    #endif
    for (int row=0;row<cache_rows;++row) {
        const int k=row%nk,slot=row/nk,m=count_value[slot],g=count_group[slot];
        if (!std::isfinite(K[k]) || !(K[k]>0)) continue;
        const size_t offset=size_t(row)*nd;
        for (int j=0;j<nd;++j) {
            const double logthreshold=K[k]-softplus(lower+j*grid_step);
            cache[offset+j]=logthreshold>=0 ? 0 :
                R::pbeta(std::exp(logthreshold),m+.5,totals[g]-m+.5,false,false);
        }
        double first=0,second=0;
        for (int j=0;j<nd-1;++j) {
            const double bin=std::max(0.0,cache[offset+j+1]-cache[offset+j]),
                value=lower+(j+.5)*grid_step;
            first+=bin*value;second+=bin*value*value;
        }
        const double mass=cache[offset+nd-1]-cache[offset];
        cache_mass[row]=mass;
        if (mass>0) {
            cache_mean[row]=first/mass;
            cache_variance[row]=std::max(0.0,second/mass-cache_mean[row]*cache_mean[row]);
        }
        // The same spike-and-beta positivity prior used by the mean channel.
        // Counts contribute a common factor to every nuisance state, so the
        // already fitted positive-part population prior remains unchanged.
        if (structural_mass>0 && m==totals[g]) {
            const double lp=std::log(structural_mass),
                lc=std::log1p(-structural_mass)+R::lbeta(m+.5,totals[g]-m+.5)-R::lbeta(.5,.5),
                peak=std::max(lp,lc),
                spike=std::exp(lp-peak)/(std::exp(lp-peak)+std::exp(lc-peak)),
                point=K[k]+std::log(-std::expm1(-K[k])),
                mixed_mass=spike+(1-spike)*cache_mass[row],
                probability=spike/mixed_mass,
                old_mean=cache_mean[row],old_variance=cache_variance[row];
            cache_mean[row]=probability*point+(1-probability)*old_mean;
            cache_variance[row]=(1-probability)*old_variance+
                probability*(1-probability)*(point-old_mean)*(point-old_mean);
            cache_mass[row]=mixed_mass;
            for (int j=0;j<nd;++j)
                cache[offset+j]=(1-spike)*cache[offset+j]+
                    spike*(lower+j*grid_step>=point);
        }
    }
    std::vector<int> count_pair_index(ng,-1);
    std::vector<std::pair<int,int>> count_pairs;
    std::vector<double> pair_cdf_lower,pair_cdf_upper;
    if (joint_pair) {
        std::map<std::pair<int,int>,int> pair_lookup;
        for (int i=0;i<ng;++i) if (usable[i] && usable[i+ng]) {
            const std::pair<int,int> key={positive[i],positive[i+ng]};
            const auto inserted=pair_lookup.emplace(key,count_pairs.size());
            if (inserted.second) count_pairs.push_back(key);
            count_pair_index[i]=inserted.first->second;
        }
        const int components=joint_weight.size();
        const size_t combinations=count_pairs.size()*size_t(components);
        pair_cdf_lower.resize(combinations);pair_cdf_upper.resize(combinations);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,8)
        #endif
        for (size_t index=0;index<combinations;++index) {
            const int component=index%components,count_pair=index/components;
            const auto counts=count_pairs[count_pair];
            const size_t offset0=size_t(count_slot[counts.first]*nk+pair0[component])*nd;
            const size_t offset1=size_t(count_slot[maximum_n+1+counts.second]*nk+pair1[component])*nd;
            double lo=0,hi=0;
            for (int j=0;j<nd-1;++j) {
                const double bin=std::max(0.0,cache[offset0+j+1]-cache[offset0+j]);
                lo+=bin*cache[offset1+j];hi+=bin*cache[offset1+j+1];
            }
            pair_cdf_lower[index]=lo;pair_cdf_upper[index]=hi;
        }
    }
    const std::vector<std::string> names={"effect","se","cdf_zero_lower","cdf_zero_upper",
        "cdf_truth_lower","cdf_truth_upper","available","undefined_mass","quadrature_tail_mass",
        "p","covered","Q0_mean","Q1_mean","log_sigma0_mean","log_sigma1_mean",
        "shared_shape_posterior","shared_shape_bayes_factor",
        "frequency_variance","positive_model_variance","log_cv2_0_mean","log_cv2_1_mean"};
    Rcpp::NumericMatrix out(ng,names.size());std::fill(out.begin(),out.end(),NA_REAL);
    Rcpp::colnames(out)=Rcpp::wrap(names);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,8)
    #endif
    for (int i=0;i<ng;++i) {
        out(i,6)=0;
        if (!(usable[i] && usable[i+ng])) continue;
        const double shared=joint_pair ? 0 : rho*bayes_factor[i]/(1-rho+rho*bayes_factor[i]);
        out(i,15)=joint_pair ? NA_REAL : shared;
        out(i,16)=joint_pair ? NA_REAL : bayes_factor[i];
        std::vector<double> pair_posterior(joint_weight.size());double pair_sum=0;
        if (joint_pair) {
            double peak=R_NegInf;
            for (int h=0;h<int(joint_weight.size());++h) {
                pair_posterior[h]=std::log(joint_weight[h])+loglike(i,pair0[h])+loglike(i+ng,pair1[h]);
                peak=std::max(peak,pair_posterior[h]);
            }
            for (double& value:pair_posterior) {value=std::exp(value-peak);pair_sum+=value;}
        }
        double f0lower=0,f0upper=0,ftlower=0,ftupper=0,represented=0,finite_mass=0;
        double first_moment=0,second_moment=0,frequency_moment=0,contrast_second_moment=0;
        double qmean[2]={0,0},smean[2]={0,0},cvmean[2]={0,0};
        const double shift=truth[i]*std::log(2.0)/grid_step;
        const int shift_lower=std::isfinite(shift) ? int(std::floor(shift)) : 0;
        const int shift_upper=std::isfinite(shift) ? int(std::ceil(shift)) : 0;
        const int components=joint_pair ? joint_weight.size() : (shared_shape ? nq+1 : 1);
        for (int component=0;component<components;++component) {
            const int q=component-1;
            const double component_weight=joint_pair ? pair_posterior[component]/pair_sum :
                (component==0 ? 1-shared :
                rho*posterior_shape[i*nq+q]*posterior_shape[(i+ng)*nq+q]/
                    (prior_shape[q]*(1-rho+rho*bayes_factor[i])));
            if (!(component_weight>0)) continue;
            if (joint_pair) {
                const int k0=pair0[component],k1=pair1[component];
                const int row0=count_slot[positive[i]]*nk+k0,
                    row1=count_slot[maximum_n+1+positive[i+ng]]*nk+k1;
                const double mass0=cache_mass[row0],mass1=cache_mass[row1],
                    mean0=cache_mean[row0],mean1=cache_mean[row1],
                    variance_sum=cache_variance[row0]+cache_variance[row1],
                    difference=mean1-mean0,component_mass=component_weight*mass0*mass1;
                qmean[0]+=component_weight*shape[k0];qmean[1]+=component_weight*shape[k1];
                smean[0]+=component_weight*log_sigma[k0];smean[1]+=component_weight*log_sigma[k1];
                represented+=component_mass;
                finite_mass+=component_weight*(std::isfinite(K[k0]) && K[k0]>0)*
                    (std::isfinite(K[k1]) && K[k1]>0);
                first_moment+=component_mass*difference;
                cvmean[0]+=component_mass*mean0;cvmean[1]+=component_mass*mean1;
                second_moment+=component_mass*(variance_sum+difference*difference);
                frequency_moment+=component_mass*variance_sum;
                contrast_second_moment+=component_mass*difference*difference;
                const size_t index=size_t(count_pair_index[i])*components+component;
                f0lower+=component_weight*pair_cdf_lower[index];
                f0upper+=component_weight*pair_cdf_upper[index];
                if (std::isfinite(truth[i])) {
                    if (shift_lower==0 && shift_upper==0) {
                        ftlower+=component_weight*pair_cdf_lower[index];
                        ftupper+=component_weight*pair_cdf_upper[index];
                    } else {
                        const size_t offset0=size_t(row0)*nd,offset1=size_t(row1)*nd;
                        for (int j=0;j<nd-1;++j) {
                            const double bin=component_weight*
                                std::max(0.0,cache[offset0+j+1]-cache[offset0+j]);
                            const int lo=j+shift_lower,hi=j+1+shift_upper;
                            ftlower+=bin*(lo<0 ? 0 : cache[offset1+std::min(lo,nd-1)]);
                            ftupper+=bin*(hi<0 ? 0 : cache[offset1+std::min(hi,nd-1)]);
                        }
                    }
                }
                continue;
            }
            std::vector<double> cdf[2]={std::vector<double>(nd,0),std::vector<double>(nd,0)};
            double finite[2]={0,0},mean[2]={0,0},variance[2]={0,0},mass[2]={0,0};
            for (int g=0;g<2;++g) {
                const int row=i+g*ng,m=positive[row];
                const int first_atom=joint_pair ? (g==0 ? pair0[component] : pair1[component]) : 0;
                const int last_atom=joint_pair ? first_atom+1 : nk;
                for (int k=first_atom;k<last_atom;++k) {
                    if (!joint_pair && component>0 && shape_index[k]!=q) continue;
                    const double p=joint_pair ? 1 :
                        posterior[row*nk+k]/(component==0 ? 1 : posterior_shape[row*nq+q]);
                    qmean[g]+=component_weight*p*shape[k];smean[g]+=component_weight*p*log_sigma[k];
                    if (!std::isfinite(K[k]) || !(K[k]>0)) continue;
                    finite[g]+=p;
                    const size_t offset=size_t(count_slot[g*(maximum_n+1)+m]*nk+k)*nd;
                    for (int j=0;j<nd;++j) cdf[g][j]+=p*cache[offset+j];
                }
                double first=0,second=0;
                for (int j=0;j<nd-1;++j) {
                    const double bin=std::max(0.0,cdf[g][j+1]-cdf[g][j]),t=lower+(j+.5)*grid_step;
                    first+=bin*t;second+=bin*t*t;
                }
                mass[g]=cdf[g].back()-cdf[g].front();
                if (mass[g]>0) {
                    mean[g]=first/mass[g];variance[g]=std::max(0.0,second/mass[g]-mean[g]*mean[g]);
                }
            }
            const double component_mass=component_weight*mass[0]*mass[1],difference=mean[1]-mean[0];
            represented+=component_mass;finite_mass+=component_weight*finite[0]*finite[1];
            first_moment+=component_mass*difference;
            cvmean[0]+=component_mass*mean[0];cvmean[1]+=component_mass*mean[1];
            second_moment+=component_mass*(variance[0]+variance[1]+difference*difference);
            frequency_moment+=component_mass*(variance[0]+variance[1]);
            contrast_second_moment+=component_mass*difference*difference;
            for (int j=0;j<nd-1;++j) {
                const double bin=component_weight*std::max(0.0,cdf[0][j+1]-cdf[0][j]);
                f0lower+=bin*cdf[1][j];f0upper+=bin*cdf[1][j+1];
                ftlower+=bin*grid_cdf(cdf[1],j+shift_lower);
                ftupper+=bin*grid_cdf(cdf[1],j+1+shift_upper);
            }
        }
        const double undefined=std::max(0.0,1-finite_mass),unknown=std::max(0.0,1-represented);
        out(i,11)=qmean[0];out(i,12)=qmean[1];out(i,13)=smean[0];out(i,14)=smean[1];
        out(i,7)=undefined;out(i,8)=std::max(0.0,unknown-undefined);
        out(i,2)=f0lower;out(i,3)=f0upper;
        out(i,4)=std::isfinite(truth[i]) ? ftlower : NA_REAL;
        out(i,5)=std::isfinite(truth[i]) ? ftupper : NA_REAL;
        if (unknown>alpha) continue;
        const double mean=first_moment/represented,variance=std::max(0.0,second_moment/represented-mean*mean);
        out(i,6)=1;out(i,0)=mean/std::log(2.0);out(i,1)=std::sqrt(variance)/std::log(2.0);
        out(i,19)=cvmean[0]/represented;out(i,20)=cvmean[1]/represented;
        if (joint_pair) {
            out(i,17)=frequency_moment/(represented*std::log(2.0)*std::log(2.0));
            out(i,18)=std::max(0.0,contrast_second_moment/represented-mean*mean)/(std::log(2.0)*std::log(2.0));
        }
        out(i,9)=2*std::max(0.0,std::min({.5,f0upper+unknown,1-f0lower}));
        out(i,10)=std::isfinite(truth[i]) && ftlower<=1-alpha/2 && ftupper+unknown>=alpha/2;
    }
    return Rcpp::List::create(Rcpp::Named("genes")=out,
        Rcpp::Named("shared_shape_probability")=(joint_pair ? NA_REAL : rho),
        Rcpp::Named("log_positive_moment_ratio")=Rcpp::wrap(K),
        Rcpp::Named("cdf_step")=grid_step,
        Rcpp::Named("cdf_lower")=lower,Rcpp::Named("cdf_upper")=lower+(nd-1)*grid_step,
        Rcpp::Named("cdf_grid_points")=nd);
}

struct MixtureIntegral {double log_likelihood,error;int nodes;double location_mean=0,location_variance=0;};
struct LocationGrid {double lower,upper,spacing,error,prune_probability;};

// Given a latent allocation, the location posterior is a translation of one
// log-gamma density (normal for Q=0). Its characteristic-function modulus is
// |Gamma(a+it)/Gamma(a)|, a=m/Q^2, t=sigma*omega/|Q|. The gamma product gives
// log|CF| <= -a*(r*atan(r)-log(1+r^2)/2), r=t/a. Convexity bounds every alias
// in Poisson summation by a geometric series. Translations and their mixtures
// preserve this bound. The resulting rule has a uniform relative-error bound
// 2*exp(-B)/(1-exp(-B))+2*tail_probability, B=log(4/tolerance).
// One quarter of the tolerance is reserved for bounded node pruning below.
static LocationGrid mixture_location_grid(double m,double q,double tail_probability,double tolerance) {
    const double pi=3.14159265358979323846,B=std::log(4/tolerance);
    LocationGrid out;
    if (q==0) {
        out.lower=R::qnorm(tail_probability,0,1,true,false)/std::sqrt(m);
        out.upper=-out.lower;out.spacing=std::sqrt(2*pi*pi/(m*B));
    } else {
        const double a=m/(q*q);
        const double quantile=R::qgamma(tail_probability,a,1/a,true,false);
        // qgamma can underflow for one positive observation and a small shape.
        // In that regime F(x)=x^a/Gamma(a+1) to double precision.
        const double log_quantile=quantile>0 ? std::log(quantile) :
            (std::log(tail_probability)+R::lgammafn(a+1))/a-std::log(a);
        const double x=-log_quantile/q;
        const double z=-std::log(R::qgamma(tail_probability,a,1/a,false,false))/q;
        out.lower=std::min(x,z);out.upper=std::max(x,z);
        double low=0,high=std::sqrt(2*B/a);
        while (a*(high*std::atan(high)-.5*std::log1p(high*high))<B) high*=2;
        for (int iteration=0;iteration<60;++iteration) {
            const double middle=(low+high)/2;
            if (a*(middle*std::atan(middle)-.5*std::log1p(middle*middle))<B) low=middle;else high=middle;
        }
        out.spacing=2*pi/(std::abs(q)*a*high);
    }
    out.prune_probability=tolerance/4;
    const double relative=2*std::exp(-B)/(1-std::exp(-B))+2*tail_probability+out.prune_probability;
    out.error=-std::log1p(-relative);
    return out;
}

// Compute the analytic GG location once for all mixture proportions and shifts.
// The location derivative and objective used by profile tests are not needed here.
static constexpr int mixture_allocation_maximum=5;
struct MixturePrepared {
    double sigma,constant,location,pure_log_likelihood;
    std::vector<double> standardized,exponential;
};
static MixturePrepared prepare_mixture(const std::vector<double>& y,
    const std::vector<double>& centered,double center,double s,double q,
    bool exact_allocations=false) {
    MixturePrepared out;out.sigma=std::exp(s);
    out.constant=y.size()*(lgg_constant(q)-s);
    double location=center;
    if (q!=0) {
        double maximum=R_NegInf,sum=0;
        for (double value:centered) maximum=std::max(maximum,q*value/out.sigma);
        if (maximum<.5) {
            for (double value:centered) sum+=std::expm1(q*value/out.sigma)/y.size();
            location+=out.sigma/q*std::log1p(sum);
        } else {
            for (double value:centered) sum+=std::exp(q*value/out.sigma-maximum)/y.size();
            location+=out.sigma/q*(maximum+std::log(sum));
        }
    }
    out.location=location;
    out.pure_log_likelihood=exact_allocations && y.size()>=2 &&
        y.size()<=mixture_allocation_maximum ? location_integral_constant(y.size(),q)+
        marginal_scale_kernel(centered,s,q).log_density : 0;
    out.standardized.resize(y.size());out.exponential.resize(y.size());
    for (int j=0;j<int(y.size());++j) {
        out.standardized[j]=(y[j]-location)/out.sigma;
        out.exponential[j]=std::exp(q*out.standardized[j]);
    }
    return out;
}

// With at most five positive observations, sum every latent allocation.
// Integrating the common GG location is analytic for each allocation.
static long double allocation_log1pmx(long double x) {
    if (std::abs(x)>.02L) return std::log1p(x)-x;
    long double term=-x*x/2,result=term;
    for (int k=3;k<=12;++k) {term*=-x*(k-1)/k;result+=term;}
    return result;
}
static long double allocation_expm1mx(long double x) {
    long double term=x*x/2,result=term;
    for (int k=3;k<=12;++k) {term*=x/k;result+=term;}
    return result;
}
static std::vector<MixtureIntegral> allocation_mixture_integrals(
    const MixturePrepared& data,double q,double relative_shift,
    const std::vector<int>& atoms,const Rcpp::NumericVector& epsilon) {
    const int m=data.standardized.size(),n=1<<m,ne=atoms.size(),full=n-1;
    const long double r=relative_shift,t=-q*r,a=q==0 ? 1 : m/(static_cast<long double>(q)*q);
    std::vector<int> count(n,0);
    std::vector<long double> sum_exp(n,0),sum_expm1(n,0),sum_z(n,0),ratio(n,0),shift(n,0);
    for (int mask=1;mask<n;++mask) {
        const int previous=mask&(mask-1),bit=mask^previous;
        int j=0;for (int value=bit;value>1;value>>=1) ++j;
        count[mask]=count[previous]+1;
        sum_exp[mask]=sum_exp[previous]+data.exponential[j];
        sum_expm1[mask]=sum_expm1[previous]+std::expm1(static_cast<long double>(q)*data.standardized[j]);
        sum_z[mask]=sum_z[previous]+data.standardized[j];
    }
    long double location=data.location,variance=data.sigma*data.sigma/m;
    if (q==0) location+=data.sigma*sum_z[full]/m;
    else {
        location+=data.sigma/q*std::log1p(sum_expm1[full]/m);
        const long double u=1/a,u2=u*u;
        const long double correction=a>=16 ?
            .5L+u*(1.L/12+u2*(-1.L/120+u2*(1.L/252+u2*(-1.L/240+
                u2*(1.L/132+u2*(-691.L/32760+u2/12)))))) :
            a*(std::log(a)-R::digamma(static_cast<double>(a)));
        location+=data.sigma*q/m*correction;
        variance*=a*R::trigamma(static_cast<double>(a));
    }
    shift[full]=-data.sigma*r;
    for (int mask=1;mask<full;++mask) {
        const long double T=count[mask],fraction=T/m;
        if (q==0) {
            ratio[mask]=r*(sum_z[mask]-fraction*sum_z[full])-r*r*T*(1-fraction)/2;
            shift[mask]=-data.sigma*r*fraction;
        } else {
            long double logR;
            if (std::abs(t)<.01L) {
                const long double p=sum_exp[mask]/sum_exp[full],v=p*std::expm1(t);
                const long double dev=(m*sum_expm1[mask]-T*sum_expm1[full])/(m+sum_expm1[full]);
                logR=std::log1p(v);
                ratio[mask]=r*dev/q-a*(allocation_log1pmx(v)+p*allocation_expm1mx(t));
            } else {
                const long double left=std::log(sum_exp[mask])+t,right=std::log(sum_exp[full^mask]);
                const long double maximum=std::max(left,right);
                logR=maximum+std::log1p(std::exp(-std::abs(left-right)))-std::log(sum_exp[full]);
                ratio[mask]=-r*T/q-a*logR;
            }
            shift[mask]=data.sigma/q*logR;
        }
    }
    std::vector<MixtureIntegral> out(ne);
    for (int e=0;e<ne;++e) {
        const long double p=epsilon[atoms[e]];
        std::vector<long double> values(n);
        long double peak=-std::numeric_limits<long double>::infinity();
        for (int mask=0;mask<n;++mask) {
            const int T=count[mask];
            const long double prior=(T==0 ? 0 : T*std::log(p))+
                (T==m ? 0 : (m-T)*std::log1p(-p));
            values[mask]=ratio[mask]+prior;peak=std::max(peak,values[mask]);
        }
        long double total=0,first=0,second=0;
        for (int mask=0;mask<n;++mask) {
            const long double w=std::exp(values[mask]-peak);
            total+=w;first+=w*shift[mask];second+=w*shift[mask]*shift[mask];
        }
        const long double mean=first/total;
        out[e]={static_cast<double>(data.pure_log_likelihood+peak+std::log(total)),0,0,
            static_cast<double>(location+mean),
            static_cast<double>(variance+std::max(0.L,second/total-mean*mean))};
    }
    return out;
}

// All proportions at the same scale, shape and relative shift use the same
// quadrature nodes and component densities. Each integral retains its own peak
// and pruning bound; grouping changes computation, not the likelihood model.
static std::vector<MixtureIntegral> positive_mixture_integrals(
    const MixturePrepared& data,double q,double relative_shift,const LocationGrid& grid,
    const std::vector<int>& atoms,const Rcpp::NumericVector& epsilon,
    bool exact_allocations=false) {
    if (exact_allocations && data.standardized.size()<=mixture_allocation_maximum)
        return allocation_mixture_integrals(data,q,relative_shift,atoms,epsilon);
    const int first=int(std::floor((grid.lower-relative_shift)/grid.spacing))-1;
    const int last=int(std::ceil(grid.upper/grid.spacing))+1,n=last-first+1,ne=atoms.size();
    const double lower=data.sigma*first*grid.spacing,upper=data.sigma*last*grid.spacing;
    const double half=(upper-lower)/2,middle=(upper+lower)/2;
    const double shift_factor=std::exp(-q*relative_shift),cutoff=std::log(n/grid.prune_probability);
    std::vector<double> probability(ne),peak(ne,R_NegInf),values(size_t(ne)*n);
    double average_epsilon=0,minimum_factor=1;
    for (int e=0;e<ne;++e) {
        probability[e]=epsilon[atoms[e]];
        average_epsilon+=epsilon[atoms[e]]/ne;
        minimum_factor=std::min({minimum_factor,probability[e],1-probability[e]});
    }
    // Accumulate at most 16 probability factors before taking their logarithm.
    // This only groups arithmetic; the smaller bound prevents underflow for
    // extreme mixture probabilities. Unfinished products are at most one,
    // preserving the same rigorous upper bound used to prune a quadrature node.
    const int accumulation_block=std::min(16,std::max(1,
        int(std::floor(std::log(std::numeric_limits<double>::min())/std::log(minimum_factor)))-1));
    const double expected_location=q==0 ? -average_epsilon*relative_shift :
        std::log1p(average_epsilon*std::expm1(-q*relative_shift))/q;
    const int pivot=int(std::lround((data.sigma*expected_location-lower)*(n-1)/(upper-lower)));
    std::vector<int> order(n);order[0]=pivot;
    for (int distance=1,position=1;position<n;++distance) {
        if (pivot-distance>=0) order[position++]=pivot-distance;
        if (pivot+distance<n) order[position++]=pivot+distance;
    }
    const int sample_count=data.standardized.size();
    std::vector<double> value(ne),product(ne),component_difference(sample_count);
    const double inverse_q=q==0 ? 0 : 1/q,shift_expm1=std::expm1(-q*relative_shift);
    const double log_shift_difference=q==0 ? 0 :
        (q>0 ? std::log(-std::expm1(-q*relative_shift)) :
            -q*relative_shift+std::log(-std::expm1(q*relative_shift)))-2*std::log(std::abs(q));
    const Eigen::Map<const Eigen::ArrayXd> proportions(probability.data(),ne);
    const Eigen::ArrayXd complements=1-proportions;
    Eigen::Map<Eigen::ArrayXd> products(product.data(),ne),log_products(value.data(),ne);
    for (int step=0;step<n;++step) {
        const int h=order[step];
        const double node=-1.0+2.0*h/(n-1);
        const double location=(middle+half*node)/data.sigma,location_factor=std::exp(-q*location);
        std::fill(value.begin(),value.end(),0);
        std::fill(product.begin(),product.end(),1);
        // h1-h0 is affine in exp(Q*z); evaluate only the smaller component penalty.
        double common=data.constant;int remaining=ne;
        const double slope=q==0 ? -relative_shift : location_factor*shift_expm1*inverse_q*inverse_q,
            intercept=q==0 ? relative_shift*(location+relative_shift/2) : relative_shift*inverse_q;
        for(int j=0;j<sample_count;++j) {
            double delta=intercept+slope*(q==0 ? data.standardized[j] : data.exponential[j]);
            if(q!=0 && !std::isfinite(delta)) {
                const double difference=std::exp(q*(data.standardized[j]-location)+log_shift_difference);
                delta=intercept+(q>0 ? -difference : difference);
            }
            const bool shifted=delta<0;
            const double z=data.standardized[j]-location-(shifted ? relative_shift : 0),t=q*z;
            double exponential=data.exponential[j]*location_factor*(shifted ? shift_factor : 1);
            if(q!=0 && (!std::isfinite(exponential) || exponential==0))exponential=std::exp(t);
            const double minimum=q==0 ? z*z/2 :
                (std::abs(t)<.005 ? std::expm1(t)-t : exponential-1-t)*inverse_q*inverse_q;
            component_difference[j]=delta;common-=std::max(minimum,0.0);
        }
        for (int e=0;e<ne;++e)
            if (common<peak[e]-cutoff) {value[e]=R_NegInf;--remaining;}
        for (int first_sample=0;first_sample<sample_count && remaining>0;first_sample+=8) {
            const int block_samples=std::min(8,sample_count-first_sample);
            Eigen::Array<double,8,1> difference;
            for (int offset=0;offset<8;++offset) {
                const int j=first_sample+std::min(offset,block_samples-1);
                difference[offset]=std::abs(component_difference[j]);
            }
            const Eigen::Array<double,8,1> ratios=(-difference).exp();
            for (int offset=0;offset<block_samples && remaining>0;++offset) {
                const int j=first_sample+offset;
                const bool base_dominates=component_difference[j]>=0;
                const double ratio=ratios[offset];
                const bool accumulate=(j+1)%accumulation_block==0 || j+1==sample_count;
                if(base_dominates)products*=1+proportions*(ratio-1);
                else products*=proportions+complements*ratio;
                if(accumulate){
                    log_products+=products.log();products.setOnes();
                    remaining=0;
                    for(int e=0;e<ne;++e){
                        if(common+value[e]<peak[e]-cutoff)value[e]=R_NegInf;
                        else ++remaining;}
                }
            }
        }
        for (int e=0;e<ne;++e) {
            value[e]+=common;
            values[size_t(e)*n+h]=value[e];peak[e]=std::max(peak[e],value[e]);
        }
    }
    std::vector<MixtureIntegral> out(ne);
    for (int e=0;e<ne;++e) {
        double sum=0,first_moment=0,second_moment=0;
        for (int h=0;h<n;++h) {
            const double weight=2.0/(n-1)*std::exp(values[size_t(e)*n+h]-peak[e]);
            const double location=(first+h)*grid.spacing;
            sum+=weight;first_moment+=weight*location;second_moment+=weight*location*location;
        }
        const double mean=first_moment/sum;
        out[e]={peak[e]+std::log(half*sum),grid.error,n,
            data.location+data.sigma*mean,
            data.sigma*data.sigma*std::max(0.0,second_moment/sum-mean*mean)};
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List population_mixture_likelihood_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,Rcpp::NumericVector epsilon,
    Rcpp::NumericVector relative_shift,int threads=1,double tolerance=1e-7,
    double tail_probability=1e-12,std::string checkpoint_directory="",int checkpoint_rows=256,
    int minimum_positive=4,Rcpp::Nullable<Rcpp::NumericMatrix> previous_loglike=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> previous_parameters=R_NilValue,
    bool previous_reversed=false,double previous_tolerance=1e-7,
    double previous_tail_probability=1e-12,bool sparse_output=false,
    SEXP connection=R_NilValue,double sparse_tolerance=1e-6,int sparse_block_rows=256) {
    const int ng=x.nrow(),nk=shape.size();int total[2]={0,0};
    const bool has_previous=previous_loglike.isNotNull();
    const Rcpp::NumericMatrix prior=has_previous ?
        Rcpp::NumericMatrix(previous_loglike) : Rcpp::NumericMatrix(0,0);
    const Rcpp::NumericMatrix prior_parameters=has_previous ?
        Rcpp::NumericMatrix(previous_parameters) : Rcpp::NumericMatrix(0,0);
    if(has_previous && (prior.nrow()!=2*ng || prior.ncol()!=prior_parameters.nrow() ||
        prior_parameters.ncol()!=4))Rcpp::stop("Cached kernel dimensions differ");
    std::vector<int> previous_column(nk,-1);
    std::map<std::tuple<double,double,double,double>,int> previous_lookup;
    if(has_previous) {
        for(int k=0;k<prior.ncol();++k)previous_lookup.emplace(std::make_tuple(
            prior_parameters(k,0),prior_parameters(k,1),prior_parameters(k,2),prior_parameters(k,3)),k);
        for(int k=0;k<nk;++k) {
            const auto found=previous_lookup.find(std::make_tuple(
                log_sigma[k],shape[k],epsilon[k],relative_shift[k]));
            if(found!=previous_lookup.end())previous_column[k]=found->second;
        }
    }
    for (int g:group) ++total[g];
    const int maximum_n=std::max(total[0],total[1]);
    std::vector<double> distinct(shape.begin(),shape.end());
    std::sort(distinct.begin(),distinct.end());
    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
    const int nq=distinct.size();std::vector<int> qindex(nk);
    for (int k=0;k<nk;++k) qindex[k]=std::lower_bound(distinct.begin(),distinct.end(),shape[k])-distinct.begin();
    std::map<std::pair<double,double>,std::map<double,std::vector<int>>> batches,missing_batches;
    for (int k=0;k<nk;++k) if (epsilon[k]>0) {
        batches[{log_sigma[k],shape[k]}][relative_shift[k]].push_back(k);
        if(previous_column[k]<0)missing_batches[{log_sigma[k],shape[k]}][relative_shift[k]].push_back(k);
    }
    std::vector<LocationGrid> grid((maximum_n+1)*nq);
    for (int m=std::max(2,minimum_positive);m<=maximum_n;++m) for (int q=0;q<nq;++q)
        grid[m*nq+q]=mixture_location_grid(m,distinct[q],tail_probability,tolerance);
    if(sparse_output && (connection==R_NilValue || sparse_block_rows<1 || !(sparse_tolerance>0)))
        Rcpp::stop("Sparse mixture output requires a binary connection and positive block/tolerance");
    Rcpp::NumericMatrix loglike(sparse_output ? std::min(2*ng,sparse_block_rows) : 2*ng,nk);
    Rcpp::IntegerVector count(2*ng),maximum_nodes(2*ng);
    Rcpp::LogicalVector usable(2*ng);Rcpp::NumericVector error(2*ng),reused(2*ng);
    std::unique_ptr<SparsePopulationBuffer> sparse;
    std::vector<int> sparse_row;
    std::vector<size_t> sparse_count;
    if(sparse_output) {
        sparse=std::make_unique<SparsePopulationBuffer>();
        sparse_row.assign(2*ng,-1);sparse_count.assign(nk,0);
        for(int row=0;row<2*ng;++row) {
            const int gene=row%ng,g=row/ng;int m=0;double center=0,ss=0;
            for(int j=0;j<x.ncol();++j)if(group[j]==g && x(gene,j)>0)++m;
            // Match the row's original summation order, including division
            // before addition when finding its log-location center.
            for(int j=0;j<x.ncol();++j)if(group[j]==g && x(gene,j)>0)
                center+=std::log(x(gene,j))/m;
            for(int j=0;j<x.ncol();++j)if(group[j]==g && x(gene,j)>0) {
                const double residual=std::log(x(gene,j))-center;ss+=residual*residual;
            }
            count[row]=m;usable[row]=m>=minimum_positive && (minimum_positive<4 || ss>0);
            if(usable[row]) {sparse_row[row]=sparse->rows.size();sparse->rows.push_back(row);}
        }
        sparse->n=sparse->rows.size();sparse->nk=nk;
        sparse->omitted_bound=sparse_tolerance/4;
        sparse->threshold=sparse->omitted_bound/sparse->n;
        sparse->row_scale.resize(sparse->n);sparse->umbrella.resize(sparse->n);
    }
    std::vector<LocationGrid> cached_grid=grid;
    if(has_previous) {
        if(mixture_location_grid(6,0,previous_tail_probability,previous_tolerance).error>
            mixture_location_grid(6,0,tail_probability,tolerance).error)
            Rcpp::stop("Cached quadrature error exceeds the requested bound");
        for(int m=std::max(6,minimum_positive);m<=maximum_n;++m)
            for(int q=0;q<nq;++q)cached_grid[m*nq+q]=
                mixture_location_grid(m,distinct[q],previous_tail_probability,previous_tolerance);
    }
    std::vector<int> cached_nodes(maximum_n+1,0);
    std::vector<double> cached_error(maximum_n+1,0);
    if(has_previous)for(int m=std::max(6,minimum_positive);m<=maximum_n;++m)
        for(int k=0;k<nk;++k)if(epsilon[k]>0 && previous_column[k]>=0) {
            const auto& at=cached_grid[m*nq+qindex[k]];
            const int first=int(std::floor((at.lower-relative_shift[k])/at.spacing))-1,
                last=int(std::ceil(at.upper/at.spacing))+1;
            cached_nodes[m]=std::max(cached_nodes[m],last-first+1);
            cached_error[m]=std::max(cached_error[m],at.error);
        }
    std::fill(loglike.begin(),loglike.end(),NA_REAL);
    Rcpp::Function file_exists=Rcpp::Environment::base_env()["file.exists"];
    Rcpp::Function read=Rcpp::Environment::base_env()["readRDS"];
    const int block_size=sparse_output ? sparse_block_rows :
        (checkpoint_directory.empty() ? 2*ng : checkpoint_rows);
    for (int begin=0;begin<2*ng;begin+=block_size) {
        const int end=std::min(begin+block_size,2*ng),block_rows=end-begin;
        if(sparse_output) {
            if(loglike.nrow()!=block_rows)loglike=Rcpp::NumericMatrix(block_rows,nk);
            std::fill(loglike.begin(),loglike.end(),NA_REAL);
        }
        const std::string path=population_checkpoint_path(checkpoint_directory,"rows_",begin);
        if (!sparse_output && !checkpoint_directory.empty() && Rcpp::as<bool>(file_exists(Rcpp::String(path,CE_UTF8)))) {
            const Rcpp::List saved=read(Rcpp::String(path,CE_UTF8));
            const Rcpp::NumericMatrix values=saved["log_likelihood"];
            const Rcpp::IntegerVector positive=saved["positive"],nodes=saved["maximum_nodes"];
            const Rcpp::LogicalVector available=saved["usable"];
            const Rcpp::NumericVector bound=saved["quadrature_error"];
            if (Rcpp::as<int>(saved["first_row"])!=begin || Rcpp::as<int>(saved["total_rows"])!=2*ng ||
                values.nrow()!=block_rows || values.ncol()!=nk) Rcpp::stop("Checkpoint likelihood mismatch");
            for (int k=0;k<nk;++k) std::copy_n(values.begin()+size_t(k)*block_rows,block_rows,
                loglike.begin()+size_t(k)*2*ng+begin);
            for (int i=0;i<block_rows;++i) {
                count[begin+i]=positive[i];maximum_nodes[begin+i]=nodes[i];
                usable[begin+i]=available[i];error[begin+i]=bound[i];
            }
            Rcpp::Rcout<<"likelihood rows="<<end<<"/"<<2*ng<<" restored\n";
            continue;
        }
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for (int row=begin;row<end;++row) {
        const int output_row=sparse_output ? row-begin : row;
        const int gene=row%ng,g=row/ng;std::vector<double> y;double center=0,ss=0;
        for (int j=0;j<x.ncol();++j) if (group[j]==g && x(gene,j)>0) y.push_back(std::log(x(gene,j)));
        const int m=y.size();count[row]=m;
        for (double value:y) center+=value/m;
        for (double value:y) ss+=(value-center)*(value-center);
        usable[row]=m>=minimum_positive && (minimum_positive<4 || ss>0);
        if (!usable[row]) continue;
        if (m==1) {
            for (int k=0;k<nk;++k) loglike(output_row,k)=0;
            continue;
        }
        std::sort(y.begin(),y.end());
        std::vector<double> ordered(m);
        for (int j=0;j<m;++j) ordered[j]=(j%2==0 ? y[j/2] : y[m-1-j/2]);
        y.swap(ordered);
        double ordered_center=0;for (double value:y) ordered_center+=value/m;
        std::vector<double> centered=y;for (double& value:centered) value-=ordered_center;
        // Small-sample allocation sums and CV quadrature are distinct rules.
        const bool reuse_row=has_previous && m>mixture_allocation_maximum;
        const int previous_row=previous_reversed ? (1-g)*ng+gene : row;
        if(reuse_row){maximum_nodes[row]=cached_nodes[m];error[row]=cached_error[m];}
        for(int k=0;k<nk;++k) {
            if(previous_column[k]>=0 && (epsilon[k]==0 || reuse_row)) {
                loglike(output_row,k)=prior(previous_row,previous_column[k]);reused[row]+=1;
            } else if(epsilon[k]==0)loglike(output_row,k)=location_integral_constant(m,shape[k])+
                marginal_scale_kernel(centered,log_sigma[k],shape[k]).log_density;
        }
        const auto& active_batches=reuse_row ? missing_batches : batches;
        for (const auto& parameter:active_batches) {
            const double s=parameter.first.first,q=parameter.first.second;
            const int qi=std::lower_bound(distinct.begin(),distinct.end(),q)-distinct.begin();
            const MixturePrepared prepared=prepare_mixture(y,centered,ordered_center,s,q,minimum_positive<4);
            for (const auto& batch:parameter.second) {
                const std::vector<MixtureIntegral> integrals=positive_mixture_integrals(
                    prepared,q,batch.first,grid[m*nq+qi],batch.second,epsilon,minimum_positive<4);
                for (int j=0;j<int(batch.second.size());++j) {
                    const int k=batch.second[j];const MixtureIntegral& at=integrals[j];
                    loglike(output_row,k)=at.log_likelihood;error[row]=std::max(error[row],at.error);
                    maximum_nodes[row]=std::max(maximum_nodes[row],at.nodes);
                }
            }
        }
    }
        if(sparse_output) {
            for(int i=0;i<block_rows;++i) {
                const int slot=sparse_row[begin+i];if(slot<0)continue;
                double largest=R_NegInf;
                for(int k=0;k<nk;++k)largest=std::max(largest,loglike(i,k));
                sparse->row_scale[slot]=largest;
            }
            std::vector<double> partial(size_t(threads)*block_rows,0);
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(threads)
            #endif
            for(int k=0;k<nk;++k) {
                int thread=0;
                #ifdef _OPENMP
                thread=omp_get_thread_num();
                #endif
                double* sum=partial.data()+size_t(thread)*block_rows;
                for(int i=0;i<block_rows;++i) {
                    const int slot=sparse_row[begin+i];if(slot<0)continue;
                    const double value=std::exp(loglike(i,k)-sparse->row_scale[slot]);
                    if(value>sparse->threshold){++sparse_count[k];sum[i]+=value/nk;}
                }
            }
            for(int i=0;i<block_rows;++i) {
                const int slot=sparse_row[begin+i];if(slot<0)continue;
                double sum=0;for(int thread=0;thread<threads;++thread)
                    sum+=partial[size_t(thread)*block_rows+i];
                sparse->umbrella[slot]=sum;
            }
            Rcpp::Function write=Rcpp::Environment::base_env()["writeBin"];
            const Rcpp::NumericVector block(loglike.begin(),loglike.end());
            write(block,connection,Rcpp::Named("size")=8,Rcpp::Named("endian")="little");
        }
        if (!sparse_output && !checkpoint_directory.empty()) {
            Rcpp::NumericMatrix values(block_rows,nk);
            for (int k=0;k<nk;++k) std::copy_n(loglike.begin()+size_t(k)*2*ng+begin,block_rows,
                values.begin()+size_t(k)*block_rows);
            const Rcpp::Range selected(begin,end-1);
            population_save_checkpoint(Rcpp::List::create(Rcpp::Named("first_row")=begin,
                Rcpp::Named("total_rows")=2*ng,Rcpp::Named("log_likelihood")=values,
                Rcpp::Named("positive")=count[selected],Rcpp::Named("usable")=usable[selected],
                Rcpp::Named("maximum_nodes")=maximum_nodes[selected],Rcpp::Named("quadrature_error")=error[selected]),path);
            Rcpp::Rcout<<"likelihood rows="<<end<<"/"<<2*ng<<" saved\n";
        }
        Rcpp::checkUserInterrupt();
    }
    if(sparse_output) {
        Rcpp::Function flush=Rcpp::Environment::base_env()["flush"],
            seek=Rcpp::Environment::base_env()["seek"],
            read_bin=Rcpp::Environment::base_env()["readBin"],
            gc=Rcpp::Environment::base_env()["gc"];
        flush(connection);loglike=Rcpp::NumericMatrix(0,0);gc(false);
        const double temporary_bytes=Rcpp::as<double>(seek(connection,
            Rcpp::Named("where")=NA_REAL,Rcpp::Named("rw")="write"));
        seek(connection,Rcpp::Named("where")=0,Rcpp::Named("origin")="start",Rcpp::Named("rw")="read");
        sparse->column_start.assign(size_t(nk)+1,0);
        for(int k=0;k<nk;++k)sparse->column_start[k+1]=sparse->column_start[k]+sparse_count[k];
        sparse->column_row.resize(sparse->column_start.back());
        sparse->column_value.resize(sparse->column_start.back());
        std::vector<size_t> at=sparse->column_start;
        for(int begin=0;begin<2*ng;begin+=block_size) {
            Rcpp::checkUserInterrupt();
            const int nr=std::min(block_size,2*ng-begin);
            const R_xlen_t length=R_xlen_t(nr)*nk;
            gc(false);
            const Rcpp::NumericVector block=read_bin(connection,Rcpp::Named("what")="double",
                Rcpp::Named("n")=double(length),Rcpp::Named("size")=8,Rcpp::Named("endian")="little");
            if(block.size()!=length)Rcpp::stop("Incomplete temporary mixture likelihood block");
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(threads)
            #endif
            for(int k=0;k<nk;++k) {
                size_t entry=at[k];
                for(int i=0;i<nr;++i) {
                    const int slot=sparse_row[begin+i];if(slot<0)continue;
                    const double value=std::exp(block[R_xlen_t(k)*nr+i]-sparse->row_scale[slot]);
                    if(value>sparse->threshold) {
                        sparse->column_row[entry]=slot;sparse->column_value[entry++]=value;
                    }
                }
                at[k]=entry;
            }
        }
        const double entries=double(sparse->column_value.size());
        Rcpp::XPtr<SparsePopulationBuffer> pointer(sparse.release(),true);
        return Rcpp::List::create(Rcpp::Named("dictionary")=pointer,
            Rcpp::Named("positive")=count,Rcpp::Named("usable")=usable,
            Rcpp::Named("quadrature_error")=error,Rcpp::Named("maximum_nodes")=maximum_nodes,
            Rcpp::Named("reused_cells")=std::accumulate(reused.begin(),reused.end(),0.0),
            Rcpp::Named("matched_columns")=nk-std::count(previous_column.begin(),previous_column.end(),-1),
            Rcpp::Named("temporary_bytes")=temporary_bytes,Rcpp::Named("sparse_entries")=entries,
            Rcpp::Named("integration_method")=minimum_positive<4 ?
                "analytic_allocation_le5_or_location_loggamma_trapezoid_shared_components" :
                "location_loggamma_trapezoid_shared_components");
    }
    return Rcpp::List::create(Rcpp::Named("log_likelihood")=loglike,Rcpp::Named("positive")=count,
        Rcpp::Named("usable")=usable,Rcpp::Named("quadrature_error")=error,
        Rcpp::Named("maximum_nodes")=maximum_nodes,
        Rcpp::Named("reused_cells")=std::accumulate(reused.begin(),reused.end(),0.0),
        Rcpp::Named("matched_columns")=nk-std::count(previous_column.begin(),previous_column.end(),-1),
        Rcpp::Named("integration_method")=minimum_positive<4 ?
            "analytic_allocation_le5_or_location_loggamma_trapezoid_shared_components" :
            "location_loggamma_trapezoid_shared_components");
}



// The temporary stream contains every input row, in row blocks. Columns use
// R's one-based indices; unavailable rows retain the original NA likelihoods.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix read_streamed_population_columns_cpp(SEXP connection,int total_rows,
    int total_columns,Rcpp::IntegerVector columns,int block_rows=256) {
    if(connection==R_NilValue || total_rows<0 || total_columns<0 || block_rows<1)
        Rcpp::stop("Invalid streamed likelihood dimensions or connection");
    for(int column:columns)if(column<1 || column>total_columns)
        Rcpp::stop("Streamed likelihood column is outside the original grid");
    Rcpp::Function seek=Rcpp::Environment::base_env()["seek"],
        read=Rcpp::Environment::base_env()["readBin"],gc=Rcpp::Environment::base_env()["gc"];
    seek(connection,Rcpp::Named("where")=0,Rcpp::Named("origin")="start",Rcpp::Named("rw")="read");
    Rcpp::NumericMatrix out(total_rows,columns.size());
    for(int begin=0;begin<total_rows;begin+=block_rows) {
        Rcpp::checkUserInterrupt();
        const int nr=std::min(block_rows,total_rows-begin);
        const R_xlen_t length=R_xlen_t(nr)*total_columns;
        gc(false);
        const Rcpp::NumericVector block=read(connection,Rcpp::Named("what")="double",
            Rcpp::Named("n")=double(length),Rcpp::Named("size")=8,Rcpp::Named("endian")="little");
        if(block.size()!=length)Rcpp::stop("Incomplete retained likelihood block");
        for(int j=0;j<columns.size();++j)
            std::copy_n(block.begin()+R_xlen_t(columns[j]-1)*nr,nr,
                out.begin()+R_xlen_t(j)*total_rows+begin);
    }
    return out;
}

// Positive-part arithmetic-mean location posterior under the same GG kernels.
// A flat log-location prior yields a proper location posterior with >=1 positive.
static double population_shift_logmean(double e,double D) {
    const double a=std::log1p(-e),b=std::log(e)+D,m=std::max(a,b);
    return m+std::log(std::exp(a-m)+std::exp(b-m));
}
// [[Rcpp::export]]
Rcpp::List population_mixture_mean_likelihood_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,Rcpp::NumericVector epsilon,
    Rcpp::NumericVector relative_shift,int threads=1,double tolerance=1e-10,
    double tail_probability=1e-13) {
    const int ng=x.nrow(),ns=shape.size();int total[2]={0,0};for (int g:group) ++total[g];
    const int maximum_n=std::max(total[0],total[1]);
    std::vector<double> distinct(shape.begin(),shape.end());std::sort(distinct.begin(),distinct.end());
    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
    const int nq=distinct.size();
    std::vector<LocationGrid> grid((maximum_n+1)*nq);
    for (int m=2;m<=maximum_n;++m) for (int j=0;j<nq;++j)
        grid[m*nq+j]=mixture_location_grid(m,distinct[j],tail_probability,tolerance);
    std::map<std::pair<double,double>,std::map<double,std::vector<int>>> batches;
    Rcpp::NumericVector moment(ns);Rcpp::LogicalVector mean_finite(ns);
    for (int s=0;s<ns;++s) {
        const double sigma=std::exp(log_sigma[s]),q=shape[s],k=q==0 ? R_PosInf : 1/(q*q);
        mean_finite[s]=q>=0 || k+sigma/q>0;
        moment[s]=mean_finite[s] ? log_moment(1,V3{0,log_sigma[s],q})+
            population_shift_logmean(epsilon[s],relative_shift[s]*sigma) : NA_REAL;
        if (epsilon[s]>0) batches[{log_sigma[s],shape[s]}][relative_shift[s]].push_back(s);
    }
    Rcpp::NumericMatrix value(2*ng,ns),alpha(2*ng,ns),variance(2*ng,ns);
    Rcpp::IntegerVector count(2*ng);Rcpp::LogicalVector usable(2*ng),cv_usable(2*ng);
    Rcpp::NumericVector error(2*ng);
    std::fill(value.begin(),value.end(),NA_REAL);std::fill(alpha.begin(),alpha.end(),NA_REAL);
    std::fill(variance.begin(),variance.end(),NA_REAL);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for (int row=0;row<2*ng;++row) {
        const int gene=row%ng,g=row/ng;std::vector<double> y;double center=0,ss=0;
        for (int j=0;j<x.ncol();++j) if (group[j]==g && x(gene,j)>0) y.push_back(std::log(x(gene,j)));
        const int m=y.size();count[row]=m;usable[row]=m>=1;
        if (m==0) continue;
        for (double v:y) center+=v/m;for (double v:y) ss+=(v-center)*(v-center);
        cv_usable[row]=m>=4 && ss>0;
        if (m==1) {
            for (int s=0;s<ns;++s) {
                const double sigma=std::exp(log_sigma[s]),q=shape[s],e=epsilon[s],D=relative_shift[s]*sigma;
                const double k=q==0 ? 1 : 1/(q*q);
                const double ez=q==0 ? 0 : (R::digamma(k)-std::log(k))/q;
                const double vz=q==0 ? 1 : R::trigamma(k)/(q*q);
                value(row,s)=0;
                if (mean_finite[s]) {alpha(row,s)=y[0]-sigma*ez-e*D+moment[s];
                    variance(row,s)=sigma*sigma*vz+e*(1-e)*D*D;}
            }
            continue;
        }
        std::sort(y.begin(),y.end());std::vector<double> ordered(m);
        for (int j=0;j<m;++j) ordered[j]=(j%2==0 ? y[j/2] : y[m-1-j/2]);
        y.swap(ordered);double ordered_center=0;for (double v:y) ordered_center+=v/m;
        std::vector<double> centered=y;for (double& v:centered) v-=ordered_center;
        for (int s=0;s<ns;++s) if (epsilon[s]==0) {
            const double q=shape[s],sigma=std::exp(log_sigma[s]);
            const MixturePrepared prepared=prepare_mixture(y,centered,ordered_center,log_sigma[s],q);
            value(row,s)=location_integral_constant(m,q)+marginal_scale_kernel(centered,log_sigma[s],q).log_density;
            if (mean_finite[s]) {
                const double k=q==0 ? 1 : 1/(q*q),h=q/sigma;
                alpha(row,s)=q==0 ? prepared.location+moment[s] :
                    prepared.location+(std::log(m*k)-R::digamma(m*k))/h+moment[s];
                variance(row,s)=q==0 ? sigma*sigma/m : R::trigamma(m*k)/(h*h);
            }
        }
        for (const auto& parameter:batches) {
            const double sig=parameter.first.first,q=parameter.first.second;
            const int qi=std::lower_bound(distinct.begin(),distinct.end(),q)-distinct.begin();
            const MixturePrepared prepared=prepare_mixture(y,centered,ordered_center,sig,q,true);
            for (const auto& batch:parameter.second) {
                const std::vector<MixtureIntegral> integrals=positive_mixture_integrals(prepared,q,
                    batch.first,grid[m*nq+qi],batch.second,epsilon,true);
                for (int j=0;j<int(batch.second.size());++j) {
                    const int s=batch.second[j];const MixtureIntegral& at=integrals[j];
                    value(row,s)=at.log_likelihood;error[row]=std::max(error[row],at.error);
                    if (mean_finite[s]) {alpha(row,s)=at.location_mean+moment[s];variance(row,s)=at.location_variance;}
                }
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("log_likelihood")=value,Rcpp::Named("log_mean")=alpha,
        Rcpp::Named("log_mean_variance")=variance,Rcpp::Named("positive")=count,
        Rcpp::Named("usable")=usable,Rcpp::Named("cv_usable")=cv_usable,
        Rcpp::Named("mean_finite")=mean_finite,Rcpp::Named("quadrature_error")=error);
}

#include <unsupported/Eigen/FFT>
#include <tuple>
#include <complex>

// Deterministic conditional location distributions. Each positive likelihood
// is evaluated on a common log-mean lattice; the beta log-frequency law has
// analytic cell masses. FFT correlations replace all nuisance-state pairs.
struct PopulationMeanGridState {
    MixturePrepared data;
    double q,epsilon,relative_shift,moment,log_marginal,weight,lower,upper;
    int index=0;
};
static double population_mean_grid_logadd(double a,double b) {
    if (a==R_NegInf) return b;if (b==R_NegInf) return a;
    const double peak=std::max(a,b);return peak+std::log1p(std::exp(std::min(a,b)-peak));
}
static double population_mean_grid_density(const PopulationMeanGridState& state,
    double alpha,int m) {
    const double mu=alpha-state.moment,d=(mu-state.data.location)/state.data.sigma;
    if (state.epsilon==0) {
        const double q=state.q/std::sqrt(double(m)),z=-d*std::sqrt(double(m));
        const double t=q*z,penalty=q==0 ? z*z/2 : (std::expm1(t)-t)/(q*q);
        return std::exp(std::log(std::sqrt(double(m))/state.data.sigma)+lgg_constant(q)-penalty);
    }
    double value=state.data.constant,product=1;
    const double q=state.q,e=state.epsilon,r=state.relative_shift,iq=q==0 ? 0 : 1/q;
    const double location_factor=std::exp(-q*d),shift_factor=std::exp(-q*r),
        slope=q==0 ? -r : location_factor*std::expm1(-q*r)*iq*iq,
        intercept=q==0 ? r*(d+r/2) : r*iq;
    const int sample_count=state.data.standardized.size(),block=std::min(16,std::max(1,
        int(std::floor(std::log(std::numeric_limits<double>::min())/std::log(std::min(e,1-e))))-1));
    for(int first=0;first<sample_count;first+=8) {
        const int count=std::min(8,sample_count-first);
        Eigen::Array<double,8,1> difference;
        double penalty[8],delta[8];
        for(int offset=0;offset<8;++offset) {
            const int j=first+std::min(offset,count-1);
            delta[offset]=intercept+slope*(q==0 ?
                state.data.standardized[j] : state.data.exponential[j]);
            const bool shifted=delta[offset]<0;
            const double z=state.data.standardized[j]-d-(shifted ? r : 0),t=q*z,
                exponential=state.data.exponential[j]*location_factor*(shifted ? shift_factor : 1);
            penalty[offset]=q==0 ? z*z/2 :
                (std::abs(t)<.005 || !std::isfinite(exponential) ?
                std::expm1(t)-t : exponential-1-t)*iq*iq;
            difference[offset]=std::abs(delta[offset]);
        }
        const Eigen::Array<double,8,1> ratio=(-difference).exp();
        for(int offset=0;offset<count;++offset) {
            const int j=first+offset;
            value-=std::max(penalty[offset],0.0);
            product*=delta[offset]>=0 ? 1+e*(ratio[offset]-1) : e+(1-e)*ratio[offset];
            if((j+1)%block==0 || j+1==sample_count){value+=std::log(product);product=1;}
        }
    }
    return std::exp(value-state.log_marginal);
}
static double population_mean_beta_cell(double low,double high,double a,double b) {
    if (!(high>low)) return 0;
    const double F=R::pbeta(high,a,b,true,false);
    return F>.5 ? std::max(0.0,R::pbeta(low,a,b,false,false)-R::pbeta(high,a,b,false,false)) :
        std::max(0.0,F-R::pbeta(low,a,b,true,false));
}
struct PopulationMeanBetaFFT {std::vector<std::complex<double>> frequency;double mass_error;};
// [[Rcpp::export]]
Rcpp::List population_mixture_mean_grid_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,Rcpp::NumericVector epsilon,
    Rcpp::NumericVector relative_shift,Rcpp::NumericVector prior,Rcpp::NumericMatrix positive_loglike,
    Rcpp::NumericVector slab_sd,double spacing=.01,double structural_mass=.5,
    double tail_probability=1e-10,int threads=1,
    Rcpp::Nullable<Rcpp::IntegerVector> pair_left=R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> pair_right=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> pair_weight=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> slab_mean=R_NilValue) {
    const int ng=x.nrow(),ns=shape.size(),nk=slab_sd.size();
    const Rcpp::NumericVector slab_center=slab_mean.isNotNull() ?
        Rcpp::NumericVector(slab_mean) : Rcpp::NumericVector(nk,0.0);
    if (slab_center.size()!=nk) Rcpp::stop("Effect means and scales must have equal lengths");
    const bool joint=pair_weight.isNotNull();
    std::vector<int> left_pair,right_pair;std::vector<double> pair_prior;
    if (joint) {
        const Rcpp::IntegerVector pl(pair_left),pr(pair_right);
        const Rcpp::NumericVector pw(pair_weight);
        left_pair.assign(pl.begin(),pl.end());right_pair.assign(pr.begin(),pr.end());
        pair_prior.assign(pw.begin(),pw.end());
        for (int& v:left_pair) --v;for (int& v:right_pair) --v;
    }
    const int np=pair_prior.size();
    std::vector<std::map<std::tuple<int,int,int>,PopulationMeanBetaFFT>> beta_cache(threads);
    std::vector<std::map<std::tuple<int,double,int>,double>> bandwidth_cache(threads);
    const double logtwo=std::log(2.0),pi=3.14159265358979323846;
    for (int s=0;s<ns;++s) if (prior[s]>0 && shape[s]<0 &&
        !(1/(shape[s]*shape[s])+std::exp(log_sigma[s])/shape[s]>0))
        Rcpp::stop("An arithmetic-mean nuisance state has no finite first moment");
    int totals[2]={0,0};for (int g:group) ++totals[g];
    Rcpp::NumericMatrix out(ng,5),ll(ng,nk),first(ng,nk),second(ng,nk),check(ng,7);
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("effect","se","cdf_zero","p","available");
    Rcpp::colnames(check)=Rcpp::CharacterVector::create("grid_points","positive_mass_error",
        "frequency_mass_error","contrast_mass_error","omitted_state_mass","negative_roundoff","beta_tail");
    std::fill(out.begin(),out.end(),NA_REAL);Rcpp::LogicalVector usable(ng);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for (int gene=0;gene<ng;++gene) {
        std::vector<double> y[2];
        for (int j=0;j<x.ncol();++j) if (x(gene,j)>0) y[group[j]].push_back(std::log(x(gene,j)));
        out(gene,3)=1;out(gene,4)=0;
        if (y[0].empty() || y[1].empty()) continue;
        std::vector<PopulationMeanGridState> states[2];
        double lower=R_PosInf,upper=R_NegInf,beta_lower=0,omitted=0;
        double spike[2],continuous[2];
        int count[2]={int(y[0].size()),int(y[1].size())};
        std::vector<double> marginal[2],pair_probability(np,0);
        marginal[0].assign(ns,0);marginal[1].assign(ns,0);
        if (joint) {
            double logz=R_NegInf;
            for (int p=0;p<np;++p) if (pair_prior[p]>0)
                logz=population_mean_grid_logadd(logz,std::log(pair_prior[p])+
                    positive_loglike(gene,left_pair[p])+positive_loglike(gene+ng,right_pair[p]));
            for (int p=0;p<np;++p) if (pair_prior[p]>0) {
                const double w=std::exp(std::log(pair_prior[p])+positive_loglike(gene,left_pair[p])+
                    positive_loglike(gene+ng,right_pair[p])-logz);
                pair_probability[p]=w;marginal[0][left_pair[p]]+=w;marginal[1][right_pair[p]]+=w;
            }
        }
        for (int g=0;g<2;++g) {
            const int row=gene+g*ng,m=count[g],n=totals[g];
            std::sort(y[g].begin(),y[g].end());double center=0;
            for (double v:y[g]) center+=v/m;
            std::vector<double> centered=y[g];for (double& v:centered) v-=center;
            double logz=R_NegInf;
            for (int s=0;s<ns;++s) if (prior[s]>0)
                logz=population_mean_grid_logadd(logz,std::log(prior[s])+positive_loglike(row,s));
            for (int s=0;s<ns;++s) if (prior[s]>0) {
                const double w=joint ? marginal[g][s] :
                    std::exp(std::log(prior[s])+positive_loglike(row,s)-logz);
                if (w<tail_probability/ns) {omitted+=w;continue;}
                const double sigma=std::exp(log_sigma[s]),q=shape[s],e=epsilon[s],r=relative_shift[s];
                const double moment=log_moment(1,V3{0,log_sigma[s],q})+population_shift_logmean(e,r*sigma);
                const MixturePrepared prepared=prepare_mixture(y[g],centered,center,log_sigma[s],q);
                const double tail=std::min(.1,tail_probability/(ns*w));
                const LocationGrid grid=mixture_location_grid(m,q,tail,1e-9);
                const double state_lower=(prepared.location+sigma*(grid.lower-(e>0 ? r : 0))+moment)/logtwo,
                    state_upper=(prepared.location+sigma*grid.upper+moment)/logtwo;
                lower=std::min(lower,state_lower);upper=std::max(upper,state_upper);
                states[g].push_back({prepared,q,e,r,moment,positive_loglike(row,s),w,state_lower,state_upper,s});
            }
            spike[g]=0;continuous[g]=1;
            if (m==n) {
                const double lpoint=std::log(structural_mass),
                    lbeta=std::log1p(-structural_mass)+R::lbeta(m+.5,n-m+.5)-R::lbeta(.5,.5);
                const double z=population_mean_grid_logadd(lpoint,lbeta);
                spike[g]=std::exp(lpoint-z);continuous[g]=std::exp(lbeta-z);
            }
            if (continuous[g]>0)
                beta_lower=std::min(beta_lower,std::log(R::qbeta(tail_probability,m+.5,n-m+.5,true,false))/logtwo);
        }
        lower+=beta_lower-2*spacing;upper+=2*spacing;
        const double origin=std::floor(lower/spacing)*spacing;
        int ngrid=2;while ((ngrid-1)*spacing<upper-origin) ngrid*=2;
        const int nfft=2*ngrid;check(gene,0)=ngrid;check(gene,4)=omitted;
        std::vector<std::vector<std::complex<double>>> transform[2];
        transform[0].resize(ns);transform[1].resize(ns);
        std::vector<int> spectral_radius[2]={std::vector<int>(ns,0),std::vector<int>(ns,0)};
        std::vector<std::complex<double>> marginal_phi[2],c(nfft,0);
        marginal_phi[0].assign(nfft,0);marginal_phi[1].assign(nfft,0);
        const std::vector<std::complex<double>>* beta_phi[2];
        Eigen::FFT<double> fft;
        int worker=0;
        #ifdef _OPENMP
        worker=omp_get_thread_num();
        #endif
        for (int g=0;g<2;++g) {
            double weighted_mass_error=0;
            for (const auto& state:states[g]) {
                // All latent allocations share the same log-gamma CF modulus.
                // Choose a state-specific sampling rate from its proven bound;
                // zero padding then restores the common contrast frequencies.
                const auto bandwidth_key=std::make_tuple(count[g],state.q,nfft);
                auto bandwidth=bandwidth_cache[worker].find(bandwidth_key);
                if (bandwidth==bandwidth_cache[worker].end()) {
                    // Reserve one quarter of the requested tail budget for the
                    // positive and negative spectral aliases over all frequencies.
                    const double step=mixture_location_grid(count[g],state.q,tail_probability/4,
                        tail_probability/(4*nfft)).spacing;
                    bandwidth=bandwidth_cache[worker].emplace(bandwidth_key,step).first;
                }
                const double allowed=state.data.sigma*bandwidth->second/(2*logtwo),period=nfft*spacing;
                int local_size=2;while (local_size<nfft && period/local_size>allowed) local_size*=2;
                const double local_step=period/local_size;
                spectral_radius[g][state.index]=local_size/2;
                std::vector<double> positive(local_size,0);
                const int left=std::max(0,int(std::floor((state.lower-origin)/local_step))-1),
                    right=std::min(local_size-1,int(std::ceil((state.upper-origin)/local_step))+1);
                long double mass=0;
                for (int j=left;j<=right;++j) {
                    positive[j]=population_mean_grid_density(state,(origin+j*local_step)*logtwo,count[g])*local_step*logtwo;
                    mass+=positive[j];
                }
                weighted_mass_error+=state.weight*std::abs(double(mass)-1);
                for (double& p:positive) p/=mass;
                std::vector<std::complex<double>> small;fft.fwd(small,positive);
                auto& frequency=transform[g][state.index];frequency.assign(nfft,0);
                if (local_size==nfft) frequency=std::move(small);
                else {
                    frequency[0]=small[0];
                    for (int j=1;j<local_size/2;++j) {
                        frequency[j]=small[j];frequency[nfft-j]=small[local_size-j];
                    }
                    const double nyquist=small[local_size/2].real()/2;
                    frequency[local_size/2]=nyquist;frequency[nfft-local_size/2]=nyquist;
                }
                if (!joint) for (int j=0;j<nfft;++j)
                    marginal_phi[g][j]+=state.weight*transform[g][state.index][j];
            }
            check(gene,1)=std::max(check(gene,1),weighted_mass_error);
            const auto cache_key=std::make_tuple(totals[g],count[g],ngrid);
            auto found=beta_cache[worker].find(cache_key);
            if (found==beta_cache[worker].end()) {
                std::vector<double> beta(nfft,0);
                beta[0]=spike[g]+continuous[g]*population_mean_beta_cell(std::exp(-.5*spacing*logtwo),1,
                    count[g]+.5,totals[g]-count[g]+.5);
                double mass=beta[0];
                for (int j=1;j<ngrid;++j) {
                    const double low=std::exp(-(j+.5)*spacing*logtwo),high=std::exp(-(j-.5)*spacing*logtwo);
                    beta[nfft-j]=continuous[g]*population_mean_beta_cell(low,high,count[g]+.5,totals[g]-count[g]+.5);
                    mass+=beta[nfft-j];
                }
                PopulationMeanBetaFFT saved;saved.mass_error=std::abs(1-mass);
                for (double& p:beta) p/=mass;
                fft.fwd(saved.frequency,beta);
                found=beta_cache[worker].emplace(cache_key,std::move(saved)).first;
            }
            beta_phi[g]=&found->second.frequency;
            check(gene,2)=std::max(check(gene,2),found->second.mass_error);
            check(gene,6)=std::max(check(gene,6),found->second.mass_error);
        }
        if (joint) {
            for (int p=0;p<np;++p) if (pair_probability[p]>0 &&
                !transform[0][left_pair[p]].empty() && !transform[1][right_pair[p]].empty()) {
                // Frequency embedding sets every value outside this intersection
                // exactly to zero. Preserve the original order inside the band.
                const int radius=std::min(spectral_radius[0][left_pair[p]],
                    spectral_radius[1][right_pair[p]]);
                for (int j=0;j<=radius;++j)
                    c[j]+=pair_probability[p]*transform[1][right_pair[p]][j]*
                        std::conj(transform[0][left_pair[p]][j]);
                for (int j=std::max(radius+1,nfft-radius);j<nfft;++j)
                    c[j]+=pair_probability[p]*transform[1][right_pair[p]][j]*
                        std::conj(transform[0][left_pair[p]][j]);
            }
        } else for (int j=0;j<nfft;++j) c[j]=marginal_phi[1][j]*std::conj(marginal_phi[0][j]);
        for (int j=0;j<nfft;++j) c[j]*=(*beta_phi[1])[j]*std::conj((*beta_phi[0])[j]);
        std::vector<double> contrast;fft.inv(contrast,c);
        long double mass=0;
        for (double& p:contrast) {check(gene,5)=std::max(check(gene,5),std::max(0.0,-p));p=std::max(0.0,p);mass+=p;}
        check(gene,3)=std::abs(double(mass)-1);
        for (double& p:contrast) p/=mass;
        long double mean=0,square=0,F=.5*contrast[0];
        for (int j=0;j<nfft;++j) {
            const double d=(j<ngrid ? j : j-nfft)*spacing,p=contrast[j];
            mean+=p*d;square+=p*d*d;if (j>=ngrid) F+=p;
        }
        out(gene,0)=mean;out(gene,1)=std::sqrt(std::max(0.0,double(square-mean*mean)));
        out(gene,2)=F;out(gene,3)=2*std::min(double(F),1-double(F));out(gene,4)=1;usable[gene]=true;
        for (int k=0;k<nk;++k) {
            const double tau=slab_sd[k],center=slab_center[k];
            if (tau==0) {
                double density=0;
                if (std::abs(center)<(ngrid-1)*spacing) {
                    const double coordinate=center/spacing;
                    const int low=int(std::floor(coordinate)),j=(low%nfft+nfft)%nfft;
                    const double fraction=coordinate-low;
                    density=((1-fraction)*contrast[j]+fraction*contrast[(j+1)%nfft])/spacing;
                }
                ll(gene,k)=std::log(density);first(gene,k)=center;second(gene,k)=center*center;
                continue;
            }
            long double norm=0,m1=0,m2=0;
            for (int j=0;j<nfft;++j) {
                const double d=(j<ngrid ? j : j-nfft)*spacing,
                    w=contrast[j]*std::exp(-.5*(d-center)*(d-center)/(tau*tau))/(std::sqrt(2*pi)*tau);
                norm+=w;m1+=w*d;m2+=w*d*d;
            }
            ll(gene,k)=std::log(double(norm));first(gene,k)=m1/norm;second(gene,k)=m2/norm;
        }
    }
    return Rcpp::List::create(Rcpp::Named("genes")=out,Rcpp::Named("log_likelihood")=ll,
        Rcpp::Named("slab_first")=first,Rcpp::Named("slab_second")=second,
        Rcpp::Named("usable")=usable,Rcpp::Named("checks")=check,
        Rcpp::Named("spacing")=spacing,
        Rcpp::Named("nuisance_structure")=joint ? "paired" : "independent");
}

// Joint mean inference integrates b(theta), location, nonzero frequency and
// both class-specific populations. Mehler/Stein moments share the likelihood.
// The fixed support is an empirical approximation; score calibration requires
// the exchangeable cross-gene population model, not uniform nuisance validity.
static std::complex<double> mean_reference_loggamma(std::complex<double> z) {
    static const double coefficient[9]={
        .99999999999980993,676.5203681218851,-1259.1392167224028,
        771.32342877765313,-176.61502916214059,12.507343278686905,
        -.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7};
    const bool recur=z.real()<.5;
    const std::complex<double> initial=z;
    if (recur) z+=1;
    z-=1;
    std::complex<double> sum=coefficient[0];
    for (int k=1;k<9;++k) sum+=coefficient[k]/(z+double(k));
    const std::complex<double> t=z+7.5;
    std::complex<double> result=.91893853320467274178+(z+.5)*std::log(t)-t+std::log(sum);
    if (recur) result-=std::log(initial);
    return result;
}



struct MeanFullGroupCurve {
    std::vector<std::complex<double>> phi,gaussian_phi;
    std::vector<double> density,coefficient,null_smoothed,point_density;
    std::vector<double> gaussian_density,gaussian_coefficient,gaussian_point_density;
    std::vector<double> conditional,conditional_first,conditional_second,gaussian_conditional;
    double maximum=0;
};
static std::vector<double> mean_full_fold_gaussian(
    const std::vector<std::complex<double>>& phi,int nfft,int stride,double spacing,double origin,
    const Rcpp::NumericVector& A,double sd,SGCBFFT& fft) {
    const int nf=nfft/stride;
    const double pi=3.14159265358979323846,period=nfft*spacing,step=stride*spacing,logtwo=std::log(2.0);
    const double translation=A[0]/logtwo-origin-std::round((A[0]/logtwo-origin)/step)*step;
    std::vector<std::complex<double>> folded(nf,0);
    for (int j=0;j<int(phi.size());++j) {
        const double omega=2*pi*j/period;
        const double attenuation=std::exp(-.5*sd*sd*omega*omega);
        if (attenuation==0) break; // The remaining monotone Gaussian multipliers are exactly floating zero.
        const auto value=phi[j]*attenuation*std::exp(std::complex<double>(0,omega*translation));
        folded[j%nf]+=value;
        if (j>0 && j<nfft/2) folded[(nf-j%nf)%nf]+=std::conj(value);
    }
    std::vector<double> out;fft.inv(out,folded);
    std::vector<double> result(A.size());
    for (int a=0;a<A.size();++a) {
        const long long at=std::llround((A[a]/logtwo-origin)/step);
        result[a]=out[int((at%nf+nf)%nf)]/step;
    }
    return result;
}
static std::vector<double> mean_full_mehler_coefficients(
    const std::vector<std::complex<double>>& phi,int nfft,int stride,double spacing,double origin,
    const Rcpp::NumericVector& A,double sd,int degree,SGCBFFT& fft) {
    const int nf=nfft/stride,nA=A.size();
    const double pi=3.14159265358979323846,period=nfft*spacing,step=stride*spacing,logtwo=std::log(2.0);
    const double translation=A[0]/logtwo-origin-std::round((A[0]/logtwo-origin)/step)*step;
    std::vector<std::complex<double>> current;current.reserve(phi.size());
    for (int j=0;j<int(phi.size());++j) {
        const double omega=2*pi*j/period;
        const double attenuation=std::exp(-.5*sd*sd*omega*omega);
        if (attenuation==0) break; // All subsequent coefficients remain exactly zero at every Hermite order.
        current.push_back(phi[j]*attenuation*std::exp(std::complex<double>(0,omega*translation)));
    }
    std::vector<double> result(size_t(degree+1)*nA);
    for (int k=0;k<=degree;++k) {
        std::vector<std::complex<double>> folded(nf,0);
        for (int j=0;j<int(current.size());++j) {
            folded[j%nf]+=current[j];
            if (j>0 && j<nfft/2) folded[(nf-j%nf)%nf]+=std::conj(current[j]);
        }
        std::vector<double> values;fft.inv(values,folded);
        for (int a=0;a<nA;++a) {
            const long long at=std::llround((A[a]/logtwo-origin)/step);
            result[a+size_t(nA)*k]=values[int((at%nf+nf)%nf)]/step;
        }
        const double denominator=std::sqrt(double(k+1));
        for (int j=0;j<int(current.size());++j)
            current[j]*=std::complex<double>(0,sd*2*pi*j/period/denominator);
    }
    return result;
}
static std::vector<double> mean_full_mehler_coefficients_reused(
    const std::vector<std::complex<double>>& phi,int nfft,int stride,double spacing,double origin,
    const Rcpp::NumericVector& A,double sd,int degree,SGCBFFT& fft) {
    const int nf=nfft/stride,nA=A.size();
    const double pi=3.14159265358979323846,period=nfft*spacing,step=stride*spacing,logtwo=std::log(2.0);
    const double translation=A[0]/logtwo-origin-std::round((A[0]/logtwo-origin)/step)*step;
    std::vector<std::complex<double>> current;current.reserve(phi.size());
    for (int j=0;j<int(phi.size());++j) {
        const double omega=2*pi*j/period;
        const double attenuation=std::exp(-.5*sd*sd*omega*omega);
        if (attenuation==0) break; // All subsequent coefficients remain exactly zero at every Hermite order.
        current.push_back(phi[j]*attenuation*std::exp(std::complex<double>(0,omega*translation)));
    }
    std::vector<double> result(size_t(degree+1)*nA);
    std::vector<std::complex<double>> folded(nf);
    std::vector<double> values;
    for (int k=0;k<=degree;++k) {
        std::fill(folded.begin(),folded.end(),std::complex<double>(0,0));
        for (int j=0;j<int(current.size());++j) {
            folded[j%nf]+=current[j];
            if (j>0 && j<nfft/2) folded[(nf-j%nf)%nf]+=std::conj(current[j]);
        }
        fft.inv(values,folded);
        for (int a=0;a<nA;++a) {
            const long long at=std::llround((A[a]/logtwo-origin)/step);
            result[a+size_t(nA)*k]=values[int((at%nf+nf)%nf)]/step;
        }
        const double denominator=std::sqrt(double(k+1));
        for (int j=0;j<int(current.size());++j)
            current[j]*=std::complex<double>(0,sd*2*pi*j/period/denominator);
    }
    return result;
}
// Conditional Gaussian factorization: U=A-m/2+T+e0,
// V=A+m/2+sign(rho)*T+e1, Var(e0)=Var(e1)=s^2*(1-|rho|).
// Each integrand is a positive mixture of Gaussians in T with common
// variance (1/Var(T)+2/Var(e))^-1. Poisson summation bounds the lattice
// error uniformly over the mixture centers; no interpolation is used.
static std::vector<double> mean_conditional_fold(
    const std::vector<std::complex<double>>& phi,int nfft,int fold,
    double spacing,double origin,const Rcpp::NumericVector& A,
    double sd,int derivative,SGCBFFT& fft,int max_frequency=std::numeric_limits<int>::max()) {
    const int nf=nfft/fold;
    const double pi=3.14159265358979323846,period=nfft*spacing,
        step=fold*spacing,logtwo=std::log(2.0);
    const long long first=std::llround((A[0]/logtwo-origin)/step);
    const double translation=A[0]/logtwo-origin-first*step;
    std::vector<std::complex<double>> folded(nf,0);
    for (int j=0;j<int(phi.size()) && j<=max_frequency;++j) {
        const double omega=2*pi*j/period;
        const double attenuation=std::exp(-.5*sd*sd*omega*omega);
        if (attenuation==0) break;
        auto value=phi[j]*attenuation*std::exp(std::complex<double>(0,omega*translation));
        if (derivative==1) value*=std::complex<double>(0,omega);
        if (derivative==2) value*=-omega*omega;
        folded[j%nf]+=value;
        if (j>0 && j<nfft/2) folded[(nf-j%nf)%nf]+=std::conj(value);
    }
    std::vector<double> values;fft.inv(values,folded);
    std::vector<double> result(nf);
    for (int j=0;j<nf;++j) {
        const int at=int(((first+j)%nf+nf)%nf);
        const double value=values[at]/step;
        result[j]=derivative==0 ? std::max(0.0,value) : value;
    }
    return result;
}
static std::vector<double> mean_full_null_curve(const std::vector<double>& a,
    const std::vector<double>& b,int nfft,int stride,double spacing,double origin,
    const Rcpp::NumericVector& A,double sd,SGCBFFT& fft) {
    std::vector<double> product(nfft,0);
    for (int j=0;j<int(a.size());++j) product[j]=a[j]*b[j]*spacing;
    std::vector<std::complex<double>> transform;fft.fwd(transform,product);
    transform.resize(nfft/2+1);
    return mean_full_fold_gaussian(transform,nfft,stride,spacing,origin,A,sd,fft);
}
static double mean_full_gaussian_pair_density(double A,double mean0,double mean1,
    double variance0,double variance1,double a2,double d2) {
    const double u=variance0+a2+d2/4,v=variance1+a2+d2/4,c=a2-d2/4;
    const double determinant=variance0*variance1+(a2+d2/4)*(variance0+variance1)+a2*d2;
    const double x=A-mean0,y=A-mean1;
    const double quadratic=variance1*x*x+variance0*y*y+a2*(x-y)*(x-y)+d2/4*(x+y)*(x+y);
    return std::exp(-.5*quadratic/determinant)/
        (2*3.14159265358979323846*std::sqrt(determinant));
}
// Private profiling dependency.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector mean_full_beta_transform_controls_cpp() {
    const double omega[7]={0,.1,1,3,10,30,100};
    double error=0;
    for (double w:omega) {
        const std::complex<double> z(1,-w);
        const auto phi=std::exp(mean_reference_loggamma(z)-mean_reference_loggamma(z+1.0));
        error=std::max(error,std::abs(phi-1.0/z));
    }
    return Rcpp::NumericVector::create(Rcpp::Named("uniform_beta_CF_identity")=error);
}
// Private profiling dependency.

// Stable log-gamma increments, DLMF 5.11.1 and 5.11(ii).
// In Re(z)>0, eight Stirling terms beyond |z|=16 leave a
// first-omitted-term bound below 3.2e-19 (including the sector factor).
using MeanCLD=std::complex<long double>;
static MeanCLD mean_spectral_log1p(MeanCLD w) {
    const long double u=w.real(),v=w.imag();
    return MeanCLD(.5L*std::log1p(2*u+u*u+v*v),std::atan2(v,1+u));
}
static MeanCLD mean_spectral_stirling(MeanCLD z) {
    const MeanCLD u=1.L/z,u2=u*u;
    return u*(1.L/12+u2*(-1.L/360+u2*(1.L/1260+u2*(-1.L/1680+
        u2*(1.L/1188+u2*(-691.L/360360+u2*(1.L/156-u2*3617.L/122400.L)))))));
}
static MeanCLD mean_spectral_gamma_increment(MeanCLD z,long double d) {
    int shift=0;while (std::abs(z+static_cast<long double>(shift))<16) ++shift;
    const MeanCLD w=z+static_cast<long double>(shift);
    MeanCLD out=d*std::log(w)+(w+d-.5L)*mean_spectral_log1p(d/w)-d+
        mean_spectral_stirling(w+d)-mean_spectral_stirling(w);
    for (int k=0;k<shift;++k) out-=mean_spectral_log1p(d/(z+static_cast<long double>(k)));
    return out;
}
static MeanCLD mean_spectral_centered_gamma(long double a,long double t) {
    const int shift=std::max(0,int(std::ceil(16-a)));
    const long double z=a+shift;
    const MeanCLD d(0,t),w(z,0);
    MeanCLD out=(w+d-.5L)*mean_spectral_log1p(d/w)-d+
        mean_spectral_stirling(w+d)-mean_spectral_stirling(w);
    if (shift) {
        out+=d*std::log1p(shift/a);
        for (int k=0;k<shift;++k) out-=mean_spectral_log1p(d/(a+k));
    }
    return out;
}
struct MeanSpectralCache {
    std::mutex location_mutex,beta_mutex;
    std::map<std::tuple<int,double,double,int,int>,
        std::shared_ptr<std::vector<std::complex<double>>>> location;
    std::map<std::tuple<int,double,double,int,int>,int> radius;
    std::map<std::tuple<int,int,int,int>,
        std::shared_ptr<std::vector<std::complex<double>>>> beta;
};
// For a translated log-gamma allocation mixture, |CF| is bounded by
// the common gamma CF. DLMF 5.7.6 implies Im psi(a+it)>=atan(t/a).
// Its monotone exponential envelope bounds the entire discrete Fourier tail.
static int mean_spectral_location_radius(MeanSpectralCache& cache,int m,double sigma,
    double q,int base_nfft,int maximum_radius,double period,double budget) {
    const auto key=std::make_tuple(m,sigma,q,base_nfft,maximum_radius);
    std::lock_guard<std::mutex> lock(cache.location_mutex);
    auto found=cache.radius.find(key);
    if (found!=cache.radius.end()) return found->second;
    const long double pi=3.141592653589793238462643383279502884L,
        delta=2*pi/(period*std::log(2.L)),target=std::log(budget);
    const auto log_tail=[&](int radius) {
        const long double j=radius+1,frequency=delta*j;
        long double log_modulus,rate;
        if (q==0) {
            const long double variance=static_cast<long double>(sigma)*sigma/m;
            log_modulus=-variance*frequency*frequency/2;
            rate=variance*delta*frequency;
        } else {
            const long double a=static_cast<long double>(m)/(static_cast<long double>(q)*q),
                scale=static_cast<long double>(sigma)/std::abs(q),t=scale*frequency;
            log_modulus=mean_spectral_centered_gamma(a,t).real();
            rate=scale*delta*std::atan(t/a);
        }
        return std::log(2.L/period)+log_modulus-std::log(-std::expm1(-rate));
    };
    int low=0,high=maximum_radius;
    while (high-low>1) {
        const int middle=low+(high-low)/2;
        if (log_tail(middle)<=target) high=middle;else low=middle;
    }
    cache.radius.emplace(key,high);
    return high;
}

static std::shared_ptr<std::vector<std::complex<double>>> mean_spectral_location_base(
    MeanSpectralCache& cache,int m,double sigma,double q,int base_nfft,int radius,double period) {
    const auto key=std::make_tuple(m,sigma,q,base_nfft,radius);
    std::lock_guard<std::mutex> lock(cache.location_mutex);
    auto found=cache.location.find(key);
    if (found!=cache.location.end()) return found->second;
    auto values=std::make_shared<std::vector<std::complex<double>>>(radius+1);
    (*values)[0]=1;
    const long double pi=3.141592653589793238462643383279502884L,logtwo=std::log(2.L);
    for (int j=1;j<=radius;++j) {
        const long double omega=2*pi*j/period;
        if (q==0) (*values)[j]=std::exp(-.5L*omega*omega*sigma*sigma/(m*logtwo*logtwo));
        else {
            const auto value=std::exp(mean_spectral_centered_gamma(
                static_cast<long double>(m)/(static_cast<long double>(q)*q),omega*sigma/(q*logtwo)));
            (*values)[j]={double(value.real()),double(value.imag())};
        }
    }
    cache.location.emplace(key,values);
    return values;
}
static std::shared_ptr<std::vector<std::complex<double>>> mean_spectral_beta(
    MeanSpectralCache& cache,int n,int m,int base_nfft,int radius,double period,
    double spike,double continuous) {
    const auto key=std::make_tuple(n,m,base_nfft,radius);
    std::lock_guard<std::mutex> lock(cache.beta_mutex);
    auto found=cache.beta.find(key);
    if (found!=cache.beta.end()) return found->second;
    auto values=std::make_shared<std::vector<std::complex<double>>>(radius+1);
    (*values)[0]=1;
    const long double a=m+.5L,d=n-m+.5L,
        pi=3.141592653589793238462643383279502884L,logtwo=std::log(2.L);
    const MeanCLD normalizer=mean_spectral_gamma_increment(MeanCLD(a,0),d);
    for (int j=1;j<=radius;++j) {
        const long double omega=2*pi*j/(period*logtwo);
        const auto ratio=std::exp(normalizer-mean_spectral_gamma_increment(MeanCLD(a,-omega),d));
        (*values)[j]=spike+continuous*std::complex<double>(double(ratio.real()),double(ratio.imag()));
    }
    cache.beta.emplace(key,values);
    return values;
}
// Pure GG locations, and every latent allocation for at most five positives,
// are translated log-gamma laws. This evaluates their CF without sampling
// the narrow positive-density peaks. No allocation probability is truncated.
static std::vector<std::complex<double>> mean_spectral_location_phi(
    const PopulationMeanGridState& state,int m,int radius,int base_nfft,
    double period,double origin,double tolerance,MeanSpectralCache& cache,double& mass_error) {
    const auto& data=state.data;
    const long double sigma=data.sigma,q=state.q,r=state.relative_shift;
    long double total_z=0,total_expm1=0,total_exp=0;
    for (int j=0;j<m;++j) {
        total_z+=data.standardized[j];
        total_expm1+=std::expm1(q*data.standardized[j]);
        total_exp+=data.exponential[j];
    }
    const long double mode=data.location+sigma*(q==0 ? total_z/m :
        std::log1p(total_expm1/m)/q);
    std::map<long double,long double> nodes;
    long double log_normalizer=data.pure_log_likelihood;
    if (state.epsilon==0) nodes[0]=1;
    else {
        const int n=1<<m,full=n-1;
        const long double a=q==0 ? 1 : m/(q*q),t=-q*r,p=state.epsilon;
        std::vector<int> count(n,0);
        std::vector<long double> sx(n,0),sm(n,0),sz(n,0),ratio(n,0),shift(n,0),score(n);
        for (int mask=1;mask<n;++mask) {
            const int previous=mask&(mask-1),bit=mask^previous;
            int j=0;for (int value=bit;value>1;value>>=1) ++j;
            count[mask]=count[previous]+1;
            sx[mask]=sx[previous]+data.exponential[j];
            sm[mask]=sm[previous]+std::expm1(q*data.standardized[j]);
            sz[mask]=sz[previous]+data.standardized[j];
        }
        shift[full]=-sigma*r;
        for (int mask=1;mask<full;++mask) {
            const long double T=count[mask],fraction=T/m;
            if (q==0) {
                ratio[mask]=r*(sz[mask]-fraction*total_z)-r*r*T*(1-fraction)/2;
                shift[mask]=-sigma*r*fraction;
            } else {
                long double logR;
                if (std::abs(t)<.01L) {
                    const long double part=sx[mask]/total_exp,v=part*std::expm1(t),
                        dev=(m*sm[mask]-T*total_expm1)/(m+total_expm1);
                    logR=std::log1p(v);
                    ratio[mask]=r*dev/q-a*(allocation_log1pmx(v)+part*allocation_expm1mx(t));
                } else {
                    const long double left=std::log(sx[mask])+t,right=std::log(sx[full^mask]),
                        peak=std::max(left,right);
                    logR=peak+std::log1p(std::exp(-std::abs(left-right)))-std::log(total_exp);
                    ratio[mask]=-r*T/q-a*logR;
                }
                shift[mask]=sigma/q*logR;
            }
        }
        long double peak=-std::numeric_limits<long double>::infinity();
        for (int mask=0;mask<n;++mask) {
            const int T=count[mask];
            score[mask]=ratio[mask]+(T==0 ? 0 : T*std::log(p))+
                (T==m ? 0 : (m-T)*std::log1p(-p));
            peak=std::max(peak,score[mask]);
        }
        long double normalizer=0;
        for (long double value:score) normalizer+=std::exp(value-peak);
        for (int mask=0;mask<n;++mask) nodes[shift[mask]]+=std::exp(score[mask]-peak)/normalizer;
        log_normalizer+=peak+std::log(normalizer);
    }
    mass_error=std::abs(std::expm1(double(log_normalizer)-state.log_marginal));
    const auto base=mean_spectral_location_base(cache,m,data.sigma,state.q,base_nfft,radius,period);
    std::vector<std::complex<double>> result(radius+1,0);
    const long double pi=3.141592653589793238462643383279502884L,logtwo=std::log(2.L);
    // Refresh phase recurrences within 1/64 of the specified error budget.
    const int phase_block=std::max(1,std::min(1024,int(tolerance/
        (64*std::numeric_limits<double>::epsilon()))));
    for (const auto& node:nodes) {
        const long double phase=-2*pi/period*((mode+node.first+state.moment)/logtwo-origin);
        const MeanCLD increment=std::exp(MeanCLD(0,phase));
        MeanCLD current(0,0);
        for (int j=0;j<=radius;++j) {
            if (j%phase_block==0) current=node.second*std::exp(MeanCLD(0,phase*j));
            result[j]+=std::complex<double>(double(current.real()),double(current.imag()));
            current*=increment;
        }
    }
    for (int j=0;j<=radius;++j) result[j]*=(*base)[j];
    result[0]=1;
    return result;
}
static int mean_spectral_gaussian_limit(double period,double sd,double tolerance,int derivative) {
    const double pi=3.14159265358979323846,delta=2*pi/period;
    const auto tail=[&](int j) {
        const double x=sd*delta*j,e=std::exp(-x*x/2),erfc=std::erfc(x/std::sqrt(2.0));
        if (derivative==0) return erfc/(sd*std::sqrt(2*pi));
        if (derivative==1) return e/(pi*sd*sd);
        return (x*e+std::sqrt(pi/2)*erfc)/(pi*sd*sd*sd);
    };
    int low=int(std::ceil(std::sqrt(double(derivative))/(sd*delta))),high=std::max(1,low);
    while (tail(high)>tolerance) high*=2;
    while (high-low>1) {
        const int middle=low+(high-low)/2;
        if (tail(middle)>tolerance) low=middle;else high=middle;
    }
    return high;
}
static std::vector<std::complex<double>> mean_spectral_prefix(
    const std::vector<std::complex<double>>& phi,int limit) {
    return std::vector<std::complex<double>>(phi.begin(),phi.begin()+std::min(int(phi.size()),limit+1));
}
// Only low common frequencies survive the Gaussian common-level prior.
// For |p|<=J, convolution phi0(q)*phi1(p-q) cannot use |q|>r1+J.
// The two input spectra can therefore be shortened separately; zero padding
// to N>r0+r1+J prevents every circular alias in those retained frequencies.
static std::vector<double> mean_spectral_null(
    const std::vector<std::complex<double>>& phi0,const std::vector<std::complex<double>>& phi1,
    int base_nfft,int base_stride,double base_spacing,double origin,
    const Rcpp::NumericVector& A,double sd,int limit,SGCBFFT& fft,
    double shift0=0,double shift1=0) {
    const int r0=std::min(int(phi0.size())-1,int(phi1.size())-1+limit),
        r1=std::min(int(phi1.size())-1,int(phi0.size())-1+limit);
    const double period=base_nfft*base_spacing,pi=3.14159265358979323846;
    int nfft=base_nfft/base_stride;while (nfft<=r0+r1+limit) nfft*=2;
    const int stride=nfft/(base_nfft/base_stride);
    const double spacing=period/nfft;
    std::vector<std::complex<double>> full(nfft,0);
    for (int j=0;j<=r0;++j) {
        const auto value=phi0[j]*std::exp(std::complex<double>(0,-2*pi*j*shift0/period));
        full[j]=value;if (j>0) full[nfft-j]=std::conj(value);
    }
    std::vector<double> a,b;fft.inv(a,full);std::fill(full.begin(),full.end(),0);
    for (int j=0;j<=r1;++j) {
        const auto value=phi1[j]*std::exp(std::complex<double>(0,-2*pi*j*shift1/period));
        full[j]=value;if (j>0) full[nfft-j]=std::conj(value);
    }
    fft.inv(b,full);
    for (int j=0;j<nfft;++j) a[j]*=b[j]/spacing;
    fft.fwd(full,a);full.resize(std::min(limit,nfft/2)+1);
    return mean_full_fold_gaussian(full,nfft,stride,spacing,origin,A,sd,fft);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List mean_full_two_group_dictionary_cpp(Rcpp::NumericMatrix x,Rcpp::IntegerVector group,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,Rcpp::NumericVector epsilon,
    Rcpp::NumericVector relative_shift,Rcpp::NumericVector bvalue,
    Rcpp::NumericMatrix positive_loglike,Rcpp::IntegerVector pair_left,Rcpp::IntegerVector pair_right,
    Rcpp::NumericVector A,double A_step,double A_sd,double slab_sd,int stride,
    double tolerance,double tail_probability,double structural_mass,int threads,
    Rcpp::Nullable<Rcpp::IntegerMatrix> posterior_labels=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> reference_likelihood=R_NilValue,double slab_mean=0,
    Rcpp::Nullable<Rcpp::NumericVector> posterior_weights=R_NilValue,
    bool selected_likelihood=false) {
    constexpr bool reuse=true;
    const int ng=x.nrow(),ns=shape.size(),np=pair_left.size(),nA=A.size(),nk=np*nA;
    const double logtwo=std::log(2.0),pi=3.14159265358979323846,
        a_sd=A_sd/logtwo,s_sd=std::sqrt(a_sd*a_sd+slab_sd*slab_sd/4),
        rho=(a_sd*a_sd-slab_sd*slab_sd/4)/(s_sd*s_sd),
        spacing=A_step/(stride*logtwo),qstep=stride*spacing;
    int totals[2]={0,0};for (int g:group) ++totals[g];
    double minimum_sigma=R_PosInf;
    for (double v:log_sigma) minimum_sigma=std::min(minimum_sigma,std::exp(v));
    const double maximum_density=logtwo*std::sqrt(double(std::max(totals[0],totals[1])))/
        (minimum_sigma*std::sqrt(2*pi));
    const bool weighted=posterior_weights.isNotNull(),selected=selected_likelihood;
    if(selected && posterior_labels.isNull())Rcpp::stop("Selected likelihoods require state labels");
    const bool effects=!selected && (posterior_labels.isNotNull() || weighted),subset=effects || selected;
    const Rcpp::NumericVector mixture_weights=weighted ? Rcpp::NumericVector(posterior_weights) : Rcpp::NumericVector(0);
    const Rcpp::IntegerMatrix labels=posterior_labels.isNotNull() ? Rcpp::IntegerMatrix(posterior_labels) : Rcpp::IntegerMatrix(0,0);
    const Rcpp::NumericMatrix reference=effects ? Rcpp::NumericMatrix(reference_likelihood) : Rcpp::NumericMatrix(0,0);
    Rcpp::NumericVector slab_A0=Rcpp::clone(A),slab_A1=Rcpp::clone(A);
    for (int a=0;a<nA;++a) {
        slab_A0[a]-=slab_mean*logtwo/2;slab_A1[a]+=slab_mean*logtwo/2;
    }
    const bool point=slab_sd==0;
    const double t2=slab_sd*slab_sd,q=std::abs(rho);
    int degree=point || q==0 ? 0 : std::max(0,int(std::ceil(
        std::log(tolerance/(maximum_density*maximum_density))/std::log(q)))-1);
    if (effects && !point && q>0) {
        while (maximum_density*maximum_density*std::pow(q,degree+1)*
            std::max({1.0,std::abs(slab_mean)+t2/s_sd*std::sqrt(double(degree+2)),
                slab_mean*slab_mean+t2+2*std::abs(slab_mean)*t2/s_sd*std::sqrt(double(degree+2))+t2*t2/(s_sd*s_sd)*(degree+3)})>tolerance) ++degree;
    }
    // Only a single observed group uses the marginal coefficients.
    const int coefficient_degree=effects && !point ? 2 : 0;
    degree=0;
    const double conditional_variance=point ? 0 : 2*std::min(a_sd*a_sd,t2/4),
        conditional_sd=std::sqrt(conditional_variance),
        latent_variance=std::abs(a_sd*a_sd-t2/4),
        latent_sd=std::sqrt(latent_variance),
        second_prior=slab_mean*slab_mean+t2,
        fourth_prior=std::pow(slab_mean,4)+6*slab_mean*slab_mean*t2+3*t2*t2,
        moment_factor=std::max({1.0,std::sqrt(second_prior),second_prior}),
        tail_factor=std::max({1.0,std::sqrt(second_prior),std::sqrt(fourth_prior)}),
        density_bound=maximum_density*maximum_density;
    const int conditional_sign=rho<0 ? -1 : 1;
    int conditional_stride=stride,conditional_half=0;
    double conditional_step=qstep,conditional_bound=0,conditional_limit=0;
    std::vector<double> conditional_weight(1,1);
    if (!point && latent_variance>0) {
        const double product_variance=1/(1/latent_variance+2/conditional_variance);
        for (;;) {
            conditional_step=conditional_stride*spacing;
            const double alias=std::exp(-2*pi*pi*product_variance/(conditional_step*conditional_step));
            conditional_bound=density_bound*moment_factor*2*alias/(1-alias);
            if (conditional_bound<=tolerance/16) break;
            conditional_stride/=2;
            if (conditional_stride==0) Rcpp::stop("Conditional lattice exceeds the positive-density resolution");
        }
        const double tail_mass=std::pow(tolerance/(16*std::max(1.0,density_bound)*tail_factor),2);
        conditional_limit=latent_sd*R::qnorm(tail_mass/2,0,1,false,false);
        conditional_half=int(std::ceil(conditional_limit/conditional_step));
        conditional_limit=conditional_half*conditional_step;
        conditional_weight.resize(2*conditional_half+1);
        for (int j=-conditional_half;j<=conditional_half;++j) {
            const double t=j*conditional_step;
            conditional_weight[j+conditional_half]=conditional_step*
                std::exp(-t*t/(2*latent_variance))/std::sqrt(2*pi*latent_variance);
        }
        conditional_bound+=tolerance/16;
    }
    // DLMF18.14.9 gives a uniform Gaussian-Hermite spatial envelope.
    const double spatial_budget=tolerance*(point ? 1 : 1-std::abs(rho))/(16*(1+maximum_density)),
        full_margin=std::max(std::abs(slab_mean)/2+s_sd*std::sqrt(4*std::log(1/spatial_budget)),
            point ? 0.0 : std::abs(slab_mean)/2+conditional_limit+conditional_sd*
                R::qnorm(tolerance/(64*(1+maximum_density)),0,1,false,false)),
        null_margin=a_sd*R::qnorm(tolerance/(16*(1+maximum_density*maximum_density)),0,1,false,false);
    const double window_origin=std::floor((A[0]/logtwo-null_margin)/qstep)*qstep;
    const int window_size=int(std::ceil((A[nA-1]/logtwo+null_margin-window_origin)/spacing))+1;
    std::vector<int> used[2]={std::vector<int>(ns,0),std::vector<int>(ns,0)};
    for (int p=0;p<np;++p) {used[0][pair_left[p]-1]=1;used[1][pair_right[p]-1]=1;}
    std::vector<int> all_pairs(np);std::iota(all_pairs.begin(),all_pairs.end(),0);
    Rcpp::NumericMatrix likelihood(subset ? 0 : ng,subset ? 0 : 2*nk),checks(ng,9),effect_moments(ng,5),
        selected_values(selected ? ng : 0,selected ? 2*labels.nrow() : 0);
    Rcpp::NumericVector offset(ng); Rcpp::IntegerVector resolution_factor(ng); Rcpp::NumericMatrix gaussian_relative(ng,2);
    Rcpp::colnames(checks)=Rcpp::CharacterVector::create("positive_mass_error","eta_CF_zero_error",
        "eta_CF_modulus_excess","conditional_quadrature_bound","negative_joint_density",
        "Gaussian_null_identity","Gaussian_slab_identity","grid_points","Gaussian_control_pairs");
    MeanSpectralCache spectral_cache;
    std::vector<std::map<std::tuple<int,double,int>,double>> bandwidth_cache(threads);
    std::vector<SGCBFFT> fft_cache(reuse ? threads : 0);
    Rcpp::NumericMatrix performance(ng,6);
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for (int gene=0;gene<ng;++gene) {
        const auto gene_started=std::chrono::steady_clock::now();
        double mehler_seconds=0,null_seconds=0;
        std::vector<double> y[2];double center[2]={0,0};
        for (int j=0;j<x.ncol();++j) if (x(gene,j)>0) y[group[j]].push_back(std::log(x(gene,j)));
        int m[2]={int(y[0].size()),int(y[1].size())},worker=0;
        #ifdef _OPENMP
        worker=omp_get_thread_num();
        #endif
        double peak=R_NegInf,frequency_logmass=0,lower=window_origin-full_margin,
            upper=window_origin+(window_size-1)*spacing+full_margin;
        for (int p=0;p<np;++p) peak=std::max(peak,
            (m[0]>0 ? positive_loglike(gene,pair_left[p]-1) : 0)+
            (m[1]>0 ? positive_loglike(gene+ng,pair_right[p]-1) : 0));
        std::vector<PopulationMeanGridState> state[2];
        std::vector<double> pair_probability(np);
        for (int p=0;p<np;++p) pair_probability[p]=std::exp(
            (m[0]>0 ? positive_loglike(gene,pair_left[p]-1) : 0)+
            (m[1]>0 ? positive_loglike(gene+ng,pair_right[p]-1) : 0)-peak);
        std::vector<double> requested(effects ? nk : 0,0);
        std::vector<int> selected_states,active_pairs=subset?std::vector<int>():all_pairs;
        std::vector<int> gene_used[2]={used[0],used[1]};
        if(subset) {
            std::fill(gene_used[0].begin(),gene_used[0].end(),0);
            std::fill(gene_used[1].begin(),gene_used[1].end(),0);
            if(selected) {
                for(int d=0;d<labels.nrow();++d) {
                    const int z=labels(d,gene)-1;
                    if(z<0)continue;
                    const int state=z%nk,p=state/nA;
                    selected_states.push_back(state);active_pairs.push_back(p);
                    gene_used[0][pair_left[p]-1]=1;gene_used[1][pair_right[p]-1]=1;
                }
                std::sort(selected_states.begin(),selected_states.end());
            } else if(weighted) {
                long double density=0;
                for(int z=0;z<2*nk;++z)density+=mixture_weights[z]*reference(gene,z);
                for(int z=nk;z<2*nk;++z)if(mixture_weights[z]>0 && reference(gene,z)>0) {
                    requested[z-nk]=mixture_weights[z]*reference(gene,z)/density;
                    const int p=(z-nk)/nA;
                    active_pairs.push_back(p);
                    gene_used[0][pair_left[p]-1]=1;gene_used[1][pair_right[p]-1]=1;
                }
            } else for(int d=0;d<labels.nrow();++d) {
                const int z=labels(d,gene)-1;
                if(z>=nk) {
                    const int state=z%nk,p=state/nA;
                    requested[state]+=1.0/labels.nrow();active_pairs.push_back(p);
                    gene_used[0][pair_left[p]-1]=1;gene_used[1][pair_right[p]-1]=1;
                }
            }
            std::sort(active_pairs.begin(),active_pairs.end());
            active_pairs.erase(std::unique(active_pairs.begin(),active_pairs.end()),active_pairs.end());
        }
        double continuous[2]={0,0},spike[2]={0,0};
        for (int g=0;g<2;++g) {
            const double lc=std::log1p(-structural_mass)+R::lbeta(m[g]+.5,totals[g]-m[g]+.5)-R::lbeta(.5,.5),
                ls=m[g]==totals[g] ? std::log(structural_mass) : R_NegInf,
                mass=population_mean_grid_logadd(lc,ls);
            frequency_logmass+=mass;continuous[g]=std::exp(lc-mass);spike[g]=std::exp(ls-mass);
            state[g].resize(ns);
            if (m[g]==0) continue;
            std::sort(y[g].begin(),y[g].end());
            for (double v:y[g]) center[g]+=v/m[g];
            std::vector<double> centered=y[g];for (double& v:centered) v-=center[g];
            const double beta_lower=std::log(R::qbeta(tail_probability/(8*ns),
                m[g]+.5,totals[g]-m[g]+.5,true,false))/logtwo;
            for (int s=0;s<ns;++s) if (gene_used[g][s]) {
                auto& st=state[g][s];
                st.data=prepare_mixture(y[g],centered,center[g],log_sigma[s],shape[s]);
                st.q=shape[s];st.epsilon=epsilon[s];st.relative_shift=relative_shift[s];
                st.moment=bvalue[s];st.log_marginal=positive_loglike(gene+g*ng,s);st.index=s;
                if (st.epsilon==0 || m[g]<=mixture_allocation_maximum)
                    st.data.pure_log_likelihood=m[g]>=2 ? location_integral_constant(m[g],st.q)+
                        marginal_scale_kernel(centered,log_sigma[s],st.q).log_density : 0;
                const auto grid=mixture_location_grid(m[g],shape[s],tail_probability/(8*ns),tail_probability/(8*ns));
                st.lower=(st.data.location+st.data.sigma*(grid.lower-(st.epsilon>0 ? st.relative_shift : 0))+st.moment)/logtwo;
                st.upper=(st.data.location+st.data.sigma*grid.upper+st.moment)/logtwo;
                lower=std::min(lower,st.lower+beta_lower-2*spacing);
                upper=std::max(upper,st.upper+2*spacing);
            }
        }
        offset[gene]=peak+frequency_logmass-(int(m[0]>0)+int(m[1]>0))*std::log(logtwo);
        const double origin=std::floor(lower/qstep)*qstep;
        int ngrid=stride;while ((ngrid-1)*spacing<upper-origin) ngrid*=2;
        const int nfft=2*ngrid; resolution_factor[gene]=1;
        checks(gene,7)=nfft;
        const double period=nfft*spacing;
        SGCBFFT local_fft;
        auto& fft=reuse ? fft_cache[worker] : local_fft;
        performance(gene,0)=std::chrono::duration<double>(std::chrono::steady_clock::now()-gene_started).count();
        const auto curve_started=std::chrono::steady_clock::now();
        std::vector<MeanFullGroupCurve> curve[2];
        curve[0].resize(ns);curve[1].resize(ns);

        const double frequency_sd=point ? s_sd : std::min({a_sd,s_sd,conditional_sd}),
            frequency_budget=tolerance/(256*(1+density_bound)*(1+std::pow(conditional_limit,4))*
                (1+std::pow(1/frequency_sd,6)));
        const int common_frequency_limit=mean_spectral_gaussian_limit(period,a_sd,frequency_budget,0),
            marginal_frequency_limit=std::max({
                mean_spectral_gaussian_limit(period,s_sd,frequency_budget,0),
                mean_spectral_gaussian_limit(period,s_sd,frequency_budget,1),
                mean_spectral_gaussian_limit(period,s_sd,frequency_budget,2)}),
            conditional_frequency_limit=point ? 0 : std::max({
                mean_spectral_gaussian_limit(period,conditional_sd,frequency_budget,0),
                mean_spectral_gaussian_limit(period,conditional_sd,frequency_budget,1),
                mean_spectral_gaussian_limit(period,conditional_sd,frequency_budget,2)});
        for (int g=0;g<2;++g) {
            if (m[g]==0) continue;
            std::vector<int> spectral_size(ns,0);
            int maximum_radius=0;
            for (int s=0;s<ns;++s) if (gene_used[g][s]) {
                int local_size=2;
                for (;;) {
                    const int reference_size=std::max(nfft,local_size);
                    const auto key=std::make_tuple(m[g],shape[s],reference_size);
                    auto width=bandwidth_cache[worker].find(key);
                    if (width==bandwidth_cache[worker].end()) {
                        const double value=mixture_location_grid(m[g],shape[s],tail_probability/4,
                            tail_probability/(4*reference_size)).spacing;
                        width=bandwidth_cache[worker].emplace(key,value).first;
                    }
                    const double allowed=state[g][s].data.sigma*width->second/(2*logtwo);
                    int required=2;while (period/required>allowed) required*=2;
                    if (required<=local_size) break;
                    local_size=required;
                }
                if (state[g][s].epsilon==0 || m[g]<=mixture_allocation_maximum) {
                    const double density_budget=tolerance/
                        (64*(1+maximum_density)*moment_factor);
                    local_size=2*mean_spectral_location_radius(spectral_cache,m[g],
                        state[g][s].data.sigma,state[g][s].q,nfft,local_size/2,
                        period,density_budget);
                }
                spectral_size[s]=local_size;
                maximum_radius=std::max(maximum_radius,local_size/2);
                resolution_factor[gene]=std::max(int(resolution_factor[gene]),std::max(1,local_size/nfft));
                checks(gene,7)=std::max(checks(gene,7),double(local_size));
            }
            const auto beta=mean_spectral_beta(spectral_cache,totals[g],m[g],nfft,
                maximum_radius,period,spike[g],continuous[g]);
            for (int s=0;s<ns;++s) if (gene_used[g][s]) {
                const auto& st=state[g][s];auto& cv=curve[g][s];
                cv.maximum=logtwo*std::sqrt(double(m[g]))/st.data.sigma*
                    std::exp(lgg_constant(st.q/std::sqrt(double(m[g]))));
                const int local_size=spectral_size[s],radius=local_size/2;
                std::vector<std::complex<double>> mu_phi;
                if (st.epsilon==0 || m[g]<=mixture_allocation_maximum) {
                    double error=0;
                    mu_phi=mean_spectral_location_phi(st,m[g],radius,nfft,period,origin,
                        tolerance,spectral_cache,error);
                    checks(gene,0)=std::max(checks(gene,0),error);
                } else {
                    const double local_step=period/local_size;
                    std::vector<double> positive(local_size,0);long double mass=0;
                    const int first=std::max(0,int(std::floor((st.lower-origin)/local_step))-1),
                        last=std::min(local_size-1,int(std::ceil((st.upper-origin)/local_step))+1);
                    for (int j=first;j<=last;++j) {
                        positive[j]=population_mean_grid_density(st,(origin+j*local_step)*logtwo,m[g])*local_step*logtwo;
                        mass+=positive[j];
                    }
                    checks(gene,0)=std::max(checks(gene,0),std::abs(double(mass)-1));
                    for (double& value:positive) value/=mass;
                    fft.fwd(mu_phi,positive);mu_phi.resize(radius+1);
                    mu_phi[radius]*=.5;
                }
                if (st.q==0 && st.epsilon==0) cv.gaussian_phi=mu_phi;
                cv.phi=std::move(mu_phi);
                for (int j=0;j<=radius;++j) cv.phi[j]*=(*beta)[j];
                checks(gene,1)=std::max(checks(gene,1),std::abs(cv.phi[0]-std::complex<double>(1,0)));
                for (const auto& value:cv.phi)
                    checks(gene,2)=std::max(checks(gene,2),std::max(0.0,std::abs(value)-1));
                const auto transform_started=std::chrono::steady_clock::now();
                if (point || m[0]==0 || m[1]==0) {
                    const auto prefix=mean_spectral_prefix(cv.phi,marginal_frequency_limit);
                    cv.coefficient=mean_full_mehler_coefficients_reused(prefix,nfft,stride,
                        spacing,origin,(g==0 ? slab_A0 : slab_A1),s_sd,coefficient_degree,fft);
                } else {
                    const auto& centers=g==0 ? slab_A0 : slab_A1;
                    cv.conditional=mean_conditional_fold(cv.phi,nfft,conditional_stride,
                        spacing,origin,centers,conditional_sd,0,fft,conditional_frequency_limit);
                    if (effects) {
                        cv.conditional_first=mean_conditional_fold(cv.phi,nfft,conditional_stride,
                            spacing,origin,centers,conditional_sd,1,fft,conditional_frequency_limit);
                        cv.conditional_second=mean_conditional_fold(cv.phi,nfft,conditional_stride,
                            spacing,origin,centers,conditional_sd,2,fft,conditional_frequency_limit);
                    }
                }
                if (m[0]==0 || m[1]==0) {
                    const auto prefix=mean_spectral_prefix(cv.phi,
                        mean_spectral_gaussian_limit(period,a_sd,frequency_budget,0));
                    cv.null_smoothed=mean_full_fold_gaussian(prefix,nfft,stride,spacing,origin,A,a_sd,fft);
                }
                if (st.q==0 && st.epsilon==0 && !point && m[0]>0 && m[1]>0)
                    cv.gaussian_conditional=mean_conditional_fold(cv.gaussian_phi,nfft,conditional_stride,
                        spacing,origin,(g==0 ? slab_A0 : slab_A1),conditional_sd,0,fft,conditional_frequency_limit);
                mehler_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-transform_started).count();
            }
        }

        performance(gene,1)=std::chrono::duration<double>(std::chrono::steady_clock::now()-curve_started).count();
        performance(gene,2)=mehler_seconds;
        const auto pairs_started=std::chrono::steady_clock::now();
        long double posterior_first=0,posterior_second=0;
        for (int p:active_pairs) {
            if (effects && std::accumulate(requested.begin()+p*nA,requested.begin()+(p+1)*nA,0.0)==0) continue;
            const int s0=pair_left[p]-1,s1=pair_right[p]-1;
            const auto& c0=curve[0][s0];const auto& c1=curve[1][s1];
            std::vector<double> k0(nA,1),gaussian0;
            std::vector<long double> k1(nA,1),n1(nA,0),n2(nA,t2);
            if (m[0]>0 && m[1]>0) {
                const auto null_started=std::chrono::steady_clock::now();
                k0=mean_spectral_null(c0.phi,c1.phi,nfft,stride,spacing,origin,A,a_sd,common_frequency_limit,fft);
                null_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-null_started).count();
                std::fill(k1.begin(),k1.end(),0);
                std::fill(n2.begin(),n2.end(),0);
                if (point) {
                    const auto exact=mean_spectral_null(c0.phi,c1.phi,nfft,stride,spacing,origin,A,a_sd,common_frequency_limit,fft,slab_mean/2,-slab_mean/2);
                    std::copy(exact.begin(),exact.end(),k1.begin());
                }
                if (!point) {
                    const int nf=c0.conditional.size(),a_stride=stride/conditional_stride;
                    for (int a=0;a<nA;++a) {
                        if (effects && requested[p*nA+a]==0) continue;
                        if (selected && !std::binary_search(selected_states.begin(),selected_states.end(),p*nA+a)) continue;
                        // Advance the exact circular indices instead of dividing at every node.
                        int i0=((a*a_stride-conditional_half)%nf+nf)%nf,
                            i1=((a*a_stride-conditional_sign*conditional_half)%nf+nf)%nf;
                        for (int j=-conditional_half;j<=conditional_half;++j) {
                            const long double u=c0.conditional[i0],v=c1.conditional[i1],
                                w=conditional_weight[j+conditional_half],uv=u*v;
                            k1[a]+=w*uv;
                            if (effects) {
                                const long double du=c0.conditional_first[i0],
                                    dv=c1.conditional_first[i1],
                                    ddu=c0.conditional_second[i0],
                                    ddv=c1.conditional_second[i1],
                                    d=(conditional_sign-1)*j*conditional_step,
                                    r=conditional_variance*(u*dv-du*v);
                                n1[a]+=w*(d*uv+r);
                                n2[a]+=w*((d*d+2*conditional_variance)*uv+2*d*r+
                                    conditional_variance*conditional_variance*
                                        (ddu*v+u*ddv-2*du*dv));
                            }
                            if (++i0==nf) i0=0;
                            i1+=conditional_sign;
                            if (i1==nf) i1=0;
                            if (i1<0) i1=nf-1;
                        }
                    }
                    checks(gene,3)=std::max(checks(gene,3),pair_probability[p]*conditional_bound);
                }
                if (shape[s0]==0 && shape[s1]==0 && epsilon[s0]==0 && epsilon[s1]==0) {
                    checks(gene,8)+=1;
                    gaussian0=mean_spectral_null(c0.gaussian_phi,c1.gaussian_phi,nfft,stride,spacing,origin,A,a_sd,common_frequency_limit,fft);
                    std::vector<double> gaussian1(nA,0);long double gp=1;
                    if (point) gaussian1=mean_spectral_null(c0.gaussian_phi,c1.gaussian_phi,nfft,stride,spacing,origin,A,a_sd,common_frequency_limit,fft,slab_mean/2,-slab_mean/2);
                    if (!point) {
                        const int nf=c0.gaussian_conditional.size(),a_stride=stride/conditional_stride;
                        for (int a=0;a<nA;++a) {
                            long double sum=0;
                            int i0=((a*a_stride-conditional_half)%nf+nf)%nf,
                                i1=((a*a_stride-conditional_sign*conditional_half)%nf+nf)%nf;
                            for (int j=-conditional_half;j<=conditional_half;++j) {
                                sum+=conditional_weight[j+conditional_half]*
                                    c0.gaussian_conditional[i0]*c1.gaussian_conditional[i1];
                                if (++i0==nf) i0=0;
                                i1+=conditional_sign;
                                if (i1==nf) i1=0;
                                if (i1<0) i1=nf-1;
                            }
                            gaussian1[a]=double(sum);
                        }
                    }
                    const double mean0=(center[0]+bvalue[s0])/logtwo,mean1=(center[1]+bvalue[s1])/logtwo,
                        v0=std::exp(2*log_sigma[s0])/(m[0]*logtwo*logtwo),
                        v1=std::exp(2*log_sigma[s1])/(m[1]*logtwo*logtwo);
                    for (int a=0;a<nA;++a) {
                        const double exact0=mean_full_gaussian_pair_density(A[a]/logtwo,mean0,mean1,v0,v1,a_sd*a_sd,0),
                            exact1=mean_full_gaussian_pair_density(A[a]/logtwo,mean0+slab_mean/2,mean1-slab_mean/2,v0,v1,a_sd*a_sd,slab_sd*slab_sd),
                            error0=std::abs(gaussian0[a]-exact0),error1=std::abs(gaussian1[a]-exact1);
                        checks(gene,5)=std::max(checks(gene,5),error0);
                        checks(gene,6)=std::max(checks(gene,6),error1);
                        gaussian_relative(gene,0)=std::max(gaussian_relative(gene,0),error0/(1+exact0));
                        gaussian_relative(gene,1)=std::max(gaussian_relative(gene,1),error1/(1+exact1));
                    }
                }
            } else if (m[0]>0 || m[1]>0) {
                const auto& available=m[0]>0 ? c0 : c1;
                k0=available.null_smoothed;
                for (int a=0;a<nA;++a) {
                    k1[a]=available.coefficient[a];
                    if (effects && !point) {
                        n1[a]=t2*(m[0]>0 ? -1 : 1)*available.coefficient[a+nA]/(2*s_sd);
                        n2[a]=t2*t2*std::sqrt(2.0)*available.coefficient[a+2*nA]/(4*s_sd*s_sd);
                    }
                }
            }
            if (effects && (m[0]>0 || m[1]>0) && (point || m[0]==0 || m[1]==0))
                for (int a=0;a<nA;++a) n2[a]+=t2*k1[a];
            if (effects) for (int a=0;a<nA;++a) {
                n2[a]+=2*slab_mean*n1[a]+slab_mean*slab_mean*k1[a];
                n1[a]+=slab_mean*k1[a];
            }
            if(selected) {
                for(int d=0;d<labels.nrow();++d) {
                    const int z=labels(d,gene)-1;
                    if(z<0 || (z%nk)/nA!=p)continue;
                    const int a=z%nA;
                    checks(gene,4)=std::max(checks(gene,4),pair_probability[p]*
                        std::max({0.0,-k0[a],-double(k1[a])}));
                    selected_values(gene,2*d)=pair_probability[p]*std::max(0.0,k0[a]);
                    selected_values(gene,2*d+1)=pair_probability[p]*std::max(0.0,double(k1[a]));
                }
                continue;
            }
            for (int a=0;a<nA;++a) {
                if (effects) {
                    const double mass=requested[p*nA+a];
                    if (mass==0) continue;
                    const long double fraction=mass,
                        mean=n1[a]/k1[a],second=n2[a]/k1[a];
                    posterior_first+=fraction*mean;posterior_second+=fraction*second;
                    const double expected=reference(gene,nk+p*nA+a),fresh=pair_probability[p]*double(k1[a]);
                    effect_moments(gene,2)+=double(fraction)*std::abs(fresh/expected-1);
                    effect_moments(gene,3)=std::max(effect_moments(gene,3),double(std::max(0.0L,mean*mean-second)));
                    effect_moments(gene,4)+=!std::isfinite(double(mean)) || !std::isfinite(double(second));
                    continue;
                }
                checks(gene,4)=std::max(checks(gene,4),pair_probability[p]*
                    std::max({0.0,-k0[a],-double(k1[a])}));
                likelihood(gene,p*nA+a)=pair_probability[p]*std::max(0.0,k0[a]);
                likelihood(gene,nk+p*nA+a)=pair_probability[p]*std::max(0.0,double(k1[a]));
            }
        }
        performance(gene,3)=null_seconds;
        performance(gene,4)=std::chrono::duration<double>(std::chrono::steady_clock::now()-pairs_started).count();
        performance(gene,5)=std::chrono::duration<double>(std::chrono::steady_clock::now()-gene_started).count();
        if (effects) {
            effect_moments(gene,0)=double(posterior_first);
            effect_moments(gene,1)=double(posterior_second-posterior_first*posterior_first);
        }
    }
    Rcpp::colnames(performance)=Rcpp::CharacterVector::create("preparation","group_curves","conditional_transforms","null_pairs","pair_total","gene_total");
    Rcpp::colnames(effect_moments)=Rcpp::CharacterVector::create("effect","variance",
        "posterior_weighted_likelihood_identity","negative_conditional_variance","nonfinite_moments");
    return Rcpp::List::create(Rcpp::Named("performance_time")=performance,Rcpp::Named("effect_moments")=effect_moments,Rcpp::Named("likelihood")=likelihood,Rcpp::Named("log_offset")=offset,
        Rcpp::Named("checks")=checks,Rcpp::Named("selected_likelihood")=selected_values,
        Rcpp::Named("resolution_factor")=resolution_factor,Rcpp::Named("Gaussian_scaled_error")=gaussian_relative,Rcpp::Named("rho")=rho,Rcpp::Named("mehler_degree")=degree,
        Rcpp::Named("conditional_nodes")=2*conditional_half+1,
        Rcpp::Named("conditional_step")=conditional_step,
        Rcpp::Named("conditional_remainder_bound")=conditional_bound,
        Rcpp::Named("spacing")=spacing,Rcpp::Named("maximum_group_density")=maximum_density,
        Rcpp::Named("spatial_margin")=full_margin,Rcpp::Named("null_margin")=null_margin);
}

// Deterministic empirical-Bayes summaries from the same complete joint likelihood.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix mean_full_mixture_posterior_cpp(Rcpp::NumericMatrix likelihood,
    Rcpp::NumericVector weights,Rcpp::NumericVector b_contrast,int threads=1) {
    const int n=likelihood.nrow(),k=likelihood.ncol(),half=k/2,
        np=b_contrast.size(),nA=half/np;
    if (k!=2*np*nA || weights.size()!=k)
        Rcpp::stop("Joint mixture dimensions differ");
    std::vector<double> partial(size_t(threads)*n*4,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for (int j=0;j<k;++j) {
        if (!(weights[j]>0)) continue;
        int worker=0;
        #ifdef _OPENMP
        worker=omp_get_thread_num();
        #endif
        double* total=partial.data()+size_t(worker)*n*4;
        const double* column=likelihood.begin()+size_t(j)*n;
        const double b=b_contrast[(j%half)/nA];
        for (int i=0;i<n;++i) {
            const double value=weights[j]*column[i];
            total[i]+=value;
            if (j<half) total[n+i]+=value;
            total[2*n+i]+=value*b;total[3*n+i]+=value*b*b;
        }
    }
    Rcpp::NumericMatrix out(n,4);
    Rcpp::colnames(out)=Rcpp::CharacterVector::create(
        "posterior_null","b_contrast","b_contrast_variance","density");
    for (int i=0;i<n;++i) {
        double total=0,null=0,first=0,second=0;
        for (int t=0;t<threads;++t) {
            const double* p=partial.data()+size_t(t)*n*4;
            total+=p[i];null+=p[n+i];first+=p[2*n+i];second+=p[3*n+i];
        }
        out(i,0)=null/total;out(i,1)=first/total;
        out(i,2)=std::max(0.0,second/total-first*first/(total*total));
        out(i,3)=total;
    }
    return out;
}


#include <Rcpp.h>
#include <chrono>
#include <mutex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <mutex>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
using MeanJointLD=long double;
static MeanJointLD mean_joint_round_gamma(int operations) {
    const MeanJointLD x=operations*std::numeric_limits<MeanJointLD>::epsilon();
    return x/(1-x);
}
// [[Rcpp::export(rng = false)]]
Rcpp::List mean_full_joint_initial_labels_cpp(Rcpp::NumericMatrix L,int threads) {
    const int n=L.nrow(),k=L.ncol(),half=k/2,tile=64;
    const double* input=L.begin();
    Rcpp::IntegerMatrix labels(n,4);Rcpp::NumericVector row_upper(n);
    long long nonfinite=0,negative=0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) reduction(+:nonfinite,negative)
    #endif
    for (int low=0;low<n;low+=tile) {
        const int high=std::min(n,low+tile);
        std::vector<double> m0(high-low,R_NegInf),m1(high-low,R_NegInf);
        std::vector<int> z0(high-low,0),z1(high-low,half);
        for (int j=0;j<k;++j) for (int i=low;i<high;++i) {
            const double value=input[i+size_t(n)*j];const int at=i-low;
            nonfinite+=!std::isfinite(value);negative+=value<0;
            if (j<half) {if (value>m0[at]) {m0[at]=value;z0[at]=j;}}
            else if (value>m1[at]) {m1[at]=value;z1[at]=j;}
        }
        for (int i=low;i<high;++i) {
            const int at=i-low;
            labels(i,0)=z0[at]+1;labels(i,1)=z1[at]+1;
            labels(i,2)=(m0[at]>=m1[at] ? z0[at] : z1[at])+1;
            row_upper[i]=std::max(m0[at],m1[at]);
        }
    }
    if (nonfinite || negative) Rcpp::stop("Nonfinite or negative original likelihood");
    return Rcpp::List::create(Rcpp::_["labels"]=labels,Rcpp::_["row_upper"]=row_upper);
}
// Exact heatbath on the occupied partition: component labels may change
// jointly for an entire cluster, including class. Other occupied labels are
// excluded, preserving the partition; single-gene Gibbs handles splits/merges.
static void mean_factorized_cluster_heatbath(
    const double* L,int n,int nA,int nh,
    const Rcpp::NumericVector& alpha,const Rcpp::NumericVector& beta,
    const std::vector<int>& prior_index,const std::vector<std::vector<double>>& rising,
    std::vector<int>& z,std::vector<int>& counts,std::vector<int>& ca,std::vector<int>& nc,
    std::vector<int>& occupied,std::vector<int>& pos,
    double& attempts,double& changes,double& class_changes,int threads) {
    const int H=nh/2;double bt[2]={0,0};
    for (int a=0;a<2*nA;++a) bt[a/nA]+=beta[a];
    std::vector<int> lookup(nh,-1);
    std::vector<std::vector<int>> partition;
    for (int i=0;i<n;++i) {
        const int h=z[i]/nA;
        if (lookup[h]<0) {lookup[h]=partition.size();partition.emplace_back();}
        partition[lookup[h]].push_back(i);
    }
    for (const auto& members:partition) {
        const int oldh=z[members[0]]/nA,oldc=oldh/H,m=members.size();
        std::vector<int> number(nA,0);
        for (int i:members) ++number[z[i]%nA];
        counts[oldh]=0;
        const int at=pos[oldh],last=occupied.back();occupied[at]=last;pos[last]=at;
        occupied.pop_back();pos[oldh]=-1;
        for (int a=0;a<nA;++a) ca[oldc*nA+a]-=number[a];
        nc[oldc]-=m;
        double class_log[2]={0,0};
        for (int c=0;c<2;++c) {
            int denominator=0;
            for (int a=0;a<nA;++a) for (int r=0;r<number[a];++r) {
                class_log[c]+=std::log((ca[c*nA+a]+beta[c*nA+a]+r)/
                    (nc[c]+bt[c]+denominator));
                ++denominator;
            }
        }
        std::vector<double> lp(nh,R_NegInf),cdf(nh);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for (int h=0;h<nh;++h) if (counts[h]==0) {
            double value=rising[prior_index[h]][m]+class_log[h/H];
            for (int i:members) {
                const double density=L[i+size_t(n)*(h*nA+z[i]%nA)];
                if (density==0) {value=R_NegInf;break;}
                value+=std::log(density);
            }
            lp[h]=value;
        }
        const double peak=*std::max_element(lp.begin(),lp.end());
        double total=0;
        for (int h=0;h<nh;++h) {total+=std::exp(lp[h]-peak);cdf[h]=total;}
        const double u=R::runif(0,1)*total;
        const int h=std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin(),c=h/H;
        ++attempts;changes+=h!=oldh;class_changes+=c!=oldc;
        counts[h]=m;pos[h]=occupied.size();occupied.push_back(h);
        for (int a=0;a<nA;++a) ca[c*nA+a]+=number[a];
        nc[c]+=m;
        for (int i:members) z[i]=h*nA+z[i]%nA;
    }
}
// [[Rcpp::export]]
Rcpp::List mean_full_factorized_collapsed_cpp(Rcpp::NumericMatrix L,Rcpp::NumericVector alpha,
    Rcpp::NumericVector beta,Rcpp::NumericVector bvalue,Rcpp::NumericVector A,
    int warmup,int retained,int batch_length,int save_every,int cluster_every,Rcpp::IntegerMatrix initial_labels,int threads) {
    const int n=L.nrow(),k=L.ncol(),nA=A.size(),nh=alpha.size(),H=nh/2,
        chains=initial_labels.ncol(),saved=retained/save_every,tile=64;
    const double* likelihood=L.begin();double at=0,bt[2]={0,0};
    for (double value:alpha) at+=value;
    for (int j=0;j<2*nA;++j) bt[j/nA]+=beta[j];
    const size_t bases=size_t(n)*2*nA;
    std::vector<std::vector<double>> base_cdf(bases);
    std::vector<std::once_flag> base_ready(bases);
    std::vector<double> base_total(bases,0);
    const auto start=std::chrono::steady_clock::now();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for (int lower=0;lower<n;lower+=tile) {
        const int upper=std::min(n,lower+tile);
        for (int h=0;h<nh;++h) for (int a=0;a<nA;++a) {
            const int c=h/H;const double* column=likelihood+size_t(n)*(h*nA+a);
            for (int i=lower;i<upper;++i) if (column[i]>0) {
                const size_t at=i+size_t(n)*(a+nA*c);
                base_total[at]+=alpha[h]*column[i];
            }
        }
    }
    // Prior CDFs are deterministic and shared across chains, materialized only
    // when the exact prior branch is selected. Zero-probability flats are retained.
    auto prior_component=[&](int i,int c,int a,double u) {
        const size_t at=i+size_t(n)*(a+nA*c);
        std::call_once(base_ready[at],[&]() {
            auto& cdf=base_cdf[at];cdf.resize(H);double sum=0;
            for (int local=0;local<H;++local) {
                const int h=c*H+local;
                sum+=alpha[h]*likelihood[i+size_t(n)*(h*nA+a)];cdf[local]=sum;
            }
        });
        const auto& cdf=base_cdf[at];
        return c*H+int(std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin());
    };
    const double prepare_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::vector<std::vector<int>> label(chains,std::vector<int>(n)),count_h(chains,std::vector<int>(nh,0)),
        count_a(chains,std::vector<int>(2*nA,0)),active(chains),position(chains,std::vector<int>(nh,-1));
    std::vector<std::vector<int>> count_class(chains,std::vector<int>(2,0));
    for (int c=0;c<chains;++c) for (int i=0;i<n;++i) {
        int state=initial_labels(i,c)-1;
        if (state<0) {
            std::vector<double> cdf(2*nA);double total=0;
            for (int g=0;g<2;++g) for (int a=0;a<nA;++a) {
                total+=base_total[i+size_t(n)*(a+nA*g)]*beta[g*nA+a]/bt[g];cdf[g*nA+a]=total;
            }
            const int ga=std::upper_bound(cdf.begin(),cdf.end(),R::runif(0,1)*total)-cdf.begin(),
                g=ga/nA,a=ga%nA;
            const size_t at=i+size_t(n)*ga;
            state=prior_component(i,g,a,R::runif(0,1)*base_total[at])*nA+a;
        }
        const int h=state/nA,g=h/H,a=state%nA;label[c][i]=state;
        if (count_h[c][h]++==0) {position[c][h]=active[c].size();active[c].push_back(h);}
        ++count_a[c][g*nA+a];++count_class[c][g];
    }
    std::vector<double> prior_values;std::vector<int> prior_index(nh);
    for (int h=0;h<nh;++h) {
        const auto found=std::find(prior_values.begin(),prior_values.end(),alpha[h]);
        if (found==prior_values.end()) prior_values.push_back(alpha[h]);
        prior_index[h]=std::find(prior_values.begin(),prior_values.end(),alpha[h])-prior_values.begin();
    }
    std::vector<std::vector<double>> rising(prior_values.size(),std::vector<double>(n+1,0));
    for (int j=0;j<int(prior_values.size());++j) {
        rising[j][1]=std::log(prior_values[j]);
        for (int m=2;m<=n;++m) rising[j][m]=rising[j][m-1]+std::log1p(prior_values[j]/(m-1));
    }
    std::vector<double> cluster_attempts(chains,0),cluster_changes(chains,0),cluster_class_changes(chains,0);
    Rcpp::IntegerMatrix labels(retained*chains,n),final_labels(n,chains);
    Rcpp::NumericMatrix b_trace(retained,n*chains),global_trace(retained,6*chains),
        b_first(n,chains),b_second(n,chains),population_draws(nh,saved*chains),mean_draws(2*nA,saved*chains);
    std::vector<double> U(size_t(n)*2*chains),minimum_total(chains,R_PosInf),count_error(chains,0);
    const auto sample_start=std::chrono::steady_clock::now();
    for (int iteration=-warmup;iteration<retained;++iteration) {
        for (double& value:U) value=R::runif(0,1);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for (int chain=0;chain<chains;++chain) {
            auto& z=label[chain];auto& counts=count_h[chain];auto& ca=count_a[chain];
            auto& nc=count_class[chain];auto& occupied=active[chain];auto& pos=position[chain];
            std::vector<double> group_cdf(2*nA),shape_cdf;
            for (int i=0;i<n;++i) {
                const int old=z[i],oldh=old/nA,oldclass=oldh/H,oldA=old%nA,local=oldh%H;
                if (--counts[oldh]==0) {
                    const int at=pos[oldh],last=occupied.back();occupied[at]=last;pos[last]=at;
                    occupied.pop_back();pos[oldh]=-1;
                }
                --ca[oldclass*nA+oldA];--nc[oldclass];
                double total=0;
                for (int g=0;g<2;++g) for (int a=0;a<nA;++a) {
                    const int h=g*H+local,state=h*nA+a;
                    total+=(counts[h]+alpha[h])*(ca[g*nA+a]+beta[g*nA+a])/(nc[g]+bt[g])*
                        likelihood[i+size_t(n)*state];
                    group_cdf[g*nA+a]=total;
                }
                minimum_total[chain]=std::min(minimum_total[chain],total);
                const double uc=U[size_t(n)*2*chain+2*i]*total;
                const int ga=std::min(int(std::upper_bound(group_cdf.begin(),group_cdf.end(),uc)-group_cdf.begin()),2*nA-1),
                    g=ga/nA,a=ga%nA;
                const size_t base=i+size_t(n)*ga;
                total=base_total[base];shape_cdf.resize(occupied.size());
                for (int j=0;j<int(occupied.size());++j) {
                    const int h=occupied[j];
                    if (h/H==g) total+=counts[h]*likelihood[i+size_t(n)*(h*nA+a)];
                    shape_cdf[j]=total;
                }
                minimum_total[chain]=std::min(minimum_total[chain],total);
                const double uh=U[size_t(n)*2*chain+2*i+1]*total;int h;
                if (uh<base_total[base]) {
                    h=prior_component(i,g,a,uh);
                } else {
                    const int j=std::min(int(std::upper_bound(shape_cdf.begin(),shape_cdf.end(),uh)-shape_cdf.begin()),int(occupied.size())-1);
                    h=occupied[j];
                }
                z[i]=h*nA+a;
                if (counts[h]++==0) {pos[h]=occupied.size();occupied.push_back(h);}
                ++ca[g*nA+a];++nc[g];
            }
        }
        if ((iteration+1)%cluster_every==0) for (int chain=0;chain<chains;++chain)
            mean_factorized_cluster_heatbath(likelihood,n,nA,nh,alpha,beta,prior_index,rising,
                label[chain],count_h[chain],count_a[chain],count_class[chain],active[chain],position[chain],
                cluster_attempts[chain],cluster_changes[chain],cluster_class_changes[chain],threads);
        if (iteration>=0) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(threads)
            #endif
            for (int c=0;c<chains;++c) {
                double ma=0,mb=0,ll=0,ac[2]={0,0};
                for (int i=0;i<n;++i) {
                    const int state=label[c][i],h=state/nA,a=state%nA;const double b=bvalue[h];
                    labels(iteration+retained*c,i)=state+1;b_trace(iteration,i+n*c)=b;
                    b_first(i,c)+=b/retained;b_second(i,c)+=b*b/retained;
                    ma+=A[a];mb+=b;ll+=std::log(likelihood[i+size_t(n)*state]);
                }
                for (int g=0;g<2;++g) for (int a=0;a<nA;++a)
                    ac[g]+=(count_a[c][g*nA+a]+beta[g*nA+a])*A[a]/(count_class[c][g]+bt[g]);
                global_trace(iteration,6*c)=ma/n;global_trace(iteration,6*c+1)=mb/n;
                global_trace(iteration,6*c+2)=ll;global_trace(iteration,6*c+3)=active[c].size();
                global_trace(iteration,6*c+4)=double(count_class[c][1])/n;
                global_trace(iteration,6*c+5)=ac[1]-ac[0];
            }
            if (iteration%save_every==0) {
                const int draw=iteration/save_every;
                for (int c=0;c<chains;++c) {
                    std::vector<double> logmass(nh);double peak=R_NegInf;
                    for (int h=0;h<nh;++h) {
                        const double shape=count_h[c][h]+alpha[h];
                        logmass[h]=shape<1 ? std::log(R::rgamma(shape+1,1))+std::log(R::runif(0,1))/shape :
                            std::log(R::rgamma(shape,1));
                        peak=std::max(peak,logmass[h]);
                    }
                    MeanJointLD norm=0;for (double value:logmass) norm+=std::exp(MeanJointLD(value-peak));
                    for (int h=0;h<nh;++h) population_draws(h,draw+saved*c)=double(std::exp(MeanJointLD(logmass[h]-peak))/norm);
                    for (int g=0;g<2;++g) {
                        std::vector<double> value(nA);MeanJointLD normA=0;
                        for (int a=0;a<nA;++a) {value[a]=R::rgamma(count_a[c][g*nA+a]+beta[g*nA+a],1);normA+=value[a];}
                        for (int a=0;a<nA;++a) mean_draws(g*nA+a,draw+saved*c)=value[a]/double(normA);
                    }
                }
            }
        }
        if ((iteration+1)%batch_length==0) {
            Rcpp::Rcout<<"factorized_collapsed "<<(iteration<0 ? "warmup_remaining=" : "retained=")
                <<(iteration<0 ? -(iteration+1) : iteration+1)<<" active/nonnull=";
            for (int c=0;c<chains;++c) {
                int h0=0,h1=0,a0=0,a1=0;
                for (int h=0;h<nh;++h) (h<H ? h0 : h1)+=count_h[c][h];
                for (int a=0;a<nA;++a) {a0+=count_a[c][a];a1+=count_a[c][nA+a];}
                count_error[c]=std::max(count_error[c],double(std::max({std::abs(h0-count_class[c][0]),
                    std::abs(h1-count_class[c][1]),std::abs(a0-h0),std::abs(a1-h1),std::abs(h0+h1-n)})));
                Rcpp::Rcout<<active[c].size()<<"/"<<count_class[c][1]<<(c+1<chains ? "," : "\n");
            }
            Rcpp::checkUserInterrupt();
        }
    }
    const double sample_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sample_start).count();
    for (int c=0;c<chains;++c) for (int i=0;i<n;++i) final_labels(i,c)=label[c][i]+1;
    return Rcpp::List::create(Rcpp::_["labels"]=labels,Rcpp::_["b_trace"]=b_trace,
        Rcpp::_["global_trace"]=global_trace,Rcpp::_["b_first"]=b_first,Rcpp::_["b_second"]=b_second,
        Rcpp::_["population_draws"]=population_draws,Rcpp::_["mean_draws"]=mean_draws,
        Rcpp::_["final_labels"]=final_labels,Rcpp::_["minimum_categorical_mass"]=Rcpp::wrap(minimum_total),
        Rcpp::_["maximum_count_marginal_error"]=Rcpp::wrap(count_error),
        Rcpp::_["prepare_seconds"]=prepare_seconds,Rcpp::_["sample_seconds"]=sample_seconds,
        Rcpp::_["cluster_attempts"]=Rcpp::wrap(cluster_attempts),
        Rcpp::_["cluster_changes"]=Rcpp::wrap(cluster_changes),
        Rcpp::_["cluster_class_changes"]=Rcpp::wrap(cluster_class_changes),
        Rcpp::_["prior_cdf_entries"]=std::accumulate(base_cdf.begin(),base_cdf.end(),0.0,
            [](double total,const std::vector<double>& value) {return total+value.size();}));
}
// [[Rcpp::export(rng = false)]]
Rcpp::List mean_factorized_collapsed_cavity_cpp(Rcpp::NumericMatrix L,
    Rcpp::IntegerMatrix labels,Rcpp::NumericVector alpha,Rcpp::NumericVector beta,
    int nA,int threads,bool verify_dense=false) {
    const int n=L.nrow(),nh=alpha.size(),H=nh/2,draws=labels.nrow(),tile=32;
    const double* likelihood=L.begin();
    std::vector<MeanJointLD> B(size_t(n)*nh,0),C(size_t(n)*2*nA,0),base(size_t(n)*2,0);
    MeanJointLD at=0,ac[2]={0,0},bt[2]={0,0};
    for (int h=0;h<nh;++h) {at+=alpha[h];ac[h/H]+=alpha[h];}
    for (int a=0;a<2*nA;++a) bt[a/nA]+=beta[a];
    const auto start=std::chrono::steady_clock::now();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for (int low=0;low<n;low+=tile) {
        const int high=std::min(n,low+tile);
        for (int h=0;h<nh;++h) for (int a=0;a<nA;++a) {
            const int ca=(h/H)*nA+a;const double* column=likelihood+size_t(n)*(h*nA+a);
            for (int i=low;i<high;++i) {
                B[i+size_t(n)*h]+=MeanJointLD(beta[ca])*column[i];
                C[i+size_t(n)*ca]+=MeanJointLD(alpha[h])*column[i];
            }
        }
        for (int ca=0;ca<2*nA;++ca) for (int i=low;i<high;++i)
            base[i+size_t(n)*(ca/nA)]+=MeanJointLD(beta[ca])*C[i+size_t(n)*ca];
    }
    const double prepare=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    Rcpp::NumericMatrix U0(n,draws),V0(n,draws),U1(n,draws),V1(n,draws),cavity(n,draws);
    Rcpp::IntegerVector states(n);
    double minimum=R_PosInf,identity=0,mass_error=0;long long maximum_terms=0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads) reduction(min:minimum) reduction(max:identity,mass_error,maximum_terms)
    #endif
    for (int d=0;d<draws;++d) {
        std::vector<int> count_h(nh,0),count_a(2*nA,0),nc(2,0),hs[2],as[2];
        for (int i=0;i<n;++i) {
            const int z=labels(d,i)-1,h=z/nA,c=h/H,a=z%nA;
            ++count_h[h];++count_a[c*nA+a];++nc[c];
        }
        for (int h=0;h<nh;++h) if (count_h[h]>0) hs[h/H].push_back(h);
        for (int ca=0;ca<2*nA;++ca) if (count_a[ca]>0) as[ca/nA].push_back(ca%nA);
        long long terms=0;
        for (int c=0;c<2;++c) terms+=1+hs[c].size()+as[c].size()+hs[c].size()*as[c].size();
        maximum_terms=std::max(maximum_terms,terms);
        for (int i=0;i<n;++i) {
            const int own=labels(d,i)-1,oh=own/nA,oc=oh/H,oa=own%nA;
            MeanJointLD density[2]={0,0},mass[2]={0,0},denominator=n-1+at;
            for (int c=0;c<2;++c) {
                MeanJointLD total=base[i+size_t(n)*c];
                for (int a:as[c]) {
                    const int number=count_a[c*nA+a]-(c==oc && a==oa);
                    total+=number*C[i+size_t(n)*(c*nA+a)];
                }
                for (int h:hs[c]) {
                    const int number_h=count_h[h]-(h==oh);
                    total+=number_h*B[i+size_t(n)*h];
                    for (int a:as[c]) {
                        const int number_a=count_a[c*nA+a]-(c==oc && a==oa);
                        total+=MeanJointLD(number_h)*number_a*likelihood[i+size_t(n)*(h*nA+a)];
                    }
                }
                density[c]=total/(denominator*(nc[c]-(c==oc)+bt[c]));
                mass[c]=(nc[c]-(c==oc)+ac[c])/denominator;
                if (verify_dense) {
                    MeanJointLD expected=0;
                    for (int h=c*H;h<(c+1)*H;++h) for (int a=0;a<nA;++a)
                        expected+=(count_h[h]-(h==oh)+MeanJointLD(alpha[h]))*
                            (count_a[c*nA+a]-(c==oc && a==oa)+MeanJointLD(beta[c*nA+a]))*
                            likelihood[i+size_t(n)*(h*nA+a)];
                    expected/=denominator*(nc[c]-(c==oc)+bt[c]);
                    identity=std::max(identity,double(std::abs(density[c]/expected-1)));
                }
            }
            const MeanJointLD total=density[0]+density[1],inverse=1/total;
            minimum=std::min(minimum,double(total));
            mass_error=std::max(mass_error,double(std::abs(mass[0]+mass[1]-1)));
            U0(i,d)=double(density[0]*inverse);U1(i,d)=double(density[1]*inverse);
            V0(i,d)=double(mass[0]*inverse);V1(i,d)=double(mass[1]*inverse);cavity(i,d)=double(inverse);
        }
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    // Positive summation and two products/divisions in extended precision;
    // final storage adds one double rounding. No adaptive omission is used.
    const MeanJointLD bound=mean_joint_round_gamma(H+nA+2*maximum_terms+12)+
        2*std::numeric_limits<double>::epsilon();
    std::fill(states.begin(),states.end(),std::min<long long>(maximum_terms,std::numeric_limits<int>::max()));
    return Rcpp::List::create(Rcpp::_["U0"]=U0,Rcpp::_["V0"]=V0,Rcpp::_["U1"]=U1,Rcpp::_["V1"]=V1,
        Rcpp::_["cavity"]=cavity,Rcpp::_["minimum_density"]=minimum,
        Rcpp::_["weight_mass_error"]=mass_error,Rcpp::_["maximum_density_relative_error_bound"]=double(bound),
        Rcpp::_["states_summed"]=states,Rcpp::_["maximum_dense_identity"]=identity,
        Rcpp::_["prepare_seconds"]=prepare,Rcpp::_["prediction_seconds"]=seconds-prepare,
        Rcpp::_["population_integration"]="exact conjugate Dirichlet integration conditional on all other gene labels");
}



#include <Eigen/Eigenvalues>

// Expand log L_positive + log L_binomial directly into the existing CSC buffer.
// Flat column = theta + n_theta * pi, matching R's expand.grid ordering.
// [[Rcpp::export(rng = false)]]
Rcpp::List product_population_dictionary_cpp(Rcpp::NumericMatrix positive_loglike,
    Rcpp::NumericMatrix frequency_loglike,int threads=1,double tolerance=1e-6) {
    const int n=positive_loglike.nrow(),nt=positive_loglike.ncol(),
        np=frequency_loglike.ncol(),nk=nt*np;
    if(frequency_loglike.nrow()!=n)Rcpp::stop("Likelihood factors have different row counts");
    std::unique_ptr<SparsePopulationBuffer> buffer(new SparsePopulationBuffer);
    buffer->n=n;buffer->nk=nk;buffer->omitted_bound=tolerance/4;
    buffer->threshold=buffer->omitted_bound/n;
    buffer->rows.resize(n);std::iota(buffer->rows.begin(),buffer->rows.end(),0);
    buffer->row_scale.resize(n);buffer->umbrella.assign(n,0);
    std::vector<double> a(size_t(n)*nt),b(size_t(n)*np);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int i=0;i<n;++i) {
        double ma=R_NegInf,mb=R_NegInf;
        for(int k=0;k<nt;++k)ma=std::max(ma,positive_loglike(i,k));
        for(int p=0;p<np;++p)mb=std::max(mb,frequency_loglike(i,p));
        buffer->row_scale[i]=ma+mb;
        for(int k=0;k<nt;++k)a[size_t(k)*n+i]=std::exp(positive_loglike(i,k)-ma);
        for(int p=0;p<np;++p)b[size_t(p)*n+i]=std::exp(frequency_loglike(i,p)-mb);
    }
    std::vector<size_t> count(nk,0);
    std::vector<double> partial(size_t(threads)*n,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int k=0;k<nk;++k) {
        int thread=0;
        #ifdef _OPENMP
        thread=omp_get_thread_num();
        #endif
        double* sum=partial.data()+size_t(thread)*n;
        const double* lhs=a.data()+size_t(k%nt)*n,*rhs=b.data()+size_t(k/nt)*n;
        for(int i=0;i<n;++i) {
            const double value=lhs[i]*rhs[i];
            if(value>buffer->threshold){++count[k];sum[i]+=value/nk;}
        }
    }
    for(int i=0;i<n;++i)for(int t=0;t<threads;++t)
        buffer->umbrella[i]+=partial[size_t(t)*n+i];
    buffer->column_start.assign(size_t(nk)+1,0);
    for(int k=0;k<nk;++k)buffer->column_start[k+1]=buffer->column_start[k]+count[k];
    buffer->column_row.resize(buffer->column_start.back());
    buffer->column_value.resize(buffer->column_start.back());
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int k=0;k<nk;++k) {
        size_t at=buffer->column_start[k];
        const double* lhs=a.data()+size_t(k%nt)*n,*rhs=b.data()+size_t(k/nt)*n;
        for(int i=0;i<n;++i) {
            const double value=lhs[i]*rhs[i];
            if(value>buffer->threshold){buffer->column_row[at]=i;buffer->column_value[at++]=value;}
        }
    }
    const double entries=buffer->column_value.size(),
        bytes=entries*(sizeof(double)+sizeof(int))+
            buffer->column_start.size()*sizeof(size_t)+
            (buffer->row_scale.size()+buffer->umbrella.size())*sizeof(double)+
            buffer->rows.size()*sizeof(int);
    Rcpp::XPtr<SparsePopulationBuffer> pointer(buffer.release(),true);
    return Rcpp::List::create(Rcpp::Named("dictionary")=pointer,Rcpp::Named("rows")=n,
        Rcpp::Named("states")=nk,Rcpp::Named("entries")=entries,
        Rcpp::Named("storage_bytes")=bytes,Rcpp::Named("sparsification_bound")=tolerance/4);
}

static Eigen::Vector2d population_cv_profile_moment_gradient(double order,double s,double q) {
    const double rs=order*std::exp(s);
    if(q==0)return Eigen::Vector2d(rs*rs,-rs/2-rs*rs*rs/6);
    const double k=1/(q*q),t=rs/q;
    if(!(k+t>0))return Eigen::Vector2d(NA_REAL,NA_REAL);
    if(std::abs(q)>=.05) {
        const double a=R::digamma(k+t)-std::log(k);
        return Eigen::Vector2d(t*a,(-2*k*(R::digamma(k+t)-R::digamma(k))-t*a+2*t)/q);
    }
    const double nodes[8]={.0198550717512319,.101666761293187,.237233795041836,
        .408282678752175,.591717321247825,.762766204958164,.898333238706813,.980144928248768};
    const double weights[8]={.0506142681451881,.111190517226687,.156853322938944,
        .181341891689181,.181341891689181,.156853322938944,.111190517226687,.0506142681451881};
    Eigen::Vector2d out=Eigen::Vector2d::Zero();
    for(int j=0;j<8;++j) {
        const double v=nodes[j]*rs*q,inv=q*q/(1+v),w=v/(1+v),
            c=-inv/2-inv*inv/12+std::pow(inv,4)/120-std::pow(inv,6)/252,
            cp=-.5-inv/6+std::pow(inv,3)/30-std::pow(inv,5)/42,
            wl=std::abs(v)<.005?v*v*(-.5+v*(2.0/3+v*(-.75+
                v*(.8+v*(-5.0/6+v*(6.0/7-v*.875)))))):w-std::log1p(v);
        out[0]+=weights[j]*t*(std::log1p(v)+c+w*(1-inv*cp));
        out[1]+=weights[j]*(t/q)*(wl-c+inv*cp*(2-w));
    }
    return out;
}
static Eigen::Vector3d population_cv_profile_base_score(double z,double q,double derivative) {
    const double t=q*z;
    double u,dq;
    if(std::abs(t)<.005) {
        u=z*(1+t*(.5+t*(1.0/6+t*(1.0/24+t/120))));
        dq=z*z*z*(1.0/6+t*(1.0/12+t*(1.0/40+t/180)));
    } else {
        const double ex=std::expm1(t);
        u=ex/q;dq=(t*ex-2*(ex-t))/(q*q*q);
    }
    return Eigen::Vector3d(u,-1+z*u,derivative-dq);
}

// Profiled positive-law information and sample-aligned target influences.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix population_joint_profile_information_cpp(
    Rcpp::NumericMatrix positive_log,Rcpp::IntegerVector group,
    Rcpp::NumericMatrix log_mean,Rcpp::NumericMatrix responsibility,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,
    Rcpp::NumericVector epsilon,Rcpp::NumericVector relative_shift,int threads,int positive_family=0) {
    const int ng=positive_log.nrow(),ns=positive_log.ncol(),nk=shape.size();
    std::vector<int> positive_count(2*ng,0);
    for(int j=0;j<ns;++j)for(int gene=0;gene<ng;++gene)
        if(std::isfinite(positive_log(gene,j)))++positive_count[gene+group[j]*ng];
    Rcpp::NumericMatrix out(nk,17);std::fill(out.begin(),out.end(),NA_REAL);
    Rcpp::NumericMatrix influence(nk,ns);std::fill(influence.begin(),influence.end(),NA_REAL);
    Rcpp::NumericMatrix base_sample(nk,ns);std::fill(base_sample.begin(),base_sample.end(),NA_REAL);
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("log_positive_cv2","log_positive_cv2_sd",
        "parameter_dimension","information_rank","gradient_unidentified_fraction",
        "weighted_observations","effective_samples","gene_group_mass",
        "log_sigma_sd","Q_sd","absolute_shift_scaled_sd","logit_epsilon_sd",
        "minimum_information_eigenvalue","maximum_information_eigenvalue",
        "log_base_cv2","log_base_cv2_sd","base_total_logcv2_covariance");
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for(int k=0;k<nk;++k) {
        const double s=log_sigma[k],q=shape[k],e=epsilon[k],r=relative_shift[k],
            sigma=std::exp(s),delta=r*sigma;
        if(!(1+2*sigma*q>0))continue;
        const bool mixture=e>0 && r>0;
        const int nb=positive_family==0?2:1,d=nb+(mixture?2:0),p=d+1;
        const double moment=log_moment(1,V3{0,s,q})+
            population_shift_logmean(e,delta),
            base=log_moment(2,V3{0,s,q})-2*log_moment(1,V3{0,s,q}),
            K=(base+population_shift_logmean(e,2*delta)-2*population_shift_logmean(e,delta)),g=(K+std::log(-std::expm1(-K))),
            derivative=std::abs(q)<.05?-q/6+std::pow(q,5)/60-std::pow(q,9)/126:
                1/q-2/(q*q*q)*(std::log(1/(q*q))-R::digamma(1/(q*q)));
        Eigen::MatrixXd bread=Eigen::MatrixXd::Zero(d,d),clusters=Eigen::MatrixXd::Zero(ns,d);
        Eigen::VectorXd exposure=Eigen::VectorXd::Zero(ns);
        double total=0,mass=0;
        for(int row=0;row<2*ng;++row) {
            const double w=responsibility(row,k);
            double mu=log_mean(row,k)-moment;
            if(!(w>0) || !std::isfinite(mu))continue;
            // Every latent allocation has the same conjugate location offset.
            // This is the exact pure-GG profile location and the allocation-
            // averaged conditional profile location for shifted mixtures.
            if(q!=0) {
                const double mk=positive_count[row]/(q*q);
                mu-=(std::log(mk)-R::digamma(mk))*sigma/q;
            }
            const int gene=row%ng,gr=row/ng;
            Eigen::MatrixXd scores=Eigen::MatrixXd::Zero(ns,p);
            Eigen::VectorXd active=Eigen::VectorXd::Zero(ns);
            int n=0;
            for(int j=0;j<ns;++j) {
                const double y=positive_log(gene,j);
                if(group[j]!=gr || !std::isfinite(y))continue;
                const double z=(y-mu)/sigma;
                double gamma=0;
                if(mixture) {
                    double difference=r*(r/2-z);
                    if(q!=0)difference=r/q+(q>0?-1:1)*
                        std::exp(q*z+std::log(std::abs(std::expm1(-q*r)))-2*std::log(std::abs(q)));
                    gamma=R::plogis(std::log(e)-std::log1p(-e)-difference,0,1,true,false);
                }
                Eigen::VectorXd score=Eigen::VectorXd::Zero(p);
                Eigen::Vector3d base_score=Eigen::Vector3d::Zero();double shifted_location=0;
                if(gamma<1)base_score=(1-gamma)*population_cv_profile_base_score(z,q,derivative);
                if(gamma>0) {
                    const Eigen::Vector3d shifted=population_cv_profile_base_score(z-r,q,derivative);
                    base_score+=gamma*shifted;shifted_location=gamma*shifted[0];
                }
                score[0]=base_score[0];
                score[1]=base_score[1]+(positive_family==1?q*base_score[2]:0);
                if(nb==2)score[2]=base_score[2];
                if(mixture){score[nb+1]=shifted_location;score[nb+2]=gamma-e;}
                if(!score.allFinite())continue;
                scores.row(j)=score.transpose();active[j]=1;++n;
            }
            if(n<2)continue;
            const Eigen::VectorXd mean=scores.colwise().sum().transpose()/n;
            const Eigen::MatrixXd centered=scores-active*mean.transpose(),
                information=centered.transpose()*centered;
            if(!(information(0,0)>0))continue;
            const Eigen::VectorXd project=information.block(1,0,d,1)/information(0,0);
            const Eigen::MatrixXd efficient=centered.rightCols(d)-centered.col(0)*project.transpose();
            bread.noalias()+=w*(information.block(1,1,d,d)-
                information.block(1,0,d,1)*information.block(0,1,1,d)/information(0,0));
            clusters.noalias()+=w*efficient;
            exposure+=w*active;total+=w*n;mass+=w;
        }
        out(k,0)=g;out(k,14)=(base+std::log(-std::expm1(-base)));
        if(!(total>0)){out(k,1)=R_PosInf;continue;}
        const double nc=total*total/exposure.squaredNorm();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(bread);
        const Eigen::VectorXd values=eig.eigenvalues();
        const double threshold=std::numeric_limits<double>::epsilon()*d*std::max(0.0,values.maxCoeff());
        Eigen::VectorXd reciprocal=Eigen::VectorXd::Zero(d);
        int rank=0;
        for(int j=0;j<d;++j)if(values[j]>threshold){reciprocal[j]=1/values[j];++rank;}
        const Eigen::MatrixXd rotated_scores=clusters*eig.eigenvectors(),
            parameter_influence=std::sqrt(nc/(nc-1))*rotated_scores*
                reciprocal.asDiagonal()*eig.eigenvectors().transpose(),
            covariance=parameter_influence.transpose()*parameter_influence;
        Eigen::VectorXd gradient=Eigen::VectorXd::Zero(d);
        const Eigen::Vector2d moment_gradient=population_cv_profile_moment_gradient(2,s,q)-
            2*population_cv_profile_moment_gradient(1,s,q);
        gradient[0]=moment_gradient[0]+(positive_family==1?q*moment_gradient[1]:0);
        if(nb==2)gradient[1]=moment_gradient[1];
        if(mixture) {
            const double eta=std::log(e)-std::log1p(-e),
                p1=R::plogis(eta+delta,0,1,true,false),p2=R::plogis(eta+2*delta,0,1,true,false);
            gradient[nb]=2*sigma*(p2-p1);gradient[nb+1]=p2-2*p1+e;
        }
        gradient/=(-std::expm1(-K));
        const Eigen::VectorXd projected=eig.eigenvectors().transpose()*gradient,
            target_influence=std::sqrt(nc/(nc-1))*rotated_scores*
                (reciprocal.array()*projected.array()).matrix();
        double missing=0;
        for(int j=0;j<d;++j)if(!(values[j]>threshold))missing+=projected[j]*projected[j];
        const double fraction=missing/std::max(std::numeric_limits<double>::min(),gradient.squaredNorm());
        out(k,0)=g;out(k,1)=fraction>std::sqrt(std::numeric_limits<double>::epsilon())?
            R_PosInf:target_influence.norm();
        out(k,2)=d;out(k,3)=rank;out(k,4)=fraction;out(k,5)=total;out(k,6)=nc;out(k,7)=mass;
        out(k,8)=std::sqrt(std::max(0.0,covariance(0,0)));
        out(k,9)=nb==2?std::sqrt(std::max(0.0,covariance(1,1))):positive_family==1?q*out(k,8):0;
        if(mixture){out(k,10)=std::sqrt(std::max(0.0,covariance(nb,nb)));
            out(k,11)=std::sqrt(std::max(0.0,covariance(nb+1,nb+1)));}
        out(k,12)=values.minCoeff();out(k,13)=values.maxCoeff();
        Eigen::VectorXd base_gradient=Eigen::VectorXd::Zero(d);
        base_gradient[0]=(moment_gradient[0]+(positive_family==1?q*moment_gradient[1]:0))/(-std::expm1(-base));
        if(nb==2)base_gradient[1]=moment_gradient[1]/(-std::expm1(-base));
        out(k,14)=(base+std::log(-std::expm1(-base)));
        const Eigen::VectorXd base_projected=eig.eigenvectors().transpose()*base_gradient,
            base_influence=std::sqrt(nc/(nc-1))*rotated_scores*
                (reciprocal.array()*base_projected.array()).matrix();
        out(k,15)=base_influence.norm();
        out(k,16)=base_influence.dot(target_influence);
        if(nc>1 && std::isfinite(out(k,15)))for(int j=0;j<ns;++j)base_sample(k,j)=base_influence[j];
        if(std::isfinite(out(k,1)) && nc>1) {
            for(int j=0;j<ns;++j)influence(k,j)=target_influence[j];
        }
    }
    out.attr("sample_influence_log_positive_cv2")=influence;
    out.attr("sample_influence_log_base_cv2")=base_sample;
    return out;
}



// [[Rcpp::export(rng = false)]]
Rcpp::List joint_state_responsibilities_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::IntegerVector theta,Rcpp::IntegerVector pi_index,int ntheta,int npi,
    Rcpp::IntegerVector left,Rcpp::IntegerVector right,Rcpp::NumericVector prior,int threads=1) {
    const int ng=loglike.nrow()/2,H=prior.size();
    Rcpp::NumericMatrix rt(2*ng,ntheta),rp(2*ng,npi);
    Rcpp::NumericVector error(ng);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,8) num_threads(threads)
    #endif
    for(int i=0;i<ng;++i) {
        std::vector<double> post(H);double peak=R_NegInf,total=0;
        for(int h=0;h<H;++h) {
            post[h]=std::log(prior[h])+loglike(i,left[h])+loglike(ng+i,right[h]);
            peak=std::max(peak,post[h]);
        }
        for(double& w:post){w=std::exp(w-peak);total+=w;}
        double mass=0;
        for(int h=0;h<H;++h) {
            const double w=post[h]/total;mass+=w;
            rt(i,theta[left[h]])+=w;rt(ng+i,theta[right[h]])+=w;
            rp(i,pi_index[left[h]])+=w;rp(ng+i,pi_index[right[h]])+=w;
        }
        error[i]=std::abs(mass-1);
    }
    return Rcpp::List::create(Rcpp::Named("positive")=rt,Rcpp::Named("frequency")=rp,
        Rcpp::Named("maximum_mass_error")=*std::max_element(error.begin(),error.end()));
}


// Studentized reference bounds respecting mixture CV2 >= its positive base CV2.
// Four scalar contrast references give Frechet bounds for max(base,total).
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix joint_logcv_reference_summary_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::IntegerVector positive,Rcpp::LogicalVector usable,Rcpp::NumericVector atom_log_cv2,
    Rcpp::LogicalVector mixed,Rcpp::IntegerVector left,Rcpp::IntegerVector right,
    Rcpp::NumericVector prior,Rcpp::NumericMatrix contrast_mean,
    Rcpp::NumericMatrix contrast_sd,Rcpp::NumericMatrix contrast_df,
    Rcpp::LogicalMatrix contrast_identified,Rcpp::NumericVector truth,
    double alpha=.05,double probability_tolerance=1e-8,double interval_step=.025,
    int threads=1,bool compute_intervals=true) {
    const int ng=loglike.nrow()/2,H=prior.size();const double logtwo=std::log(2.0);
    std::vector<double> point(H,NA_REAL);
    std::vector<char> defined(H,0),known(H,0);
    for(int h=0;h<H;++h) {
        const double a=atom_log_cv2[left[h]],b=atom_log_cv2[right[h]];
        defined[h]=std::isfinite(a)&&std::isfinite(b);
        if(!defined[h])continue;
        point[h]=b-a;
        if(std::abs(point[h])<=64*std::numeric_limits<double>::epsilon()*(1+std::abs(a)+std::abs(b)))point[h]=0;
        known[h]=contrast_identified(h,0)&&std::isfinite(contrast_sd(h,0))&&
            contrast_sd(h,0)>=0&&contrast_df(h,0)>0;
    }
    auto pair_cdf=[&](int h,double value,bool strict,bool reverse) {
        if(left[h]==right[h]) {
            const double p=double(strict?0<value:0<=value);return std::make_pair(p,p);
        }
        double lo[4],hi[4];const int order[4]={0,2,1,3};
        for(int c=0;c<4;++c) {
            const int j=reverse?order[c]:c;
            if(!contrast_identified(h,j)||!std::isfinite(contrast_mean(h,j))||
                !std::isfinite(contrast_sd(h,j))||!(contrast_df(h,j)>0)) {
                lo[c]=0;hi[c]=1;continue;
            }
            const double mean=(reverse?-1:1)*contrast_mean(h,j),sd=contrast_sd(h,j);
            const double p=sd>0?R::pt((value-mean)/sd,contrast_df(h,j),true,false):
                double(strict?mean<value:mean<=value);
            lo[c]=p;hi[c]=p;
        }
        const bool p0=!mixed[reverse?right[h]:left[h]],p1=!mixed[reverse?left[h]:right[h]];
        if(p0&&p1)return std::make_pair(lo[0],hi[0]);
        if(p0)return std::make_pair(std::max(0.0,lo[0]+lo[1]-1),std::min(hi[0],hi[1]));
        if(p1)return std::make_pair(std::max(lo[0],lo[2]),std::min(1.0,hi[0]+hi[2]));
        return std::make_pair(std::max({0.0,lo[0]+lo[1]-1,lo[2]+lo[3]-1}),
            std::min(1.0,std::min(hi[0],hi[1])+std::min(hi[2],hi[3])));
    };
    const std::vector<std::string> names={"effect","se","p","cdf_zero_lower","cdf_zero_upper",
        "cdf_truth_lower","cdf_truth_upper","available","undefined_mass",
        "parameter_unidentified_mass","numerical_cdf_error","lower","upper","width",
        "interval_covered","log_cv2_0_mean","log_cv2_1_mean"};
    Rcpp::NumericMatrix out(ng,names.size());std::fill(out.begin(),out.end(),NA_REAL);
    Rcpp::colnames(out)=Rcpp::wrap(names);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,8) num_threads(threads)
    #endif
    for(int i=0;i<ng;++i) {
        out(i,7)=0;
        if(!(usable[i]&&usable[i+ng])||positive[i]<4||positive[i+ng]<4)continue;
        std::vector<double> post(H);double peak=R_NegInf,total=0;
        for(int h=0;h<H;++h) {
            post[h]=std::log(prior[h])+loglike(i,left[h])+loglike(i+ng,right[h]);
            peak=std::max(peak,post[h]);
        }
        for(double& w:post){w=std::exp(w-peak);total+=w;}
        std::vector<int> active;
        double finite_mass=0,first=0,second=0,mu0=0,mu1=0,undefined=0,unidentified=0,
            omitted=0,lo_tail=0,hi_tail=0,flo0=0,fhi0=0,lower=R_PosInf,upper=R_NegInf;
        for(int h=0;h<H;++h) {
            const double w=post[h]/total;post[h]=w;if(!(w>0))continue;
            if(!defined[h]){undefined+=w;lo_tail+=w;hi_tail+=w;continue;}
            finite_mass+=w;first+=w*point[h];
            mu0+=w*atom_log_cv2[left[h]];mu1+=w*atom_log_cv2[right[h]];
            if(!known[h]){unidentified+=w;lo_tail+=w;hi_tail+=w;continue;}
            second+=w*(point[h]*point[h]+contrast_sd(h,0)*contrast_sd(h,0));
            const auto fl=pair_cdf(h,0,true,false),fu=pair_cdf(h,0,false,false),
                right_cdf=pair_cdf(h,0,false,true);
            flo0+=w*fl.first;fhi0+=w*fu.second;
            lo_tail+=w*fu.second;hi_tail+=w*right_cdf.second;
            if(w<=probability_tolerance/H){omitted+=w;continue;}
            active.push_back(h);
            for(int c=0;c<4;++c)if(contrast_identified(h,c)&&
                std::isfinite(contrast_mean(h,c))&&std::isfinite(contrast_sd(h,c))&&contrast_df(h,c)>0) {
                const double z=R::qt(probability_tolerance/32,contrast_df(h,c),false,false);
                lower=std::min(lower,contrast_mean(h,c)-z*contrast_sd(h,c));
                upper=std::max(upper,contrast_mean(h,c)+z*contrast_sd(h,c));
            }
        }
        const double unknown=undefined+unidentified;
        out(i,8)=undefined;out(i,9)=unidentified;out(i,10)=omitted;
        out(i,3)=flo0;out(i,4)=std::min(1.0,fhi0+unknown);
        if(finite_mass>0){out(i,0)=first/finite_mass/logtwo;out(i,15)=mu0/finite_mass;out(i,16)=mu1/finite_mass;}
        if(unknown>alpha||!(finite_mass>0))continue;
        out(i,7)=1;out(i,2)=2*std::min({.5,lo_tail,hi_tail});
        if(unknown<=probability_tolerance) {
            const double m=first/finite_mass;
            out(i,1)=std::sqrt(std::max(0.0,second/finite_mass-m*m))/logtwo;
        }
        auto cdf=[&](double value,bool strict) {
            double lo=0,hi=unknown+omitted;
            for(int h:active){const auto p=pair_cdf(h,value,strict,false);lo+=post[h]*p.first;hi+=post[h]*p.second;}
            return std::make_pair(std::max(0.0,lo),std::min(1.0,hi));
        };
        if(std::isfinite(truth[i])) {
            out(i,5)=cdf(truth[i]*logtwo,true).first;
            out(i,6)=cdf(truth[i]*logtwo,false).second;
        }
        if(!compute_intervals)continue;
        double qlo=R_NegInf,qhi=R_PosInf;
        if(unknown+omitted<alpha/2&&!active.empty()&&std::isfinite(lower)&&std::isfinite(upper)) {
            lower-=interval_step;upper+=interval_step;
            if(cdf(lower,false).second<alpha/2) {
                double lo=lower,hi=upper;
                while(hi-lo>interval_step){const double m=lo+(hi-lo)/2;if(!(m>lo&&m<hi))break;
                    if(cdf(m,false).second<alpha/2)lo=m;else hi=m;}
                qlo=lo/logtwo;
            }
            if(cdf(upper,false).first>=1-alpha/2) {
                double lo=lower,hi=upper;
                while(hi-lo>interval_step){const double m=lo+(hi-lo)/2;if(!(m>lo&&m<hi))break;
                    if(cdf(m,false).first<1-alpha/2)lo=m;else hi=m;}
                qhi=hi/logtwo;
            }
        }
        out(i,11)=qlo;out(i,12)=qhi;out(i,13)=qhi-qlo;
        if(std::isfinite(truth[i]))out(i,14)=double(qlo<=truth[i]&&qhi>=truth[i]);
    }
    return out;
}
