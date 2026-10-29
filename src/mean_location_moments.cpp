// Conditional GG location approximation with exact first four pure-GG cumulants.
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <numeric>
#include <chrono>
#include "sgcb_sparse_mixture.h"
// [[Rcpp::depends(RcppEigen)]]
#ifdef _OPENMP
#include <omp.h>
#endif
struct LocationNormalPair {
    double weight[2]={1,0},shift[2]={0,0},point_variance=0;
    int count=1;
};
// [[Rcpp::export(rng = false)]]
Rcpp::List gaussian_location_dictionary_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::NumericMatrix mean,Rcpp::NumericMatrix variance,Rcpp::IntegerVector positive,
    Rcpp::IntegerVector left,Rcpp::IntegerVector right,Rcpp::NumericVector A,
    double A_sd,double slab_mean,double slab_sd,int threads=1,
    Rcpp::Nullable<Rcpp::NumericVector> weights=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> log_sigma=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> shape=R_NilValue) {
    const int n=loglike.nrow()/2,ns=mean.ncol(),np=left.size(),na=A.size(),half=np*na,nk=2*half;
    const bool summarize=weights.isNotNull(),skew=shape.isNotNull();
    const Rcpp::NumericVector w=summarize?Rcpp::NumericVector(weights):Rcpp::NumericVector(0),
        sig=skew?Rcpp::NumericVector(log_sigma):Rcpp::NumericVector(0),
        Q=skew?Rcpp::NumericVector(shape):Rcpp::NumericVector(0);
    const int maximum=*std::max_element(positive.begin(),positive.end());
    const double l2=std::log(2.0),av=A_sd*A_sd/(l2*l2),tv=slab_sd*slab_sd,
        pi=3.14159265358979323846;
    std::vector<LocationNormalPair> component((maximum+1)*ns);
    double moment_error=0,minimum_variance=R_PosInf;
    if(skew) for(int m=1;m<=maximum;++m) for(int s=0;s<ns;++s) if(Q[s]!=0) {
        const double a=m/(Q[s]*Q[s]),p1=R::trigamma(a),
            g1=-std::copysign(1.0,Q[s])*R::psigamma(a,2)/std::pow(p1,1.5),
            g2=R::psigamma(a,3)/(p1*p1),base_variance=std::exp(2*sig[s])*p1/(Q[s]*Q[s]);
        double lower=0,upper=1;
        // A two-point law convolved with a Gaussian matches all four cumulants:
        // 2*u^3 + excess_kurtosis*u - skewness^2 = 0.
        for(int j=0;j<80;++j) {
            const double u=(lower+upper)/2;
            if(2*u*u*u+g2*u<g1*g1)lower=u;else upper=u;
        }
        const double u=(lower+upper)/2,d=std::sqrt(4*u+(g1/u)*(g1/u)),
            probability=(1-g1/(u*d))/2,
            z0=-probability*d,z1=(1-probability)*d,residual=1-u;
        auto& c=component[m*ns+s];
        c.count=2;c.weight[0]=1-probability;c.weight[1]=probability;
        c.shift[0]=z0*std::sqrt(base_variance);
        c.shift[1]=z1*std::sqrt(base_variance);c.point_variance=u*base_variance;
        const double first=c.weight[0]*z0+c.weight[1]*z1,
            second=c.weight[0]*(z0*z0+residual)+c.weight[1]*(z1*z1+residual),
            third=c.weight[0]*(z0*z0*z0+3*z0*residual)+
                c.weight[1]*(z1*z1*z1+3*z1*residual),
            fourth=c.weight[0]*(std::pow(z0,4)+6*z0*z0*residual+3*residual*residual)+
                c.weight[1]*(std::pow(z1,4)+6*z1*z1*residual+3*residual*residual);
        moment_error=std::max({moment_error,std::abs(first),std::abs(second-1),
            std::abs(third-g1)/(1+std::abs(g1)),std::abs(fourth-3-g2)/(1+g2)});
    }
    for(int row=0;row<2*n;++row) if(positive[row]>0) for(int s=0;s<ns;++s)
        minimum_variance=std::min(minimum_variance,
            variance(row,s)-component[positive[row]*ns+s].point_variance);
    Rcpp::NumericMatrix L(summarize?0:n,summarize?0:nk),pair(n,np),out(n,3);
    std::vector<double> peak(n,R_NegInf);
    std::vector<int> active;
    for(int k=0;k<nk;++k) if(!summarize || w[k]>0) active.push_back(k);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int i=0;i<n;++i) {
        for(int p=0;p<np;++p) peak[i]=std::max(peak[i],
            (positive[i]>0?loglike(i,left[p]-1):0)+
            (positive[n+i]>0?loglike(n+i,right[p]-1):0));
        for(int p=0;p<np;++p) pair(i,p)=std::exp(
            (positive[i]>0?loglike(i,left[p]-1):0)+
            (positive[n+i]>0?loglike(n+i,right[p]-1):0)-peak[i]);
    }
    std::vector<double> partial(summarize?size_t(threads)*n*4:0,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int at=0;at<int(active.size());++at) {
        const int k=active[at],h=k/half,p=(k%half)/na,a=k%na,s0=left[p]-1,s1=right[p]-1;
        const double t=h?tv:0,center=h?slab_mean:0,common=A[a]/l2;
        int worker=0;
        #ifdef _OPENMP
        worker=omp_get_thread_num();
        #endif
        for(int i=0;i<n;++i) {
            const bool observed0=positive[i]>0,observed1=positive[n+i]>0;
            const auto& c0=component[positive[i]*ns+s0];
            const auto& c1=component[positive[n+i]*ns+s1];
            double density_sum=0,first_sum=0,second_sum=0;
            for(int j=0;j<c0.count;++j) for(int l=0;l<c1.count;++l) {
                double density=1,pm=center,pv=t;
                if(observed0 && observed1) {
                    const double v0=(variance(i,s0)-c0.point_variance)/(l2*l2),
                        v1=(variance(n+i,s1)-c1.point_variance)/(l2*l2),
                        y0=(mean(i,s0)+c0.shift[j])/l2-common+center/2,
                        y1=(mean(n+i,s1)+c1.shift[l])/l2-common-center/2,
                        s00=v0+av+t/4,s11=v1+av+t/4,s01=av-t/4,
                        det=v0*v1+(av+t/4)*(v0+v1)+av*t,
                        q=v1*y0*y0+v0*y1*y1+av*(y0-y1)*(y0-y1)+t/4*(y0+y1)*(y0+y1);
                    density=std::exp(-q/(2*det))/(2*pi*std::sqrt(det));
                    if(h) {
                        pm=center-t/2*(s11*y0-s01*y1)/det+t/2*(-s01*y0+s00*y1)/det;
                        pv=t*(v0*v1+av*(v0+v1))/det;
                    }
                } else if(observed0 || observed1) {
                    const int row=observed0?i:n+i,s=observed0?s0:s1,sign=observed0?-1:1;
                    const auto& c=observed0?c0:c1;
                    const int z=observed0?j:l;
                    const double v=(variance(row,s)-c.point_variance)/(l2*l2),
                        y=(mean(row,s)+c.shift[z])/l2-common-sign*center/2,
                        total=v+av+t/4;
                    density=std::exp(-y*y/(2*total))/std::sqrt(2*pi*total);
                    if(h){pm=center+sign*t*y/(2*total);pv=t*(v+av)/total;}
                }
                const double value=c0.weight[j]*c1.weight[l]*density;
                density_sum+=value;first_sum+=value*pm;second_sum+=value*(pv+pm*pm);
            }
            if(!summarize){L(i,k)=pair(i,p)*density_sum;continue;}
            double* q=partial.data()+size_t(worker)*n*4;
            const double factor=w[k]*pair(i,p);
            q[i]+=factor*density_sum;
            if(!h)q[n+i]+=factor*density_sum;
            q[2*n+i]+=factor*first_sum;q[3*n+i]+=factor*second_sum;
        }
    }
    if(summarize)for(int i=0;i<n;++i){
        double d=0,q=0,m=0,s=0;
        for(int t=0;t<threads;++t){const double* x=partial.data()+size_t(t)*n*4;
            d+=x[i];q+=x[n+i];m+=x[2*n+i];s+=x[3*n+i];}
        out(i,0)=m/d;out(i,1)=std::max(0.0,s/d-m*m/(d*d));out(i,2)=q/d;
    }
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("effect","variance","posterior_null");
    return Rcpp::List::create(Rcpp::Named("likelihood")=L,Rcpp::Named("genes")=out,
        Rcpp::Named("maximum_base_moment_error")=moment_error,
        Rcpp::Named("minimum_residual_variance")=minimum_variance);
}


// The flat contrast mean is a linear expectation. Reuse the already integrated
// group-wise location means instead of reconstructing a full contrast density.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector paired_mean_location_expectation_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::NumericMatrix mean,Rcpp::IntegerVector positive,Rcpp::IntegerVector left,
    Rcpp::IntegerVector right,Rcpp::NumericVector prior,int threads=1) {
    const int n=loglike.nrow()/2,np=left.size();const double logtwo=std::log(2.0);
    Rcpp::NumericVector out(n,NA_REAL);
    std::vector<double> logprior(np);
    for(int p=0;p<np;++p)logprior[p]=std::log(prior[p]);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int i=0;i<n;++i)if(positive[i]>0 && positive[n+i]>0){
        double peak=R_NegInf,total=0,first=0;
        for(int p=0;p<np;++p)if(prior[p]>0)peak=std::max(peak,
            logprior[p]+loglike(i,left[p]-1)+loglike(n+i,right[p]-1));
        for(int p=0;p<np;++p)if(prior[p]>0){
            const int s0=left[p]-1,s1=right[p]-1;
            const double mass=std::exp(logprior[p]+loglike(i,s0)+loglike(n+i,s1)-peak);
            total+=mass;first+=mass*(mean(n+i,s1)-mean(i,s0))/logtwo;
        }
        out[i]=first/total;
    }
    return out;
}

// The same Gaussian-mixture dictionary stored only where its likelihood or
// first two contrast moments exceed a row-wise, analytic error budget.
// State indices remain global; no gene, nuisance pair or A grid is removed.
struct LocationCurve {
    double log_amplitude,center,variance,intercept,slope,posterior_variance;
};
struct LocationSparseRow {
    std::vector<int> state;
    std::vector<double> value;
};
struct SparseLocationDictionary {
    int n,ns,np,na,nk,half,threads;
    int omitted_row=-1;
    double anchor=-.5; bool symmetric_anchor=false;
    double av,tv,slab_mean,l2=std::log(2.0),pi=3.14159265358979323846;
    double epsilon,moment_scale,moment_error=0,minimum_variance=R_PosInf;
    Rcpp::NumericMatrix loglike,mean,variance;
    Rcpp::IntegerVector positive,left,right;
    Rcpp::NumericVector A,log_sigma,shape;
    std::vector<LocationNormalPair> component;
    std::vector<double> peak,envelope,umbrella;
    std::vector<size_t> column_start;
    std::vector<int> column_row,constant_rows;
    std::vector<double> column_value;
    size_t entries=0;
    SparseLocationDictionary(Rcpp::NumericMatrix ll,Rcpp::NumericMatrix m,
        Rcpp::NumericMatrix v,Rcpp::IntegerVector count,Rcpp::IntegerVector lhs,
        Rcpp::IntegerVector rhs,Rcpp::NumericVector grid,double bandwidth,
        double center,double sd,int nt,Rcpp::NumericVector sig,Rcpp::NumericVector Q,
        double budget):n(ll.nrow()/2),ns(m.ncol()),np(lhs.size()),na(grid.size()),
        threads(nt),slab_mean(center),loglike(ll),mean(m),variance(v),positive(count),
        left(lhs),right(rhs),A(Rcpp::clone(grid)),log_sigma(sig),shape(Q),peak(n,R_NegInf),
        envelope(n,R_NegInf),umbrella(n,0) {
        const uint64_t states=uint64_t(2)*np*na;
        if(states>uint64_t(std::numeric_limits<int>::max()))
            Rcpp::stop("Joint location dictionary exceeds 32-bit state indexing");
        nk=int(states);half=nk/2;
        av=bandwidth*bandwidth/(l2*l2);tv=sd*sd;
        epsilon=budget/(16.0*n);
        moment_scale=std::max(1.0,std::sqrt(center*center+tv));
        for(double& x:A)x/=l2;
        const int maximum=*std::max_element(positive.begin(),positive.end());
        component.resize(size_t(maximum+1)*ns);
        for(int c=1;c<=maximum;++c)for(int s=0;s<ns;++s)if(Q[s]!=0) {
            const double a=c/(Q[s]*Q[s]),p1=R::trigamma(a),
                g1=-std::copysign(1.0,Q[s])*R::psigamma(a,2)/std::pow(p1,1.5),
                g2=R::psigamma(a,3)/(p1*p1),
                base=std::exp(2*sig[s])*p1/(Q[s]*Q[s]);
            double low=0,high=1;
            for(int j=0;j<80;++j) {
                const double u=(low+high)/2;
                if(2*u*u*u+g2*u<g1*g1)low=u;else high=u;
            }
            const double u=(low+high)/2,d=std::sqrt(4*u+(g1/u)*(g1/u)),
                p=(1-g1/(u*d))/2,z0=-p*d,z1=(1-p)*d,r=1-u;
            auto& cpt=component[size_t(c)*ns+s];
            cpt.count=2;cpt.weight[0]=1-p;cpt.weight[1]=p;
            cpt.shift[0]=z0*std::sqrt(base);cpt.shift[1]=z1*std::sqrt(base);
            cpt.point_variance=u*base;
            const double first=(1-p)*z0+p*z1,
                second=(1-p)*(z0*z0+r)+p*(z1*z1+r),
                third=(1-p)*(z0*z0*z0+3*z0*r)+p*(z1*z1*z1+3*z1*r),
                fourth=(1-p)*(std::pow(z0,4)+6*z0*z0*r+3*r*r)+
                    p*(std::pow(z1,4)+6*z1*z1*r+3*r*r);
            moment_error=std::max({moment_error,std::abs(first),std::abs(second-1),
                std::abs(third-g1)/(1+std::abs(g1)),std::abs(fourth-3-g2)/(1+g2)});
        }
        for(int i=0;i<2*n;++i)if(positive[i]>0)for(int s=0;s<ns;++s)
            minimum_variance=std::min(minimum_variance,
                variance(i,s)-component[size_t(positive[i])*ns+s].point_variance);
        if(!(minimum_variance>0))Rcpp::stop("Nonpositive conditional location variance");
    }
    int curves(int i,int p,int h,LocationCurve* out) const {
        return curves_pair(i,left[p]-1,right[p]-1,h,out);
    }
    int curves_pair(int i,int s0,int s1,int h,LocationCurve* out) const {
        if(!symmetric_anchor || !h)return curves_pair_anchor(i,s0,s1,h,out,anchor);
        const int first=curves_pair_anchor(i,s0,s1,h,out,0),
            second=curves_pair_anchor(i,s0,s1,h,out+first,-1);
        for(int j=0;j<first+second;++j)out[j].log_amplitude-=std::log(2.0);
        return first+second;
    }
    int curves_pair_anchor(int i,int s0,int s1,int h,LocationCurve* out,double first_coefficient) const {
        const bool o0=positive[i]>0,o1=positive[n+i]>0;
        const auto& c0=component[size_t(positive[i])*ns+s0];
        const auto& c1=component[size_t(positive[n+i])*ns+s1];
        const double t=h?tv:0,mu=h?slab_mean:0,a0=first_coefficient,a1=first_coefficient+1,
            lp=(o0?loglike(i,s0):0)+(o1?loglike(n+i,s1):0)-peak[i];
        int at=0;
        for(int j=0;j<c0.count;++j)for(int k=0;k<c1.count;++k) {
            const double lc=lp+std::log(c0.weight[j])+std::log(c1.weight[k]);
            LocationCurve c;
            if(o0 && o1) {
                const double v0=(variance(i,s0)-c0.point_variance)/(l2*l2),
                    v1=(variance(n+i,s1)-c1.point_variance)/(l2*l2),
                    u0=(mean(i,s0)+c0.shift[j])/l2,
                    u1=(mean(n+i,s1)+c1.shift[k])/l2,
                    vd=v0+v1+t,
                    information=a0*a0*v1+a1*a1*v0+av,
                    det=v0*v1+av*(v0+v1)+t*information,
                    contrast=u1-u0-mu;
                c.variance=det/vd;
                c.center=((v1+a1*t)*u0+(v0-a0*t)*u1-
                    (a0*v1+a1*v0)*mu)/vd;
                c.log_amplitude=lc-.5*std::log(2*pi*vd)-contrast*contrast/(2*vd)-
                    .5*std::log(2*pi*c.variance);
                c.intercept=mu+t/det*((a0*v1-av)*u0+
                    (a1*v0+av)*u1-information*mu);
                c.slope=-t*(a0*v1+a1*v0)/det;
                c.posterior_variance=t*(v0*v1+av*(v0+v1))/det;
            } else {
                const int row=o0?i:n+i,s=o0?s0:s1,z=o0?j:k;
                const double coefficient=o0?a0:a1;
                const auto& cc=o0?c0:c1;
                const double v=(variance(row,s)-cc.point_variance)/(l2*l2),
                    u=(mean(row,s)+cc.shift[z])/l2;
                c.variance=v+av+coefficient*coefficient*t;
                c.center=u-coefficient*mu;
                c.log_amplitude=lc-.5*std::log(2*pi*c.variance);
                c.intercept=mu+coefficient*t*c.center/c.variance;
                c.slope=-coefficient*t/c.variance;
                c.posterior_variance=t*(v+av)/c.variance;
            }
            out[at++]=c;
        }
        return at;
    }
    double grid_peak(const LocationCurve& c) const {
        const auto a=std::lower_bound(A.begin(),A.end(),c.center);
        int index=std::min(na-1,int(a-A.begin()));
        double distance=std::abs(A[index]-c.center);
        if(index>0)distance=std::min(distance,std::abs(A[index-1]-c.center));
        return c.log_amplitude-distance*distance/(2*c.variance);
    }
    // Component omission <= epsilon/(2*C) for each scaled moment.
    // The aggregate state omission contributes the other epsilon/2.
    void block(int i,int p,int h,std::vector<double>& v,std::vector<double>& m1,
        std::vector<double>& m2,int& begin,int& end) const {
        block_pair(i,left[p]-1,right[p]-1,h,v,m1,m2,begin,end);
    }
    void block_pair(int i,int s0,int s1,int h,std::vector<double>& v,std::vector<double>& m1,
        std::vector<double>& m2,int& begin,int& end) const {
        std::fill(v.begin(),v.end(),0);std::fill(m1.begin(),m1.end(),0);
        std::fill(m2.begin(),m2.end(),0);begin=na;end=0;
        LocationCurve cc[8];const int nc=curves_pair(i,s0,s1,h,cc);
        for(int j=0;j<nc;++j) {
            const auto& c=cc[j];
            const double max_mean=std::max(std::abs(c.intercept+c.slope*A[0]),
                std::abs(c.intercept+c.slope*A[na-1])),
                bound=std::max({1.0,max_mean/moment_scale,
                    (c.posterior_variance+max_mean*max_mean)/(moment_scale*moment_scale)}),
                threshold=std::log(epsilon/(2*nc*bound)),
                radius2=2*c.variance*(c.log_amplitude-envelope[i]-threshold);
            if(!(radius2>=0))continue;
            const double radius=std::sqrt(radius2);
            int low=std::max(0,int(std::lower_bound(A.begin(),A.end(),c.center-radius)-A.begin())-1);
            int high=std::min(na,int(std::upper_bound(A.begin(),A.end(),c.center+radius)-A.begin())+1);
            begin=std::min(begin,low);end=std::max(end,high);
            for(int a=low;a<high;++a) {
                const double z=A[a]-c.center,value=std::exp(
                    c.log_amplitude-envelope[i]-z*z/(2*c.variance)),
                    pm=c.intercept+c.slope*A[a];
                v[a]+=value;m1[a]+=value*pm;
                m2[a]+=value*(c.posterior_variance+pm*pm);
            }
        }
    }
    bool retained(double v,double m1,double m2) const {
        return std::max({v,std::abs(m1)/moment_scale,m2/(moment_scale*moment_scale)})>epsilon/2;
    }
    void build(bool fixed_scaling=false) {
        // Count and then fill column postings directly. Virtual row partitions
        // keep every column sorted without storing a second row-oriented copy.
        std::vector<size_t> cursor(size_t(threads)*nk,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int partition=0;partition<threads;++partition) {
            const int first=partition*(n/threads)+std::min(partition,n%threads),
                last=(partition+1)*(n/threads)+std::min(partition+1,n%threads);
            size_t* count=cursor.data()+size_t(partition)*nk;
            std::vector<double> v(na),m1(na),m2(na);int low,high;
            for(int i=first;i<last;++i) {
                if(positive[i]==0 && positive[n+i]==0){envelope[i]=0;umbrella[i]=1;continue;}
                if(!fixed_scaling) {
                for(int p=0;p<np;++p)peak[i]=std::max(peak[i],
                    (positive[i]>0?loglike(i,left[p]-1):0)+
                    (positive[n+i]>0?loglike(n+i,right[p]-1):0));
                for(int h=0;h<2;++h)for(int p=0;p<np;++p) {
                    LocationCurve cc[8];const int nc=curves(i,p,h,cc);
                    double top=R_NegInf,values[8];
                    for(int j=0;j<nc;++j){values[j]=grid_peak(cc[j]);top=std::max(top,values[j]);}
                    double sum=0;for(int j=0;j<nc;++j)sum+=std::exp(values[j]-top);
                    envelope[i]=std::max(envelope[i],top+std::log(sum));
                }
                }
                for(int h=0;h<2;++h)for(int p=0;p<np;++p) {
                    block(i,p,h,v,m1,m2,low,high);
                    for(int a=low;a<high;++a)if(retained(v[a],m1[a],m2[a])) {
                        ++count[h*half+p*na+a];umbrella[i]+=v[a]/nk;
                    }
                }
            }
        }
        column_start.assign(size_t(nk)+1,0);
        for(int k=0;k<nk;++k) {
            size_t at=column_start[k];
            for(int partition=0;partition<threads;++partition) {
                const size_t index=size_t(partition)*nk+k,count=cursor[index];
                cursor[index]=at;at+=count;
            }
            column_start[k+1]=at;
        }
        entries=column_start[nk];column_row.resize(entries);column_value.resize(entries);
        for(int i=0;i<n;++i)if(positive[i]==0 && positive[n+i]==0)constant_rows.push_back(i);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int partition=0;partition<threads;++partition) {
            const int first=partition*(n/threads)+std::min(partition,n%threads),
                last=(partition+1)*(n/threads)+std::min(partition+1,n%threads);
            size_t* next=cursor.data()+size_t(partition)*nk;
            std::vector<double> v(na),m1(na),m2(na);int low,high;
            for(int i=first;i<last;++i) {
                if(positive[i]==0 && positive[n+i]==0)continue;
                for(int h=0;h<2;++h)for(int p=0;p<np;++p) {
                    block(i,p,h,v,m1,m2,low,high);
                    for(int a=low;a<high;++a)if(retained(v[a],m1[a],m2[a])) {
                        const int state=h*half+p*na+a;const size_t at=next[state]++;
                        column_row[at]=i;column_value[at]=v[a];
                    }
                }
            }
        }
    }
    // Column postings give contiguous gradient scans and direct active-column
    // assembly. Row order within every column follows the original gene order.
    double value(int i,int k) const {
        if(i==omitted_row || (positive[i]==0 && positive[n+i]==0))return 1;
        const auto first=column_row.begin()+column_start[k],last=column_row.begin()+column_start[k+1];
        const auto at=std::lower_bound(first,last,i);
        return at!=last && *at==i ? column_value[at-column_row.begin()] : 0;
    }
    void fill_column(int k,double* output) const {
        std::fill(output,output+n,0);
        for(size_t j=column_start[k];j<column_start[k+1];++j)
            output[column_row[j]]=column_value[j];
        for(int i:constant_rows)output[i]=1;
    }
    std::vector<double> gradient(const std::vector<double>& density,std::vector<int>* dominant_row=nullptr) const {
        std::vector<double> factor(n),constant(threads,0),out(nk,0);
        std::vector<int> ends(threads);
        for(int i=0;i<n;++i)factor[i]=i==omitted_row?0:1/(n*density[i]);
        for(int t=0;t<threads;++t) {
            const int first=t*(n/threads)+std::min(t,n%threads);
            ends[t]=(t+1)*(n/threads)+std::min(t+1,n%threads);
            for(int i=first;i<ends[t];++i)
                if(i==omitted_row || (positive[i]==0 && positive[n+i]==0))
                    constant[t]+=1/(n*density[i]);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            size_t at=column_start[k];const size_t end=column_start[k+1];
            double largest=-1;int dominant=0;
            for(int t=0;t<threads;++t) {
                double score=0;
                for(;at<end && column_row[at]<ends[t];++at) {
                    const double contribution=column_value[at]*factor[column_row[at]];
                    score+=contribution;
                    if(contribution>largest){largest=contribution;dominant=column_row[at];}
                }
                out[k]+=score+constant[t];
            }
            if(dominant_row)(*dominant_row)[k]=dominant;
        }
        return out;
    }
};

// A constant likelihood row has zero log-likelihood on the simplex. Keep the
// original CSC postings and restore this view when the solve returns.
struct LocationOmissionScope {
    SparseLocationDictionary& op;int target;double previous_umbrella=0;bool inserted=false;
    LocationOmissionScope(SparseLocationDictionary& dictionary,int one_based):
        op(dictionary),target(one_based-1) {
        if(one_based==0)return;
        if(target<0 || target>=op.n)Rcpp::stop("Omitted gene index is outside the dictionary");
        previous_umbrella=op.umbrella[target];op.umbrella[target]=1;op.omitted_row=target;
        const auto at=std::lower_bound(op.constant_rows.begin(),op.constant_rows.end(),target);
        inserted=at==op.constant_rows.end() || *at!=target;
        if(inserted)op.constant_rows.insert(at,target);
    }
    ~LocationOmissionScope() {
        if(target<0)return;
        if(inserted)op.constant_rows.erase(std::lower_bound(op.constant_rows.begin(),op.constant_rows.end(),target));
        op.umbrella[target]=previous_umbrella;op.omitted_row=-1;
    }
};

struct LocationClassScope {
    SparseLocationDictionary& op;int fitted_class;
    std::vector<double> previous;
    LocationClassScope(SparseLocationDictionary& dictionary,int selected):op(dictionary),fitted_class(selected) {
        if(selected<0)return;
        if(selected>1)Rcpp::stop("Location class must be zero, one, or unrestricted");
        previous=op.umbrella;
        std::vector<double> partial(size_t(op.threads)*op.n,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(op.threads)
        #endif
        for(int j=selected*op.half;j<(selected+1)*op.half;++j) {
            int thread=0;
            #ifdef _OPENMP
            thread=omp_get_thread_num();
            #endif
            double* total=partial.data()+size_t(thread)*op.n;
            for(size_t at=op.column_start[j];at<op.column_start[j+1];++at)
                total[op.column_row[at]]+=op.column_value[at]/op.half;
        }
        std::fill(op.umbrella.begin(),op.umbrella.end(),0);
        for(int i=0;i<op.n;++i)for(int t=0;t<op.threads;++t)op.umbrella[i]+=partial[size_t(t)*op.n+i];
        for(int i:op.constant_rows)op.umbrella[i]=1;
    }
    ~LocationClassScope(){if(fitted_class>=0)op.umbrella.swap(previous);}
};

// [[Rcpp::export(rng = false)]]
Rcpp::List sparse_location_dictionary_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::NumericMatrix mean,Rcpp::NumericMatrix variance,Rcpp::IntegerVector positive,
    Rcpp::IntegerVector left,Rcpp::IntegerVector right,Rcpp::NumericVector A,
    double A_sd,double slab_mean,double slab_sd,int threads,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,double error_budget=1e-12) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(new SparseLocationDictionary(loglike,mean,
        variance,positive,left,right,A,A_sd,slab_mean,slab_sd,threads,log_sigma,shape,error_budget),true);
    ptr->build();
    return Rcpp::List::create(Rcpp::Named("dictionary")=ptr,
        Rcpp::Named("states")=ptr->nk,Rcpp::Named("entries")=double(ptr->entries),
        Rcpp::Named("storage_bytes")=double(ptr->entries)*(sizeof(int)+sizeof(double))+
            double(ptr->column_start.size())*sizeof(size_t),
        Rcpp::Named("dense_bytes")=double(ptr->n)*ptr->nk*sizeof(double),
        Rcpp::Named("row_error_budget")=ptr->epsilon,
        Rcpp::Named("log_row_envelope")=Rcpp::wrap(ptr->envelope),
        Rcpp::Named("maximum_base_moment_error")=ptr->moment_error,
        Rcpp::Named("minimum_residual_variance")=ptr->minimum_variance);
}

// Reconstruct only the saved nonzero support at the original row scale.
// A fixed-support leave-one-out diagnostic need not rebuild the full dictionary.
// [[Rcpp::export(rng = false)]]
Rcpp::List sparse_location_saved_support_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::NumericMatrix mean,Rcpp::NumericMatrix variance,Rcpp::IntegerVector positive,
    Rcpp::IntegerVector left,Rcpp::IntegerVector right,Rcpp::NumericVector A,
    double A_sd,double slab_mean,double slab_sd,int threads,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,
    Rcpp::NumericVector log_row_envelope,Rcpp::NumericVector saved_weights,
    double error_budget=1e-12) {
    SparseLocationDictionary op(loglike,mean,variance,positive,left,right,A,A_sd,
        slab_mean,slab_sd,threads,log_sigma,shape,error_budget);
    std::copy(log_row_envelope.begin(),log_row_envelope.end(),op.envelope.begin());
    std::vector<int> selected;
    for(int k=0;k<op.nk;++k)if(saved_weights[k]>0)selected.push_back(k);
    const int count=selected.size();
    std::vector<std::vector<std::pair<int,int>>> groups(size_t(2)*op.np);
    std::vector<int> active;
    Rcpp::NumericVector weights(count);Rcpp::IntegerVector state(count),classes(count);
    for(int j=0;j<count;++j){const int k=selected[j],h=k/op.half,p=(k%op.half)/op.na,a=k%op.na;
        groups[h*op.np+p].push_back({a,j});weights[j]=saved_weights[k];
        state[j]=k+1;classes[j]=h;}
    for(int g=0;g<2*op.np;++g)if(!groups[g].empty())active.push_back(g);
    Rcpp::NumericMatrix likelihood(op.n,count);Rcpp::NumericVector density(op.n),null(op.n);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int i=0;i<op.n;++i){
        const bool constant=positive[i]==0 && positive[op.n+i]==0;
        if(constant){double total=0,pn=0;
            for(int j=0;j<count;++j){likelihood(i,j)=1;total+=weights[j];
                if(classes[j]==0)pn+=weights[j];}density[i]=total;null[i]=pn/total;continue;}
        for(int p=0;p<op.np;++p)op.peak[i]=std::max(op.peak[i],
            (positive[i]>0?loglike(i,left[p]-1):0)+
            (positive[op.n+i]>0?loglike(op.n+i,right[p]-1):0));
        double total=0,pn=0;
        for(int g:active){const int h=g/op.np,p=g%op.np;
            const auto& knots=groups[g];const int size=knots.size();
            std::vector<double> v(size,0),m1(size,0),m2(size,0);
            LocationCurve cc[8];const int nc=op.curves(i,p,h,cc);
            for(int c=0;c<nc;++c){const auto& curve=cc[c];
                const double maximum_mean=std::max(std::abs(curve.intercept+curve.slope*op.A[0]),
                    std::abs(curve.intercept+curve.slope*op.A[op.na-1])),
                    bound=std::max({1.0,maximum_mean/op.moment_scale,
                        (curve.posterior_variance+maximum_mean*maximum_mean)/(op.moment_scale*op.moment_scale)}),
                    threshold=std::log(op.epsilon/(2*nc*bound)),
                    radius2=2*curve.variance*(curve.log_amplitude-op.envelope[i]-threshold);
                if(!(radius2>=0))continue;
                const double radius=std::sqrt(radius2);
                const int low=std::max(0,int(std::lower_bound(op.A.begin(),op.A.end(),curve.center-radius)-op.A.begin())-1),
                    high=std::min(op.na,int(std::upper_bound(op.A.begin(),op.A.end(),curve.center+radius)-op.A.begin())+1);
                for(int j=0;j<size;++j){const int a=knots[j].first;if(a<low || a>=high)continue;
                    const double z=op.A[a]-curve.center,value=std::exp(
                        curve.log_amplitude-op.envelope[i]-z*z/(2*curve.variance)),
                        pm=curve.intercept+curve.slope*op.A[a];
                    v[j]+=value;m1[j]+=value*pm;m2[j]+=value*(curve.posterior_variance+pm*pm);}
            }
            for(int j=0;j<size;++j)if(op.retained(v[j],m1[j],m2[j])){
                const int column=knots[j].second;likelihood(i,column)=v[j];
                const double mass=v[j]*weights[column];total+=mass;if(h==0)pn+=mass;}
        }
        density[i]=total;null[i]=pn/total;
    }
    return Rcpp::List::create(Rcpp::Named("likelihood")=likelihood,
        Rcpp::Named("state")=state,Rcpp::Named("class")=classes,
        Rcpp::Named("weights")=weights,Rcpp::Named("density")=density,
        Rcpp::Named("posterior_null")=null,
        Rcpp::Named("original_states")=op.nk,
        Rcpp::Named("reconstructed_states")=count,
        Rcpp::Named("minimum_residual_variance")=op.minimum_variance);
}

// Price omitted pairs from the existing atom Cartesian product. The Gaussian
// peak bound screens every pair without constructing its gene-by-mean matrix.
// [[Rcpp::export(rng = false)]]
Rcpp::List sparse_location_missing_pair_prices_cpp(SEXP dictionary,
    Rcpp::NumericVector density,int threads=1,int maximum_pairs=64,
    double tolerance=1e-6,double maximum_seconds=60,int priced_class=-1,int start_rank=0) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);const auto& op=*ptr;
    const int n=op.n,ns=op.ns,na=op.na;
    const int first_class=priced_class<0?0:priced_class,last_class=priced_class<0?2:priced_class+1;
    const auto started=std::chrono::steady_clock::now();
    Rcpp::NumericMatrix left(n,ns),right(n,ns),weighted(n,ns),upper(ns,ns);
    std::vector<double> factor(n);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int i=0;i<n;++i) {
        const bool o0=op.positive[i]>0,o1=op.positive[n+i]>0;
        double peak0=o0?R_NegInf:0,peak1=o1?R_NegInf:0;
        double minimum0=R_PosInf,minimum1=R_PosInf;
        for(int s=0;s<ns;++s) {
            if(o0) {
                peak0=std::max(peak0,op.loglike(i,s));
                const double v=(op.variance(i,s)-op.component[size_t(op.positive[i])*ns+s].point_variance)/(op.l2*op.l2);
                minimum0=std::min(minimum0,v);left(i,s)=v;
            }
            if(o1) {
                peak1=std::max(peak1,op.loglike(n+i,s));
                const double v=(op.variance(n+i,s)-op.component[size_t(op.positive[n+i])*ns+s].point_variance)/(op.l2*op.l2);
                minimum1=std::min(minimum1,v);right(i,s)=v;
            }
        }
        const double dependence=o0 && o1?std::sqrt((op.av+minimum0)*(op.av+minimum1)/
            (op.av*(minimum0+minimum1)+minimum0*minimum1)):1;
        const double scale=o0 || o1?std::exp(peak0+peak1-op.peak[i]-op.envelope[i])*dependence:1;
        factor[i]=1/(n*density[i]);
        for(int s=0;s<ns;++s) {
            left(i,s)=o0?std::exp(op.loglike(i,s)-peak0)/std::sqrt(2*op.pi*(op.av+left(i,s))):1;
            right(i,s)=o1?std::exp(op.loglike(n+i,s)-peak1)/std::sqrt(2*op.pi*(op.av+right(i,s))):1;
            weighted(i,s)=left(i,s)*scale*factor[i];
        }
    }
    const Eigen::Map<const Eigen::MatrixXd> lhs(weighted.begin(),n,ns),rhs(right.begin(),n,ns);
    Eigen::Map<Eigen::MatrixXd> out(upper.begin(),ns,ns);
    Eigen::setNbThreads(threads);out.noalias()=lhs.transpose()*rhs;
    const double machine=std::numeric_limits<double>::epsilon(),gamma=(2.0*n+32)*machine/(1-(2.0*n+32)*machine);
    for(double& value:upper)value/=(1-gamma);
    std::vector<char> present(size_t(ns)*ns,0);
    for(int p=0;p<op.np;++p)present[size_t(op.left[p]-1)*ns+op.right[p]-1]=1;
    std::vector<int> order;
    for(int a=0;a<ns;++a)for(int b=0;b<ns;++b)
        if(!present[size_t(a)*ns+b])order.push_back(a*ns+b);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return upper(a/ns,a%ns)>upper(b/ns,b%ns);});
    const double row_cutoff=tolerance*tolerance/n,step=(op.A[na-1]-op.A[0])/(na-1),
        rounding_terms=8.0*na*na+4.0*n+64,
        recurrence_rounding=rounding_terms*machine/(1-rounding_terms*machine);
    auto price=[&](int a,int b,double* omitted=nullptr) {
        std::vector<double> partial(size_t(threads)*2*na,0),score(2*na,0);
        std::vector<double> dropped(threads,0);
        #ifdef _OPENMP
        #pragma omp parallel num_threads(threads)
        #endif
        {
            int thread=0;
            #ifdef _OPENMP
            thread=omp_get_thread_num();
            #endif
            double* sum=partial.data()+size_t(thread)*2*na;
            std::vector<double> v(na),m1(na),m2(na);int low,high;
            #ifdef _OPENMP
            #pragma omp for schedule(static)
            #endif
            for(int i=0;i<n;++i) {
                if(op.positive[i]==0 && op.positive[n+i]==0) {
                    for(int j=0;j<2*na;++j)sum[j]+=factor[i];
                    continue;
                }
                const double row_bound=weighted(i,a)*right(i,b);
                if(row_bound<=row_cutoff){dropped[thread]+=row_bound;continue;}
                for(int h=first_class;h<last_class;++h) {
                    std::fill(v.begin(),v.end(),0);std::fill(m1.begin(),m1.end(),0);std::fill(m2.begin(),m2.end(),0);
                    low=na;high=0;LocationCurve curves[8];const int nc=op.curves_pair(i,a,b,h,curves);
                    for(int c=0;c<nc;++c) {
                        const auto& curve=curves[c];
                        const double maximum_mean=std::max(std::abs(curve.intercept+curve.slope*op.A[0]),
                            std::abs(curve.intercept+curve.slope*op.A[na-1])),
                            moment_bound=std::max({1.0,maximum_mean/op.moment_scale,
                                (curve.posterior_variance+maximum_mean*maximum_mean)/(op.moment_scale*op.moment_scale)}),
                            threshold=std::log(op.epsilon/(2*nc*moment_bound)),
                            radius2=2*curve.variance*(curve.log_amplitude-op.envelope[i]-threshold);
                        if(!(radius2>=0))continue;
                        const double radius=std::sqrt(radius2);
                        const int lo=std::max(0,int(std::lower_bound(op.A.begin(),op.A.end(),curve.center-radius)-op.A.begin())-1),
                            hi=std::min(na,int(std::upper_bound(op.A.begin(),op.A.end(),curve.center+radius)-op.A.begin())+1);
                        low=std::min(low,lo);high=std::max(high,hi);
                        const bool recurrence=step*step<=curve.variance;
                        const double first=lo<hi?op.A[lo]-curve.center:0,
                            decrease=recurrence?std::exp(-step*step/curve.variance):1;
                        double value=recurrence?std::exp(curve.log_amplitude-op.envelope[i]-first*first/(2*curve.variance)):0,
                            ratio=recurrence?std::exp(-step*(first+step/2)/curve.variance):1;
                        for(int j=lo;j<hi;++j) {
                            if(!recurrence){const double z=op.A[j]-curve.center;value=std::exp(curve.log_amplitude-op.envelope[i]-z*z/(2*curve.variance));}
                            const double pm=curve.intercept+curve.slope*op.A[j];
                            v[j]+=value;m1[j]+=value*pm;m2[j]+=value*(curve.posterior_variance+pm*pm);
                            if(recurrence){value*=ratio;ratio*=decrease;}
                        }
                    }
                    for(int j=low;j<high;++j)if(op.retained(v[j],m1[j],m2[j]))sum[h*na+j]+=v[j]*factor[i];
                }
            }
        }
        for(int t=0;t<threads;++t)for(int j=0;j<2*na;++j)score[j]+=partial[size_t(t)*2*na+j];
        if(omitted)*omitted=std::accumulate(dropped.begin(),dropped.end(),0.0)/(1-gamma);
        return score;
    };
    const auto known=price(op.left[0]-1,op.right[0]-1);
    const std::vector<double> fitted_density(density.begin(),density.end());
    const auto existing=op.gradient(fitted_density);
    double price_error=0;
    for(int h=first_class;h<last_class;++h)for(int a=0;a<na;++a)
        price_error=std::max(price_error,std::abs(known[h*na+a]-existing[h*op.half+a]));
    std::vector<int> lhs_index,rhs_index,knot0,knot1;
    std::vector<double> bounds,null_score,alternative_score,error_bound;
    int candidate_pairs=0;
    for(int state:order)if(upper(state/ns,state%ns)>1+tolerance)++candidate_pairs;
    int at=std::min(std::max(0,start_rank),int(order.size()));
    const int stop=std::min(int(order.size()),at+maximum_pairs);
    for(;at<stop;++at) {
        const int state=order[at],a=state/ns,b=state%ns;
        if(upper(a,b)<=1+tolerance)break;
        if(maximum_seconds>0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>=maximum_seconds)break;
        double omitted=0;Rcpp::checkUserInterrupt();const auto score=price(a,b,&omitted);
        const int j0=std::max_element(score.begin(),score.begin()+na)-score.begin(),
            j1=std::max_element(score.begin()+na,score.end())-score.begin();
        lhs_index.push_back(a+1);rhs_index.push_back(b+1);bounds.push_back(upper(a,b));
        null_score.push_back(priced_class==1?NA_REAL:score[j0]);alternative_score.push_back(priced_class==0?NA_REAL:score[j1]);
        knot0.push_back(priced_class==1?NA_INTEGER:j0+1);knot1.push_back(priced_class==0?NA_INTEGER:j1-na+1);
        const double maximum=priced_class==0?score[j0]:(priced_class==1?score[j1]:std::max(score[j0],score[j1]));
        error_bound.push_back(omitted+recurrence_rounding*maximum);
    }
    const double remaining=at<int(order.size())?upper(order[at]/ns,order[at]%ns):0;
    return Rcpp::List::create(Rcpp::Named("prices")=Rcpp::DataFrame::create(
        Rcpp::Named("left")=lhs_index,Rcpp::Named("right")=rhs_index,Rcpp::Named("upper_bound")=bounds,
        Rcpp::Named("null_score")=null_score,Rcpp::Named("alternative_score")=alternative_score,
        Rcpp::Named("null_knot")=knot0,Rcpp::Named("alternative_knot")=knot1,
        Rcpp::Named("gradient_error_bound")=error_bound),
        Rcpp::Named("remaining_upper_bound")=remaining,Rcpp::Named("complete")=(remaining<=1+tolerance),
        Rcpp::Named("retained_pair_price_error")=price_error,
        Rcpp::Named("retained_dictionary_KKT")=*std::max_element(existing.begin()+first_class*op.half,existing.begin()+last_class*op.half)-1,
        Rcpp::Named("atom_pairs")=double(ns)*ns,Rcpp::Named("retained_pairs")=op.np,
        Rcpp::Named("candidate_pairs")=candidate_pairs,Rcpp::Named("next_rank")=at,
        Rcpp::Named("priced_class")=priced_class);
}

// Append only newly priced pairs. Keep the original row scale and CSC values,
// and remap the old alternative-class indices to the enlarged half dictionary.
// [[Rcpp::export(rng = false)]]
Rcpp::List append_sparse_location_pairs_cpp(SEXP dictionary,Rcpp::IntegerVector added_left,
    Rcpp::IntegerVector added_right,Rcpp::NumericVector weights) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);auto& op=*ptr;
    const int old_half=op.half,old_states=op.nk,added_pairs=added_left.size(),
        added_half=added_pairs*op.na,new_half=old_half+added_half;
    SparseLocationDictionary extra(op.loglike,op.mean,op.variance,op.positive,
        added_left,added_right,op.A,std::sqrt(op.av)*op.l2,op.slab_mean,std::sqrt(op.tv),
        op.threads,op.log_sigma,op.shape,op.epsilon*16*op.n);
    extra.anchor=op.anchor;extra.symmetric_anchor=op.symmetric_anchor;
    extra.A=Rcpp::clone(op.A);extra.av=op.av;extra.tv=op.tv;
    extra.epsilon=op.epsilon;extra.moment_scale=op.moment_scale;
    extra.peak=op.peak;extra.envelope=op.envelope;extra.build(true);
    const auto previous=op.column_start;
    const size_t old_entries=op.entries,old_split=previous[old_half],
        extra_split=extra.column_start[added_half],total_entries=old_entries+extra.entries;
    op.column_row.reserve(total_entries);op.column_value.reserve(total_entries);
    op.column_row.resize(total_entries);op.column_value.resize(total_entries);
    std::copy_backward(op.column_row.begin()+old_split,op.column_row.begin()+old_entries,
        op.column_row.begin()+old_entries+extra_split);
    std::copy_backward(op.column_value.begin()+old_split,op.column_value.begin()+old_entries,
        op.column_value.begin()+old_entries+extra_split);
    std::copy(extra.column_row.begin(),extra.column_row.begin()+extra_split,op.column_row.begin()+old_split);
    std::copy(extra.column_value.begin(),extra.column_value.begin()+extra_split,op.column_value.begin()+old_split);
    std::copy(extra.column_row.begin()+extra_split,extra.column_row.end(),op.column_row.begin()+old_entries+extra_split);
    std::copy(extra.column_value.begin()+extra_split,extra.column_value.end(),op.column_value.begin()+old_entries+extra_split);
    std::vector<size_t> starts(size_t(2)*new_half+1);
    for(int j=0;j<old_half;++j)starts[j]=previous[j];
    for(int j=0;j<added_half;++j)starts[old_half+j]=old_split+extra.column_start[j];
    for(int j=0;j<old_half;++j)starts[new_half+j]=previous[old_half+j]+extra_split;
    for(int j=0;j<added_half;++j)starts[new_half+old_half+j]=old_entries+extra.column_start[added_half+j];
    starts.back()=total_entries;op.column_start.swap(starts);
    Rcpp::IntegerVector left(op.np+added_pairs),right(op.np+added_pairs);
    std::copy(op.left.begin(),op.left.end(),left.begin());std::copy(op.right.begin(),op.right.end(),right.begin());
    std::copy(added_left.begin(),added_left.end(),left.begin()+op.np);
    std::copy(added_right.begin(),added_right.end(),right.begin()+op.np);
    Rcpp::NumericVector initial(2*new_half);
    std::copy(weights.begin(),weights.begin()+old_half,initial.begin());
    std::copy(weights.begin()+old_half,weights.end(),initial.begin()+new_half);
    for(int i=0;i<op.n;++i)op.umbrella[i]=(old_states*op.umbrella[i]+extra.nk*extra.umbrella[i])/(2*new_half);
    op.left=left;op.right=right;op.np+=added_pairs;op.half=new_half;op.nk=2*new_half;op.entries=total_entries;
    std::vector<double> partial(size_t(op.threads)*op.n*2,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(op.threads)
    #endif
    for(int j=0;j<op.nk;++j)if(initial[j]>0) {
        int thread=0;
        #ifdef _OPENMP
        thread=omp_get_thread_num();
        #endif
        double* total=partial.data()+size_t(thread)*op.n*2;
        for(size_t at=op.column_start[j];at<op.column_start[j+1];++at) {
            const int i=op.column_row[at];const double mass=initial[j]*op.column_value[at];
            total[i]+=mass;if(j<op.half)total[op.n+i]+=mass;
        }
    }
    Rcpp::NumericVector density(op.n),null(op.n);
    const double total_mass=std::accumulate(initial.begin(),initial.end(),0.0),
        null_mass=std::accumulate(initial.begin(),initial.begin()+op.half,0.0);
    for(int i=0;i<op.n;++i) {
        double d=0,q=0;
        for(int t=0;t<op.threads;++t){const double* total=partial.data()+size_t(t)*op.n*2;d+=total[i];q+=total[op.n+i];}
        if(op.positive[i]==0 && op.positive[op.n+i]==0){d=total_mass;q=null_mass;}
        density[i]=d;null[i]=q/d;
    }
    return Rcpp::List::create(Rcpp::Named("weights")=initial,Rcpp::Named("left")=left,
        Rcpp::Named("right")=right,Rcpp::Named("added_pairs")=added_pairs,
        Rcpp::Named("density")=density,Rcpp::Named("posterior_null")=null,
        Rcpp::Named("added_entries")=double(extra.entries),Rcpp::Named("states")=op.nk,
        Rcpp::Named("storage_bytes")=double(op.entries)*(sizeof(int)+sizeof(double))+
            double(op.column_start.size())*sizeof(size_t));
}

// Conditional MAP states identify where the fitted observation operator places
// its mass. Return both class kernels at each selected physical pair and knot.
// [[Rcpp::export(rng = false)]]
Rcpp::List sparse_location_class_map_cpp(SEXP dictionary,Rcpp::NumericVector weights,
    Rcpp::NumericVector density,int threads=1) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);const auto& op=*ptr;
    std::vector<double> best(size_t(threads)*2*op.n,0);
    std::vector<int> index(size_t(threads)*2*op.n,-1);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int j=0;j<op.nk;++j)if(weights[j]>0) {
        int thread=0;
        #ifdef _OPENMP
        thread=omp_get_thread_num();
        #endif
        const int h=j/op.half;const size_t offset=(size_t(thread)*2+h)*op.n;
        for(size_t at=op.column_start[j];at<op.column_start[j+1];++at) {
            const int i=op.column_row[at];const double value=weights[j]*op.column_value[at];
            if(value>best[offset+i]){best[offset+i]=value;index[offset+i]=j;}
        }
    }
    const int maximum0=std::max_element(weights.begin(),weights.begin()+op.half)-weights.begin(),
        maximum1=std::max_element(weights.begin()+op.half,weights.end())-weights.begin();
    Rcpp::IntegerMatrix labels(2,op.n);Rcpp::NumericMatrix approximate(op.n,4),selected_weights(op.n,4),mass(op.n,2);
    for(int i=0;i<op.n;++i)for(int h=0;h<2;++h) {
        double value=0;int state=-1;
        for(int t=0;t<threads;++t) {
            const size_t at=(size_t(t)*2+h)*op.n+i;
            if(best[at]>value || (best[at]==value && index[at]>=0 && (state<0 || index[at]<state))) {
                value=best[at];state=index[at];
            }
        }
        if(op.positive[i]==0 && op.positive[op.n+i]==0) {
            state=h?maximum1:maximum0;value=weights[state];
        }
        labels(h,i)=state+1;mass(i,h)=value/density[i];
        if(state<0)continue;
        const int physical=state%op.half;
        for(int c=0;c<2;++c) {
            const int j=physical+c*op.half;
            approximate(i,2*h+c)=op.value(i,j);selected_weights(i,2*h+c)=weights[j];
        }
    }
    return Rcpp::List::create(Rcpp::Named("labels")=labels,Rcpp::Named("approximate")=approximate,
        Rcpp::Named("selected_weights")=selected_weights,Rcpp::Named("posterior_mass")=mass);
}

Rcpp::List population_weights_cpp(Rcpp::NumericMatrix likelihood,Rcpp::NumericVector start,
    double restricted_tolerance,int threads);

static void location_exchange(const Rcpp::NumericMatrix& L,std::vector<double>& w,
    std::vector<double>& d,int enter,int leave) {
    const int n=L.nrow(),k=L.ncol();const double mass=w[leave];
    std::vector<double> end(n,0);
    for(int j=0;j<k;++j) {
        const double z=j==leave?0:(w[j]+(j==enter?mass:0));
        if(z>0)for(int i=0;i<n;++i)end[i]+=z*L(i,j);
    }
    double low=0,high=1,derivative=0;
    for(int i=0;i<n;++i)derivative+=(L(i,enter)-L(i,leave))/end[i];
    if(derivative<=0)for(int search=0;search<50;++search) {
        const double step=(low+high)/2;derivative=0;
        for(int i=0;i<n;++i)derivative+=(L(i,enter)-L(i,leave))/
            ((1-step)*d[i]+step*end[i]);
        if(derivative>0)low=step;else high=step;
    }
    const double step=low==0 && high==1 && derivative>0?1:(low+high)/2;
    w[enter]+=step*mass;w[leave]-=step*mass;
    for(int i=0;i<n;++i)d[i]=(1-step)*d[i]+step*end[i];
}

// [[Rcpp::export(rng = false)]]
Rcpp::List fit_sparse_location_cpp(SEXP dictionary,int threads=1,
    int maximum_columns=2000,double tolerance=1e-6,int column_batch=64,
    Rcpp::Nullable<Rcpp::NumericVector> start_weights=R_NilValue,
    int omitted_row=0,double maximum_seconds=0,int fitted_class=-1) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);auto& op=*ptr;
    LocationOmissionScope omission(op,omitted_row);
    LocationClassScope class_scope(op,fitted_class);
    const auto start_time=std::chrono::steady_clock::now();
    const int n=op.n,k=op.nk;
    const int first_state=fitted_class<0?0:fitted_class*op.half,
        last_state=fitted_class<0?k:(fitted_class+1)*op.half;
    const double objective_scale=omitted_row>0?double(n)/(n-1):1;
    tolerance/=objective_scale;
    std::vector<int> support;std::vector<double> weights(1,1),density=op.umbrella;
    if(start_weights.isNotNull()) {
        const Rcpp::NumericVector initial(start_weights);
        if(initial.size()!=k)Rcpp::stop("Initial weights do not match the fixed dictionary");
        const double floor=*std::min_element(initial.begin()+first_state,initial.begin()+last_state);
        weights[0]=floor*(last_state-first_state);double total=weights[0];
        for(int j=first_state;j<last_state;++j)if(initial[j]>floor) {
            const double mass=initial[j]-floor;
            support.push_back(j);weights.push_back(mass);total+=mass;
        }
        if(!(total>0))Rcpp::stop("Initial weights have no positive mass");
        for(double& mass:weights)mass/=total;
        const LocationSparseActive active(op,support,threads);density=active.multiply(weights);
    }
    const int starting_columns=support.size();
    std::vector<char> selected(k,0);std::vector<int> pricing_row(k);
    int added=0,iteration=0;double gap=R_PosInf,initial_gap=R_PosInf;std::string status;
    for(;iteration<maximum_columns;++iteration) {
        Rcpp::checkUserInterrupt();
        auto gradient=op.gradient(density,&pricing_row);
        gap=*std::max_element(gradient.begin()+first_state,gradient.begin()+last_state)-1;
        if(iteration==0)initial_gap=gap*objective_scale;
        Rcpp::Rcout<<"sparse joint added="<<added<<", active="<<support.size()
            <<", KKT gap="<<gap*objective_scale<<"\n";
        if(gap<tolerance || added>=maximum_columns)break;
        if(maximum_seconds>0 && std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start_time).count()>=maximum_seconds)break;
        std::fill(selected.begin(),selected.end(),0);for(int x:support)selected[x]=1;
        // A block contains at most one candidate dominated by each observation.
        // Adjacent mean knots otherwise fill the block with nearly identical
        // columns for one gene. Every dictionary state remains in the KKT scan.
        std::vector<int> entering,best(n,-1);
        for(int j=first_state;j<last_state;++j)if(gradient[j]>1+tolerance && !selected[j]) {
            const int row=pricing_row[j],previous=best[row];
            if(previous<0 || gradient[j]>gradient[previous])best[row]=j;
        }
        for(int j:best)if(j>=0)entering.push_back(j);
        const int count=std::min({int(entering.size()),column_batch,maximum_columns-added});
        std::partial_sort(entering.begin(),entering.begin()+count,entering.end(),
            [&](int a,int b){return gradient[a]!=gradient[b]?gradient[a]>gradient[b]:a<b;});
        entering.resize(count);
        if(count>0) {
            std::vector<double> candidate(n,0);
            for(int j:entering)for(size_t at=op.column_start[j];at<op.column_start[j+1];++at)
                candidate[op.column_row[at]]+=op.column_value[at]/count;
            for(int i:op.constant_rows)candidate[i]=1;
            double low=0,high=1;
            for(int search=0;search<50;++search) {
                const double step=(low+high)/2;double derivative=0;
                for(int i=0;i<n;++i){const double diff=candidate[i]-density[i];
                    derivative+=diff/(density[i]+step*diff);}
                if(derivative>0)low=step;else high=step;
            }
            const double step=(low+high)/2;
            for(double& x:weights)x*=1-step;
            for(int j:entering){support.push_back(j);weights.push_back(step/count);}
            added+=count;
        }
        const int na=support.size();const LocationSparseActive active(op,support,threads);
        const double restricted=std::max(tolerance*.1,std::min(.1,gap*.1));
        const double remaining=maximum_seconds>0?std::max(1e-6,maximum_seconds-
            std::chrono::duration<double>(std::chrono::steady_clock::now()-start_time).count()):0;
        const auto fitted=location_sparse_weights(active,weights,restricted,remaining);
        weights=fitted.weights;
        status=fitted.gap<restricted?"converged":"restricted solve incomplete";
        density=active.multiply(weights);
        for(int exchange=0;exchange<2000;++exchange) {
            if(maximum_seconds>0 && std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start_time).count()>=maximum_seconds)break;
            const auto score=active.scores(density);
            int enter=0,leave=-1;
            for(int j=0;j<=na;++j) {
                if(score[j]>score[enter])enter=j;
                if(weights[j]>0 && (leave<0 || score[j]<score[leave]))leave=j;
            }
            // The maximum score minus one is the simplex duality gap.
            if(score[enter]-1<restricted)break;
            location_sparse_exchange(active,weights,density,enter,leave);
        }
        location_prune(support,weights);
    }
    const auto final_gradient=op.gradient(density);
    gap=*std::max_element(final_gradient.begin()+first_state,final_gradient.begin()+last_state)-1;
    const double unrestricted_gap=*std::max_element(final_gradient.begin(),final_gradient.end())-1,
        opposing_gap=fitted_class<0?NA_REAL:*std::max_element(
            final_gradient.begin()+(1-fitted_class)*op.half,
            final_gradient.begin()+(2-fitted_class)*op.half)-1;
    Rcpp::NumericVector fitted(k);
    for(int j=first_state;j<last_state;++j)fitted[j]=weights[0]/(last_state-first_state);
    for(int j=0;j<int(support.size());++j)fitted[support[j]]+=weights[j+1];
    Rcpp::NumericVector null(n);double density_error=0,total_mass=0,null_mass=0;
    for(int j=0;j<k;++j){total_mass+=fitted[j];if(j<op.half)null_mass+=fitted[j];}
    std::vector<double> partial(size_t(threads)*n*2,0),check_density(n);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int j=0;j<k;++j)if(fitted[j]>0) {
        int thread=0;
        #ifdef _OPENMP
        thread=omp_get_thread_num();
        #endif
        double* total=partial.data()+size_t(thread)*n*2;
        for(size_t at=op.column_start[j];at<op.column_start[j+1];++at) {
            const int i=op.column_row[at];
            const double mass=fitted[j]*op.column_value[at];
            total[i]+=mass;if(j<op.half)total[n+i]+=mass;
        }
    }
    for(int i=0;i<n;++i) {
        double d=0,q=0;
        for(int t=0;t<threads;++t){const double* total=partial.data()+size_t(t)*n*2;
            d+=total[i];q+=total[n+i];}
        if(i==op.omitted_row || (op.positive[i]==0 && op.positive[n+i]==0)){d=total_mass;q=null_mass;}
        check_density[i]=d;null[i]=q/d;
    }
    for(int i=0;i<n;++i)density_error=std::max(density_error,std::abs(check_density[i]-density[i])/density[i]);
    const double minimum=*std::min_element(density.begin(),density.end()),
        bound=op.epsilon/minimum;
    return Rcpp::List::create(Rcpp::Named("weights")=fitted,
        Rcpp::Named("kkt_gap")=gap*objective_scale,Rcpp::Named("converged")=(gap<tolerance),
        Rcpp::Named("initial_kkt_gap")=initial_gap,
        Rcpp::Named("starting_columns")=starting_columns,Rcpp::Named("omitted_row")=omitted_row,
        Rcpp::Named("fitted_class")=fitted_class,
        Rcpp::Named("unconstrained_kkt_gap")=unrestricted_gap*objective_scale,
        Rcpp::Named("opposing_class_gap")=opposing_gap*objective_scale,
        Rcpp::Named("columns")=added,Rcpp::Named("restricted_solves")=iteration,
        Rcpp::Named("retained_columns")=support.size(),
        Rcpp::Named("uniform_component_weight")=weights[0],
        Rcpp::Named("initial_component_weight")=weights[0],
        Rcpp::Named("last_solver_status")=status,Rcpp::Named("density")=Rcpp::wrap(density),
        Rcpp::Named("posterior_null")=null,
        Rcpp::Named("mixture_density_identity_error")=density_error,
        Rcpp::Named("maximum_relative_density_omission")=bound,
        Rcpp::Named("maximum_posterior_null_omission")=2*bound);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix sparse_location_posterior_cpp(SEXP dictionary,
    Rcpp::NumericVector weights,Rcpp::NumericVector b_contrast,int threads=1) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);auto& op=*ptr;
    const int n=op.n;Rcpp::NumericMatrix out(n,6);
    std::vector<char> active(2*op.np,0);
    for(int k=0;k<op.nk;++k)if(weights[k]>0)active[k/op.na]=1;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for(int i=0;i<n;++i) {
        double density=0,null=0,first=0,second=0,bfirst=0,bsecond=0;
        if(op.positive[i]==0 && op.positive[n+i]==0) {
            for(int k=0;k<op.nk;++k)if(weights[k]>0) {
                const int h=k/op.half,p=(k%op.half)/op.na;
                const double w=weights[k],b=b_contrast[p];
                density+=w;if(!h)null+=w;
                if(h){first+=w*op.slab_mean;second+=w*(op.tv+op.slab_mean*op.slab_mean);}
                bfirst+=w*b;bsecond+=w*b*b;
            }
        } else {
            std::vector<double> v(op.na),m1(op.na),m2(op.na);int low,high;
            for(int h=0;h<2;++h)for(int p=0;p<op.np;++p)if(active[h*op.np+p]) {
                op.block(i,p,h,v,m1,m2,low,high);
                const double b=b_contrast[p];
                for(int a=low;a<high;++a) {
                    const double w=weights[h*op.half+p*op.na+a];
                    if(w<=0 || !op.retained(v[a],m1[a],m2[a]))continue;
                    const double mass=w*v[a];
                    density+=mass;if(!h)null+=mass;
                    first+=w*m1[a];second+=w*m2[a];bfirst+=mass*b;bsecond+=mass*b*b;
                }
            }
        }
        out(i,0)=first/density;out(i,1)=std::max(0.0,second/density-first*first/(density*density));
        out(i,2)=null/density;out(i,3)=bfirst/density;
        out(i,4)=std::max(0.0,bsecond/density-bfirst*bfirst/(density*density));out(i,5)=density;
    }
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("effect","variance","posterior_null",
        "b_contrast","b_contrast_variance","density");
    return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List anchored_location_dictionary_cpp(Rcpp::NumericMatrix loglike,
    Rcpp::NumericMatrix mean,Rcpp::NumericMatrix variance,Rcpp::IntegerVector positive,
    Rcpp::IntegerVector left,Rcpp::IntegerVector right,Rcpp::NumericVector A,
    double A_sd,double slab_mean,double slab_sd,int threads,
    Rcpp::NumericVector log_sigma,Rcpp::NumericVector shape,double error_budget=1e-12,double anchor=-.5,bool symmetric=false) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(new SparseLocationDictionary(loglike,mean,
        variance,positive,left,right,A,A_sd,slab_mean,slab_sd,threads,log_sigma,shape,error_budget),true);
    ptr->anchor=anchor;ptr->symmetric_anchor=symmetric;
    ptr->build();
    return Rcpp::List::create(Rcpp::Named("dictionary")=ptr,
        Rcpp::Named("states")=ptr->nk,Rcpp::Named("entries")=double(ptr->entries),
        Rcpp::Named("storage_bytes")=double(ptr->entries)*(sizeof(int)+sizeof(double))+
            double(ptr->column_start.size())*sizeof(size_t),
        Rcpp::Named("dense_bytes")=double(ptr->n)*ptr->nk*sizeof(double),
        Rcpp::Named("row_error_budget")=ptr->epsilon,
        Rcpp::Named("log_row_envelope")=Rcpp::wrap(ptr->envelope),
        Rcpp::Named("maximum_base_moment_error")=ptr->moment_error,
        Rcpp::Named("minimum_residual_variance")=ptr->minimum_variance);
}



// One global null fraction and one shared nuisance distribution. Every physical
// nuisance state remains in the pricing scan; no unrestricted fit is required.
// [[Rcpp::export(rng = false)]]
Rcpp::List fit_shared_location_cpp(SEXP dictionary,int threads=1,
    int maximum_columns=2000,double tolerance=1e-6,int column_batch=64,
    Rcpp::Nullable<Rcpp::NumericVector> start_weights=R_NilValue,
    double initial_alternative=.5,int maximum_iterations=1000,
    double maximum_seconds=0) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);const auto& op=*ptr;
    const int n=op.n,k=op.half;
    const auto started=std::chrono::steady_clock::now();
    std::vector<double> base0(n,0),base1(n,0),partial(size_t(threads)*n*2,0);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
    #endif
    for(int j=0;j<2*k;++j) {
        int thread=0;
        #ifdef _OPENMP
        thread=omp_get_thread_num();
        #endif
        double* total=partial.data()+size_t(thread)*n*2+(j>=k?n:0);
        for(size_t at=op.column_start[j];at<op.column_start[j+1];++at)
            total[op.column_row[at]]+=op.column_value[at]/k;
    }
    for(int i=0;i<n;++i)for(int t=0;t<threads;++t) {
        base0[i]+=partial[size_t(t)*n*2+i];
        base1[i]+=partial[size_t(t)*n*2+n+i];
    }
    for(int i:op.constant_rows)base0[i]=base1[i]=1;
    std::vector<int> support;std::vector<double> weights(1,1);
    double alternative=initial_alternative;
    if(start_weights.isNotNull()) {
        const Rcpp::NumericVector initial(start_weights);
        if(initial.size()!=2*k)Rcpp::stop("Shared initial weights have the wrong size");
        std::vector<double> common(k);
        for(int j=0;j<k;++j)common[j]=initial[j]+initial[k+j];
        const double floor=*std::min_element(common.begin(),common.end());
        weights[0]=floor*k;double total=weights[0];
        for(int j=0;j<k;++j)if(common[j]>floor) {
            support.push_back(j);weights.push_back(common[j]-floor);
            total+=common[j]-floor;
        }
        for(double& value:weights)value/=total;
    }
    auto component_density=[&]() {
        std::fill(partial.begin(),partial.end(),0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int a=0;a<int(support.size());++a) {
            int thread=0;
            #ifdef _OPENMP
            thread=omp_get_thread_num();
            #endif
            double* total=partial.data()+size_t(thread)*n*2;
            for(int h=0;h<2;++h) {
                const int j=support[a]+h*k;
                for(size_t at=op.column_start[j];at<op.column_start[j+1];++at)
                    total[h*n+op.column_row[at]]+=weights[a+1]*op.column_value[at];
            }
        }
        std::pair<std::vector<double>,std::vector<double>> out{base0,base1};
        for(int i=0;i<n;++i) {
            out.first[i]*=weights[0];out.second[i]*=weights[0];
            for(int t=0;t<threads;++t) {
                out.first[i]+=partial[size_t(t)*n*2+i];
                out.second[i]+=partial[size_t(t)*n*2+n+i];
            }
        }
        for(int i:op.constant_rows)out.first[i]=out.second[i]=1;
        return out;
    };
    auto component=component_density();
    std::vector<double> density(n),gradient(k),candidate(n),zero(n),one(n);
    std::vector<int> pricing_row(2*k),best(n);
    std::vector<char> selected(k,0);
    for(int i=0;i<n;++i)density[i]=(1-alternative)*component.first[i]+
        alternative*component.second[i];
    double nuisance_gap=R_PosInf,class_gap=R_PosInf,initial_gap=R_PosInf;
    int added=0,iteration=0,starting_columns=support.size();
    for(;iteration<maximum_iterations;++iteration) {
        Rcpp::checkUserInterrupt();
        const auto score=op.gradient(density,&pricing_row);
        for(int j=0;j<k;++j)gradient[j]=(1-alternative)*score[j]+alternative*score[k+j];
        nuisance_gap=*std::max_element(gradient.begin(),gradient.end())-1;
        double c0=0,c1=0;
        for(int i=0;i<n;++i) {
            c0+=component.first[i]/density[i]/n;
            c1+=component.second[i]/density[i]/n;
        }
        class_gap=std::max(c0,c1)-1;
        if(iteration==0)initial_gap=std::max(nuisance_gap,class_gap);
        Rcpp::Rcout<<"shared location iteration="<<iteration<<", active="<<support.size()
            <<", alternative="<<alternative<<", gap="<<std::max(nuisance_gap,class_gap)<<"\n";
        if(std::max(nuisance_gap,class_gap)<tolerance)break;
        const double elapsed=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        if(maximum_seconds>0 && elapsed>=maximum_seconds)break;
        std::fill(selected.begin(),selected.end(),0);
        for(int state:support)selected[state]=1;
        std::fill(best.begin(),best.end(),-1);
        for(int j=0;j<k;++j)if(!selected[j] && gradient[j]>1+tolerance) {
            const int h=alternative*score[k+j]>(1-alternative)*score[j]?1:0,
                row=pricing_row[h*k+j],previous=best[row];
            if(previous<0 || gradient[j]>gradient[previous])best[row]=j;
        }
        std::vector<int> entering;
        for(int j:best)if(j>=0)entering.push_back(j);
        const int count=std::min({int(entering.size()),column_batch,maximum_columns-added});
        std::partial_sort(entering.begin(),entering.begin()+count,entering.end(),
            [&](int a,int b){return gradient[a]!=gradient[b]?gradient[a]>gradient[b]:a<b;});
        entering.resize(count);
        if(count>0) {
            std::fill(candidate.begin(),candidate.end(),0);
            for(int j:entering) {
                op.fill_column(j,zero.data());op.fill_column(k+j,one.data());
                for(int i=0;i<n;++i)candidate[i]+=
                    ((1-alternative)*zero[i]+alternative*one[i])/count;
            }
            double low=0,high=1;
            for(int step=0;step<50;++step) {
                const double fraction=(low+high)/2;double derivative=0;
                for(int i=0;i<n;++i) {
                    const double difference=candidate[i]-density[i];
                    derivative+=difference/(density[i]+fraction*difference);
                }
                if(derivative>0)low=fraction;else high=fraction;
            }
            const double fraction=(low+high)/2;
            for(double& value:weights)value*=1-fraction;
            for(int j:entering){support.push_back(j);weights.push_back(fraction/count);}
            added+=count;
        }
        {
            const LocationSparseActive active(op,support,threads,alternative,base0,base1);
            const double restricted=std::max(tolerance*.1,std::min(.1,nuisance_gap*.1)),
                remaining=maximum_seconds>0?std::max(1e-6,maximum_seconds-
                    std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()):0;
            const auto fitted=location_sparse_weights(active,weights,restricted,remaining);
            weights=fitted.weights;density=active.multiply(weights);
            for(int exchange=0;exchange<2000;++exchange) {
                const auto scores=active.scores(density);int enter=0,leave=-1;
                for(int j=0;j<active.k;++j) {
                    if(scores[j]>scores[enter])enter=j;
                    if(weights[j]>0 && (leave<0 || scores[j]<scores[leave]))leave=j;
                }
                if(scores[enter]-1<restricted)break;
                location_sparse_exchange(active,weights,density,enter,leave);
                if(maximum_seconds>0 && std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-started).count()>=maximum_seconds)break;
            }
        }
        location_prune(support,weights);
        component=component_density();
        double low=0,high=1,d0=0,d1=0;
        for(int i=0;i<n;++i) {
            const double difference=component.second[i]-component.first[i];
            d0+=difference/component.first[i];d1+=difference/component.second[i];
        }
        if(d0<=0)alternative=0;
        else if(d1>=0)alternative=1;
        else {
            for(int step=0;step<60;++step) {
                const double middle=(low+high)/2;double derivative=0;
                for(int i=0;i<n;++i) {
                    const double difference=component.second[i]-component.first[i];
                    derivative+=difference/(component.first[i]+middle*difference);
                }
                if(derivative>0)low=middle;else high=middle;
            }
            alternative=(low+high)/2;
        }
        for(int i=0;i<n;++i)density[i]=(1-alternative)*component.first[i]+
            alternative*component.second[i];
    }
    const auto score=op.gradient(density);
    double c0=0,c1=0;
    for(int j=0;j<k;++j)gradient[j]=(1-alternative)*score[j]+alternative*score[k+j];
    nuisance_gap=*std::max_element(gradient.begin(),gradient.end())-1;
    Rcpp::NumericVector fitted(2*k),null(n);
    for(int j=0;j<k;++j){fitted[j]=(1-alternative)*weights[0]/k;fitted[k+j]=alternative*weights[0]/k;}
    for(int a=0;a<int(support.size());++a) {
        fitted[support[a]]+=(1-alternative)*weights[a+1];
        fitted[k+support[a]]+=alternative*weights[a+1];
    }
    for(int i=0;i<n;++i) {
        c0+=component.first[i]/density[i]/n;c1+=component.second[i]/density[i]/n;
        null[i]=(1-alternative)*component.first[i]/density[i];
    }
    class_gap=std::max(c0,c1)-1;
    const double gap=std::max(nuisance_gap,class_gap),
        bound=op.epsilon/(*std::min_element(density.begin(),density.end()));
    return Rcpp::List::create(Rcpp::Named("weights")=fitted,
        Rcpp::Named("posterior_null")=null,Rcpp::Named("density")=Rcpp::wrap(density),
        Rcpp::Named("prior_null")=1-alternative,Rcpp::Named("kkt_gap")=gap,
        Rcpp::Named("nuisance_gap")=nuisance_gap,Rcpp::Named("class_gap")=class_gap,
        Rcpp::Named("initial_kkt_gap")=initial_gap,Rcpp::Named("converged")=(gap<tolerance),
        Rcpp::Named("starting_columns")=starting_columns,Rcpp::Named("columns")=added,
        Rcpp::Named("retained_columns")=support.size(),Rcpp::Named("iterations")=iteration,
        Rcpp::Named("uniform_component_weight")=weights[0],
        Rcpp::Named("maximum_relative_density_omission")=bound,
        Rcpp::Named("maximum_posterior_null_omission")=2*bound);
}

// Query the fitted continuous contrast distribution, including its point null.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix anchored_location_cdf_cpp(SEXP dictionary,
    Rcpp::NumericVector weights,Rcpp::NumericVector query,int threads=1) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);auto& op=*ptr;
    Rcpp::NumericMatrix out(op.n,3);
    std::vector<char> active(2*op.np,0);
    for(int k=0;k<op.nk;++k)if(weights[k]>0)active[k/op.na]=1;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for(int i=0;i<op.n;++i) {
        double density=0,lower=0,upper=0;
        if(op.positive[i]==0 && op.positive[op.n+i]==0) {
            for(int k=0;k<op.nk;++k)if(weights[k]>0) {
                const int h=k/op.half;
                const double w=weights[k];
                density+=w;
                if(h) {
                    const double cdf=R::pnorm5(query[i],op.slab_mean,std::sqrt(op.tv),1,0);
                    lower+=w*cdf;upper+=w*cdf;
                } else {lower+=w*(query[i]>0);upper+=w*(query[i]>=0);}
            }
        } else {
            std::vector<double> v(op.na),m1(op.na),m2(op.na);int low,high;
            for(int h=0;h<2;++h)for(int p=0;p<op.np;++p)if(active[h*op.np+p]) {
                op.block(i,p,h,v,m1,m2,low,high);
                LocationCurve curves[8];
                const int nc=op.curves(i,p,h,curves);
                for(int a=low;a<high;++a) {
                    const double w=weights[h*op.half+p*op.na+a];
                    if(w<=0 || !op.retained(v[a],m1[a],m2[a]))continue;
                    for(int c=0;c<nc;++c) {
                        const auto& curve=curves[c];
                        const double z=op.A[a]-curve.center;
                        const double mass=w*std::exp(curve.log_amplitude-op.envelope[i]-
                            z*z/(2*curve.variance));
                        density+=mass;
                        if(h) {
                            const double pm=curve.intercept+curve.slope*op.A[a],
                                cdf=R::pnorm5(query[i],pm,std::sqrt(curve.posterior_variance),1,0);
                            lower+=mass*cdf;upper+=mass*cdf;
                        } else {lower+=mass*(query[i]>0);upper+=mass*(query[i]>=0);}
                    }
                }
            }
        }
        out(i,0)=lower/density;out(i,1)=upper/density;out(i,2)=density;
    }
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("cdf_lower","cdf_upper","density");
    return out;
}

// Posterior intervals from the same normal-mixture likelihood and point null.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix anchored_location_summary_cpp(SEXP dictionary,
    Rcpp::NumericVector weights,double alpha,int threads=1,double tolerance=1e-10) {
    Rcpp::XPtr<SparseLocationDictionary> ptr(dictionary);auto& op=*ptr;
    Rcpp::NumericMatrix out(op.n,6);
    std::vector<char> active(2*op.np,0);
    for(int k=0;k<op.nk;++k)if(weights[k]>0)active[k/op.na]=1;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for(int i=0;i<op.n;++i) {
        double density=0,null=0,first=0,second=0;
        std::vector<double> mass,center,sd;
        if(op.positive[i]==0 && op.positive[op.n+i]==0) {
            double alternative=0;
            for(int k=0;k<op.nk;++k) {
                if(k<op.half)null+=weights[k];else alternative+=weights[k];
            }
            density=null+alternative;
            first=alternative*op.slab_mean;
            second=alternative*(op.tv+op.slab_mean*op.slab_mean);
            mass.push_back(alternative);center.push_back(op.slab_mean);sd.push_back(std::sqrt(op.tv));
        } else {
            std::vector<double> v(op.na),m1(op.na),m2(op.na);int low,high;
            for(int h=0;h<2;++h)for(int p=0;p<op.np;++p)if(active[h*op.np+p]) {
                op.block(i,p,h,v,m1,m2,low,high);
                LocationCurve curves[8];const int nc=op.curves(i,p,h,curves);
                for(int a=low;a<high;++a) {
                    const double w=weights[h*op.half+p*op.na+a];
                    if(w<=0 || !op.retained(v[a],m1[a],m2[a]))continue;
                    for(int c=0;c<nc;++c) {
                        const auto& curve=curves[c];
                        const double z=op.A[a]-curve.center,
                            value=w*std::exp(curve.log_amplitude-op.envelope[i]-
                                z*z/(2*curve.variance));
                        density+=value;
                        if(!h)null+=value;
                        else {
                            const double m=curve.intercept+curve.slope*op.A[a];
                            first+=value*m;second+=value*(curve.posterior_variance+m*m);
                            mass.push_back(value);center.push_back(m);sd.push_back(std::sqrt(curve.posterior_variance));
                        }
                    }
                }
            }
        }
        out(i,0)=first/density;
        out(i,1)=std::max(0.0,second/density-out(i,0)*out(i,0));
        out(i,2)=null/density;out(i,3)=density;
        double at_zero=0;
        for(size_t k=0;k<mass.size();++k)at_zero+=mass[k]*R::pnorm5(0,center[k],sd[k],1,0);
        for(int tail=0;tail<2;++tail) {
            const double probability=tail?1-alpha/2:alpha/2,target=probability*density;
            if(target>=at_zero && target<=at_zero+null) {out(i,4+tail)=0;continue;}
            const double quantile=R::qnorm5(probability,0,1,1,0);
            double low=0,high=0;
            for(size_t k=0;k<mass.size();++k)if(mass[k]>0) {
                const double q=center[k]+quantile*sd[k];
                low=std::min(low,q);high=std::max(high,q);
            }
            for(int iteration=0;iteration<64;++iteration) {
                const double point=(low+high)/2;
                double cdf=null*(point>=0);
                for(size_t k=0;k<mass.size();++k)cdf+=mass[k]*R::pnorm5(point,center[k],sd[k],1,0);
                if(cdf>=target)high=point;else low=point;
                if(high-low<=tolerance*(1+std::max(std::abs(low),std::abs(high))))break;
            }
            out(i,4+tail)=(low+high)/2;
        }
    }
    Rcpp::colnames(out)=Rcpp::CharacterVector::create("effect","variance","posterior_null","density","lower","upper");
    return out;
}
