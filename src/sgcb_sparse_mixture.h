#ifndef SGCB_SPARSE_MIXTURE_H
#define SGCB_SPARSE_MIXTURE_H
#include <Rcpp.h>
#include <Eigen/SparseCore>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <chrono>
#include <utility>

// Likelihood storage shared by the streamed generator and the existing solver.
struct SparsePopulationBuffer {
    int n=0,nk=0;
    bool consumed=false;
    double omitted_bound=0,threshold=0;
    std::vector<int> rows,column_row;
    std::vector<size_t> column_start;
    std::vector<double> row_scale,column_value,umbrella;
};
static void location_prune(std::vector<int>& support,std::vector<double>& weights) {
    std::vector<int> kept;std::vector<double> mass(1,weights[0]);
    for(int j=0;j<int(support.size());++j)if(weights[j+1]!=0) {
        kept.push_back(support[j]);mass.push_back(weights[j+1]);
    }
    support.swap(kept);weights.swap(mass);
}
// Restricted Newton products retain the dictionary's exact sparse entries.
// Both storage orders avoid dense n-by-active matrices and thread-local copies.
struct LocationSparseActive {
    int n,k,threads;
    Eigen::SparseMatrix<double> column;
    Eigen::SparseMatrix<double,Eigen::RowMajor> row;
    template<class Dictionary>
    LocationSparseActive(const Dictionary& op,
        const std::vector<int>& support,int nt):n(op.n),k(support.size()+1),
        threads(nt),column(n,k) {
        size_t entries=n;
        for(int state:support)entries+=op.column_start[state+1]-op.column_start[state]+
            op.constant_rows.size();
        column.reserve(entries);
        column.startVec(0);
        for(int i=0;i<n;++i)column.insertBack(i,0)=op.umbrella[i];
        for(int j=1;j<k;++j) {
            column.startVec(j);
            const int state=support[j-1];
            size_t at=op.column_start[state],end=op.column_start[state+1],constant=0;
            while(at<end || constant<op.constant_rows.size()) {
                if(constant<op.constant_rows.size() &&
                    (at==end || op.constant_rows[constant]<=op.column_row[at])) {
                    const int i=op.constant_rows[constant++];
                    column.insertBack(i,j)=1;
                    if(at<end && op.column_row[at]==i)++at;
                } else {
                    column.insertBack(op.column_row[at],j)=op.column_value[at];++at;
                }
            }
        }
        column.finalize();row=column;
    }
    // The null and alternative share one nuisance law. Merge their CSC
    // postings only for the active columns used by the existing Newton solver.
    template<class Dictionary>
    LocationSparseActive(const Dictionary& op,const std::vector<int>& support,
        int nt,double alternative,const std::vector<double>& base0,
        const std::vector<double>& base1):n(op.n),k(support.size()+1),
        threads(nt),column(n,k) {
        size_t entries=n;
        for(int state:support)entries+=
            op.column_start[state+1]-op.column_start[state]+
            op.column_start[op.half+state+1]-op.column_start[op.half+state]+
            op.constant_rows.size();
        column.reserve(entries);column.startVec(0);
        for(int i=0;i<n;++i)column.insertBack(i,0)=
            (1-alternative)*base0[i]+alternative*base1[i];
        const int sentinel=std::numeric_limits<int>::max();
        for(int j=1;j<k;++j) {
            const int state=support[j-1];
            size_t at0=op.column_start[state],end0=op.column_start[state+1],
                at1=op.column_start[op.half+state],
                end1=op.column_start[op.half+state+1],constant=0;
            column.startVec(j);
            while(at0<end0 || at1<end1 || constant<op.constant_rows.size()) {
                const int row0=at0<end0?op.column_row[at0]:sentinel,
                    row1=at1<end1?op.column_row[at1]:sentinel,
                    rowc=constant<op.constant_rows.size()?op.constant_rows[constant]:sentinel,
                    i=std::min({row0,row1,rowc});
                double value=0;
                if(row0==i)value+=(1-alternative)*op.column_value[at0++];
                if(row1==i)value+=alternative*op.column_value[at1++];
                if(rowc==i){value=1;++constant;}
                if(value>0)column.insertBack(i,j)=value;
            }
        }
        column.finalize();row=column;
    }

    std::vector<double> multiply(const std::vector<double>& x) const {
        std::vector<double> out(n,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int i=0;i<n;++i)
            for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator at(row,i);at;++at)
                out[i]+=at.value()*x[at.col()];
        return out;
    }
    std::vector<double> transpose(const std::vector<double>& factor) const {
        std::vector<double> out(k,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int j=0;j<k;++j)
            for(Eigen::SparseMatrix<double>::InnerIterator at(column,j);at;++at)
                out[j]+=at.value()*factor[at.row()];
        return out;
    }
    std::vector<double> scores(const std::vector<double>& density) const {
        std::vector<double> factor(n);
        for(int i=0;i<n;++i)factor[i]=1/(n*density[i]);
        return transpose(factor);
    }
};

struct LocationSparseSolve {
    std::vector<double> weights;
    int iterations=0,cg_iterations=0;
    double gap=R_PosInf;
};

static long double location_sparse_objective_difference(const LocationSparseActive& op,
    const std::vector<double>& density,const std::vector<double>& next,
    const std::vector<double>& gradient,const std::vector<double>& difference) {
    const auto projected=op.multiply(difference);
    long double value=0;
    for(int j=0;j<op.k;++j)value+=static_cast<long double>(gradient[j])*difference[j];
    for(int i=0;i<op.n;++i) {
        const long double r=static_cast<long double>(projected[i])/density[i];
        long double remainder;
        if(std::abs(r)<.125L) {
            long double term=r*r;remainder=term/2;
            for(int j=3;j<128;++j) {
                term*=-r;const long double updated=remainder+term/j;
                if(updated==remainder)break;
                remainder=updated;
            }
        } else remainder=r-std::log(static_cast<long double>(next[i])/density[i]);
        value+=remainder/op.n;
    }
    return value;
}

static LocationSparseSolve location_sparse_weights(const LocationSparseActive& op,
    const std::vector<double>& initial,double tolerance,double maximum_seconds=0) {
    const auto start_time=std::chrono::steady_clock::now();
    Rcpp::Environment ns=Rcpp::Environment::namespace_env("mixsqp");
    Rcpp::Function defaults=ns["mixsqp_control_default"];
    const Rcpp::List control=defaults();
    const int maximum=Rcpp::as<int>(control["maxiter.sqp"]);
    const double sufficient=Rcpp::as<double>(control["suffdecr.linesearch"]),
        reduction=Rcpp::as<double>(control["stepsizereduce"]),
        minimum_step=Rcpp::as<double>(control["minstepsize"]),
        machine=std::numeric_limits<double>::epsilon(),zero_threshold=1e-12;
    const int n=op.n,k=op.k;
    LocationSparseSolve fit;fit.weights=initial;
    auto& x=fit.weights;
    double total=std::accumulate(x.begin(),x.end(),0.0);
    for(double& value:x)value/=total;
    auto density=op.multiply(x);
    std::vector<double> gradient(k),diagonal(k),free(k),residual(k),z(k),direction(k),
        step_direction(k),factor(n),product(k),candidate(k),difference(k);
    for(;fit.iterations<maximum;++fit.iterations) {
        if(maximum_seconds>0 && std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start_time).count()>=maximum_seconds)break;
        const auto score=op.scores(density);
        for(int j=0;j<k;++j)gradient[j]=1-score[j];
        fit.gap=std::max(0.0,-*std::min_element(gradient.begin(),gradient.end()));
        if(fit.gap<tolerance)break;
        Rcpp::checkUserInterrupt();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(op.threads)
        #endif
        for(int j=0;j<k;++j) {
            free[j]=x[j]>0 || gradient[j]<0;
            diagonal[j]=0;
            for(Eigen::SparseMatrix<double>::InnerIterator at(op.column,j);at;++at) {
                const double ratio=at.value()/density[at.row()];
                diagonal[j]+=ratio*ratio/n;
            }
        }
        const double ridge=std::sqrt(machine) * (*std::max_element(diagonal.begin(),diagonal.end()));
        double residual2=0,rz=0;
        for(int j=0;j<k;++j) {
            diagonal[j]+=ridge;residual[j]=-gradient[j]*free[j];
            residual2+=residual[j]*residual[j];z[j]=residual[j]/diagonal[j];
            direction[j]=z[j];step_direction[j]=0;rz+=residual[j]*z[j];
        }
        const double target=std::min(.5,std::sqrt(std::max(fit.gap,machine)))*std::sqrt(residual2);
        for(int cg=0;cg<k && std::sqrt(residual2)>target;++cg) {
            if(maximum_seconds>0 && std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start_time).count()>=maximum_seconds)break;
            const auto projected=op.multiply(direction);
            for(int i=0;i<n;++i)factor[i]=projected[i]/density[i]/density[i]/n;
            product=op.transpose(factor);
            double curvature=0;
            for(int j=0;j<k;++j) {
                product[j]=(product[j]+ridge*direction[j])*free[j];
                curvature+=direction[j]*product[j];
            }
            if(!(curvature>0) || !(rz>0))break;
            double scale=rz/curvature;
            int boundary=-1;
            // Keep the quadratic step feasible while solving it. Clipping an
            // unconstrained Newton solution afterwards discards its curvature
            // calculation and can cause many small projected iterations.
            for(int j=0;j<k;++j)if(free[j] && direction[j]<0) {
                const double limit=std::max(0.0,x[j]+step_direction[j])/(-direction[j]);
                if(limit<scale){scale=limit;boundary=j;}
            }
            if(boundary>=0)free[boundary]=0;
            residual2=0;
            for(int j=0;j<k;++j) {
                step_direction[j]+=scale*direction[j];residual[j]-=scale*product[j];
                residual[j]*=free[j];
                residual2+=residual[j]*residual[j];
            }
            if(boundary>=0)step_direction[boundary]=-x[boundary];
            ++fit.cg_iterations;
            if(std::sqrt(residual2)<=target)break;
            double next_rz=0;
            for(int j=0;j<k;++j){z[j]=residual[j]/diagonal[j];next_rz+=residual[j]*z[j];}
            for(int j=0;j<k;++j)direction[j]=boundary>=0?z[j]:z[j]+next_rz/rz*direction[j];
            rz=next_rz;
        }
        double slope=0;
        for(int j=0;j<k;++j) {
            if(x[j]==0 && step_direction[j]<0)step_direction[j]=0;
            slope+=gradient[j]*step_direction[j];
        }
        if(!(slope<0))for(int j=0;j<k;++j)
            step_direction[j]=-gradient[j]*free[j]/diagonal[j];
        bool accepted=false;
        for(double step=1;step>=minimum_step;step*=reduction) {
            if(maximum_seconds>0 && std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start_time).count()>=maximum_seconds)break;
            total=0;
            for(int j=0;j<k;++j) {
                candidate[j]=std::max(0.0,x[j]+step*step_direction[j]);
                if(candidate[j]<zero_threshold)candidate[j]=0;
                total+=candidate[j];
            }
            if(!(total>0))continue;
            for(double& value:candidate)value/=total;
            const auto next=op.multiply(candidate);
            if(!(*std::min_element(next.begin(),next.end())>0))continue;
            slope=0;
            for(int j=0;j<k;++j){difference[j]=candidate[j]-x[j];slope+=gradient[j]*difference[j];}
            const long double change=location_sparse_objective_difference(op,density,next,
                gradient,difference);
            if(slope<0 && change<=sufficient*slope) {
                x=candidate;density=next;accepted=true;break;
            }
        }
        if(!accepted)break;
    }
    const auto score=op.scores(density);
    fit.gap=std::max(0.0,*std::max_element(score.begin(),score.end())-1);
    return fit;
}

static void location_sparse_exchange(const LocationSparseActive& op,
    std::vector<double>& weights,std::vector<double>& density,int enter,int leave) {
    const double mass=weights[leave];auto endpoint_weights=weights;
    endpoint_weights[enter]+=mass;endpoint_weights[leave]=0;
    const auto end=op.multiply(endpoint_weights);
    std::vector<double> difference(op.n,0);
    for(Eigen::SparseMatrix<double>::InnerIterator at(op.column,enter);at;++at)
        difference[at.row()]+=at.value();
    for(Eigen::SparseMatrix<double>::InnerIterator at(op.column,leave);at;++at)
        difference[at.row()]-=at.value();
    double low=0,high=1,boundary=0;
    for(int i=0;i<op.n;++i)if(difference[i]!=0)boundary+=difference[i]/end[i];
    if(boundary<=0)for(int search=0;search<50;++search) {
        const double step=(low+high)/2;double derivative=0;
        for(int i=0;i<op.n;++i)if(difference[i]!=0)
            derivative+=difference[i]/((1-step)*density[i]+step*end[i]);
        if(derivative>0)low=step;else high=step;
    }
    const double step=boundary>0?1:(low+high)/2;
    weights[enter]+=step*mass;weights[leave]-=step*mass;
    for(int i=0;i<op.n;++i)density[i]=(1-step)*density[i]+step*end[i];
}

#endif
