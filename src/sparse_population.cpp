// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include "sgcb_sparse_mixture.h"
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
// Every normalized row retains a maximum entry of one. A truncated KKT gap g
// gives density_i >= 1/(n*(1+g)). Dropping entries <= tolerance/(4*n) therefore
// increases any complete score by at most tolerance*(1+g)/4. Solve the truncated
// problem to tolerance/2. Full scans or the omitted-score bound certify all columns.
struct PopulationSparseDictionary {
    Rcpp::NumericMatrix input;
    int n,nk,threads;bool log_input,sparse_only;
    double omitted_bound,threshold;
    std::vector<int> rows,column_row,constant_rows;
    std::vector<size_t> column_start;
    std::vector<double> row_scale,column_value,umbrella;
    PopulationSparseDictionary(Rcpp::NumericMatrix L,Rcpp::Nullable<Rcpp::LogicalVector> eligible,
        bool logs,int nt,double tolerance,SEXP prepared=R_NilValue):
        input(L),nk(L.ncol()),threads(nt),log_input(logs),sparse_only(prepared!=R_NilValue) {
        if(sparse_only) {
            Rcpp::XPtr<SparsePopulationBuffer> buffer(prepared);
            if(buffer->consumed)Rcpp::stop("Sparse likelihood buffer has already been consumed");
            if(buffer->threshold>tolerance/(4.0*buffer->n))
                Rcpp::stop("Sparse likelihood buffer does not meet the requested tolerance");
            n=buffer->n;nk=buffer->nk;
            omitted_bound=buffer->omitted_bound;threshold=buffer->threshold;
            rows=std::move(buffer->rows);row_scale=std::move(buffer->row_scale);
            column_start=std::move(buffer->column_start);
            column_row=std::move(buffer->column_row);
            column_value=std::move(buffer->column_value);
            umbrella=std::move(buffer->umbrella);
            buffer->consumed=true;
            return;
        }
        if(eligible.isNotNull()) {
            Rcpp::LogicalVector included(eligible);
            for(int i=0;i<L.nrow();++i)if(included[i])rows.push_back(i);
        } else for(int i=0;i<L.nrow();++i)rows.push_back(i);
        n=rows.size();omitted_bound=tolerance/4;threshold=omitted_bound/n;
        row_scale.assign(n,log_input?R_NegInf:0);
        std::vector<double> partial(size_t(threads)*n,log_input?R_NegInf:0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            int t=0;
            #ifdef _OPENMP
            t=omp_get_thread_num();
            #endif
            double* largest=partial.data()+size_t(t)*n;
            for(int i=0;i<n;++i)largest[i]=std::max(largest[i],input(rows[i],k));
        }
        for(int i=0;i<n;++i)for(int t=0;t<threads;++t)
            row_scale[i]=std::max(row_scale[i],partial[size_t(t)*n+i]);
        std::fill(partial.begin(),partial.end(),0);
        std::vector<size_t> count(nk,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            int t=0;
            #ifdef _OPENMP
            t=omp_get_thread_num();
            #endif
            double* sum=partial.data()+size_t(t)*n;
            for(int i=0;i<n;++i) {
                const double value=full_value(i,k);
                if(value>threshold){++count[k];sum[i]+=value/nk;}
            }
        }
        umbrella.assign(n,0);
        for(int i=0;i<n;++i)for(int t=0;t<threads;++t)umbrella[i]+=partial[size_t(t)*n+i];
        column_start.assign(size_t(nk)+1,0);
        for(int k=0;k<nk;++k)column_start[k+1]=column_start[k]+count[k];
        column_row.resize(column_start.back());column_value.resize(column_start.back());
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            size_t at=column_start[k];
            for(int i=0;i<n;++i) {
                const double value=full_value(i,k);
                if(value>threshold){column_row[at]=i;column_value[at++]=value;}
            }
        }
    }
    double full_value(int i,int k) const {
        return log_input?std::exp(input(rows[i],k)-row_scale[i]):input(rows[i],k)/row_scale[i];
    }
    std::vector<double> gradient(const std::vector<double>& density,
        std::vector<int>* dominant_row=nullptr) const {
        std::vector<double> score(nk,0),factor(n);
        for(int i=0;i<n;++i)factor[i]=1/(n*density[i]);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k) {
            double maximum=-1;int dominant=0;
            for(size_t at=column_start[k];at<column_start[k+1];++at) {
                const int i=column_row[at];const double contribution=column_value[at]*factor[i];
                score[k]+=contribution;
                if(contribution>maximum){maximum=contribution;dominant=i;}
            }
            if(dominant_row)(*dominant_row)[k]=dominant;
        }
        return score;
    }
    std::vector<double> full_density(const Rcpp::NumericVector& weights) const {
        std::vector<double> partial(size_t(threads)*n,0),out(n,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k)if(weights[k]>0) {
            int t=0;
            #ifdef _OPENMP
            t=omp_get_thread_num();
            #endif
            double* sum=partial.data()+size_t(t)*n;
            for(int i=0;i<n;++i)sum[i]+=weights[k]*full_value(i,k);
        }
        for(int i=0;i<n;++i)for(int t=0;t<threads;++t)out[i]+=partial[size_t(t)*n+i];
        return out;
    }
    double full_gap(const std::vector<double>& density) const {
        std::vector<double> score(nk,0);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
        #endif
        for(int k=0;k<nk;++k)for(int i=0;i<n;++i)score[k]+=full_value(i,k)/(n*density[i]);
        return *std::max_element(score.begin(),score.end())-1;
    }
};

// [[Rcpp::export(rng = false)]]
Rcpp::List fit_sparse_population_dictionary_cpp(Rcpp::NumericMatrix loglike,Rcpp::Function solver,
    int threads=1,int maximum_columns=200,double tolerance=1e-6,
    Rcpp::Nullable<Rcpp::NumericVector> start_weights=R_NilValue,bool dense_start=false,
    Rcpp::Nullable<Rcpp::LogicalVector> eligible=R_NilValue,bool reuse_buffer=false,
    std::string checkpoint_directory="",int column_batch=64,bool input_is_log=true,
    double maximum_seconds=0,SEXP prepared_dictionary=R_NilValue) {
    if(!checkpoint_directory.empty())Rcpp::stop("Sparse solver does not consume legacy checkpoints");
    const double requested_tolerance=tolerance;
    const auto started=std::chrono::steady_clock::now();
    PopulationSparseDictionary op(loglike,eligible,input_is_log,threads,tolerance,prepared_dictionary);
    tolerance*=.5;    const int n=op.n,k=op.nk;
    std::vector<int> support;std::vector<double> weights(1,1),density=op.umbrella;
    if(start_weights.isNotNull()) {
        const Rcpp::NumericVector initial(start_weights);weights[0]=0;
        for(int j=0;j<k;++j)if(initial[j]>0){support.push_back(j);weights.push_back(initial[j]);}
        const double total=std::accumulate(weights.begin(),weights.end(),0.0);
        for(double& mass:weights)mass/=total;
        const LocationSparseActive active(op,support,threads);density=active.multiply(weights);
        double derivative=0;
        for(int i=0;i<n;++i)derivative+=(op.umbrella[i]-density[i])/density[i];
        if(derivative>0) {
            double low=0,high=1;
            for(int search=0;search<50;++search) {
                const double step=(low+high)/2;derivative=0;
                for(int i=0;i<n;++i) {
                    const double difference=op.umbrella[i]-density[i];
                    derivative+=difference/(density[i]+step*difference);
                }
                if(derivative>0)low=step;else high=step;
            }
            const double step=(low+high)/2;
            for(double& mass:weights)mass*=1-step;
            weights[0]=step;
            for(int i=0;i<n;++i)density[i]=(1-step)*density[i]+step*op.umbrella[i];
        }
    }
    std::vector<char> selected(k,0);std::vector<int> pricing_row(k);
    int added=0,iteration=0;double gap=R_PosInf,initial_gap=R_PosInf;std::string status;
    for(;iteration<maximum_columns;++iteration) {
        Rcpp::checkUserInterrupt();
        auto gradient=op.gradient(density,&pricing_row);
        gap=*std::max_element(gradient.begin(),gradient.end())-1;
        if(iteration==0)initial_gap=gap;
        if(gap<tolerance || added>=maximum_columns)break;
        if(maximum_seconds>0 && std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count()>=maximum_seconds)break;
        std::fill(selected.begin(),selected.end(),0);for(int x:support)selected[x]=1;
        // A block contains at most one candidate dominated by each observation.
        // Adjacent mean knots otherwise fill the block with nearly identical
        // columns for one gene. Every dictionary state remains in the KKT scan.
        std::vector<int> entering,best(n,-1);
        for(int j=0;j<k;++j)if(gradient[j]>1+tolerance && !selected[j]) {
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
            std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()):0;
        const auto fitted=location_sparse_weights(active,weights,restricted,remaining);
        weights=fitted.weights;
        status=fitted.gap<restricted?"converged":"restricted solve incomplete";
        density=active.multiply(weights);
        for(int exchange=0;exchange<2000;++exchange) {
            if(maximum_seconds>0 && std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count()>=maximum_seconds)break;
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
        if(added>0 && added%40<column_batch)
            Rcpp::Rcout<<"sparse joint columns="<<added<<", KKT gap="<<gap<<"\n";
    }
    Rcpp::NumericVector fitted(k,weights[0]/k);
    for(int j=0;j<int(support.size());++j)fitted[support[j]]+=weights[j+1];
    double exact_gap=R_PosInf,density_error=NA_REAL,density_error_bound=NA_REAL;
    if(op.sparse_only) {
        const auto score=op.gradient(density);
        double inverse_sum=0,minimum=R_PosInf;
        for(double d:density){inverse_sum+=1/d;minimum=std::min(minimum,d);}
        // Full density is no smaller than sparse density, and each omitted
        // normalized entry is at most threshold. This certifies every column.
        exact_gap=*std::max_element(score.begin(),score.end())-1+
            op.threshold*inverse_sum/n;
        density_error_bound=op.threshold*std::accumulate(fitted.begin(),fitted.end(),0.0)/minimum;
    } else {
        const auto exact_density=op.full_density(fitted);
        exact_gap=op.full_gap(exact_density);
        density_error=0;
        for(int i=0;i<n;++i)density_error=std::max(density_error,
            std::abs(exact_density[i]-density[i])/exact_density[i]);
        density_error_bound=density_error;
    }
    return Rcpp::List::create(Rcpp::Named("weights")=fitted,
        Rcpp::Named("kkt_gap")=exact_gap,Rcpp::Named("converged")=(exact_gap<requested_tolerance),
        Rcpp::Named("initial_kkt_gap")=initial_gap,
        Rcpp::Named("truncated_kkt_gap")=gap,Rcpp::Named("columns")=added,
        Rcpp::Named("retained_columns")=support.size(),Rcpp::Named("iterations")=iteration,
        Rcpp::Named("uniform_component_weight")=weights[0],
        Rcpp::Named("last_solver_status")=status,
        Rcpp::Named("maximum_density_difference")=density_error,
        Rcpp::Named("maximum_density_difference_bound")=density_error_bound,
        Rcpp::Named("kkt_gap_type")=op.sparse_only ? "certified_full_dictionary_upper_bound" : "evaluated_full_dictionary",
        Rcpp::Named("sparsification_bound")=op.omitted_bound,
        Rcpp::Named("sparse_entries")=double(op.column_value.size()),
        Rcpp::Named("dictionary_rows")=n,Rcpp::Named("dictionary_states")=k);
}
