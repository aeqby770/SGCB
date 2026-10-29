// Prior replica exchange for the same finite-dictionary collapsed posterior.
// van Havre et al. (2015), doi:10.1371/journal.pone.0131739.
// Only the concentration of auxiliary replicas changes; cold draws retain
// the original concentration. Swaps integrate out all Dirichlet weights.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <mutex>
#include <cstdint>
#include <chrono>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

// Conditional on the partition by (class,pair), every empty destination has
// the same uniform-Dirichlet factor: the counts at each level are unchanged.
// The analogous level/class move preserves pair counts. Alternating these
// exact heatbaths moves coherent groups while retaining the full joint prior.
static void mean_joint_coarse_heatbath(const double* L,int n,int V,int nA,
    std::vector<int>& z,std::vector<int>& count,std::vector<int>& active,
    std::vector<int>& position,std::vector<int>& nc,double& attempts,double& changes,int threads) {
    const int P=V/nA;
    for (int axis=0;axis<2;++axis) {
        const int units=axis==0 ? P : nA,G=2*units;
        std::vector<std::vector<int>> members(G);
        std::vector<int> number(G,0),partition;
        for (int i=0;i<n;++i) {
            const int state=z[i],v=state%V,
                group=(state/V)*units+(axis==0 ? v/nA : v%nA);
            if (number[group]++==0) partition.push_back(group);
            members[group].push_back(i);
        }
        std::vector<double> logp(G),cdf(G);
        for (int id:partition) {
            const auto& genes=members[id];const int first=z[genes[0]],
                old=(first/V)*units+(axis==0 ? (first%V)/nA : first%nA);
            number[old]-=genes.size();
            for (int i:genes) {
                const int state=z[i];
                if (--count[state]==0) {
                    const int at=position[state],last=active.back();
                    active[at]=last;position[last]=at;active.pop_back();position[state]=-1;
                }
                --nc[state/V];
            }
            auto candidate_log_likelihood=[&](int group) {
                if (number[group]>0) return R_NegInf;
                const int c=group/units,unit=group%units;double value=0;
                for (int i:genes) {
                    const int v=z[i]%V,state=c*V+(axis==0 ? unit*nA+v%nA : (v/nA)*nA+unit);
                    const double probability=L[i+size_t(n)*state];
                    if (probability<=0) return R_NegInf;
                    value+=std::log(probability);
                }
                return value;
            };
            // Use the unchanged worker pool only for substantial partitions;
            // small partitions remain serial without resizing OpenMP teams.
            if (G>=threads && genes.size()>=size_t(threads)) {
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
                #endif
                for (int group=0;group<G;++group) logp[group]=candidate_log_likelihood(group);
            } else {
                for (int group=0;group<G;++group) logp[group]=candidate_log_likelihood(group);
            }
            const double peak=*std::max_element(logp.begin(),logp.end());
            double total=0;
            for (int group=0;group<G;++group) {total+=std::exp(logp[group]-peak);cdf[group]=total;}
            const int next=std::min(G-1,int(std::upper_bound(cdf.begin(),cdf.end(),R::runif(0,1)*total)-cdf.begin())),
                c=next/units,unit=next%units;
            for (int i:genes) {
                const int v=z[i]%V,state=c*V+(axis==0 ? unit*nA+v%nA : (v/nA)*nA+unit);
                if (count[state]++==0) {position[state]=active.size();active.push_back(state);}
                ++nc[c];z[i]=state;
            }
            number[next]+=genes.size();++attempts;changes+=next!=old;
        }
    }
}
static double mean_joint_logadd(double a,double b) {
    if (a==R_NegInf) return b;
    if (b==R_NegInf) return a;
    const double peak=std::max(a,b);
    return peak+std::log1p(std::exp(std::min(a,b)-peak));
}

// Resample every non-anchor allocation in two complete cells. Integrating
// their Dirichlet weights gives a beta-weighted Poisson-binomial count law.
// The polynomial coefficients and backward conditional draws are exact.
// Fixed anchors keep state-pair selection reversible while permitting the
// other members to split, merge or exchange their occupied cells.
static void mean_joint_allocation_block(const double* L,int n,int V,double alpha,
    std::vector<int>& z,std::vector<int>& count,std::vector<int>& active,
    std::vector<int>& position,std::vector<int>& nc,int anchor,int other,int a,int b,
    double& attempts,double& changes,double& class_changes,double& births,double& deaths,
    double selection_total,double& rejections) {
    std::vector<int> genes;
    for (int i=0;i<n;++i) if (i!=anchor && i!=other && (z[i]==a || z[i]==b)) genes.push_back(i);
    const int m=genes.size();
    if (m==0) return;
    ++attempts;
    const int fixed_a=1,fixed_b=other>=0 ? 1 : 0;
    std::vector<double> loga(m),logb(m),coefficient(size_t(m+1)*(m+2)/2,R_NegInf);
    coefficient[0]=0;
    for (int j=1;j<=m;++j) {
        const int i=genes[j-1];
        const double u=L[i+size_t(n)*a],v=L[i+size_t(n)*b],normalizer=std::log(u+v);
        loga[j-1]=std::log(u)-normalizer;logb[j-1]=std::log(v)-normalizer;
        const size_t row=size_t(j)*(j+1)/2,previous=size_t(j-1)*j/2;
        for (int k=0;k<=j;++k) coefficient[row+k]=mean_joint_logadd(
            k<j ? coefficient[previous+k]+logb[j-1] : R_NegInf,
            k>0 ? coefficient[previous+k-1]+loga[j-1] : R_NegInf);
    }
    const size_t last=size_t(m)*(m+1)/2;
    std::vector<double> mass(m+1),cdf(m+1);
    double peak=R_NegInf;
    for (int k=0;k<=m;++k) {
        mass[k]=coefficient[last+k]+std::lgamma(alpha+fixed_a+k)+std::lgamma(alpha+fixed_b+m-k);
        peak=std::max(peak,mass[k]);
    }
    double total=0;
    for (int k=0;k<=m;++k) {total+=std::exp(mass[k]-peak);cdf[k]=total;}
    const double draw=std::min(std::nextafter(total,0.0),R::runif(0,1)*total);
    int k=std::upper_bound(cdf.begin(),cdf.end(),draw)-cdf.begin();
    const int next_a=fixed_a+k,next_b=fixed_b+m-k,old_a=count[a],old_b=count[b];
    std::vector<int> proposed(m);
    for (int j=m;j>0;--j) {
        const size_t row=size_t(j)*(j+1)/2,previous=size_t(j-1)*j/2;
        const bool use_a=k==j || (k>0 && std::log(R::runif(0,1))<
            coefficient[previous+k-1]+loga[j-1]-coefficient[row+k]);
        proposed[j-1]=use_a ? a : b;k-=use_a;
    }
    // The second anchor is selected with weight L[first anchor, its cell].
    // The anchors stay fixed, but the opposite-class selection normalizer
    // changes when other members move. Its reverse/forward ratio is Z/Z'.
    if (selection_total>0) {
        const double next_total=selection_total+(next_b-old_b)*L[anchor+size_t(n)*b];
        if (std::log(R::runif(0,1))>=std::log(selection_total)-std::log(next_total)) {
            ++rejections;return;
        }
    }
    for (int j=0;j<m;++j) {
        const int i=genes[j],next=proposed[j];
        changes+=z[i]!=next;class_changes+=z[i]/V!=next/V;
        z[i]=next;
    }
    count[a]=next_a;count[b]=next_b;
    nc[a/V]+=next_a-old_a;nc[b/V]+=next_b-old_b;
    if (old_b==0 && next_b>0) {
        position[b]=active.size();active.push_back(b);++births;
    }
    if (old_b>0 && next_b==0) {
        const int at=position[b],last_state=active.back();
        active[at]=last_state;position[last_state]=at;active.pop_back();position[b]=-1;++deaths;
    }
}
template<int RTYPE>
static Rcpp::Matrix<RTYPE> mean_joint_prefix_columns(const Rcpp::Matrix<RTYPE>& x,int columns) {
    Rcpp::Matrix<RTYPE> out(x.nrow(),columns);
    std::copy(x.begin(),x.begin()+size_t(x.nrow())*columns,out.begin());
    return out;
}
template<int RTYPE>
static Rcpp::Matrix<RTYPE> mean_joint_prefix_rows(const Rcpp::Matrix<RTYPE>& x,int rows) {
    Rcpp::Matrix<RTYPE> out(rows,x.ncol());
    for (int j=0;j<x.ncol();++j)
        std::copy(x.begin()+size_t(x.nrow())*j,x.begin()+size_t(x.nrow())*j+rows,
            out.begin()+size_t(rows)*j);
    return out;
}
static long double mean_joint_concentration_change(const std::vector<int>& count,
    const std::vector<int>& occupied,double from,double to) {
    long double out=0;
    for (int k:occupied) out+=std::lgamma(to+count[k])-std::lgamma(to)
        -std::lgamma(from+count[k])+std::lgamma(from);
    return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector mean_joint_prior_concentration_changes_cpp(Rcpp::IntegerMatrix counts,
    double strength_from,double strength_to) {
    const int K=counts.ncol();
    Rcpp::NumericVector out(counts.nrow());
    for (int i=0;i<counts.nrow();++i) {
        std::vector<int> count(K),occupied;
        for (int k=0;k<K;++k) {count[k]=counts(i,k);if (count[k]>0) occupied.push_back(k);}
        out[i]=double(mean_joint_concentration_change(count,occupied,strength_from/K,strength_to/K));
    }
    return out;
}

// The nuisance state is the complete (shape pair, common mean) cell.
// A shared prior keeps this state independent of the mean hypothesis;
// an unshared prior learns the complete joint class/nuisance distribution.
// [[Rcpp::export]]
Rcpp::List mean_joint_nuisance_prior_exchange_cpp(Rcpp::NumericMatrix L,
    double strength,Rcpp::NumericVector A,Rcpp::NumericVector bvalue,
    int warmup,int retained,int batch_length,int save_every,int cluster_every,
    Rcpp::IntegerMatrix initial_labels,int threads,bool shared,bool joint_moves,bool allocation_moves,bool informed_moves,Rcpp::NumericVector concentration_factors,int cold_chains,int exchange_every=16,int tempered_refresh=1) {
    const int n=L.nrow(),K=L.ncol(),V=K/2,nA=A.size(),chains=cold_chains*concentration_factors.size(),
        nh=shared ? V : K,saved=retained/save_every;
    if (K%2 || V%nA || bvalue.size()!=V/nA || strength<=0 ||
        retained%save_every || initial_labels.nrow()!=n)
        Rcpp::stop("Incompatible joint nuisance dimensions");
    if (shared || cold_chains<1 || exchange_every<1 || tempered_refresh<1 || concentration_factors.size()<1 || concentration_factors[0]!=1 ||
        (initial_labels.ncol()!=cold_chains && initial_labels.ncol()!=chains))
        Rcpp::stop("Invalid prior-exchange replica configuration");
    for (int j=0;j<concentration_factors.size();++j)
        if (!std::isfinite(concentration_factors[j]) || concentration_factors[j]<=0 ||
            (j>0 && concentration_factors[j]<=concentration_factors[j-1]))
            Rcpp::stop("Concentration factors must increase from one");
    if (initial_labels.ncol()==cold_chains && chains>cold_chains) {
        Rcpp::IntegerMatrix expanded(n,chains);
        for (int c=0;c<chains;++c) for (int i=0;i<n;++i)
            expanded(i,c)=initial_labels(i,c%cold_chains);
        initial_labels=expanded;
    }
    const int levels=concentration_factors.size();
    Rcpp::NumericMatrix exchange_attempts(std::max(0,levels-1),cold_chains),
        exchange_accepts(std::max(0,levels-1),cold_chains);
    const double alpha=strength/nh,ac=strength/2;
    const double* likelihood=L.begin();
    std::vector<double> base(size_t(n)*2,0);
    std::vector<std::vector<double>> prior_cdf(size_t(n)*2);
    std::vector<std::vector<uint32_t>> prior_index(size_t(n)*2);
    std::vector<std::once_flag> prior_ready(size_t(n)*2);
    const auto started=std::chrono::steady_clock::now();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    #endif
    for (int i=0;i<n;++i) for (int c=0;c<2;++c) {
        long double total=0;
        for (int v=0;v<V;++v) total+=alpha*likelihood[i+size_t(n)*(c*V+v)];
        base[i+size_t(n)*c]=double(total);
    }
    auto prior_state=[&](int i,int c,double u) {
        const size_t row=i+size_t(n)*c;
        std::call_once(prior_ready[row],[&]() {
            auto& cdf=prior_cdf[row];cdf.resize(V);long double total=0;
            size_t nodes=0;double previous=0;
            for (int v=0;v<V;++v) {
                total+=alpha*likelihood[i+size_t(n)*(c*V+v)];cdf[v]=double(total);
                nodes+=cdf[v]>previous;previous=cdf[v];
            }
            if (nodes*(sizeof(double)+sizeof(uint32_t))<size_t(V)*sizeof(double)) {
                std::vector<double> compact;compact.reserve(nodes);
                auto& indices=prior_index[row];indices.reserve(nodes);previous=0;
                for (int v=0;v<V;++v) {
                    if (cdf[v]>previous) {compact.push_back(cdf[v]);indices.push_back(v);}
                    previous=cdf[v];
                }
                cdf.swap(compact);
            }
        });
        const auto& cdf=prior_cdf[row];
        const size_t at=std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin();
        return at>=cdf.size() ? V-1 : (prior_index[row].empty() ? int(at) : int(prior_index[row][at]));
    };
    std::vector<std::vector<int>> z(chains,std::vector<int>(n)),
        count(chains,std::vector<int>(nh,0)),active(chains),position(chains,std::vector<int>(nh,-1)),
        nc(chains,std::vector<int>(2,0));
    for (int c=0;c<chains;++c) for (int i=0;i<n;++i) {
        int state=initial_labels(i,c)-1;
        if (state<0 || state>=K || likelihood[i+size_t(n)*state]<=0) {
            const double p0=base[i],p1=base[i+size_t(n)];
            const int g=R::runif(0,1)*(p0+p1)>=p0;
            state=g*V+prior_state(i,g,R::runif(0,1)*base[i+size_t(n)*g]);
        }
        z[c][i]=state;const int h=shared ? state%V : state;
        if (count[c][h]++==0) {position[c][h]=active[c].size();active[c].push_back(h);}
        ++nc[c][state/V];
    }
    Rcpp::IntegerMatrix labels(saved*cold_chains,n),final_labels(n,chains);
    Rcpp::NumericMatrix U0(n,saved*cold_chains),qmean(n,cold_chains),weights(nh,cold_chains),
        global_trace(retained,6*cold_chains),b_first(n,cold_chains),b_second(n,cold_chains);
    std::vector<double> uniforms(size_t(n)*2*chains),minimum(chains,R_PosInf),count_error(chains,0),
        cluster_attempts(chains,0),cluster_changes(chains,0),
        joint_attempts(chains,0),joint_changes(chains,0),joint_class_changes(chains,0),
        allocation_attempts(chains,0),allocation_changes(chains,0),allocation_class_changes(chains,0),
        allocation_births(chains,0),allocation_deaths(chains,0),allocation_rejections(chains,0);
    // Concentration changes depend only on an integer occupancy count.
    // Retain the original arithmetic and occupied-state summation order.
    std::vector<std::vector<double>> concentration_delta(2*std::max(0,levels-1),
        std::vector<double>(n+1));
    for (int level=0;level+1<levels;++level) for (int direction=0;direction<2;++direction) {
        const double from=alpha*concentration_factors[level+direction],
            to=alpha*concentration_factors[level+1-direction];
        for (int m=0;m<=n;++m)
            concentration_delta[2*level+direction][m]=std::lgamma(to+m)-std::lgamma(to)
                -std::lgamma(from+m)+std::lgamma(from);
    }
    auto concentration_change=[&](int chain,int index) {
        long double result=0;
        for (int h:active[chain]) result+=concentration_delta[index][count[chain][h]];
        return result;
    };
    const double prepare=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    for (int iteration=-warmup;iteration<retained;++iteration) {
        const bool refresh_auxiliary=iteration%tempered_refresh==0;
        const size_t draws=size_t(n)*2*(refresh_auxiliary ? chains : cold_chains);
        for (size_t draw=0;draw<draws;++draw) uniforms[draw]=R::runif(0,1);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static,1) num_threads(threads)
        #endif
        for (int c=0;c<chains;++c) {
            if (c>=cold_chains && !refresh_auxiliary) continue;
            auto& counts=count[c];auto& occupied=active[c];auto& pos=position[c];
            const double factor=concentration_factors[c/cold_chains];
            std::vector<double> cdf0,cdf1;
            for (int i=0;i<n;++i) {
                const int old=z[c][i],oldh=shared ? old%V : old;
                if (--counts[oldh]==0) {
                    const int at=pos[oldh],last=occupied.back();occupied[at]=last;pos[last]=at;
                    occupied.pop_back();pos[oldh]=-1;
                }
                --nc[c][old/V];
                double s0=factor*base[i],s1=factor*base[i+size_t(n)];
                cdf0.resize(occupied.size());cdf1.resize(occupied.size());
                for (int j=0;j<int(occupied.size());++j) {
                    const int h=occupied[j],v=shared ? h : h%V;
                    if (shared || h<V) s0+=counts[h]*likelihood[i+size_t(n)*v];
                    if (shared || h>=V) s1+=counts[h]*likelihood[i+size_t(n)*(V+v)];
                    cdf0[j]=s0;cdf1[j]=s1;
                }
                const double d0=s0*(shared ? nc[c][0]+ac : 1),
                    d1=s1*(shared ? nc[c][1]+ac : 1),total=d0+d1;
                minimum[c]=std::min(minimum[c],total);
                const double q=d0/total;
                if (iteration>=0 && c<cold_chains) {
                    qmean(i,c)+=q/retained;
                    if (iteration%save_every==0) U0(i,iteration/save_every+saved*c)=q;
                }
                const size_t at=size_t(n)*2*c+2*i;
                const int g=uniforms[at]*total>=d0;
                const double u=uniforms[at+1]*(g ? s1 : s0),prior=factor*base[i+size_t(n)*g];
                int v;
                if (u<prior) v=prior_state(i,g,u/factor);
                else {
                    const auto& cdf=g ? cdf1 : cdf0;
                    const int j=std::min(int(occupied.size())-1,
                        int(std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin()));
                    v=occupied[j]%V;
                }
                const int state=g*V+v,h=shared ? v : state;
                z[c][i]=state;
                if (counts[h]++==0) {pos[h]=occupied.size();occupied.push_back(h);}
                ++nc[c][g];
            }
        }
        // Conditional on the occupied partition, a move to an empty full cell
        // preserves all Dirichlet factors. Propose from one member's complete
        // likelihood; its term cancels in the Metropolis ratio. This moves
        // shape, level and class together without scanning the full dictionary.
        if (joint_moves && !shared && (iteration+1)%save_every==0)
            for (int c=0;c<chains;++c) {
                std::vector<std::vector<int>> members(active[c].size());
                for (int i=0;i<n;++i) members[position[c][z[c][i]]].push_back(i);
                // Fix the scan order within this unchanged gene partition.
                std::sort(members.begin(),members.end(),[](const auto& a,const auto& b) {
                    return a.front()<b.front();
                });
                for (const auto& genes:members) {
                    const int m=genes.size(),old=z[c][genes[0]],
                        representative=genes[std::min(m-1,int(R::runif(0,1)*m))];
                    ++joint_attempts[c];
                    const double p0=base[representative],p1=base[representative+n],
                        u=R::runif(0,1)*(p0+p1);
                    const int g=u>=p0,next=g*V+prior_state(representative,g,u-(g ? p0 : 0));
                    if (next==old || count[c][next]>0) continue;
                    long double logratio=0;
                    for (int i:genes) if (i!=representative) {
                        const double value=likelihood[i+size_t(n)*next];
                        if (value<=0) {logratio=R_NegInf;break;}
                        logratio+=std::log(value)-std::log(likelihood[i+size_t(n)*old]);
                    }
                    if (std::log(R::runif(0,1))>=logratio) continue;
                    const int at=position[c][old];
                    active[c][at]=next;position[c][next]=at;position[c][old]=-1;
                    count[c][old]=0;count[c][next]=m;
                    nc[c][old/V]-=m;nc[c][g]+=m;
                    for (int i:genes) z[c][i]=next;
                    ++joint_changes[c];joint_class_changes[c]+=g!=old/V;
                }
            }
        if (allocation_moves && !shared && (iteration+1)%cluster_every==0)
            for (int c=0;c<chains;++c) {
                const int anchor=std::min(n-1,int(R::runif(0,1)*n)),a=z[c][anchor],opposite=1-a/V;
                const double total=base[anchor+size_t(n)*opposite];
                if (total>0) {
                    const int b=opposite*V+prior_state(anchor,opposite,R::runif(0,1)*total);
                    mean_joint_allocation_block(likelihood,n,V,alpha*concentration_factors[c/cold_chains],z[c],count[c],active[c],position[c],nc[c],
                        anchor,-1,a,b,allocation_attempts[c],allocation_changes[c],allocation_class_changes[c],
                        allocation_births[c],allocation_deaths[c],0,allocation_rejections[c]);
                }
                int other;double selection_total=0;
                if (informed_moves) {
                    std::vector<double> cdf(n);
                    for (int j=0;j<n;++j) {
                        if (z[c][j]/V==opposite) selection_total+=likelihood[anchor+size_t(n)*z[c][j]];
                        cdf[j]=selection_total;
                    }
                    if (selection_total<=0) continue;
                    const double u=std::min(std::nextafter(selection_total,0.0),R::runif(0,1)*selection_total);
                    other=std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin();
                } else other=std::min(n-1,int(R::runif(0,1)*n));
                const int b=z[c][other];
                if (b/V!=a/V) mean_joint_allocation_block(likelihood,n,V,alpha*concentration_factors[c/cold_chains],
                    z[c],count[c],active[c],position[c],nc[c],anchor,other,a,b,
                    allocation_attempts[c],allocation_changes[c],allocation_class_changes[c],
                    allocation_births[c],allocation_deaths[c],selection_total,allocation_rejections[c]);
            }
        if ((iteration+1)%cluster_every==0 && !shared) for (int c=0;c<chains;++c)
            mean_joint_coarse_heatbath(likelihood,n,V,nA,z[c],count[c],active[c],position[c],nc[c],
                cluster_attempts[c],cluster_changes[c],threads);
        if ((iteration+1)%cluster_every==0 && shared) for (int c=0;c<chains;++c) {
            auto& counts=count[c];auto& occupied=active[c];auto& pos=position[c];
            const double factor=concentration_factors[c/cold_chains];
            std::vector<std::vector<int>> members(K);std::vector<int> cells(K,0),partition;
            for (int i=0;i<n;++i) {
                const int state=z[c][i];
                if (cells[state]++==0) partition.push_back(state);
                members[state].push_back(i);
            }
            std::vector<double> lp(K),cdf(K);
            for (int id:partition) {
                const auto& genes=members[id];const int old=z[c][genes[0]],m=genes.size(),
                    oh=shared ? old%V : old,oc=old/V;
                cells[old]-=m;counts[oh]-=m;nc[c][oc]-=m;
                if (counts[oh]==0) {
                    const int at=pos[oh],last=occupied.back();occupied[at]=last;pos[last]=at;
                    occupied.pop_back();pos[oh]=-1;
                }
                const double empty=std::lgamma(alpha+m)-std::lgamma(alpha);
                double class_term[2]={0,0};
                if (shared) for (int g=0;g<2;++g)
                    class_term[g]=std::lgamma(nc[c][g]+ac+m)-std::lgamma(nc[c][g]+ac);
                double peak=R_NegInf;
                for (int k=0;k<K;++k) {
                    if (cells[k]>0) {lp[k]=R_NegInf;continue;}
                    double value=shared ? class_term[k/V] : 0;
                    if (shared) {
                        const int number=counts[k%V];
                        value+=number==0 ? empty : std::lgamma(number+alpha+m)-std::lgamma(number+alpha);
                    }
                    for (int i:genes) {
                        const double p=likelihood[i+size_t(n)*k];
                        if (p<=0) {value=R_NegInf;break;}
                        value+=std::log(p);
                    }
                    lp[k]=value;peak=std::max(peak,value);
                }
                double total=0;
                for (int k=0;k<K;++k) {total+=std::exp(lp[k]-peak);cdf[k]=total;}
                const int next=std::min(K-1,int(std::upper_bound(cdf.begin(),cdf.end(),R::runif(0,1)*total)-cdf.begin())),
                    h=shared ? next%V : next;
                if (counts[h]==0) {pos[h]=occupied.size();occupied.push_back(h);}
                counts[h]+=m;nc[c][next/V]+=m;cells[next]+=m;
                for (int i:genes) z[c][i]=next;
                ++cluster_attempts[c];cluster_changes[c]+=next!=old;
            }
        }
        // Alternating adjacent swaps are symmetric on the replica product.
        // The likelihood and all concentration normalizers cancel. Empty
        // cells contribute zero to the integrated Dirichlet ratio.
        if ((iteration+1)%exchange_every==0) {
            const int phase=((iteration+warmup)/exchange_every)%2;
            for (int level=phase;level+1<levels;level+=2) for (int c=0;c<cold_chains;++c) {
                const int left=level*cold_chains+c,right=left+cold_chains;
                const double a=alpha*concentration_factors[level],
                    b=alpha*concentration_factors[level+1];
                const long double ratio=
                    concentration_change(right,2*level+1)+
                    concentration_change(left,2*level);
                ++exchange_attempts(level,c);
                if (std::log(R::runif(0,1))<ratio) {
                    std::swap(z[left],z[right]);std::swap(count[left],count[right]);
                    std::swap(active[left],active[right]);std::swap(position[left],position[right]);
                    std::swap(nc[left],nc[right]);++exchange_accepts(level,c);
                }
            }
        }
        if (iteration>=0) for (int c=0;c<cold_chains;++c) {
            double ma=0,mb=0,mab=0,loglike=0;
            for (int i=0;i<n;++i) {
                const int state=z[c][i],v=state%V;const double a=A[v%nA],b=bvalue[v/nA];
                ma+=a;mb+=b;mab+=a*b;loglike+=std::log(likelihood[i+size_t(n)*state]);
                b_first(i,c)+=b/retained;b_second(i,c)+=b*b/retained;
                if (iteration%save_every==0) labels(iteration/save_every+saved*c,i)=state+1;
            }
            global_trace(iteration,6*c)=ma/n;global_trace(iteration,6*c+1)=mb/n;
            global_trace(iteration,6*c+2)=loglike;global_trace(iteration,6*c+3)=active[c].size();
            global_trace(iteration,6*c+4)=double(nc[c][1])/n;
            global_trace(iteration,6*c+5)=mab/n-ma*mb/(double(n)*n);
            if (iteration%save_every==0) for (int h=0;h<nh;++h)
                weights(h,c)+=(count[c][h]+alpha*concentration_factors[c/cold_chains])/((n+strength*concentration_factors[c/cold_chains])*saved);
        }
        if ((iteration+1)%batch_length==0) {
            Rcpp::Rcout<<"joint_nuisance "<<(shared ? "shared " : "class_specific ")
                <<(iteration<0 ? "warmup_remaining=" : "retained=")
                <<(iteration<0 ? -(iteration+1) : iteration+1)<<" occupied/alternative=";
            for (int c=0;c<chains;++c) {
                std::vector<int> expected(nh,0);int alternatives=0;
                for (int state:z[c]) {++expected[shared ? state%V : state];alternatives+=state>=V;}
                for (int h=0;h<nh;++h) count_error[c]=std::max(count_error[c],double(std::abs(expected[h]-count[c][h])));
                count_error[c]=std::max(count_error[c],double(std::abs(alternatives-nc[c][1])));
                Rcpp::Rcout<<active[c].size()<<"/"<<alternatives<<(c+1<chains ? "," : "\n");
            }
            Rcpp::checkUserInterrupt();
        }
    }
    for (int c=0;c<chains;++c) for (int i=0;i<n;++i) final_labels(i,c)=z[c][i]+1;
    Rcpp::List out=Rcpp::List::create(Rcpp::_["labels"]=labels,Rcpp::_["U0"]=U0,
        Rcpp::_["null_probability"]=qmean,Rcpp::_["weights"]=weights,
        Rcpp::_["global_trace"]=global_trace,Rcpp::_["b_first"]=b_first,Rcpp::_["b_second"]=b_second,
        Rcpp::_["final_labels"]=final_labels,Rcpp::_["maximum_count_error"]=Rcpp::wrap(count_error),
        Rcpp::_["minimum_categorical_mass"]=Rcpp::wrap(minimum),
        Rcpp::_["cluster_attempts"]=Rcpp::wrap(cluster_attempts),Rcpp::_["cluster_changes"]=Rcpp::wrap(cluster_changes),
        Rcpp::_["prepare_seconds"]=prepare,
        Rcpp::_["total_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count(),
        Rcpp::_["prior_cdf_entries"]=std::accumulate(prior_cdf.begin(),prior_cdf.end(),0.0,
            [](double total,const std::vector<double>& v){return total+v.size();}),
        Rcpp::_["joint_attempts"]=Rcpp::wrap(joint_attempts),
        Rcpp::_["joint_changes"]=Rcpp::wrap(joint_changes),
        Rcpp::_["joint_class_changes"]=Rcpp::wrap(joint_class_changes),
        Rcpp::_["allocation_attempts"]=Rcpp::wrap(allocation_attempts),
        Rcpp::_["allocation_changes"]=Rcpp::wrap(allocation_changes),
        Rcpp::_["allocation_class_changes"]=Rcpp::wrap(allocation_class_changes),
        Rcpp::_["allocation_births"]=Rcpp::wrap(allocation_births),
        Rcpp::_["allocation_deaths"]=Rcpp::wrap(allocation_deaths),
        Rcpp::_["allocation_rejections"]=Rcpp::wrap(allocation_rejections),
        Rcpp::_["shared_nuisance"]=shared);
    out["auxiliary_metropolis"]=false;
    out["cdf_storage_bytes"]=std::accumulate(prior_cdf.begin(),prior_cdf.end(),0.0,
        [](double total,const std::vector<double>& x){return total+sizeof(double)*x.size();})+
        std::accumulate(prior_index.begin(),prior_index.end(),0.0,
        [](double total,const std::vector<uint32_t>& x){return total+sizeof(uint32_t)*x.size();});
    out["replica_final_labels"]=final_labels;
    out["concentration_factors"]=concentration_factors;
    out["exchange_every"]=exchange_every;
    out["tempered_refresh"]=tempered_refresh;
    out["exchange_attempts"]=exchange_attempts;
    out["exchange_accepts"]=exchange_accepts;
    out["labels"]=mean_joint_prefix_rows(labels,saved*cold_chains);
    out["U0"]=mean_joint_prefix_columns(U0,saved*cold_chains);
    out["null_probability"]=mean_joint_prefix_columns(qmean,cold_chains);
    out["weights"]=mean_joint_prefix_columns(weights,cold_chains);
    out["global_trace"]=mean_joint_prefix_columns(global_trace,6*cold_chains);
    out["b_first"]=mean_joint_prefix_columns(b_first,cold_chains);
    out["b_second"]=mean_joint_prefix_columns(b_second,cold_chains);
    out["final_labels"]=mean_joint_prefix_columns(final_labels,cold_chains);
    return out;
}

