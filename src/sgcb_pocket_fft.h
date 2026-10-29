#ifndef SGCB_POCKET_FFT_H
#define SGCB_POCKET_FFT_H
// PocketFFT, Martin Reinecke's contribution mirror, commit
// 4d865cba24359b05b37da4e75aee24498c3b78e7. BSD-3-Clause notices are
// retained in pocketfft_hdronly.h. Only the real one-dimensional interface
// used by the existing mean integrals is adapted here.
#define POCKETFFT_NAMESPACE sgcb_pocketfft
#define POCKETFFT_NO_MULTITHREADING
// Share immutable plans for the power-of-two spectral grids. Each outer
// OpenMP worker keeps its own plan lookup and mutable data buffer.
#define POCKETFFT_CACHE_SIZE 16
#include "pocketfft_hdronly.h"
#include <map>
#include <memory>
#include <complex>
#include <vector>
#include <algorithm>

class SGCBFFT {
    using Plan=sgcb_pocketfft::detail::pocketfft_r<double>;
    std::map<size_t,std::shared_ptr<Plan>> plans;
    std::vector<double> work;
    const Plan& plan(size_t n) {
        auto& value=plans[n];
        if(!value) value=sgcb_pocketfft::detail::get_plan<Plan>(n);
        return *value;
    }
public:
    void fwd(std::vector<std::complex<double>>& out,const std::vector<double>& input) {
        const size_t n=input.size();work.assign(input.begin(),input.end());
        plan(n).exec(work.data(),1.0,true);out.resize(n);out[0]={work[0],0};
        size_t i=1,k=1;
        for(;i<n-1;i+=2,++k) out[k]={work[i],work[i+1]};
        if(i<n) out[k]={work[i],0};
        for(size_t j=n/2+1;j<n;++j) out[j]=std::conj(out[n-j]);
    }
    void inv(std::vector<double>& out,const std::vector<std::complex<double>>& input) {
        const size_t n=input.size();work.resize(n);work[0]=input[0].real();
        size_t i=1,k=1;
        for(;i<n-1;i+=2,++k) {work[i]=input[k].real();work[i+1]=input[k].imag();}
        if(i<n) work[i]=input[k].real();
        plan(n).exec(work.data(),1.0,false);out.resize(n);
        const double inverse=1.0/n;
        for(size_t j=0;j<n;++j) out[j]=work[j]*inverse;
    }
};
#endif
