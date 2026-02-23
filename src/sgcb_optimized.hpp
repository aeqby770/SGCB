// =============================================================================
// SGCB optimization header - modern C++17 performance optimizations
// =============================================================================
#ifndef SGCB_OPTIMIZED_HPP
#define SGCB_OPTIMIZED_HPP

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

// =============================================================================
// Compiler optimization hint macros
// =============================================================================

#if defined(__GNUC__) || defined(__clang__)
    #define SGCB_LIKELY(x)   __builtin_expect(!!(x), 1)
    #define SGCB_UNLIKELY(x) __builtin_expect(!!(x), 0)
    #define SGCB_RESTRICT    __restrict__
    #define SGCB_ALWAYS_INLINE __attribute__((always_inline)) inline
    #define SGCB_HOT         __attribute__((hot))
    #define SGCB_PURE        __attribute__((pure))
    #define SGCB_CONST       __attribute__((const))
#else
    #define SGCB_LIKELY(x)   (x)
    #define SGCB_UNLIKELY(x) (x)
    #define SGCB_RESTRICT
    #define SGCB_ALWAYS_INLINE inline
    #define SGCB_HOT
    #define SGCB_PURE
    #define SGCB_CONST
#endif

// =============================================================================
// Constexpr math functions
// =============================================================================

namespace sgcb {

constexpr double EPS = 1e-8;
constexpr double LOG2_E = 1.4426950408889634;  // 1/ln(2)
constexpr double LN_2 = 0.6931471805599453;

// Compile-time constants
template<typename T>
SGCB_CONST constexpr T sqr(T x) noexcept { return x * x; }

template<typename T>
SGCB_CONST constexpr T cube(T x) noexcept { return x * x * x; }

// Fast log2 (for LFC computation)
SGCB_ALWAYS_INLINE SGCB_PURE
double fast_log2(double x) noexcept {
    return std::log(x) * LOG2_E;
}

// Safe log (with epsilon)
SGCB_ALWAYS_INLINE SGCB_PURE
double safe_log(double x) noexcept {
    return std::log(x > EPS ? x : EPS);
}

// Safe log2 (with epsilon)
SGCB_ALWAYS_INLINE SGCB_PURE
double safe_log2(double x) noexcept {
    return fast_log2(x > EPS ? x : EPS);
}

// Fast softplus: log(1 + exp(x))
SGCB_ALWAYS_INLINE SGCB_PURE
double softplus(double x) noexcept {
    return SGCB_LIKELY(x > 20.0) ? x : std::log1p(std::exp(x));
}

// Fast sigmoid: 1 / (1 + exp(-x))
SGCB_ALWAYS_INLINE SGCB_PURE
double sigmoid(double x) noexcept {
    return SGCB_LIKELY(x > 20.0) ? 1.0 : 1.0 / (1.0 + std::exp(-x));
}

// Clipping function
SGCB_ALWAYS_INLINE SGCB_CONST
double clip(double x, double lo, double hi) noexcept {
    return x < lo ? lo : (x > hi ? hi : x);
}

// =============================================================================
// GG distribution core functions (inline optimized)
// =============================================================================

// GG log-likelihood (single point)
SGCB_ALWAYS_INLINE SGCB_HOT
double gg_loglik_scalar(double x, double alpha, double beta, double gamma) noexcept {
    const double xi = x > EPS ? x : EPS;
    const double log_xi = std::log(xi);
    const double log_beta = std::log(beta);
    const double ratio_g = std::pow(xi / beta, gamma);
    
    return std::log(gamma) - gamma * alpha * log_beta - std::lgamma(alpha) +
           (gamma * alpha - 1.0) * log_xi - ratio_g;
}

// GG log-likelihood (pre-computed version, avoids redundant computation)
struct GGParams {
    double alpha;
    double beta;
    double gamma;
    double log_gamma;
    double log_beta;
    double lgamma_alpha;
    double coef;  // gamma * alpha - 1
    
    SGCB_ALWAYS_INLINE
    void precompute(double a, double b, double g) noexcept {
        alpha = a;
        beta = b;
        gamma = g;
        log_gamma = std::log(g);
        log_beta = std::log(b);
        lgamma_alpha = std::lgamma(a);
        coef = g * a - 1.0;
    }
    
    SGCB_ALWAYS_INLINE SGCB_HOT
    double loglik(double x) const noexcept {
        const double xi = x > EPS ? x : EPS;
        const double log_xi = std::log(xi);
        const double ratio_g = std::pow(xi / beta, gamma);
        return log_gamma - gamma * alpha * log_beta - lgamma_alpha + coef * log_xi - ratio_g;
    }
};

// =============================================================================
// Vectorization helper functions
// =============================================================================

// Parallel summation (Kahan summation for improved accuracy)
template<typename Iter>
SGCB_HOT
double parallel_sum(Iter begin, Iter end) {
    double sum = 0.0;
    double c = 0.0;  // compensation term
    
    for (auto it = begin; it != end; ++it) {
        double y = *it - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

// Batch statistics computation (mean + variance, single pass)
struct MeanVar {
    double mean;
    double var;
};

SGCB_ALWAYS_INLINE
MeanVar compute_mean_var(const double* SGCB_RESTRICT data, int n) noexcept {
    double sum = 0.0;
    double sum_sq = 0.0;
    
    // Manual loop unrolling (4x)
    int i = 0;
    for (; i + 3 < n; i += 4) {
        sum += data[i] + data[i+1] + data[i+2] + data[i+3];
        sum_sq += sqr(data[i]) + sqr(data[i+1]) + sqr(data[i+2]) + sqr(data[i+3]);
    }
    for (; i < n; ++i) {
        sum += data[i];
        sum_sq += sqr(data[i]);
    }
    
    double mean = sum / n;
    double var = sum_sq / n - sqr(mean);
    return {mean, var};
}

// =============================================================================
// Memory pool (avoid frequent allocations)
// =============================================================================

template<typename T>
class MemoryPool {
private:
    std::vector<T> buffer_;
    size_t capacity_;
    
public:
    explicit MemoryPool(size_t capacity = 0) : capacity_(capacity) {
        if (capacity > 0) buffer_.reserve(capacity);
    }
    
    void reserve(size_t n) {
        if (n > capacity_) {
            buffer_.reserve(n);
            capacity_ = n;
        }
    }
    
    T* get(size_t n) {
        buffer_.resize(n);
        return buffer_.data();
    }
    
    void clear() { buffer_.clear(); }
};

// =============================================================================
// Fast median (nth_element)
// =============================================================================

SGCB_ALWAYS_INLINE
double fast_median(std::vector<double>& v) noexcept {
    const size_t n = v.size();
    const size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return v[mid];
}

// =============================================================================
// Fast quantile (partial sort)
// =============================================================================

SGCB_ALWAYS_INLINE
double fast_quantile(std::vector<double>& v, double p) noexcept {
    const size_t n = v.size();
    const double idx_d = p * (n - 1);
    const size_t idx_lo = static_cast<size_t>(idx_d);
    const size_t idx_hi = std::min(idx_lo + 1, n - 1);
    
    std::nth_element(v.begin(), v.begin() + idx_hi, v.end());
    if (idx_lo < idx_hi) {
        std::nth_element(v.begin(), v.begin() + idx_lo, v.begin() + idx_hi);
    }
    
    const double frac = idx_d - idx_lo;
    return v[idx_lo] * (1.0 - frac) + v[idx_hi] * frac;
}

} // namespace sgcb

#endif // SGCB_OPTIMIZED_HPP
