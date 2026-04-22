#ifndef SGCB_FAST_SPECIAL_H
#define SGCB_FAST_SPECIAL_H

#include <cmath>

// =============================================================================
// Fast special function approximations (avoid frequent expensive R API calls)
// Chebyshev polynomial approximations for Digamma/Trigamma
// =============================================================================

namespace fast_special {

// Digamma approximation: psi(x) ~ ln(x) - 1/(2x) - 1/(12x^2) + 1/(120x^4)
// Accuracy: |error| < 1e-8 for x > 0.5
inline double digamma(double x) {
    if (x < 0.5) {
        // Reflection formula: psi(1-x) - psi(x) = pi*cot(pi*x)
        return digamma(1.0 - x) - 3.14159265358979323846 / std::tan(3.14159265358979323846 * x);
    }
    if (x < 6.0) {
        // Recurrence to large-value region: psi(x) = psi(x+1) - 1/x
        double result = 0;
        while (x < 6.0) {
            result -= 1.0 / x;
            x += 1.0;
        }
        return result + digamma(x);
    }
    // Asymptotic expansion for x >= 6
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    return std::log(x) - 0.5 * inv_x 
           - inv_x2 * (1.0/12.0 - inv_x2 * (1.0/120.0 - inv_x2 / 252.0));
}

// Trigamma approximation: psi'(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5)
// Accuracy: |error| < 1e-8 for x > 0.5
inline double trigamma(double x) {
    if (x < 0.5) {
        // Reflection formula
        double pi = 3.14159265358979323846;
        double s = std::sin(pi * x);
        return (pi * pi) / (s * s) - trigamma(1.0 - x);
    }
    if (x < 6.0) {
        double result = 0;
        while (x < 6.0) {
            result += 1.0 / (x * x);
            x += 1.0;
        }
        return result + trigamma(x);
    }
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    return inv_x + 0.5 * inv_x2 
           + inv_x2 * inv_x * (1.0/6.0 - inv_x2 * (1.0/30.0 - inv_x2 / 42.0));
}

// lgamma approximation: Stirling's formula + recurrence
// Accuracy: |error| < 1e-7 for x > 0.5
inline double lgamma_stirling(double x) {
    // Stirling approximation for x >= 6
    double v = 1.0 / x;
    double s = v * v;
    return (x - 0.5) * std::log(x) - x + 0.91893853320467274178  // 0.5*log(2π)
           + v * (1.0/12.0 - s * (1.0/360.0 - s * (1.0/1260.0 - s/1680.0)));
}

inline double lgamma(double x) {
    if (x < 0.5) {
        // Reflection formula: Gamma(x)*Gamma(1-x) = pi/sin(pi*x)
        // lgamma(x) = log(pi) - log(sin(pi*x)) - lgamma(1-x)
        double pi = 3.14159265358979323846;
        return std::log(pi) - std::log(std::sin(pi * x)) - lgamma(1.0 - x);
    }
    if (x >= 6.0) {
        return lgamma_stirling(x);
    }
    // Recurrence formula: lgamma(x) = lgamma(x+1) - log(x)
    double result = 0;
    while (x < 6.0) {
        result -= std::log(x);
        x += 1.0;
    }
    return result + lgamma_stirling(x);
}

// Tetragamma: psi''(x) = d/dx psi'(x), used for Newton iteration of trigamma inverse
// Asymptotic expansion + recurrence, accuracy ~1e-8
inline double tetragamma(double x) {
    // ψ''(x) < 0 for all x > 0
    // Recurrence: ψ''(x) = ψ''(x+1) - 2/x³
    double result = 0.0;
    while (x < 8.0) {
        result -= 2.0 / (x * x * x);
        x += 1.0;
    }
    // Asymptotic expansion: ψ''(x) = -1/x² - 1/x³ - 1/(2x⁴) + 1/(6x⁶)
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x2 * x2;
    double x6 = x3 * x3;
    result += -1.0/x2 - 1.0/x3 - 0.5/x4 + 1.0/(6.0*x6);
    return result;
}

// Trigamma inverse: solve psi'(y) = x (Newton-Raphson)
// Reference: limma::trigammaInverse (Smyth 2004)
inline double trigamma_inverse(double x) {
    // Solve ψ'(y) = x for y, using Newton on 1/ψ'(y) = 1/x
    // Newton step: Δy = ψ'(y)·(1 - ψ'(y)/x) / (-ψ''(y))
    // Reference: limma::trigammaInverse (Smyth 2004)
    const double eps = 1e-8;
    if (x > 1e7) return 1.0 / std::sqrt(x);
    if (x < 1e-6) return 1.0 / x;
    
    double y = 0.5 + 1.0 / x;
    
    for (int iter = 0; iter < 50; iter++) {
        double tri = trigamma(y);
        double tet = tetragamma(y);  // tet < 0
        if (std::abs(tet) < eps) break;
        // Newton on g(y) = 1/ψ'(y) - 1/x:
        // Δy = -g/g' = (1/ψ'-1/x)·ψ'²/ψ'' = ψ'·(1-ψ'/x)/ψ''
        // ψ'' < 0 naturally provides correct direction
        double dif = tri * (1.0 - tri / x) / tet;
        y += dif;
        y = std::max(y, eps);
        if (std::abs(dif) < 1e-8 * y) break;
    }
    return y;
}

} // namespace fast_special

#endif // SGCB_FAST_SPECIAL_H
