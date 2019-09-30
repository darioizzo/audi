#ifndef AUDI_FUNCTIONS_HPP
#define AUDI_FUNCTIONS_HPP

#include <audi/config.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <cmath>
#include <stdexcept>

#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>
#include <audi/detail/overloads_vectorized.hpp>
#include <audi/gdual.hpp>
#include <audi/type_traits.hpp>
#include <audi/vectorized.hpp>

namespace audi
{
/// Overload for the exponential
/**
 * Implements the exponential of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ...
 * \right) \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the exponential of \p d
 */
template <typename T>
inline gdual<T> exp(const gdual<T> &d)
{
    gdual<T> retval(1.);
    gdual<T> fact(1.);
    auto p0 = d.constant_cf();
    auto phat = d - p0;
    gdual<T> tmp(phat);

    retval += phat;
    for (auto i = 2u; i <= d.get_order(); ++i) {
        phat *= tmp;
        fact *= i;
        retval += phat / fact;
    }
    return retval * audi::exp(p0);
}

/// Overload for the logarithm
/**
 * Implements the logarithm of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat
 * f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ... \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the logarithm of \p d
 *
 */
template <typename T>
inline gdual<T> log(const gdual<T> &d)
{
    gdual<T> retval(0.);
    gdual<T> fatt(1.);
    auto p0 = d.constant_cf();
    auto log_p0 = audi::log(p0);

    auto phat = (d - p0);
    phat = phat / p0;
    gdual<T> tmp(phat);

    retval = log_p0 + phat;
    for (auto i = 2u; i <= d.get_order(); ++i) {
        fatt *= -1.;
        phat *= tmp;
        retval = retval + fatt * phat / i;
    }
    return retval;
}

/// Overload for the exponentiation to a gdual power
/**
 * Computes the exponentiation to the power of an audi::gdual.
 * If the exponent is an integer constant, it calls the std::pow overload. Otherwise
 * it converts \f$a^{T_f}\f$ to \f$\exp(T_g*\log(a))\f$ and computes this
 * last expression instead.
 *
 * @param base the base for the exponent
 * @param d audi::gdual argument
 *
 */
template <typename T, typename U, enable_if_t<audi::is_arithmetic<U>::value, int> = 0>
inline gdual<T> pow(U base, const gdual<T> &d)
{
    // checks wether the exponent is a constant in which
    // case it calls for a different overload
    if (d.degree() == 0) {
        auto p0 = d.constant_cf();
        return gdual<T>(audi::pow(base, p0));
    }
    return exp(audi::log(T(base)) * d);
}

template <typename T>
T binomial(T x, unsigned y)
{
    T retval(x);
    for (auto i = 1u; i < y; ++i) {
        retval = retval * (x - i) / (i + 1);
    }
    return retval;
}

/// Overload for the exponentiation
/**
 * Implements the exponentiation of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(f^\alpha)} = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\hat f / f_0\right)^k
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 * @param alpha exponent
 *
 * @return an audi:gdual containing the Taylor expansion of \p d elevated to the power \p alpha
 *
 */
template <typename T, typename U, enable_if_t<audi::is_arithmetic<U>::value, int> = 0>
inline gdual<T> pow(const gdual<T> &d, U alpha)
{
    // We check if the exponent is representable as a positive integer,
    // in which case we just do d*d*d*d*... etc.
    int n = static_cast<int>(alpha);
    if (n == alpha && alpha > 0.) {
        gdual<T> retval(d);
        for (auto i = 1; i < n; ++i) {
            retval *= d;
        }
        return retval;
    } else {
        auto p0 = d.constant_cf();
        auto phat = d - p0;
        gdual<T> retval(audi::pow(p0, alpha));
        gdual<T> tmp(phat);
        retval += alpha * phat * audi::pow(p0, alpha - 1);
        for (auto i = 2u; i <= d.get_order(); ++i) {
            phat *= tmp;
            auto achoosei = binomial(alpha, i);
            retval += achoosei * phat * audi::pow(p0, alpha - i);
        }
        return retval;
    }
}

/// Overload for the integer exponentiation
/**
 * Implements the integer exponentiation of a audi::gdual. Essentially,
 * it uses the \f$\mathcal P_{n,m}\f$ multiplication on \p d \p n times
 *
 * @param d audi::gdual argument
 * @param n integer exponent
 */
template <typename T>
inline gdual<T> pow(const gdual<T> &d, int n)
{
    if (n < 0) {
        gdual<T> retval(T(1.) / d);
        for (auto i = 1; i < -n; ++i) {
            retval /= d;
        }
        return retval;
    } else if (n > 0) {
        gdual<T> retval(d);
        for (auto i = 1; i < n; ++i) {
            retval *= d;
        }
        return retval;
    } else {              //==0 case
        return d - d + 1; // TODO: can this be made more efficient? (preserving symbol set and order)
    }
}

/// Overload for the exponentiation of a gdual to a gdual power
/**
 * Computes the exponentiation of an audi::gdual to the power of an audi::gdual.
 * If the exponent is an integer constant, it calls a different overload. Otherwise
 * it converts \f$T_f^{T_g}\f$ to \f$\exp(T_g*\log(T_f))\f$ and computes this
 * last expression instead.
 *
 * @param d1 audi::gdual argument
 * @param d2 audi::gdual argument
 *
 * @throw std::domain_error if std::log(\f$f_0\f$) is not finite (uses std::isfinite)
 */
template <typename T>
inline gdual<T> pow(const gdual<T> &d1, const gdual<T> &d2)
{
    return exp(d2 * log(d1));
}

/// Overload for the square root
/**
 * Implements the square root of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{\sqrt{f}} = \sqrt{f_0} \sum_{k=0}^m {\frac 12 \choose k} \left(\hat f / f_0\right)^k
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the square root of \p d
 *
 */
template <typename T>
inline gdual<T> sqrt(const gdual<T> &d)
{
    T alpha(0.5);
    gdual<T> retval(1.);
    T p0 = d.constant_cf();
    T sqrt_p0 = audi::sqrt(p0);

    gdual<T> phat = d - p0;
    phat = phat / p0;
    gdual<T> tmp(phat);

    retval += alpha * phat;
    for (decltype(d.get_order()) i = 2u; i <= d.get_order(); ++i) {
        phat *= tmp;
        retval += binomial(alpha, i) * phat;
    }
    retval *= sqrt_p0;
    return retval;
}

/// Overload for the cubic root
/**
 * Implements the cubic root of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{\sqrt[3]{f}} = \sqrt[3]{f_0} \sum_{k=0}^m {\frac 13 \choose k} \left(\hat f / f_0\right)^k
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the square root of \p d
 *
 */
template <typename T>
inline gdual<T> cbrt(const gdual<T> &d)
{
    T alpha(T(1.) / T(3.));
    gdual<T> retval(1.);
    T p0 = d.constant_cf();
    T cbrt_p0 = audi::cbrt(p0);

    gdual<T> phat = d - p0;
    phat = phat / p0;
    gdual<T> tmp(phat);

    retval += alpha * phat;
    for (decltype(d.get_order()) i = 2u; i <= d.get_order(); ++i) {
        phat *= tmp;
        retval += binomial(alpha, i) * phat;
    }
    retval *= cbrt_p0;
    return retval;
}

/// Overload for the sine
/**
 * Implements the sine of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0
 * \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\ \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the sine of \p d
 */
template <typename T>
inline gdual<T> sin(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sin_p0 = audi::sin(p0);
    auto cos_p0 = audi::cos(p0);

    double factorial(1.);
    int coeff(1.);
    gdual<T> cos_taylor(1.);
    gdual<T> tmp(cos_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1); // 2!, 4!, 6!, ...
        cos_taylor += (coeff * tmp) / factorial;
    }

    factorial = 1.;
    coeff = 1;
    gdual<T> sin_taylor(phat);
    tmp = sin_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1); // 3!, 5!, 7!, ...
        sin_taylor += (coeff * tmp) / factorial;
    }
    return (sin_p0 * cos_taylor + cos_p0 * sin_taylor);
}

/// Overload for the cosine
/**
 * Implements the cosine of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\cos f)} = \cos f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0
 * \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the cosine of \p d
 */
template <typename T>
inline gdual<T> cos(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sin_p0 = audi::sin(p0);
    auto cos_p0 = audi::cos(p0);

    double factorial(1.);
    int coeff(1.);
    gdual<T> cos_taylor(1.);
    gdual<T> tmp(cos_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1); // 2!, 4!, 6!, ...
        cos_taylor += (coeff * tmp) / factorial;
    }

    factorial = 1.;
    coeff = 1;
    gdual<T> sin_taylor(phat);
    tmp = sin_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1); // 3!, 5!, 7!, ...
        sin_taylor += (coeff * tmp) / factorial;
    }
    return (cos_p0 * cos_taylor - sin_p0 * sin_taylor);
}

/// Computes both sine and cosine
/**
 * As most of the computations for the sine and cosine is the same, it is twice as fast
 * to get both sine and cosine at once rather than computing them in sequence.
 * Use this function when both sine and cosine are needed.
 *
 * @param d audi::gdual argument
 *
 * @return an std::array containing the Taylor expansions of sine and the cosine (first element, second element)
 *
 */
template <typename T>
std::array<gdual<T>, 2> sin_and_cos(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sin_p0 = audi::sin(p0);
    auto cos_p0 = audi::cos(p0);

    double factorial(1.);
    int coeff(1.);
    gdual<T> cos_taylor(1.);
    gdual<T> tmp(cos_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1); // 2!, 4!, 6!, ...
        cos_taylor += (coeff * tmp) / factorial;
    }

    factorial = 1.;
    coeff = 1;
    gdual<T> sin_taylor(phat);
    tmp = sin_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        coeff = -coeff;           // -1, 1, -1, 1, ...
        tmp *= phat2;             // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1); // 3!, 5!, 7!, ...
        sin_taylor += (coeff * tmp) / factorial;
    }
    auto sine = sin_p0 * cos_taylor + cos_p0 * sin_taylor;
    auto cosine = cos_p0 * cos_taylor - sin_p0 * sin_taylor;
    return std::array<gdual<T>, 2>{{std::move(sine), std::move(cosine)}};
}

/// Overload for the tangent
/**
 * Implements the tangent of a audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\tan f)} = \frac{\tan f_0 + \sum_{k=1}^{k \le 2k+1} B_{2k} \frac{(-4)^k(1-4^k)}{2k!}x^{2k - 1}}{1 - \tan f_0
 * \sum_{k=1}^{k \le 2k+1} \frac{B_{2k}(-4)^k(1-4^k)}{2k!}x^{2k - 1} } \f]
 *
 * where \f$T_f = f_0 + \hat f\f$ and \f$ B_{2k}\f$ are the Bernoulli numbers.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the tangent of \p d
 *
 */
template <typename T>
inline gdual<T> tan(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;
    auto tan_p0 = audi::tan(p0);

    // Pre-compute Bernoulli numbers.
    std::vector<double> bn;
    boost::math::bernoulli_b2n<double>(0, (d.get_order() + 1) / 2 + 1,
                                       std::back_inserter(bn)); // Fill vector with even Bernoulli numbers.

    gdual<T> tan_taylor = phat;
    // Factors
    double factorial(24);
    double four_k(16);
    for (auto k = 2u; 2 * k - 1 <= d.get_order(); ++k) {
        phat *= phat2;
        tan_taylor += bn[k] * four_k * (1 - audi::abs(four_k)) / factorial * phat;
        four_k *= -4.;
        factorial *= (2. * k + 1.) * (2. * k + 2.);
    }
    return (tan_p0 + tan_taylor) / (1. - tan_p0 * tan_taylor);
}

/// Overload for the hyperbolic sine
/**
 * Implements the hyperbolic sine of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\sin f)} = \sinh f_0 \left(\sum_{i=0}^{2i\le m} \frac{\hat f^{2i}}{(2i)!}\right) + \cosh f_0
 * \left(\sum_{i=0}^{(2i+1)\le m} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\ \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the hyperbolic sine of \p d
 */
template <typename T>
inline gdual<T> sinh(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sinh_p0 = audi::sinh(p0);
    auto cosh_p0 = audi::cosh(p0);

    double factorial(1.);
    gdual<T> cosh_taylor(1.);
    gdual<T> tmp(cosh_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1.); // 2!, 4!, 6!, ...
        cosh_taylor += tmp / factorial;
    }

    factorial = 1.;
    gdual<T> sinh_taylor(phat);
    tmp = sinh_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1.); // 3!, 5!, 7!, ...
        sinh_taylor += tmp / factorial;
    }
    return (sinh_p0 * cosh_taylor + cosh_p0 * sinh_taylor);
}

/// Overload for the hyperbolic cosine
/**
 * Implements the hyperbolic cosine of a audi::gdual.
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\sin f)} = \cosh f_0 \left(\sum_{i=0}^{2i\le m} \frac{\hat f^{2i}}{(2i)!}\right) + \sinh f_0
 * \left(\sum_{i=0}^{(2i+1)\le m} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\ \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the hyperbolic cosine of \p d
 */
template <typename T>
inline gdual<T> cosh(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sinh_p0 = audi::sinh(p0);
    auto cosh_p0 = audi::cosh(p0);

    double factorial(1.);
    gdual<T> cosh_taylor(1.);
    gdual<T> tmp(cosh_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1.); // 2!, 4!, 6!, ...
        cosh_taylor += tmp / factorial;
    }

    factorial = 1.;
    gdual<T> sinh_taylor(phat);
    tmp = sinh_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1.); // 3!, 5!, 7!, ...
        sinh_taylor += tmp / factorial;
    }
    return (cosh_p0 * cosh_taylor + sinh_p0 * sinh_taylor);
}

/// Computes both the hyperbolic sine and the hyperbolic cosine
/**
 * As most of the computations for the hyperbolic sine and hyperbolic cosine is the same, it is twice as fast
 * to get them both at once rather than computing them in sequence.
 * Use this function when both the hyperbolic sine and the hyperbolic cosine are needed.
 *
 * @param d audi::gdual argument
 *
 * @return an std::array containing the Taylor expansions of hyperbolic sine and cosine (first element, second element)
 *
 */
template <typename T>
std::array<gdual<T>, 2> sinh_and_cosh(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sinh_p0 = audi::sinh(p0);
    auto cosh_p0 = audi::cosh(p0);

    double factorial(1.);
    gdual<T> cosh_taylor(1.);
    gdual<T> tmp(cosh_taylor);
    for (auto i = 2u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^2, phat^4, phat^6 ...
        factorial *= i * (i - 1.); // 2!, 4!, 6!, ...
        cosh_taylor += tmp / factorial;
    }

    factorial = 1.;
    gdual<T> sinh_taylor(phat);
    tmp = sinh_taylor;
    for (auto i = 3u; i <= d.get_order(); i += 2) {
        tmp *= phat2;              // phat^3, phat^5, phat^7 ...
        factorial *= i * (i - 1.); // 3!, 5!, 7!, ...
        sinh_taylor += tmp / factorial;
    }
    auto sineh = sinh_p0 * cosh_taylor + cosh_p0 * sinh_taylor;
    auto cosineh = cosh_p0 * cosh_taylor + sinh_p0 * sinh_taylor;
    return std::array<gdual<T>, 2>{{std::move(sineh), std::move(cosineh)}};
}

/// Overload for the hyperbolic tangent
/**
 * Implements the hyperbolic tangent of a audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\tan f)} = \frac{\tanh f_0 + \sum_{k=1}^{k \le 2k+1} B_{2k} \frac{4^k(4^k-1)}{2k!}x^{2k - 1}}{1 + \tanh f_0
 * \sum_{k=1}^{k \le 2k+1} \frac{B_{2k}4^k(4^k-1)}{2k!}x^{2k - 1} } \f]
 *
 * where \f$T_f = f_0 + \hat f\f$ and \f$ B_{2k}\f$ are the Bernoulli numbers.
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the hyperbolic tangent of \p d
 *
 */
template <typename T>
inline gdual<T> tanh(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;
    auto tanh_p0 = audi::tanh(p0);

    // Pre-compute Bernoulli numbers.
    std::vector<double> bn;
    boost::math::bernoulli_b2n<double>(0, (d.get_order() + 1) / 2 + 1,
                                       std::back_inserter(bn)); // Fill vector with even Bernoulli numbers.

    gdual<T> tanh_taylor = phat;
    // Factors
    double factorial(24.);
    double four_k(16.);
    for (auto k = 2u; 2 * k - 1 <= d.get_order(); ++k) {
        phat *= phat2;
        tanh_taylor += bn[k] * four_k * (four_k - 1.) / factorial * phat;
        four_k *= 4.;
        factorial *= (2 * k + 1.) * (2 * k + 2.);
    }
    return (tanh_p0 + tanh_taylor) / (1. + tanh_p0 * tanh_taylor);
}

/// Overload for the inverse hyperbolic tangent
/**
 * Implements the inverse hyperbolic tangent of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{atanh} f)} =  \mbox{atanh} f_0 +\frac 12 \sum_{k=1}^m \left(\frac{1}{(1-f_0)^k} +
 * \frac{(-1)^{k+1}}{(1+f_0)^k}\right) \frac {\hat f^k}{k} \f]
 *
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse hyperbolic tangent of \p d
 *
 */
template <typename T>
inline gdual<T> atanh(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto powphat(phat);
    auto atanh_p0 = audi::atanh(p0);

    gdual<T> retval(0.);
    double coeff = 1.;

    for (auto k = 1u; k <= d.get_order(); ++k) {
        auto add = (1. / audi::pow(1. - p0, k) + coeff / audi::pow(1. + p0, k)) / static_cast<T>(k);
        retval += add * powphat;
        coeff *= -1;
        powphat *= phat;
    }
    return atanh_p0 + 0.5 * retval;
}

/// Overload for the inverse tangent
/**
 * Implements the inverse tangent of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{atan} f)} =  \mbox{atan} f_0 + \sum_{k=1}^{2k-1\le m} \left(\frac{1 + \sum_{j=1}^{2j\le 2k-1} {2k-1 \choose
 * 2j} f_0^{2j}(-1)^j}{(1+f_0^2)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \sum_{k=1}^{2k\le m}
 * \left(\frac{\sum_{j=1}^{2j-1\le 2k} {2k \choose 2j-1} f_0^{2j-1}(-1)^{j+1}}{(1+f_0^2)^{2k}}\right) \frac {\hat
 * f^{2k}}{2k}(-1)^k \f]
 *
 * This formula derives directly from the formula for audi::atanh noting that: \f$ \mbox{atan}(z) = i
 * \mbox{atanh}(-iz)\f$
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse tangent of \p d
 *
 */
template <typename T>
inline gdual<T> atan(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0) / (1. + p0 * p0);
    auto powphat(phat);

    auto retval = gdual<T>(audi::atan(p0));
    T coeff1(1.);
    T coeff2(-1.);

    for (auto k = 1u; k <= d.get_order(); ++k) {
        if (k % 2u) { // This is for odd powers 1..3..5
            T binom(1.);
            auto f0 = p0 * p0;
            T cf_i(-1.);
            for (auto j = 1u; 2 * j <= k; ++j) {
                binom += (binomial(k, 2u * j) * cf_i) * f0;
                f0 *= p0 * p0;
                cf_i = -cf_i;
            }
            retval += binom * powphat * coeff1 / k;
            coeff1 = -coeff1;
        } else { // This is for even powers 2..4..6
            T binom(0.);
            auto f0 = p0;
            T cf_i(1.);
            for (auto j = 1u; 2 * j - 1 <= k; ++j) {
                binom += T(binomial(k, 2u * j - 1u) * cf_i) * f0;
                f0 *= p0 * p0;
                cf_i = -cf_i;
            }
            retval += binom * powphat * coeff2 / k;
            coeff2 = -coeff2;
        }
        powphat *= phat;
    }
    return retval;
}

/// Overload for the inverse hyperbolic sine
/**
 * Implements the inverse inverse hyperbolic sine of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{asinh} f)} = T_{\left(\log\left(f + \sqrt{1 + f^2}\right)\right)}
 * \f]
 *
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse hyperbolic sine of \p d
 *
 */
template <typename T>
inline gdual<T> asinh(const gdual<T> &d)
{
    return log(d + sqrt(1. + d * d));
}

/// Overload for the inverse hyperbolic cosine
/**
 * Implements the inverse inverse hyperbolic cosine of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{acosh} f)} = T_{\left(\log\left(f + \sqrt{f^2 - 1}\right)\right)}
 * \f]
 *
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse hyperbolic cosine of \p d
 *
 */
template <typename T>
inline gdual<T> acosh(const gdual<T> &d)
{
    return log(d + sqrt(d * d - 1.));
}

/// Overload for the inverse sine
/**
 * Implements the inverse inverse sine of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{asin} f)} = T_{\left(\mbox{atan} \left(f / \sqrt{1 - f^2}\right)\right)}
 * \f]
 *
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse sine of \p d
 *
 */
template <typename T>
inline gdual<T> asin(const gdual<T> &d)
{
    return atan(d / sqrt(1. - d * d));
}

/// Overload for the inverse cosine
/**
 * Implements the inverse inverse cosine of an audi::gdual.
 * Essentially, it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\mbox{acos} f)} = T_{\left(\mbox{atan} \left(\sqrt{1 - f^2} / f\right)\right)}
 * \f]
 *
 *
 * @param d audi::gdual argument
 *
 * @return an audi::gdual containing the Taylor expansion of the inverse cosine of \p d
 *
 */
template <typename T>
inline gdual<T> acos(const gdual<T> &d)
{
    return 0.5 * boost::math::constants::pi<double>() - asin(d);
}

/// Overload for the absolute value
/**
 * Implements the absolute value of a audi::gdual.
 * Essentially, it inverts the sign of \f$T_f\f$ as follows:
 *
 * \f[
 * T_{(\mbox{abs} f)} = \left\{ \begin{array}{ll} T_f & f_0 \ge 0 \\ -T_f & f_0 < 0 \end{array} \right.
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * \note If \f$f_0\f$ is zero, the right Taylor expansion will be returned rather than nans.
 *
 * \note This operation is not availiable for std::complex types.
 *
 * @param d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the absoute value of \p d
 *
 */
template <typename T>
inline gdual<T> abs(const gdual<T> &d)
{
    auto p0 = d.constant_cf();
    if (p0 >= 0) {
        return d;
    }
    return -d;
}
// Vectorized abs specialization
template <typename T>
inline gdual<vectorized<T>> abs(const gdual<vectorized<T>> &d)
{
    gdual<vectorized<T>> retval(d);
    auto p0 = retval.constant_cf();
    for (auto it = retval._container().begin(); it != retval._container().end(); ++it) {
        // if the coefficient has one only element copy it over to the whole length of p0
        if (it->m_cf.size() == 1u) {
            it->m_cf.resize(p0.size(), it->m_cf[0]);
        }
        for (auto i = 0u; i < p0.size(); ++i) {
            if (p0[i] < 0) {
                it->m_cf.set_value(i, -it->m_cf[i]);
            }
        }
    }
    return retval;
}

#if defined(AUDI_WITH_MPPP)
// template specialization for the mppp::real128 type.  boost::math::bernoulli_b2n<mppp::real128> can also work, but
// requires to specialize std::numeric_limits for vectorized and mppp::real128
template <>
inline gdual<mppp::real128> tan(const gdual<mppp::real128> &d)
{
    auto sincos = sin_and_cos(d);
    return sincos[0] / sincos[1];
}

// template specialization for the mppp::real128 type. boost::math::constants::pi<mppp::real128> can also work, but
// requires to specialize std::numeric_limits for vectorized and mppp::real128
template <>
inline gdual<mppp::real128> acos(const gdual<mppp::real128> &d)
{
    return 0.5 * mppp::real128_pi() - asin(d);
}

// template specialization for the mppp::real128 type.  boost::math::bernoulli_b2n<mppp::real128> can also work, but
// requires to specialize std::numeric_limits for vectorized and mppp::real128
template <>
inline gdual<mppp::real128> tanh(const gdual<mppp::real128> &d)
{
    auto sinhcosh = sinh_and_cosh(d);
    return sinhcosh[0] / sinhcosh[1];
}
#endif

} // end of namespace audi

#endif
