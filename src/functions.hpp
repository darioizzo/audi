#ifndef AUDI_FUNCTIONS_HPP
#define AUDI_FUNCTIONS_HPP

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <piranha/binomial.hpp>
#include <stdexcept>

#include "gdual.hpp"



namespace audi
{

/// Overload for the exponential
/**
 * Implements the exponential of a audi::gdual. 
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ... \right)
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the exponential of \p d
*/
inline gdual exp(const gdual &d)
{
    gdual retval(1., d.get_order());
    double fact=1;
    auto p0 = d.constant_cf();
    auto phat = d - p0;
    gdual tmp(phat);

    retval+=phat;
    for (int i = 2; i <= d.get_order(); ++i) {
        phat*=tmp;
        fact*=i;
        retval+=phat / fact;
    }
    return retval * std::exp(p0);
}

/// Overload for the logarithm
/**
 * Implements the logarithm of a audi::gdual. 
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ...
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the logarithm of \p d
 *
 * @throw std::domain_error if std::log(\f$f_0\f$) is not finite (uses std::isfinite)
 */
inline gdual log(const gdual &d)
{
    gdual retval(0., d.get_order());
    double fatt = 1;
    auto p0 = d.constant_cf();
    auto log_p0 = std::log(p0);
    if (!std::isfinite(log_p0)) {
        throw std::domain_error("std::log(" + boost::lexical_cast<std::string>(p0) + ") returns " + boost::lexical_cast<std::string>(log_p0) + " in log(const gdual& d)");
    }
    auto phat = (d - p0);
    phat = phat/p0;
    gdual tmp(phat);

    retval = log_p0 + phat;
    for (auto i = 2; i <= d.get_order(); ++i) {
        fatt *= -1;
        phat*=tmp;
        retval =  retval + fatt * phat / i;
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
 * @param[in] base the base for the exponent
 * @param[in] d audi::gdual argument
 *
 * @throw std::domain_error if std::log(\p base) is not finite (uses std::isfinite)
*/
inline gdual pow(double base, const gdual &d)
{
    double int_part;
    // checks wether the exponent is an integer, in which case
    // it calls a different overload
    if (d.degree() == 0) {
        auto p0 = d.constant_cf();
        double float_part = std::modf(p0, &int_part);
        if (float_part == 0.) {
            return gdual(std::pow(base, p0), d.get_order()); //nan is possible here
        }
    }
    return exp(std::log(base) * d);
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
 * @param[in] d audi::gdual argument
 * @param[in] alpha exponent
 *
 * @return an audi:gdual containing the Taylor expansion of \p d elevated to the power \p alpha
 *
 * @throw std::domain_error if std::pow(\f$f_0, \alpha \f$) is not finite (uses std::isfinite)
 * or if f_0 is 0.
 */
inline gdual pow(const gdual &d, double alpha)
{
    gdual retval(1., d.get_order());
    auto p0 = d.constant_cf();
    double pow_p0 = std::pow(p0,alpha);
    if (!std::isfinite(pow_p0)) {
        throw std::domain_error("std::pow(" + boost::lexical_cast<std::string>(p0)+ ", " + boost::lexical_cast<std::string>(alpha) + ") returned a non finite number: " + boost::lexical_cast<std::string>(pow_p0));
    }
    if (p0 == 0) {
        throw std::domain_error("exponentiation not defined for polynomial having zero constant term");
    }
    auto phat = d - p0;
    phat = phat/p0;
    gdual tmp(phat);

    retval+=alpha * phat;
    for (auto i = 2; i <= d.get_order(); ++i) {
        phat*=tmp;
        retval+=piranha::math::binomial(alpha,i) * phat;
    }
    retval*=pow_p0;
    return retval;
}

// Its important this comes after the pow(gdual, double) overload
/// Overload for the integer exponentiation
/**
 * Implements the integer exponentiation of a audi::gdual. Essentially,
 * it uses the \f$\mathcal P_{n,m}\f$ multiplication on \p d \p n times
 *
 * @param[in] d audi::gdual argument
 * @param[in] n integer exponent
*/
inline gdual pow(const gdual& d, int n)
{
    if (n <= 0) { 
        return pow(d, (double)n);
    }
    gdual retval(d);
    for (auto i = 1; i < n; ++i) {
        retval*=d;
    }
    return retval;
}

/// Overload for the exponentiation of a gdual to a gdual power
/**
 * Computes the exponentiation of an audi::gdual to the power of an audi::gdual. 
 * If the exponent is an integer constant, it calls a different overload. Otherwise
 * it converts \f$T_f^{T_g}\f$ to \f$\exp(T_g*\log(T_f))\f$ and computes this
 * last expression instead.
 *
 * @param[in] d1 audi::gdual argument
 * @param[in] d2 audi::gdual argument
 *
 * @throw std::domain_error if std::log(\f$f_0\f$) is not finite (uses std::isfinite)
*/
inline gdual pow(const gdual &d1, const gdual &d2)
{
    double int_part;
    // checks wether the exponent is an integer, in which case it calls 
    // a different overload
    if (d2.degree() == 0) {
        auto p0 = d2.constant_cf();
        double float_part = std::modf(p0, &int_part);
        if (float_part == 0.) {
            return pow(d1, p0);
        }
    }
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
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the square root of \p d
 *
 * @throw std::domain_error if std::pow(\f$f_0, 0.5\f$) is not finite (uses std::isfinite)
 * or if f_0 is 0.
 */
inline gdual sqrt(const gdual &d)
{
    return pow(d, 0.5); // TODO: subsitute this by similar code to cbrt?
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
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the square root of \p d
 *
 * @throw std::domain_error if std::cbrt(\f$f_0\f$) is not finite (uses std::isfinite)
 */
inline gdual cbrt(const gdual& d)
{
    double alpha = 1/3.;
    gdual retval(1., d.get_order());
    auto p0 = d.constant_cf();
    double cbrt_p0 = std::cbrt(p0);
    if (!std::isfinite(cbrt_p0)) {
        throw std::domain_error("std::cbrt_p0(" + boost::lexical_cast<std::string>(p0)+ ", " + boost::lexical_cast<std::string>(alpha) + ") returned a non finite number: " + boost::lexical_cast<std::string>(cbrt_p0));
    }
    auto phat = d - p0;
    phat = phat/p0;
    gdual tmp(phat);

    retval+=alpha * phat;
    for (auto i = 2; i <= d.get_order(); ++i) {
        phat*=tmp;
        retval+=piranha::math::binomial(alpha,i) * phat;
    }
    retval*=cbrt_p0;
    return retval;
}

/// Overload for the sine
/**
 * Implements the sine of a audi::gdual. 
 * Essentially it performs the following computations in the \f$\mathcal P_{n,m}\f$
 * algebra:
 *
 * \f[
 * T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the sine of \p d
*/
inline gdual sin(const gdual& d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    double sin_p0 = std::sin(p0);
    double cos_p0 = std::cos(p0);

    double factorial=1.;
    double coeff=1.;
    gdual cos_taylor(1., d.get_order());
    gdual tmp(cos_taylor);
    for (auto i=2; i<=d.get_order(); i+=2) {
        coeff*=-1.;                             // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^2, phat^4, phat^6 ...
        factorial*=i * (i-1);                   // 2!, 4!, 6!, ...
        cos_taylor += (coeff/factorial) * tmp;
    }

    factorial=1.;
    coeff=1.;
    gdual sin_taylor(phat);
    tmp = sin_taylor;
    for (auto i=3; i<=d.get_order(); i+=2) {
        coeff*=-1.;                             // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^3, phat^5, phat^7 ...
        factorial*=i * (i-1);                   // 3!, 5!, 7!, ...
        sin_taylor += (coeff/factorial) * tmp;
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
 * T_{(\cos f)} = \cos f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
 * \f]
 *
 * where \f$T_f = f_0 + \hat f\f$.
 *
 * @param[in] d audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the cosine of \p d
*/
inline gdual cos(const gdual& d)
{
    auto p0 = d.constant_cf();
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    double sin_p0 = std::sin(p0);
    double cos_p0 = std::cos(p0);

    double factorial=1.;
    double coeff=1.;
    gdual cos_taylor(1., d.get_order());
    gdual tmp(cos_taylor);
    for (auto i=2; i<=d.get_order(); i+=2) {
        coeff*=-1.;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                              // phat^2, phat^4, phat^6 ...
        factorial*=i * (i-1);                    // 2!, 4!, 6!, ...
        cos_taylor += (coeff/factorial) * tmp;
    }

    factorial=1.;
    coeff=1.;
    gdual sin_taylor(phat);
    tmp = sin_taylor;
    for (auto i=3; i<=d.get_order(); i+=2) {
        coeff*=-1.;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                              // phat^3, phat^5, phat^7 ...
        factorial*=i * (i-1);                    // 3!, 5!, 7!, ...
        sin_taylor += (coeff/factorial) * tmp;
    }
    return (cos_p0 * cos_taylor - sin_p0 * sin_taylor);
}

} // end of namespace audi 

#endif
