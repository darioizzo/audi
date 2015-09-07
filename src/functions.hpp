#ifndef AUDI_FUNCTIONS_HPP
#define AUDI_FUNCTIONS_HPP

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <piranha/binomial.hpp>
#include <stdexcept>

#include "gdual.hpp"

namespace audi
{


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

inline gdual pow(const gdual& d, int n)
{
    gdual retval(d);
    for (auto i = 1; i < n; ++i) {
        retval*=d;
    }
    return retval;
}


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

inline gdual pow(const gdual &d1, const gdual &d2)
{
    double int_part;
    // checks wether the exponent is an integer, in which case it calls a different overload
    if (d2.degree() == 0) {
        auto p0 = d2.constant_cf();
        double float_part = std::modf(p0, &int_part);
        if (float_part == 0.) {
            return pow(d1, p0);
        }
    }
    return exp(d2 * log(d1));
}

inline gdual sqrt(const gdual &d)
{
    return pow(d, 0.5);
}

inline  gdual cbrt(const gdual& d)
{
    double alpha = 0.33333333333333333333333; // bello eh!!!
    gdual retval(1., d.get_order());
    auto p0 = d.find_cf(std::vector<int>(d.get_n_variables(),0));
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

inline gdual sin(const gdual& d)
{
    auto p0 = d.find_cf(std::vector<int>(d.get_n_variables(),0));
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sin_p0 = std::sin(p0);
    auto cos_p0 = std::cos(p0);

    double factorial=1;
    double coeff=1;
    gdual cos_taylor(1., d.get_order());
    gdual tmp(cos_taylor);
    for (auto i=2; i<=d.get_order(); i+=2)
    {
        coeff*=-1;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^2, phat^4, phat^6 ...
        factorial*=i * (i-1);                   // 2!, 4!, 6!, ...
std::cout << "i: " << i << std::endl; 
std::cout << "coeff: " << coeff << std::endl; 
std::cout << "factorial: " << factorial << std::endl; 
        cos_taylor += coeff * tmp / factorial;
    }

    factorial=1;
    coeff=1;
    gdual sin_taylor(phat);
    tmp = gdual(1, d.get_order());
    for (auto i=3; i<=d.get_order(); i+=2)
    {
        coeff*=-1;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^3, phat^5, phat^7 ...
        factorial*=i * (i-1);                   // 3!, 5!, 7!, ...
std::cout << "i: " << i << std::endl; 
std::cout << "coeff: " << coeff << std::endl; 
std::cout << "factorial: " << factorial << std::endl; 
        sin_taylor += coeff * tmp / factorial;
    }
    return (sin_p0 * cos_taylor + cos_p0 * sin_taylor);
}

inline gdual cos(const gdual& d)
{
    auto p0 = d.find_cf(std::vector<int>(d.get_n_variables(),0));
    auto phat = (d - p0);
    auto phat2 = phat * phat;

    auto sin_p0 = std::sin(p0);
    auto cos_p0 = std::cos(p0);

    double factorial=1;
    double coeff=1;
    gdual cos_taylor(1., d.get_order());
    gdual tmp(cos_taylor);
    for (auto i=2; i<d.get_order(); i+=2)
    {
        coeff*=-1;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^2, phat^4, phat^6 ...
        factorial*=i * (i-1);                   // 2!, 4!, 6!, ...
        cos_taylor += coeff * tmp / factorial;
    }

    factorial=2;
    coeff=1;
    gdual sin_taylor(phat);
    tmp = gdual(1, d.get_order());
    for (auto i=3; i<d.get_order(); i+=2)
    {
        coeff*=-1;                              // -1, 1, -1, 1, ...
        tmp*=phat2;                             // phat^3, phat^5, phat^7 ...
        factorial*=i * (i-1);                   // 3!, 5!, 7!, ...
        sin_taylor += coeff * tmp / factorial;
    }
    return (cos_p0 * cos_taylor - sin_p0 * sin_taylor);
}

} // end of namespace audi 

#endif
