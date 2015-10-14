#ifndef AUDI_FUNCTIONS_FROM_D_HPP
#define AUDI_FUNCTIONS_FROM_D_HPP

#include <cmath>
#include "gdual.hpp"

namespace audi
{

/// Computes the gdual of a function composition when the first derivative is known 
/** 
 * Finds the Taylor expansion for \f$ g(f(\mathbf x))\f$ when the 
 * Taylor expansion of \f$ f(\mathbf x)\f$ is known as well as
 * the derivative \f$ \frac{dg}{df} \f$ and the value \f$ g_0 = g(f(\mathbf x_0))\f$
 *
 * Essentially it exploits the identity:
 *
 * \f[
 * g(f(\mathbf x)) = \int \frac{dg}{df} \frac{\partial f}{\partial x_1}dx_1 + h(x_2..x_n)    
 * \f]
 * 
 * and computes \f$ h(x_2..x_n)\f$ applying the same formula. Note that \f$ \mathbf x = [x_1, x_2, ..., x_n]\f$
 *
 * \note This way of computing derivatives is slower with respect to the corresponding
 * methods based on the nilpotency of \f$ \hat f\f$ and should thus be used only when 
 * necessary (that is when no formula is found to exploit nilpotency.)
 *
 * @param[in] f Taylor expansion of the inner function
 * @param[in] dg Taylor expansion of the derivative of the outer function
 * @param[in] g0 Value of the outer function at the expansion point
*/
inline gdual _compose_from_derivative(gdual f, gdual dg, double g0)
{
    auto ss = f.get_symbol_set();
    if (ss.size() == 0) {
        return gdual(g0);
    }
    auto retval = (dg * f.partial(ss[0])).integrate(ss[0]);
    for (auto i = 1u; i < ss.size(); ++i) {
        f = f.subs(ss[i-1], 0);
        dg = dg.subs(ss[i-1], 0);
        retval += (dg * f.partial(ss[i])).integrate(ss[i]);
    }
    return g0 + retval;
}

/// Possible overload for the inverse hyperbolic tangent
/**
 * Implements the inverse hyperbolic tangent of an audi::gdual. The audi:atanh
 * overload is actually faster and this is only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param[in] f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse hyperbolic tangent of \p d
 *
*/
inline gdual atanh_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::atanh(f0);
    auto dg = 1. / (1 - f*f);
    return _compose_from_derivative(f, dg, g0);
}

inline gdual atan_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::atan(f0);
    auto dg = 1. / (1 + f*f);
    return _compose_from_derivative(f, dg, g0);
}

inline gdual asin_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::asin(f0);
    auto dg = 1. / sqrt(1 - f*f);
    return _compose_from_derivative(f, dg, g0);
}

inline gdual asinh_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::asin(f0);
    auto dg = 1. / sqrt(1 + f*f);
    return _compose_from_derivative(f, dg, g0);
}

inline gdual acos_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::asin(f0);
    auto dg = - 1. / sqrt(1 - f*f);
    return _compose_from_derivative(f, dg, g0);
}

inline gdual acosh_d(const gdual& f)
{
    auto f0 = f.constant_cf();
    double g0 = std::asin(f0);
    auto dg = 1. / sqrt(f*f-1);
    return _compose_from_derivative(f, dg, g0);
}


} // end of namespace audi 

#endif


