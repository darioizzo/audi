// Copyright © 2018–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com),
// Sean Cowan (lambertarc@icloud.com)
//
// This file is part of the audi library.
//
// The audi library is free software: you can redistribute it and/or modify
// it under the terms of either:
//   - the GNU General Public License as published by the Free Software
//     Foundation, either version 3 of the License, or (at your option)
//     any later version, or
//   - the GNU Lesser General Public License as published by the Free
//     Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
// The audi library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License and the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser General Public License along with the audi library.
// If not, see <https://www.gnu.org/licenses/>.

#ifndef AUDI_FUNCTIONS_FROM_D_HPP
#define AUDI_FUNCTIONS_FROM_D_HPP

#include <boost/math/constants/constants.hpp>
#include <cmath>

#include <audi/back_compatibility.hpp>
#include <audi/config.hpp>
#include <audi/detail/overloads.hpp>
#include <audi/gdual.hpp>

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
 * @param f Taylor expansion of the inner function
 * @param dg Taylor expansion of the derivative of the outer function
 * @param g0 Value of the outer function at the expansion point
 */
template <typename T, typename M>
inline gdual<T, M> _compose_from_derivative(gdual<T, M> f, gdual<T, M> dg, T g0)
{
    auto ss = f.get_symbol_set();
    if (ss.size() == 0) {
        return gdual<T, M>(g0);
    }
    auto retval = (dg * f.partial(ss[0])).integrate(ss[0]);
    for (auto i = 1u; i < ss.size(); ++i) {
        f = f.subs("d" + ss[i - 1], 0);
        dg = dg.subs("d" + ss[i - 1], 0);
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
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse hyperbolic tangent of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> atanh_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::atanh(f0);
    auto dg = 1. / (1. - f * f);
    return _compose_from_derivative(f, dg, g0);
}

/// Possible overload for the inverse tangent
/**
 * Implements the inverse tangent of an audi::gdual. The audi:atan
 * overload is actually faster and this is only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse tangent of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> atan_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::atan(f0);
    auto dg = 1. / (1. + f * f);
    return _compose_from_derivative(f, dg, g0);
}

/// Possible overload for the inverse sine
/**
 * Implements the inverse sine of an audi::gdual. The audi:asin
 * overload may be slower, but this is currently only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse sine of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> asin_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::asin(f0);
    auto dg = 1. / audi::sqrt(1. - f * f);
    return _compose_from_derivative(f, dg, g0);
}

/// Possible overload for the inverse hyperbolic sine
/**
 * Implements the inverse hyperbolic sine of an audi::gdual. The audi:asinh
 * overload may be slower, but this is currently only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse hyperbolic sine of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> asinh_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::asinh(f0);
    auto dg = 1. / audi::sqrt(1. + f * f);
    return _compose_from_derivative(f, dg, g0);
}

/// Possible overload for the inverse cosine
/**
 * Implements the inverse cosine of an audi::gdual. The audi::acos
 * overload may be slower, but this is currently only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse cosine of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> acos_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::acos(f0);
    auto dg = -1. / audi::sqrt(1. - f * f);
    return _compose_from_derivative(f, dg, g0);
}

/// Possible overload for the inverse hyperbolic cosine
/**
 * Implements the inverse hyperbolic cosine of an audi::gdual. The audi:cosh
 * is faster, and this is currently only provided as to benchmark the
 * performances of audi::_compose_from_derivative
 *
 * @param f audi::gdual argument
 *
 * @return an audi:gdual containing the Taylor expansion of the inverse hyperbolic cosine of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> acosh_d(const gdual<T, M> &f)
{
    auto f0 = f.constant_cf();
    auto g0 = audi::acosh(f0);
    auto dg = 1. / audi::sqrt((f - 1.) * (f + 1.));
    return _compose_from_derivative(f, dg, g0);
}

/// Overload for the error function
/**
 * Implements the error function of an audi::gdual. Essentially it
 * makes use of the fact that \f$ \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}}\exp(-x^2)\f$
 * to then use audi::_compose_from_derivative
 *
 * @param d audi::gdual type
 *
 * @return an audi:gdual containing the Taylor expansion of the error function of \p d
 *
 */
template <typename T, typename M>
inline gdual<T, M> erf(const gdual<T, M> &d)
{
    auto f0 = d.constant_cf();
    auto g0 = audi::erf(f0);
    auto dg = (2. / std::sqrt(boost::math::constants::pi<double>())) * exp(-d * d);
    return _compose_from_derivative(d, dg, g0);
}

#if defined(AUDI_WITH_QUADMATH)
template <typename M>
inline gdual<mppp::real128, M> erf(const gdual<mppp::real128, M> &d)
{
    auto f0 = d.constant_cf();
    auto g0 = audi::erf(f0);
    auto dg = (2. / audi::sqrt(mppp::real128_pi())) * exp(-d * d);
    return _compose_from_derivative(d, dg, g0);
}
#endif

} // end of namespace audi

#endif
