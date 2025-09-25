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

#ifndef AUDI_OVERLOADS_HPP
#define AUDI_OVERLOADS_HPP

#include <cmath>

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)
#include <audi/real128.hpp>
#endif

#include <audi/back_compatibility.hpp>
#include <boost/type_traits/is_complex.hpp>

// This is repeated here instead of including type_traits as to avoid a circular dependency
template <class T>
struct is_arithmetic_or_complex
    : std::integral_constant<bool, std::is_arithmetic<T>::value || boost::is_complex<T>::value> {
};

// This macro writes the overload for std::fun_name af an arithmetic or complex type. It simply calls std::fun_name
// It is used to allow calls such as audi::cos(T) [T = double] in templated functions.
#define ARITH_OR_COMPLEX_OVERLOAD(fun_name)                                                                            \
    template <typename T, enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>                                    \
    inline auto fun_name(T in)->decltype(std::fun_name(T(0.)))                                                         \
    {                                                                                                                  \
        return std::fun_name(in);                                                                                      \
    }

#if defined(AUDI_WITH_QUADMATH)
// This macro writes the overload for std::fun_name af a mppp::real128. It simply calls mppp::fun_name
// It is used to allow calls such as audi::cos(T) [T = double] in templated functions.
#define REAL128_OVERLOAD(fun_name)                                                                                     \
    inline mppp::real128 fun_name(mppp::real128 in)                                                                    \
    {                                                                                                                  \
        return mppp::fun_name(in);                                                                                     \
    }
#endif

namespace audi
{
#if defined(AUDI_WITH_QUADMATH)

REAL128_OVERLOAD(exp)
REAL128_OVERLOAD(erf)
REAL128_OVERLOAD(lgamma)
REAL128_OVERLOAD(log)
REAL128_OVERLOAD(sin)
REAL128_OVERLOAD(cos)
REAL128_OVERLOAD(sinh)
REAL128_OVERLOAD(tan)
REAL128_OVERLOAD(cosh)
REAL128_OVERLOAD(tanh)
REAL128_OVERLOAD(asin)
REAL128_OVERLOAD(acos)
REAL128_OVERLOAD(atan)
REAL128_OVERLOAD(asinh)
REAL128_OVERLOAD(acosh)
REAL128_OVERLOAD(atanh)
REAL128_OVERLOAD(sqrt)
REAL128_OVERLOAD(abs)
REAL128_OVERLOAD(cbrt)
template <typename T, typename U,
          enable_if_t<std::is_same<U, mppp::real128>::value || std::is_same<T, mppp::real128>::value, int> = 0>
inline mppp::real128 pow(const T &base, const U &d)
{
    return mppp::pow(base, d);
}
#endif

ARITH_OR_COMPLEX_OVERLOAD(exp)
ARITH_OR_COMPLEX_OVERLOAD(erf)
ARITH_OR_COMPLEX_OVERLOAD(lgamma)
ARITH_OR_COMPLEX_OVERLOAD(log)
ARITH_OR_COMPLEX_OVERLOAD(sin)
ARITH_OR_COMPLEX_OVERLOAD(cos)
ARITH_OR_COMPLEX_OVERLOAD(tan)
ARITH_OR_COMPLEX_OVERLOAD(sinh)
ARITH_OR_COMPLEX_OVERLOAD(cosh)
ARITH_OR_COMPLEX_OVERLOAD(tanh)
ARITH_OR_COMPLEX_OVERLOAD(asin)
ARITH_OR_COMPLEX_OVERLOAD(acos)
ARITH_OR_COMPLEX_OVERLOAD(atan)
ARITH_OR_COMPLEX_OVERLOAD(asinh)
ARITH_OR_COMPLEX_OVERLOAD(acosh)
ARITH_OR_COMPLEX_OVERLOAD(atanh)
ARITH_OR_COMPLEX_OVERLOAD(sqrt)
ARITH_OR_COMPLEX_OVERLOAD(abs)

template <typename T, enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T cbrt(T in)
{
    return std::cbrt(in);
}
template <typename T, enable_if_t<boost::is_complex<T>::value, int> = 0>
inline T cbrt(T in)
{
    return std::pow(in, 1. / 3.); // needs a separate template as cbrt does not exist for complex types
}

// NON UNARY FUNCTIONS -------------------------------------------------------------------------------//
template <typename T, typename U, enable_if_t<std::is_arithmetic<U>::value && std::is_arithmetic<T>::value, int> = 0>
inline double pow(const U &base, const T &d)
{
    return static_cast<double>(std::pow(base, d));
}

template <typename T, typename U, enable_if_t<boost::is_complex<U>::value && std::is_arithmetic<T>::value, int> = 0>
inline U pow(const U &base, const T &d)
{
    return std::pow(base, d);
}

template <typename T, typename U, enable_if_t<boost::is_complex<T>::value && std::is_arithmetic<U>::value, int> = 0>
inline T pow(const U &base, const T &d)
{
    return std::pow(base, d);
}

} // namespace audi
#endif
