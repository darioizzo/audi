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

#ifndef AUDI_OVERLOADS_VECTORIZED_HPP
#define AUDI_OVERLOADS_VECTORIZED_HPP

#include <cmath>

#include <audi/vectorized.hpp>
#include <audi/detail/overloads.hpp> // to find audi::fun_name(el)


// This macro writes the overload for std::fun_name performing elementwise evaluations on a vectorized<T>
#define VECTORIZED_OVERLOAD(fun_name)                                                                                  \
    template<typename T>                                                                                               \
    inline vectorized<T> fun_name(vectorized<T> in)                                                                      \
    {                                                                                                                  \
        for (auto &el : in) {                                                                                          \
            el = audi::fun_name(el);                                                                                    \
        }                                                                                                              \
        return in;                                                                                                     \
    }

namespace audi
{
VECTORIZED_OVERLOAD(exp)
VECTORIZED_OVERLOAD(erf)
VECTORIZED_OVERLOAD(lgamma)
VECTORIZED_OVERLOAD(log)
VECTORIZED_OVERLOAD(sin)
VECTORIZED_OVERLOAD(cos)
VECTORIZED_OVERLOAD(tan)
VECTORIZED_OVERLOAD(sinh)
VECTORIZED_OVERLOAD(cosh)
VECTORIZED_OVERLOAD(tanh)
VECTORIZED_OVERLOAD(asin)
VECTORIZED_OVERLOAD(acos)
VECTORIZED_OVERLOAD(atan)
VECTORIZED_OVERLOAD(asinh)
VECTORIZED_OVERLOAD(acosh)
VECTORIZED_OVERLOAD(atanh)
VECTORIZED_OVERLOAD(sqrt)
VECTORIZED_OVERLOAD(abs)
VECTORIZED_OVERLOAD(cbrt)

template<typename T>
inline vectorized<T> pow(double base, vectorized<T> in)
{
    for (auto &el : in) {
        el = std::pow(base, el);
    }
    return in;
}

template<typename T>
inline vectorized<T> pow(vectorized<T> in, double exponent)
{
    for (auto &el : in) {
        el = std::pow(el, exponent);
    }
    return in;
}

} // namespace audi
#endif
