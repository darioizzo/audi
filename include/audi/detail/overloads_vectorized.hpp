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
