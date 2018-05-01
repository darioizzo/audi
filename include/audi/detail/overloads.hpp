#ifndef AUDI_OVERLOADS_HPP
#define AUDI_OVERLOADS_HPP

#include <cmath>

#include <audi/back_compatibility.hpp>
#include <audi/vectorized_double.hpp>
#include <audi/real128.hpp>
#include <boost/type_traits/is_complex.hpp>

// This is repeated here instead of including type_traits as to avoid a circular dependency
template <class T>
struct is_arithmetic_or_complex
    : std::integral_constant<bool, std::is_arithmetic<T>::value || boost::is_complex<T>::value> {
};

// This macro writes the overload for std::fun_name performing elementwise evaluations on a vectorized_double
#define VECTORIZED_OVERLOAD(fun_name)                                                                                  \
    inline vectorized_double fun_name(vectorized_double in)                                                            \
    {                                                                                                                  \
        for (auto &el : in) {                                                                                          \
            el = std::fun_name(el);                                                                                    \
        }                                                                                                              \
        return in;                                                                                                     \
    }

// This macro writes the overload for std::fun_name af an arithmetic or complex type. It simply calls std::fun_name
// It is used to allow calls such as audi::cos(T) [T = double] in templated functions.
#define ARITH_OR_COMPLEX_OVERLOAD(fun_name)                                                                            \
    template <typename T, enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>                                    \
    inline auto fun_name(T in) -> decltype(std::fun_name(T(0.)))                                                       \
    {                                                                                                                  \
        return std::fun_name(in);                                                                                      \
    }

// This macro writes the overload for std::fun_name af a mppp::real128. It simply calls mppp::fun_name
// It is used to allow calls such as audi::cos(T) [T = double] in templated functions.
#define REAL128_OVERLOAD(fun_name)                                                                                     \
    inline mppp::real128 fun_name(mppp::real128 in)                                                                    \
    {                                                                                                                  \
        return mppp::fun_name(in);                                                                                     \
    }

namespace audi
{
VECTORIZED_OVERLOAD(exp)
ARITH_OR_COMPLEX_OVERLOAD(exp)
REAL128_OVERLOAD(exp)

VECTORIZED_OVERLOAD(erf)
ARITH_OR_COMPLEX_OVERLOAD(erf)
REAL128_OVERLOAD(erf)

VECTORIZED_OVERLOAD(lgamma)
ARITH_OR_COMPLEX_OVERLOAD(lgamma)
REAL128_OVERLOAD(lgamma)

VECTORIZED_OVERLOAD(log)
ARITH_OR_COMPLEX_OVERLOAD(log)
REAL128_OVERLOAD(log)

VECTORIZED_OVERLOAD(sin)
ARITH_OR_COMPLEX_OVERLOAD(sin)
REAL128_OVERLOAD(sin)

VECTORIZED_OVERLOAD(cos)
ARITH_OR_COMPLEX_OVERLOAD(cos)
REAL128_OVERLOAD(cos)

VECTORIZED_OVERLOAD(tan)
ARITH_OR_COMPLEX_OVERLOAD(tan)
REAL128_OVERLOAD(tan)

VECTORIZED_OVERLOAD(sinh)
ARITH_OR_COMPLEX_OVERLOAD(sinh)
REAL128_OVERLOAD(sinh)

VECTORIZED_OVERLOAD(cosh)
ARITH_OR_COMPLEX_OVERLOAD(cosh)
REAL128_OVERLOAD(cosh)

VECTORIZED_OVERLOAD(tanh)
ARITH_OR_COMPLEX_OVERLOAD(tanh)
REAL128_OVERLOAD(tanh)

VECTORIZED_OVERLOAD(asin)
ARITH_OR_COMPLEX_OVERLOAD(asin)
REAL128_OVERLOAD(asin)

VECTORIZED_OVERLOAD(acos)
ARITH_OR_COMPLEX_OVERLOAD(acos)
REAL128_OVERLOAD(acos)

VECTORIZED_OVERLOAD(atan)
ARITH_OR_COMPLEX_OVERLOAD(atan)
REAL128_OVERLOAD(atan)

VECTORIZED_OVERLOAD(asinh)
ARITH_OR_COMPLEX_OVERLOAD(asinh)
REAL128_OVERLOAD(asinh)

VECTORIZED_OVERLOAD(acosh)
ARITH_OR_COMPLEX_OVERLOAD(acosh)
REAL128_OVERLOAD(acosh)

VECTORIZED_OVERLOAD(atanh)
ARITH_OR_COMPLEX_OVERLOAD(atanh)
REAL128_OVERLOAD(atanh)

ARITH_OR_COMPLEX_OVERLOAD(abs)

VECTORIZED_OVERLOAD(cbrt)
REAL128_OVERLOAD(cbrt)

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
template <typename T, typename U, enable_if_t<std::is_same<U, mppp::real128>::value || std::is_same<T, mppp::real128>::value, int> = 0>
inline mppp::real128 pow(const T &base, const U &d)
{
    return mppp::pow(base, d);
}

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

inline vectorized_double pow(double base, vectorized_double in)
{
    for (auto &el : in) {
        el = std::pow(base, el);
    }
    return in;
}

inline vectorized_double pow(vectorized_double in, double exponent)
{
    for (auto &el : in) {
        el = std::pow(el, exponent);
    }
    return in;
}


}
#endif
