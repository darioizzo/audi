#ifndef AUDI_OVERLOADS_HPP
#define AUDI_OVERLOADS_HPP

#include <boost/type_traits/is_complex.hpp>
#include "../vectorized_double.hpp"

namespace audi {
// type traits
template <typename T>
struct is_gdual: std::false_type {};
template <typename T>
struct is_gdual<gdual<T>>: std::true_type {};

template< class T >
struct is_arithmetic_or_complex : std::integral_constant<bool, std::is_arithmetic<T>::value ||
                                                               boost::is_complex<T>::value> {};
// overloads for the coefficients
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T exp(T in) {
   return std::exp(in);
}
inline vectorized_double exp(vectorized_double in)
{
   for (auto &el : in)
   {
       el = std::exp(el);
   }
   return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T erf(T in) {
   return std::erf(in);
}
inline vectorized_double erf(vectorized_double in)
{
   for (auto &el : in)
   {
       el = std::erf(el);
   }
   return in;
}

template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T log(T in) {
    return std::log(in);
}
inline vectorized_double log(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::log(el);
    }
    return in;
}

template<typename T, typename U, std::enable_if_t<is_arithmetic_or_complex<T>::value &&
                                                  is_arithmetic_or_complex<U>::value, int> = 0>
inline T pow(const U& base, const T &d)
{
    return std::pow(base, d);
}
inline vectorized_double pow(double base, vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::pow(base, el);
    }
    return in;
}
inline vectorized_double pow(vectorized_double in, double exponent)
{
    for (auto &el : in)
    {
        el = std::pow(el, exponent);
    }
    return in;
}
template<typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
inline T cbrt(T in) {
    return std::cbrt(in);
}
template<typename T, std::enable_if_t<boost::is_complex<T>::value, int> = 0>
inline T cbrt(T in) {
    return std::pow(in, 1./3.); // needs a separate template as cbrt does not exist for complex types
}
inline vectorized_double cbrt(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::cbrt(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T sin(T in) {
    return std::sin(in);
}
inline vectorized_double sin(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::sin(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T cos(T in) {
    return std::cos(in);
}
inline vectorized_double cos(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::cos(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T tan(T in) {
    return std::tan(in);
}
inline vectorized_double tan(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::tan(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T sinh(T in) {
    return std::sinh(in);
}
inline vectorized_double sinh(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::sinh(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T cosh(T in) {
    return std::cosh(in);
}
inline vectorized_double cosh(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::cosh(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T tanh(T in) {
    return std::tanh(in);
}
inline vectorized_double tanh(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::tanh(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T atanh(T in) {
    return std::atanh(in);
}
inline vectorized_double atanh(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::atanh(el);
    }
    return in;
}
template<typename T, std::enable_if_t<is_arithmetic_or_complex<T>::value, int> = 0>
inline T atan(T in) {
    return std::atan(in);
}
inline vectorized_double atan(vectorized_double in)
{
    for (auto &el : in)
    {
        el = std::atan(el);
    }
    return in;
}
}
#endif
