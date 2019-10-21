#ifndef AUDI_TYPE_TRAITS_HPP
#define AUDI_TYPE_TRAITS_HPP

#include <boost/type_traits/is_complex.hpp>
#include <type_traits>

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)
#include <mp++/real128.hpp>
#endif

namespace audi
{
/// Type is arithmetic (includes multiple precision floats)
/**
 * Checks whether T is an arithmetic type (that is, an integral type or a floating-point type that includes multiple
 * precision floats if available)
 *
 * \tparam T a type to check
 */
#if defined(AUDI_WITH_QUADMATH)
template <class T>
struct is_arithmetic
    : std::integral_constant<bool, std::is_arithmetic<T>::value || std::is_same<T, mppp::real128>::value> {
};
#else
template <class T>
struct is_arithmetic : std::integral_constant<bool, std::is_arithmetic<T>::value> {
};
#endif
/// Type is arithmetic or complex
/**
 * Checks whether T is an arithmetic type (that is, an integral type or a floating-point type) or a complex
 * that is a type std::complex<U> for any U.
 *
 * \tparam T a type to check
 */
template <class T>
struct is_arithmetic_or_complex
    : std::integral_constant<bool, audi::is_arithmetic<T>::value || boost::is_complex<T>::value> {
};

} // namespace audi

#endif
