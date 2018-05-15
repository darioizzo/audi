#ifndef AUDI_TYPE_TRAITS_HPP
#define AUDI_TYPE_TRAITS_HPP

#include <audi/config.hpp>

#if defined(AUDI_WITH_MPPP)
#include <mp++/real128.hpp>
#endif

#include <type_traits>

#include <audi/gdual.hpp>
#include <audi/vectorized_double.hpp>

namespace audi
{
/// Type is a gdual
/**
 * Checks whether T is a gdual type. Provides the member constant value which is equal to true,
 * if T is the type gdual<U> for any U.
 *
 * \tparam T a type to check
 */

template <typename T>
struct is_gdual : std::false_type {
};
template <typename T>
struct is_gdual<gdual<T>> : std::true_type {
};

/// Type is arithmetic (includes multiple precision floats)
/**
 * Checks whether T is an arithmetic type (that is, an integral type or a floating-point type that includes multiple
 * precision floats if available)
 *
 * \tparam T a type to check
 */
#if defined(AUDI_WITH_MPPP)
template <class T>
struct is_arithmetic
    : std::integral_constant<bool, std::is_arithmetic<T>::value || std::is_same<T, mppp::real128>::value> {
};
#else
template <class T>
struct is_arithmetic
    : std::integral_constant<bool, std::is_arithmetic<T>::value> {
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

/// Type is allowed to be a coefficient for gduals
/**
 * Checks whether T is a type from which gdual construction is allowed
 *
 * \tparam T a type to check
 */
template <class T>
struct is_gdual_cf : std::integral_constant<bool, audi::is_arithmetic_or_complex<T>::value
                                                      || std::is_same<T, vectorized_double>::value> {
};
} // namespace audi

#endif
