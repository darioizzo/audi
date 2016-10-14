#ifndef AUDI_TYPE_TRAITS_HPP
#define AUDI_TYPE_TRAITS_HPP

#include <type_traits>
#include "gdual.hpp"

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
struct is_gdual: std::false_type {};
template <typename T>
struct is_gdual<gdual<T> >: std::true_type {};

/// Type is arithmetic or complex
/**
 * Checks whether T is an arithmetic type (that is, an integral type or a floating-point type) or a complex
 * that is a type std::complex<U> for any U.
 *
 * \tparam T a type to check
*/
template< class T >
struct is_arithmetic_or_complex : std::integral_constant<bool, std::is_arithmetic<T>::value || boost::is_complex<T>::value> {};

}

#endif
