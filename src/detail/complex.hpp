#ifndef AUDI_COMPLEX_HPP
#define AUDI_COMPLEX_HPP

#include<boost/type_traits/is_complex.hpp>
#include<piranha/math.hpp>
#include<type_traits>
#include<vector>

namespace piranha { namespace math {

// This overload is needed for T = std::complex<double> "meet" the rquirement: piranha::is_differentiable<T>
// necessary for the methods partial, integrate and evaluate to be defined.
template <typename T>
struct partial_impl<T, typename std::enable_if<boost::is_complex<T>::value>::type>
{
   T operator()(const T &, const std::string &) const
   {
       return T(0);
   }
};

template <typename T, typename U>
struct pow_impl<T,U,typename std::enable_if<boost::is_complex<T>::value>::type>
{
  T operator()(const T &c, const U &exp) const
  {
    return piranha::math::pow(c,exp);
  };
};
}} // end of namespace piranha math

// These operators are not in the C++ standard for the type std::complex<double> and most
// compilers do not "offer" them, but they are necessary for the methods partial and integrate
// as Cf *,/ Key must be defined and, in this case that is std::complex<double> * char(0)
// They basically allow multiplication and decision between std::complex<double> and any
// integer type
namespace std{
template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
std::complex<double> operator*(const std::complex<double> &d1, T d2)
{
	return d1 * std::complex<double>(d2);
};

template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
std::complex<double> operator*(T d2, const std::complex<double> &d1)
{
	return d1 * std::complex<double>(d2);
};

template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
std::complex<double> operator/(const std::complex<double> &d1, T d2)
{
	return d1 / std::complex<double>(d2);
};

template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
std::complex<double> operator/(T d2, const std::complex<double> &d1)
{
	return d1 / std::complex<double>(d2);
};
} // end namespace std
#endif
