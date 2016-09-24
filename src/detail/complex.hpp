#ifndef AUDI_COMPLEX_HPP
#define AUDI_COMPLEX_HPP

#include<boost/type_traits/is_complex.hpp>
#include<piranha/math.hpp>
#include<type_traits>
#include<vector>

namespace piranha { namespace math {

// These overloads are needed by piranha to define the partial and integral operations on
// the polynomial having coefficient in std::complex

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
    return std::pow(c,exp);
  };
};
}}

// We need the operator *,/ defined so that std::complex<double>(0.) * char(0) compiles
// as that is needed in partial. IF THE COMPILER DECIDES TO IMPLEMENT THESE, REMOVE.
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
}
#endif
