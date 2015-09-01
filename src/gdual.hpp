#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

#include <stdexcept>
#include <string>
#include <iostream>
#include <piranha.hpp>

#include "functions.hpp"

namespace audi
{

class gdual
{
	// The type for the truncated Taylor polynomial
	using cf_type = double;
	using p_type = piranha::polynomial<cf_type,piranha::k_monomial>;

	// We disable the overloads of the +,-,*,/ operators if none of the operands is a gdual
	template <typename T, typename U>
	using magic_type = typename std::enable_if<std::is_same<T,gdual>::value || std::is_same<U,gdual>::value,gdual>::type;

public:
	explicit gdual(const std::string &str, int order):m_p(str),m_order(order) {}
	explicit gdual(double x, int order):m_p(x),m_order(order) {}
	int get_order() const
	{
		return m_order;
	}

	friend std::ostream& operator<<(std::ostream& os, const gdual& d)
	{
	    os << d.m_p;
	    return os;
	}

	template <typename T, typename U>
	friend magic_type<T, U> operator+(const T &d1, const U &d2)
	{
		return binary_add(d1,d2);
	}

private:
	static gdual binary_add(const gdual  &d1, const gdual  &d2) {
		if (d1.get_order() != d2.get_order()) {
			throw std::invalid_argument("different truncation limit");
		}
		gdual retval(d1);
		retval.m_p += d2.m_p;
		return retval;
	}

	template <typename T>
	static gdual binary_add(const T &d1, const gdual  &d2) {
		gdual retval(d2);
		retval.m_p += d1;
		return retval;
	}

	template <typename T>
	static gdual binary_add(const gdual  &d1, const T &d2) {
		return binary_add(d2, d1);
	}

	p_type		m_p;
	const int	m_order;
};
} // end of namespace audi 
#endif