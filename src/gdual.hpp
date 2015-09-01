#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

#include <stdexcept>
#include <string>
#include <iostream>
#include <piranha/polynomial.hpp>
#include <algorithm> 

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
	using gdual_if_enabled = typename std::enable_if<std::is_same<T,gdual>::value || std::is_same<U,gdual>::value,gdual>::type;

public:
	explicit gdual(const std::string &str, int order):m_p(str),m_order(order) {
		if (order < 1) {
			throw std::invalid_argument("polynomial truncation order must be > 1");
		}
	}

	explicit gdual(double x, int order):m_p(x),m_order(order) {
		if (order < 1) {
			throw std::invalid_argument("polynomial truncation order must be > 1");
		}
	}

	int get_order() const
	{
		return m_order;
	}

	template <typename T>
	cf_type find_cf(const T &c) const
	{
		if (!std::all_of(c.begin(), c.end(), [&](int i){return i<=get_order();}) ||
		   !std::all_of(c.begin(), c.end(), [](int i){return i>=0;}))
		{
			throw std::invalid_argument("invalid order(s) requested for derivative ");
		}
		return m_p.find_cf(c);
	}

	template <typename T>
	cf_type find_cf(std::initializer_list<T> l) const
	{
		if (!std::all_of(l.begin(), l.end(), [&](T i){return i<=get_order();}) || 
		    !std::all_of(l.begin(), l.end(), [](T i){return i>=0;}))
		{
			throw std::invalid_argument("invalid order(s) requested for derivative ");
		}
		return m_p.find_cf(l);
	}

	friend std::ostream& operator<<(std::ostream& os, const gdual& d)
	{
	    os << d.m_p;
	    return os;
	}

	friend bool operator==(const gdual  &d1, const gdual &d2) 
	{
		return (d1.m_p == d2.m_p) && (d1.m_order == d2.m_order);
	}

	template <typename T, typename U>
	friend gdual_if_enabled<T, U> operator+(const T &d1, const U &d2)
	{
		return add(d1,d2);
	}

	template <typename T, typename U>
	friend gdual_if_enabled<T, U> operator-(const T &d1, const U &d2)
	{
		return sub(d1,d2);
	}

private:
	// Basic overloads for the addition
	static gdual add(const gdual  &d1, const gdual  &d2) {
		if (d1.get_order() != d2.get_order()) {
			throw std::invalid_argument("different truncation limit");
		}
		gdual retval(d1);
		retval.m_p += d2.m_p;
		return retval;
	}

	template <typename T>
	static gdual add(const T &d1, const gdual  &d2) {
		gdual retval(d2);
		retval.m_p += d1;
		return retval;
	}

	template <typename T>
	static gdual add(const gdual  &d1, const T &d2) {
		return add(d2, d1);
	}

	// Basic overloads for the subtraction
	static gdual sub(const gdual  &d1, const gdual  &d2) {
		if (d1.get_order() != d2.get_order()) {
			throw std::invalid_argument("different truncation limit");
		}
		gdual retval(d1);
		retval.m_p -= d2.m_p;
		return retval;
	}

	template <typename T>
	static gdual sub(const T &d1, const gdual  &d2) {
		gdual retval(d2);
		retval.m_p = d1 - retval.m_p;
		return retval;
	}

	template <typename T>
	static gdual sub(const gdual  &d1, const T &d2) {
		gdual retval(d1);
		retval.m_p -= d2;
		return retval;
	}

	p_type		m_p;
	const int	m_order;
};
} // end of namespace audi 
#endif