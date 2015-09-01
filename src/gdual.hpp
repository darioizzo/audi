#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

#include <stdexcept>
#include <string>
#include <piranha.hpp>

#include "operators.hpp"

namespace audi
{

class gdual
{
		friend gdual operator+(const gdual &d1, const gdual &d2);
		friend gdual operator-(const gdual &d1, const gdual &d2);
		friend gdual operator*(const gdual &d1, const gdual &d2);
		friend gdual operator/(const gdual &d1, const gdual &d2);
		using cf_type = double;
		using p_type = piranha::polynomial<cf_type,piranha::k_monomial>;
	public:
		explicit gdual(const std::string &str, int limit):m_p(str),m_limit(limit) {}
		explicit gdual(double x, int limit):m_p(x),m_limit(limit) {}
		int get_order() const
		{
			return m_order;
		}

		template <typename T>
		cf_type find_cf(const T &c) const
		{
			return m_p.find_cf(c);
		}
		
		template <typename T>
		cf_type find_cf(std::initializer_list<T> l) const
		{
			return m_p.find_cf(l);
		}

	private:
		p_type		m_p;
		const int	m_order;
};
} // end of namespace audi 
#endif