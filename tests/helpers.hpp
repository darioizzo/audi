#ifndef AUDI_TEST_HELPER_HPP
#define AUDI_TEST_HELPER_HPP

#include "../src/gdual.hpp"

/// Boost 1.59 introduces a new floating point comparison method. In which case the below macro can substitute EPSILON_COMPARE
//#define BOOST_EQUAL_GDUALS_TOL(d1, d2) \
//{gdual zero = d2 - d1;                 \
//gdual copy(d1);                        \
//double max_cf = (std::max_element(copy._container().begin(), copy._container().end(), [](const detail::p_type::term_type& a, const detail::p_type::term_type& b){return (std::fabs(a.m_cf) < std::fabs(b.m_cf));}))->m_cf; \
//max_cf = std::fabs(max_cf);	           \
//for (auto it = zero._container().begin(); it != zero._container().end(); ++it) \
//{ \
//BOOST_TEST(std::fabs(it->m_cf) / max_cf == 0, boost::test_tools::tolerance(10 * std::numeric_limits<double>::epsilon())); \
//} \
//}

namespace audi
{
	using p_type = detail::p_type;

	// Compares two gduals allowing for a small epsilon tolerance of 10 * std::numeric_limits<double>::epsilon()
	bool EPSILON_COMPARE(const gdual& d1, const gdual& d2)
	{
		gdual zero = d2 - d1;
		gdual copy(d1);
		// extracts the maximum coefficient (in absolute value) from d1
		double max_cf = (std::max_element(copy._container().begin(), copy._container().end(), [](const p_type::term_type& a, const p_type::term_type& b){return (std::fabs(a.m_cf) < std::fabs(b.m_cf));}))->m_cf;
		max_cf = std::fabs(max_cf);
		// checks that the coefficients of d2 - d1 are all within a reltol <= eps
		for (auto it = zero._container().begin(); it != zero._container().end(); ++it)
		{
			if (std::fabs(it->m_cf) / max_cf > 10 * std::numeric_limits<double>::epsilon()) 
				{
					std::cout << "Failing to be within epsilon from zero: " << zero << std::endl;
					std::cout << "Coefficient: " << std::fabs(it->m_cf) << std::endl;
					std::cout << "Normalized coefficient: " << std::fabs(it->m_cf) / max_cf << std::endl;
					std::cout << "Allowed epsilon: " << 10 * std::numeric_limits<double>::epsilon() << std::endl;
					return false;
				}
		}
		return true;
	}

} // end of namespace audi 
#endif