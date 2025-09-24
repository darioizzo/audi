#ifndef AUDI_TEST_HELPER_HPP
#define AUDI_TEST_HELPER_HPP

#include <algorithm>
#include <cmath>
#include <iostream>

#include <audi/gdual.hpp>
#include <audi/taylor_model.hpp>

/*
// Boost 1.59 introduces a new floating point comparison method. In which case the below macro can substitute
EPSILON_COMPARE
#define BOOST_EQUAL_GDUALS_TOL(d1, d2) \
{using p_type = detail::p_type;        \
gdual zero = d2 - d1;                  \
gdual copy(d1);                        \
double max_cf = (std::max_element(copy._container().begin(), copy._container().end(), [](const
detail::p_type::term_type& a, const detail::p_type::term_type& b){return (std::abs(a.m_cf) <
std::abs(b.m_cf));}))->m_cf; \
max_cf = std::abs(max_cf);	           \
for (auto it = zero._container().begin(); it != zero._container().end(); ++it) \
{ \
BOOST_TEST(std::abs(it->m_cf) / max_cf == 0, boost::test_tools::tolerance(10 * std::numeric_limits<double>::epsilon()));
\
} \
}*/

namespace audi
{

// Compares two taylor_models allowing for a small epsilon tolerance of 10 * std::numeric_limits<double>::epsilon()
inline bool EPSILON_COMPARE(const audi::taylor_model &d1, const audi::taylor_model &d2, const double epsilon)
{
    audi::taylor_model zero = d2 - d1;

    // checks that the coefficients of d2 - d1 are all within a reltol <= eps
    for (const auto &t : zero.get_tpol()._poly()) {
        if (std::abs(t.second) > epsilon) {
            std::cout << "Failing to be within epsilon from zero: " << zero << std::endl;
            std::cout << "Coefficient: " << std::abs(t.second) << std::endl;
            std::cout << "Allowed epsilon: " << epsilon << std::endl;
            return false;
        }
    }
    audi::taylor_model::map_interval_equal(d1.get_dom(), d2.get_dom(), epsilon);
    audi::taylor_model::map_equal(d1.get_exp(), d2.get_exp(), epsilon);
    
    return true;
}

// Compares two gduals allowing for a small epsilon tolerance of 10 * std::numeric_limits<double>::epsilon()
template <typename T, typename M>
inline bool EPSILON_COMPARE(const audi::gdual<T, M> &d1, const audi::gdual<T, M> &d2, const double epsilon)
{
    audi::gdual<T, M> zero = d2 - d1;

    // checks that the coefficients of d2 - d1 are all within a reltol <= eps
    for (const auto &t : zero._poly()) {
        if (std::abs(t.second) > epsilon) {
            std::cout << "Failing to be within epsilon from zero: " << zero << std::endl;
            std::cout << "Coefficient: " << std::abs(t.second) << std::endl;
            std::cout << "Allowed epsilon: " << epsilon << std::endl;
            return false;
        }
    }
    return true;
}

inline bool EPSILON_COMPARE(double d1, double d2, const double epsilon)
{
    double zero = std::abs(d2 - d1);
    if (zero > epsilon) {
        std::cout << "Failing to be within epsilon from zero: " << zero << std::endl;
        return false;
    }
    return true;
}

} // end of namespace audi

#endif

