#include "helpers.hpp"
#include <audi/gdual.hpp>
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>

#define BOOST_TEST_MODULE audi_functions_from_d_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;
using gdual_d = gdual<double>;

BOOST_AUTO_TEST_CASE(functions_from_derivative)
{
    // Testing atanh. Function domain is (-1, 1)
    {
    unsigned int order = 5;
    gdual_d x(1.1, "x",order);
    gdual_d y(1.2, "y",order);

    auto p1 = 1. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual_d x(-1.1, "x",order);
    gdual_d y(-1.2, "y",order);
    gdual_d z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    // Testing atan. Function domain is (-\infty, \infty)
    {
    unsigned int order = 5;
    gdual_d x(1.1, "x",order);
    gdual_d y(1.2, "y",order);

    auto p1 = 1. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(atan(p1), atan_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual_d x(-1.1, "x",order);
    gdual_d y(-1.2, "y",order);
    gdual_d z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(atan(p1), atan_d(p1), 1e-12) == true);
    }

    // Testing asinh. Function domain is (-\infty, \infty)
    {
    unsigned int order = 5;
    gdual_d x(1.1, "x",order);
    gdual_d y(1.2, "y",order);

    auto p1 = 1. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual_d x(-1.1, "x",order);
    gdual_d y(-1.2, "y",order);
    gdual_d z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    // Testing asin. Function domain is (-1, 1)
    {
    unsigned int order = 5;
    gdual_d x(1.1, "x",order);
    gdual_d y(1.2, "y",order);

    auto p1 = 1. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual_d x(-1.1, "x",order);
    gdual_d y(-1.2, "y",order);
    gdual_d z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    // Testing acosh. Function domain is (1, \infty). The 10. is there to avoid numerical floating point problems
    {
    unsigned int order = 5;
    gdual_d x(0.1, "x",order);
    gdual_d y(0.2, "y",order);

    auto p1 = 10. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(acosh(p1), acosh_d(p1), 1e-10) == true);
    }

    {
    unsigned int order = 5;
    gdual_d x(0.1, "x",order);
    gdual_d y(0.2, "y",order);
    gdual_d z(0.3, "z",order);

    auto p1 = 10. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(acosh(p1), acosh_d(p1), 1e-10) == true);
    }

    // Testing acos. Function domain is (-1, 1).
    {
    unsigned int order = 5;
    gdual_d x(1.1, "x",order);
    gdual_d y(1.2, "y",order);

    auto p1 = 1. / (x + y);
    BOOST_CHECK(EPSILON_COMPARE(acos(p1), acos_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual_d x(-1.1, "x",order);
    gdual_d y(-1.2, "y",order);
    gdual_d z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);
    BOOST_CHECK(EPSILON_COMPARE(acos(p1), acos_d(p1), 1e-12) == true);
    }

}

BOOST_AUTO_TEST_CASE(error_function)
{
    // We check that all the maclaurin series coefficients (one variable) are correct up to order
    unsigned int order = 9;
    gdual_d x(0., "x", order);
    auto res = erf(x);
    double c = 2.;
    c /= std::sqrt(boost::math::constants::pi<double>());
    // odd powers
    for (auto i = 0u; 2 * i + 1 <= order; ++i) {
        double mac_term = std::pow(-1, i) / boost::math::factorial<double>(i) / (2 * i + 1);
        BOOST_CHECK(EPSILON_COMPARE(res.get_derivative({2 * i + 1}) / boost::math::factorial<double>(2 * i + 1), c * mac_term, 1e-14) == true);
    }
    // even powers
    for (auto i = 0u; 2 * i  <= order; ++i) {
        BOOST_CHECK_EQUAL(res.get_derivative({2 * i}), 0.);
    }
}
