#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"
#include "../src/functions_from_d.hpp"

#define BOOST_TEST_MODULE audi_functions_from_d_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(functions_from_derivative)
{
    // Testing atanh. Function domain is (-1, 1)
    {
    unsigned int order = 5;
    gdual x(1.1, "x",order);
    gdual y(1.2, "y",order);

    auto p1 = 1. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(-1.1, "x",order);
    gdual y(-1.2, "y",order);
    gdual z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    // Testing atan. Function domain is (-\infty, \infty)
    {
    unsigned int order = 5;
    gdual x(1.1, "x",order);
    gdual y(1.2, "y",order);

    auto p1 = 1. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(atan(p1), atan_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(-1.1, "x",order);
    gdual y(-1.2, "y",order);
    gdual z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(atan(p1), atan_d(p1), 1e-12) == true);
    }

    // Testing asinh. Function domain is (-\infty, \infty)
    {
    unsigned int order = 5;
    gdual x(1.1, "x",order);
    gdual y(1.2, "y",order);

    auto p1 = 1. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(-1.1, "x",order);
    gdual y(-1.2, "y",order);
    gdual z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    // Testing asin. Function domain is (-1, 1)
    {
    unsigned int order = 5;
    gdual x(1.1, "x",order);
    gdual y(1.2, "y",order);

    auto p1 = 1. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(-1.1, "x",order);
    gdual y(-1.2, "y",order);
    gdual z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(asinh(p1), asinh_d(p1), 1e-12) == true);
    }

    // Testing acosh. Function domain is (1, \infty). The 10. is there to avoid numerical floating point problems
    {
    unsigned int order = 5;
    gdual x(0.1, "x",order);
    gdual y(0.2, "y",order);

    auto p1 = 10. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(acosh(p1), acosh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 9;
    gdual x(0.1, "x",order);
    gdual y(0.2, "y",order);
    gdual z(0.3, "z",order);

    auto p1 = 10. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(acosh(p1), acosh_d(p1), 1e-12) == true);
    }

    // Testing acos. Function domain is (-1, 1). 
    {
    unsigned int order = 5;
    gdual x(1.1, "x",order);
    gdual y(1.2, "y",order);

    auto p1 = 1. / (x + y);  
    BOOST_CHECK(EPSILON_COMPARE(acos(p1), acos_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(-1.1, "x",order);
    gdual y(-1.2, "y",order);
    gdual z(-1.3, "z",order);

    auto p1 = 1. / (x + y + z);  
    BOOST_CHECK(EPSILON_COMPARE(acos(p1), acos_d(p1), 1e-12) == true);
    }
   
}