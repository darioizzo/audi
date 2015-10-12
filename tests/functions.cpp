#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_functions_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(exponentiation)
{
    {
        gdual x("x",3);
        gdual y("y",3);

        auto p1 = x*x*y + x*y*x*x*x - 3*y*y*y*y*x*y*x + 3.2;
        BOOST_CHECK_EQUAL(pow(p1, 3), p1*p1*p1);                                    // calls pow(gdual, int)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, 3.), p1*p1*p1, 1e-12) == true);                // calls pow(gdual, double)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual(3,3)), p1*p1*p1, 1e-12) == true);        // calls pow(gdual, gdual) with an integer exponent
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual(3.1,3)), pow(p1,3.1), 1e-12) == true);   // calls pow(gdual, gdual) with an real exponent
    }

    gdual x("x",3);
    gdual y("y",3);
    gdual p1 = x+y-3*x*y+y*y;
    gdual p2 = p1 - 3.5;

    BOOST_CHECK(EPSILON_COMPARE(pow(3, gdual(3.2, 3)), std::pow(3, 3.2) * gdual(1, 3), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(pow(p2, 3), pow(p2, 3.), 1e-12) == true);

    BOOST_CHECK(EPSILON_COMPARE(pow(p2, -1), 1 / p2, 1e-12) == true);                                      // negative exponent (gdual, int)
    BOOST_CHECK(EPSILON_COMPARE(pow(p2, -1.), 1 / p2, 1e-12) == true);                                     // negative exponent (gdual, double)
    BOOST_CHECK(EPSILON_COMPARE(pow(p1 + 3.5, gdual(-1.1, 3)), pow(p1 + 3.5, -1.1), 1e-12) == true);       // negative exponent (gdual, gdual)
}

BOOST_AUTO_TEST_CASE(square_root)
{
    gdual x(2.3, "x",3);
    gdual y(1.5, "y",3);

    auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;  // positive p0
    auto p2 = x*x*y - x*y*x*x*x - 3*y*y*y*y*x*y*x;  // negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(sqrt(p1) * sqrt(p1), p1, 1e-12) == true);
}

BOOST_AUTO_TEST_CASE(cubic_root)
{
    gdual x(2.3, "x",3);
    gdual y(1.5, "y",3);

    auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;  // positive p0
    auto p2 = x*x*y - x*y*x*x*x - 3*y*y*y*y*x*y*x;  // negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p1) * cbrt(p1) * cbrt(p1), p1, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p2) * cbrt(p2) * cbrt(p2), p2, 1e-12) == true);
}

BOOST_AUTO_TEST_CASE(exponential)
{
    // we test the result of f = exp((1+dx)*dy)
    gdual x(1, "x",2);
    gdual y(0, "y",2);

    auto p1 = exp(x*y);
    BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1,0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0,1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 0.5);
}

BOOST_AUTO_TEST_CASE(logarithm)
{
    {
        gdual x(0.1, "x",4);
        gdual y(0.13, "y",4);

        auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;
        BOOST_CHECK(EPSILON_COMPARE(exp(log(p1)), p1, 1e-12) == true);
    }

    gdual x("x",3);
    gdual y("y",3);
    gdual p1 = x+y-3*x*y+y*y;
    gdual p2 = p1 - 3.5;

}


BOOST_AUTO_TEST_CASE(sine_and_cosine)
{
    unsigned int order = 8;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  

    BOOST_CHECK(EPSILON_COMPARE(sin(2. * p1), 2. * sin(p1) * cos(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), 1. - 2. * sin(p1) * sin(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), 2. * cos(p1) * cos(p1) - 1., 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), cos(p1) * cos(p1) - sin(p1) * sin(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(sin(p1) * sin(p1) + cos(p1) * cos(p1), gdual(1.), 1e-12) == true);

    gdual sine(p1);
    gdual cosine(p1);
    sin_and_cos(p1, sine, cosine);
    BOOST_CHECK_EQUAL(sine, sin(p1));
    BOOST_CHECK_EQUAL(cosine, cos(p1));
}

BOOST_AUTO_TEST_CASE(tangent)
{
    {
    unsigned int order = 5;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }

    {
    unsigned int order = 10;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }
    
    {
    unsigned int order = 11; // tolerance decreases here to 1e-11 as high order derivatives can loose precision
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-11) == true);
    }
}

BOOST_AUTO_TEST_CASE(hyperbolic_sine_and_cosine)
{
    unsigned int order = 8;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);

    auto p1 = x + y;  

    // Checking some trivial identities
    BOOST_CHECK(EPSILON_COMPARE(sinh(2. * p1), 2. * sinh(p1) * cosh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), 1. + 2. * sinh(p1) * sinh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), 2. * cosh(p1) * cosh(p1) - 1., 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), cosh(p1) * cosh(p1) + sinh(p1) * sinh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(p1) * cosh(p1) - sinh(p1) * sinh(p1), gdual(1), 1e-12) == true);

    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    BOOST_CHECK(EPSILON_COMPARE(sinh(p1), (exp(p1) - exp(-p1)) / 2, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(p1), (exp(p1) + exp(-p1)) / 2, 1e-12) == true);

    gdual sineh(p1);
    gdual cosineh(p1);
    sinh_and_cosh(p1, sineh, cosineh);
    BOOST_CHECK_EQUAL(sineh, sinh(p1));
    BOOST_CHECK_EQUAL(cosineh, cosh(p1));
}

BOOST_AUTO_TEST_CASE(hyperbolic_tangent)
{
    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    {
    unsigned int order = 5;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tanh(p1), (exp(p1) - exp(-p1)) / (exp(p1) + exp(-p1)), 1e-12) == true);
    }

    // Checking the validity and precision of the identity tanh = sinh/cosh
    {
    unsigned int order = 5;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }

    {
    unsigned int order = 10;
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }
    
    {
    unsigned int order = 11; // tolerance decreases here to 1e-11 as high order derivatives can loose precision
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-11) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_tangent)
{
    // Checking the atanh is the inverse of tanh
    {
    unsigned int order = 5;
    gdual x(0.1, "x",order);
    gdual y(0.12, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(atanh(tanh(p1)), p1, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(tanh(atanh(p1)), p1, 1e-12) == true);
    }
    {
    unsigned int order = 6;
    gdual x(0.1, "x",order);
    gdual y(-0.12, "y",order);
    auto p1 = x + y + x * y;  
    BOOST_CHECK(EPSILON_COMPARE(atanh(tanh(p1)), p1, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(tanh(atanh(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(absolute_value)
{
    unsigned int order = 5;
    // Case f_0 > 0
    {
    gdual x(2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK_EQUAL(p1, abs(p1));
    }

    // Case f_0 < 0
    {
    gdual x(-2.3, "x",order);
    gdual y(1.5, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK_EQUAL(-p1, abs(p1));
    }

    // Case f_0 == 0
    {
    gdual x(-2.3, "x",order);
    gdual y(2.3, "y",order);
    auto p1 = x + y;  
    BOOST_CHECK_EQUAL(p1, abs(p1));
    }
}