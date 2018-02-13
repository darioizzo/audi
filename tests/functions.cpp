#define BOOST_TEST_MODULE audi_functions_test
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/functions.hpp>
#include <audi/gdual.hpp>

using namespace audi;
using gdual_d = gdual<double>;

BOOST_AUTO_TEST_CASE(exponentiation)
{
    {
        gdual_d x(0., "x", 3);
        gdual_d y(0., "y", 3);

        auto p1 = x * x * y + x * y * x * x * x - 3 * y * y * y * y * x * y * x + 3.2;
        BOOST_CHECK_EQUAL(pow(p1, 3), p1 * p1 * p1);  // calls pow(gdual_d, int)
        BOOST_CHECK_EQUAL(pow(p1, 3.), p1 * p1 * p1); // calls pow(gdual_d, double) (with a positive integer exponent)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual_d(3)), p1 * p1 * p1, 1e-12)
                    == true); // calls pow(gdual_d, gdual_d) with an integer exponent
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual_d(3.1)), pow(p1, 3.1), 1e-12)
                    == true); // calls pow(gdual_d, gdual_d) with an real exponent
    }

    {
        gdual_d x(0., "x", 3);
        gdual_d y(0., "y", 3);
        gdual_d p1 = x + y - 3 * x * y + y * y;
        gdual_d p2 = p1 - 3.5;

        BOOST_CHECK(EPSILON_COMPARE(pow(3, gdual_d(3.2)), std::pow(3, 3.2) * gdual_d(1), 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(pow(p2, 3), pow(p2, 3.), 1e-12) == true);

        BOOST_CHECK(EPSILON_COMPARE(pow(p2, -1), 1 / p2, 1e-12) == true);  // negative exponent (gdual_d, int)
        BOOST_CHECK(EPSILON_COMPARE(pow(p2, -1.), 1 / p2, 1e-12) == true); // negative exponent (gdual_d, double)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1 + 3.5, gdual_d(-1.1)), pow(p1 + 3.5, -1.1), 1e-12)
                    == true); // negative exponent (gdual_d, gdual_d)
    }

    // We check the implementation of pow(gdual_d, double) with respect to the behaviour in 0.
    // We compute the Taylor expansion of f = x^3.1 around 0. Which is T_f = 0. + 0.dx + 0.dx^2 + 0.dx^3 + inf dx^4-inf
    // dx^5 ...

    gdual_d x(0., "x", 7);
    auto f = pow(x, 3.1);
    BOOST_CHECK_EQUAL(f.find_cf({0}), 0.);
    BOOST_CHECK_EQUAL(f.find_cf({1}), 0.);
    BOOST_CHECK_EQUAL(f.find_cf({2}), 0.);
    BOOST_CHECK_EQUAL(f.find_cf({3}), 0.);
    BOOST_CHECK_EQUAL(f.find_cf({4}), 1. / 0.);
    BOOST_CHECK_EQUAL(f.find_cf({5}), -1. / 0.);
    BOOST_CHECK_EQUAL(f.find_cf({6}), 1. / 0.);
    BOOST_CHECK_EQUAL(f.find_cf({7}), -1. / 0.);
}

BOOST_AUTO_TEST_CASE(square_root)
{
    gdual_d x(2.3, "x", 3);
    gdual_d y(1.5, "y", 3);

    auto p1 = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x; // positive p0
    auto p2 = x * x * y - x * y * x * x * x - 3 * y * y * y * y * x * y * x; // negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(sqrt(p1) * sqrt(p1), p1, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(sqrt(p2) * sqrt(p2), p2, 1e-12) == true);
}

BOOST_AUTO_TEST_CASE(cubic_root)
{
    gdual_d x(2.3, "x", 3);
    gdual_d y(1.5, "y", 3);

    auto p1 = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x; // positive p0
    auto p2 = x * x * y - x * y * x * x * x - 3 * y * y * y * y * x * y * x; // negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p1) * cbrt(p1) * cbrt(p1), p1, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p2) * cbrt(p2) * cbrt(p2), p2, 1e-12) == true);
}

BOOST_AUTO_TEST_CASE(exponential)
{
    // we test the result of f = exp((1+dx)*dy)
    gdual_d x(1, "x", 2);
    gdual_d y(0, "y", 2);

    auto p1 = exp(x * y);
    BOOST_CHECK_EQUAL(p1.find_cf({0, 0}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1, 0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0, 1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1, 1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({2, 0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0, 2}), 0.5);
}

BOOST_AUTO_TEST_CASE(logarithm)
{
    {
        gdual_d x(0.1, "x", 4);
        gdual_d y(0.13, "y", 4);

        auto p1 = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x;
        BOOST_CHECK(EPSILON_COMPARE(exp(log(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(sine_and_cosine)
{
    unsigned int order = 8;
    gdual_d x(2.3, "x", order);
    gdual_d y(1.5, "y", order);

    auto p1 = x + y;

    BOOST_CHECK(EPSILON_COMPARE(sin(2. * p1), 2. * sin(p1) * cos(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), 1. - 2. * sin(p1) * sin(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), 2. * cos(p1) * cos(p1) - 1., 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cos(2. * p1), cos(p1) * cos(p1) - sin(p1) * sin(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(sin(p1) * sin(p1) + cos(p1) * cos(p1), gdual_d(1.), 1e-12) == true);

    auto res = sin_and_cos(p1);
    BOOST_CHECK_EQUAL(res[0], sin(p1));
    BOOST_CHECK_EQUAL(res[1], cos(p1));
}

BOOST_AUTO_TEST_CASE(tangent)
{
    {
        unsigned int order = 5;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);

        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }

    {
        unsigned int order = 6;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);

        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }

    {
        unsigned int order = 10;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);

        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-12) == true);
    }

    {
        unsigned int order = 11; // tolerance decreases here to 1e-11 as high order derivatives can loose precision
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);

        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tan(p1), sin(p1) / cos(p1), 1e-11) == true);
    }
}

BOOST_AUTO_TEST_CASE(hyperbolic_sine_and_cosine)
{
    unsigned int order = 8;
    gdual_d x(2.3, "x", order);
    gdual_d y(1.5, "y", order);

    auto p1 = x + y;

    // Checking some trivial identities
    BOOST_CHECK(EPSILON_COMPARE(sinh(2. * p1), 2. * sinh(p1) * cosh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), 1. + 2. * sinh(p1) * sinh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), 2. * cosh(p1) * cosh(p1) - 1., 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(2. * p1), cosh(p1) * cosh(p1) + sinh(p1) * sinh(p1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(p1) * cosh(p1) - sinh(p1) * sinh(p1), gdual_d(1), 1e-12) == true);

    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    BOOST_CHECK(EPSILON_COMPARE(sinh(p1), (exp(p1) - exp(-p1)) / 2, 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE(cosh(p1), (exp(p1) + exp(-p1)) / 2, 1e-12) == true);

    gdual_d sineh(p1);
    gdual_d cosineh(p1);
    auto res = sinh_and_cosh(p1);
    BOOST_CHECK_EQUAL(res[0], sinh(p1));
    BOOST_CHECK_EQUAL(res[1], cosh(p1));
}

BOOST_AUTO_TEST_CASE(hyperbolic_tangent)
{
    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    {
        unsigned int order = 5;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tanh(p1), (exp(p1) - exp(-p1)) / (exp(p1) + exp(-p1)), 1e-12) == true);
    }

    // Checking the validity and precision of the identity tanh = sinh/cosh
    {
        unsigned int order = 5;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }

    {
        unsigned int order = 6;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }

    {
        unsigned int order = 8;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }

    {
        unsigned int order = 9;
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK(EPSILON_COMPARE(tanh(p1), sinh(p1) / cosh(p1), 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_tangent)
{
    // Checking that atanh is the inverse of tanh
    {
        unsigned int order = 5;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(atanh(tanh(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(tanh(atanh(p1)), p1, 1e-12) == true);
    }
    {
        unsigned int order = 6;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(atanh(tanh(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(tanh(atanh(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_tangent)
{
    // Checking that atan is the inverse of tan
    {
        unsigned int order = 5;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(atan(tan(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(tan(atan(p1)), p1, 1e-12) == true);
    }
    {
        unsigned int order = 6;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(atan(tan(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(tan(atan(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_sine)
{
    // Checking that asinh is the inverse of sinh
    {
        unsigned int order = 5;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(asinh(sinh(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(sinh(asinh(p1)), p1, 1e-12) == true);
    }
    {
        unsigned int order = 6;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(asinh(sinh(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(sinh(asinh(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_sine)
{
    // Checking that asin is the inverse of sin
    {
        unsigned int order = 5;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(asin(sin(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(sin(asin(p1)), p1, 1e-12) == true);
    }
    {
        unsigned int order = 6;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(asin(sin(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(sin(asin(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_cosine)
{
    // Checking that acosh is the inverse of acos (PRECISION IS A PROBLEM HERE ALREADY AT LOW ORDERS!! RELATED TO LOG?
    // OR TO SQRT? OR DIV?)
    {
        unsigned int order = 4;
        gdual_d x(0.1, "x", order);
        gdual_d y(0.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(acosh(cosh(p1)), p1, 1e-8) == true);
        BOOST_CHECK(EPSILON_COMPARE(cosh(acosh(p1)), p1, 1e-8) == true);
    }
    {
        unsigned int order = 5;
        gdual_d x(0.1, "x", order);
        gdual_d y(0.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(acosh(cosh(p1)), p1, 1e-8) == true);
        BOOST_CHECK(EPSILON_COMPARE(cosh(acosh(p1)), p1, 1e-8) == true);
    }
}

BOOST_AUTO_TEST_CASE(inverse_cosine)
{
    // Checking that acos is the inverse of cos
    {
        unsigned int order = 5;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(acos(cos(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(cos(acos(p1)), p1, 1e-12) == true);
    }
    {
        unsigned int order = 6;
        gdual_d x(1.1, "x", order);
        gdual_d y(1.2, "y", order);
        auto p1 = 1. / (x + y);
        BOOST_CHECK(EPSILON_COMPARE(acos(cos(p1)), p1, 1e-12) == true);
        BOOST_CHECK(EPSILON_COMPARE(cos(acos(p1)), p1, 1e-12) == true);
    }
}

BOOST_AUTO_TEST_CASE(absolute_value)
{
    unsigned int order = 5;
    // Case f_0 > 0
    {
        gdual_d x(2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK_EQUAL(p1, abs(p1));
    }

    // Case f_0 < 0
    {
        gdual_d x(-2.3, "x", order);
        gdual_d y(1.5, "y", order);
        auto p1 = x + y;
        BOOST_CHECK_EQUAL(-p1, abs(p1));
    }

    // Case f_0 == 0
    {
        gdual_d x(-2.3, "x", order);
        gdual_d y(2.3, "y", order);
        auto p1 = x + y;
        BOOST_CHECK_EQUAL(p1, abs(p1));
    }
}
