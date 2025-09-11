#define BOOST_TEST_MODULE taylor_model_functions_test
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/taylor_model.hpp>
#include <audi/taylor_model_functions.hpp>

#include <audi/io.hpp>

BOOST_AUTO_TEST_CASE(exponentiation)
{
    {

        // Constructing a multivariate Taylor model
        uint order = 3;
        int_d rem(0.0, 0.0);
        var_map_d exp = {{"x", 0.0}, {"y", 0.0}};
        double domain_size = 0.2;
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        audi::gdual<double> f_xy = x * x * y + x * y * x * x * x - 3 * y * y * y * y * x * y * x + 3.2;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK_EQUAL(audi::pow(tm_xy, 3), tm_xy * tm_xy * tm_xy); // calls pow(audi::gdual<double>, int)
        BOOST_CHECK_EQUAL(audi::pow(tm_xy, 3.),
                          tm_xy * tm_xy
                              * tm_xy); // calls pow(audi::gdual<double>, double) (with a positive integer exponent)
    }

    {
        uint order = 3;
        int_d rem(0.0, 0.0);
        var_map_d exp = {{"x", 0.0}, {"y", 0.0}};
        double domain_size = 0.2;
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        audi::gdual<double> f_xy = x + y - 3 * x * y + y * y;
        audi::taylor_model tm_fxy(f_xy, rem, exp, dom);
        audi::gdual<double> g_xy = f_xy - 3.5;
        audi::taylor_model tm_gxy(g_xy, rem, exp, dom);

        BOOST_CHECK(audi::pow(tm_gxy, 3) == audi::pow(tm_gxy, 3.));
        BOOST_CHECK(audi::pow(tm_gxy, -1) == 1 / tm_gxy);  // negative exponent (audi::gdual<double>, int)
        BOOST_CHECK(audi::pow(tm_gxy, -1.) == 1 / tm_gxy); // negative exponent (audi::gdual<double>, double)
    }
}

BOOST_AUTO_TEST_CASE(square_root)
{
    uint order = 3;
    int_d rem(0.0, 0.0);
    var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
    double domain_size = 0.2;
    var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                     {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
    audi::gdual<double> x(exp.find("x")->second, "x", order);
    audi::gdual<double> y(exp.find("y")->second, "y", order);

    // p2 from the gdual tests (with negative p0) did not work because the bounds become negative
    // and one cannot take the sqrt() of a negative interval
    auto f_xy = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x; // positive p0
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    // Epsilon check necessary because the extra operations introduce remainders that are non-zero
    BOOST_CHECK(audi::EPSILON_COMPARE(audi::sqrt(tm_xy) * audi::sqrt(tm_xy), tm_xy, 1e-12));
}

BOOST_AUTO_TEST_CASE(exponential)
{
    // we test the result of f = exp((1+dx)*dy)
    uint order = 2;
    int_d rem(0.0, 0.0);
    var_map_d exp = {{"x", 1}, {"y", 0}};
    double domain_size = 0.2;
    var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                     {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
    audi::gdual<double> x(exp.find("x")->second, "x", order);
    audi::gdual<double> y(exp.find("y")->second, "y", order);

    auto f_xy = audi::exp(x * y);
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({0, 0}), 1);
    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({1, 0}), 0);
    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({0, 1}), 1);
    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({1, 1}), 1);
    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({2, 0}), 0);
    BOOST_CHECK_EQUAL(tm_xy.get_tpol().find_cf({0, 2}), 0.5);
}

BOOST_AUTO_TEST_CASE(logarithm)
{
    uint order = 4;
    int_d rem(0.0, 0.0);
    var_map_d exp = {{"x", 0.1}, {"y", 0.13}};
    double domain_size = 0.2;
    var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                     {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
    audi::gdual<double> x(exp.find("x")->second, "x", order);
    audi::gdual<double> y(exp.find("y")->second, "y", order);

    // Added 5 w.r.t. gdual test because the remainder bound could be negative and the logarithm is
    // therefore invalid
    auto f_xy = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x + 5;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    // Epsilon check necessary because the extra operations introduce remainders that are non-zero
    BOOST_CHECK(EPSILON_COMPARE(audi::exp(audi::log(tm_xy)), tm_xy, 1e-12));
}

BOOST_AUTO_TEST_CASE(sine_and_cosine)
{
    uint order = 8;
    int_d rem(0.0, 0.0);
    var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
    double domain_size = 0.2;
    var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                     {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
    audi::gdual<double> x(exp.find("x")->second, "x", order);
    audi::gdual<double> y(exp.find("y")->second, "y", order);

    auto f_xy = x + y;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    BOOST_CHECK(EPSILON_COMPARE(audi::sin(2. * tm_xy), 2. * audi::sin(tm_xy) * audi::cos(tm_xy), 1e-10));
    BOOST_CHECK(EPSILON_COMPARE(audi::cos(2. * tm_xy), 1. - 2. * audi::sin(tm_xy) * audi::sin(tm_xy), 1e-10));
    BOOST_CHECK(EPSILON_COMPARE(audi::cos(2. * tm_xy), 2. * audi::cos(tm_xy) * audi::cos(tm_xy) - 1., 1e-10));
    BOOST_CHECK(EPSILON_COMPARE(audi::cos(2. * tm_xy),
                                audi::cos(tm_xy) * audi::cos(tm_xy) - audi::sin(tm_xy) * audi::sin(tm_xy), 1e-10));
    BOOST_CHECK(EPSILON_COMPARE(audi::sin(tm_xy) * audi::sin(tm_xy) + audi::cos(tm_xy) * audi::cos(tm_xy),
                                audi::taylor_model(1.), 1e-10));

    auto res = audi::sin_and_cos(tm_xy);
    BOOST_CHECK(res[0] == audi::sin(tm_xy));
    BOOST_CHECK(res[1] == audi::cos(tm_xy));
}

BOOST_AUTO_TEST_CASE(tangent)
{
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        double domain_size = 0.2;
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK(EPSILON_COMPARE(audi::tan(tm_xy), audi::sin(tm_xy) / audi::cos(tm_xy), 1e-12));
    }

    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK(EPSILON_COMPARE(audi::tan(tm_xy), audi::sin(tm_xy) / audi::cos(tm_xy), 1e-12));
    }

    {
        uint order = 10;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK(EPSILON_COMPARE(audi::tan(tm_xy), audi::sin(tm_xy) / audi::cos(tm_xy), 1e-12));
    }

    {
        uint order = 11;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::tan(tm_xy), audi::sin(tm_xy) / audi::cos(tm_xy), 1e-12));
    }
}

BOOST_AUTO_TEST_CASE(hyperbolic_sine_and_cosine)
{
    uint order = 8;
    int_d rem(0.0, 0.0);
    double domain_size = 0.2;
    var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
    var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                     {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
    audi::gdual<double> x(exp.find("x")->second, "x", order);
    audi::gdual<double> y(exp.find("y")->second, "y", order);

    auto f_xy = x + y;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    // Checking some trivial identities
    BOOST_CHECK(EPSILON_COMPARE(audi::sinh(2. * tm_xy), 2. * audi::sinh(tm_xy) * audi::cosh(tm_xy), 1e-12));
    BOOST_CHECK(EPSILON_COMPARE(audi::cosh(2. * tm_xy), 1. + 2. * audi::sinh(tm_xy) * audi::sinh(tm_xy), 1e-12));
    BOOST_CHECK(EPSILON_COMPARE(audi::cosh(2. * tm_xy), 2. * audi::cosh(tm_xy) * audi::cosh(tm_xy) - 1., 1e-12));
    BOOST_CHECK(EPSILON_COMPARE(audi::cosh(2. * tm_xy),
                                audi::cosh(tm_xy) * audi::cosh(tm_xy) + audi::sinh(tm_xy) * audi::sinh(tm_xy), 1e-12));
    BOOST_CHECK(EPSILON_COMPARE(audi::cosh(tm_xy) * audi::cosh(tm_xy) - audi::sinh(tm_xy) * audi::sinh(tm_xy),
                                audi::taylor_model(1.0), 1e-12));

    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    BOOST_CHECK(EPSILON_COMPARE(audi::sinh(tm_xy), (audi::exp(tm_xy) - audi::exp(-tm_xy)) / 2, 1e-12));
    BOOST_CHECK(EPSILON_COMPARE(audi::cosh(tm_xy), (audi::exp(tm_xy) + audi::exp(-tm_xy)) / 2, 1e-12));

    auto res = audi::sinh_and_cosh(tm_xy);
    BOOST_CHECK(res[0] == audi::sinh(tm_xy));
    BOOST_CHECK(res[1] == audi::cosh(tm_xy));
}

BOOST_AUTO_TEST_CASE(hyperbolic_tangent)
{
    // Checking the validity of the definitions of hyperbolic finctions in terms of exponentials
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK(EPSILON_COMPARE(
            audi::tanh(tm_xy), (audi::exp(tm_xy) - audi::exp(-tm_xy)) / (audi::exp(tm_xy) + audi::exp(-tm_xy)), 1e-10));
    }

    // Checking the validity and precision of the identity audi::tanh = audi::sinh/audi::cosh
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);

        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(tm_xy), audi::sinh(tm_xy) / audi::cosh(tm_xy), 1e-12));
    }

    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(tm_xy), audi::sinh(tm_xy) / audi::cosh(tm_xy), 1e-12));
    }

    {
        uint order = 8;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(tm_xy), audi::sinh(tm_xy) / audi::cosh(tm_xy), 1e-12));
    }

    {
        uint order = 9;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(tm_xy), audi::sinh(tm_xy) / audi::cosh(tm_xy), 1e-12));
    }
}
BOOST_AUTO_TEST_CASE(inverse_hyperbolic_tangent)
{
    // Checking that atanh is the inverse of tanh
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.1;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::atanh(audi::tanh(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(audi::atanh(tm_xy)), tm_xy, 1e-10));
    }
    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.1;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::atanh(audi::tanh(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::tanh(audi::atanh(tm_xy)), tm_xy, 1e-10));
    }
}

BOOST_AUTO_TEST_CASE(inverse_tangent)
{
    // Checking that atan is the inverse of tan
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::atan(audi::tan(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::tan(audi::atan(tm_xy)), tm_xy, 1e-10));
    }
    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::atan(audi::tan(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::tan(audi::atan(tm_xy)), tm_xy, 1e-10));
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_sine)
{
    // Checking that asinh is the inverse of sinh
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::asinh(audi::sinh(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::sinh(audi::asinh(tm_xy)), tm_xy, 1e-10));
    }
    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::asinh(audi::sinh(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::sinh(audi::asinh(tm_xy)), tm_xy, 1e-10));
    }
}

BOOST_AUTO_TEST_CASE(arcsine_derivative)
{
    int_d a = int_d(0.0, 0.5);
    BOOST_CHECK_EQUAL(
        audi::taylor_model::interval_equal(audi::asin_derivative(a, 0), boost::numeric::asin(a)), true);
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(audi::asin_derivative(a, 1), int_d(1.0) / boost::numeric::sqrt(int_d(1.0) - a * a)), true);
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(audi::asin_derivative(a, 2),
                a * boost::numeric::pow(int_d(1.0) - a * a, 2) / boost::numeric::pow(int_d(1.0) - a * a, 3)), true);
}

BOOST_AUTO_TEST_CASE(inverse_sine)
{
    // Checking that asin is the inverse of sin
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::asin(audi::sin(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::sin(audi::asin(tm_xy)), tm_xy, 1e-10));
    }
    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::asin(audi::sin(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::sin(audi::asin(tm_xy)), tm_xy, 1e-10));
    }
}

BOOST_AUTO_TEST_CASE(inverse_hyperbolic_cosine)
{
    // Comment from gdual tests:
    // Checking that acosh is the inverse of acos (PRECISION IS A PROBLEM HERE ALREADY AT LOW ORDERS!! RELATED TO LOG?
    // OR TO SQRT? OR DIV?)
    {
        uint order = 4;
        int_d rem(0.0, 0.0);
        double domain_size = 0.01;
        var_map_d exp = {{"x", 0.1}, {"y", 0.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::acosh(audi::cosh(tm_xy)), tm_xy, 1e-8));
        BOOST_CHECK(EPSILON_COMPARE(audi::cosh(audi::acosh(tm_xy)), tm_xy, 1e-8));
    }
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.01;
        var_map_d exp = {{"x", 0.1}, {"y", 0.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::acosh(audi::cosh(tm_xy)), tm_xy, 1e-8));
        BOOST_CHECK(EPSILON_COMPARE(audi::cosh(audi::acosh(tm_xy)), tm_xy, 1e-8));
    }
}

BOOST_AUTO_TEST_CASE(inverse_cosine)
{
    // Checking that acos is the inverse of cos
    {
        uint order = 5;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::acos(audi::cos(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::cos(audi::acos(tm_xy)), tm_xy, 1e-10));
    }
    {
        uint order = 6;
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 1.1}, {"y", 1.2}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = 1 / (x + y);
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK(EPSILON_COMPARE(audi::acos(audi::cos(tm_xy)), tm_xy, 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(audi::cos(audi::acos(tm_xy)), tm_xy, 1e-10));
    }
}

BOOST_AUTO_TEST_CASE(absolute_value)
{
    uint order = 5;
    // Case f_0 > 0
    {
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", 2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK_EQUAL(tm_xy, audi::abs(tm_xy));
    }

    // Case f_0 < 0
    {
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", -2.3}, {"y", 1.5}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK_EQUAL(-tm_xy, audi::abs(tm_xy));
    }

    // Case f_0 == 0
    {
        int_d rem(0.0, 0.0);
        double domain_size = 0.2;
        var_map_d exp = {{"x", -2.3}, {"y", 2.3}};
        var_map_i dom = {{"x", int_d(exp.find("x")->second - domain_size, exp.find("x")->second + domain_size)},
                         {"y", int_d(exp.find("y")->second - domain_size, exp.find("y")->second + domain_size)}};
        audi::gdual<double> x(exp.find("x")->second, "x", order);
        audi::gdual<double> y(exp.find("y")->second, "y", order);

        auto f_xy = x + y;
        audi::taylor_model tm_xy(f_xy, rem, exp, dom);
        BOOST_CHECK_EQUAL(tm_xy, audi::abs(tm_xy));
    }
}
