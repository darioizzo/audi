#define BOOST_TEST_MODULE audi_gdualld_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>

using namespace audi;
/*
BOOST_AUTO_TEST_CASE(construction)
{
    mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
    gdual_mp x{onethird, "x", 10};
    BOOST_CHECK(abs(onethird - x.constant_cf()) == 0.);

}

BOOST_AUTO_TEST_CASE(exp_and_log)
{
    mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
    gdual_mp x{onethird, "x", 3};
    auto res = exp(x);
    BOOST_CHECK(abs(res.find_cf({0}) - audi::exp(onethird)) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({1}) - audi::exp(onethird)) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({2}) - audi::exp(onethird) / 2.) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({3}) - audi::exp(onethird) / 6.) < 1e-33);
    BOOST_CHECK((exp(log(x)) - x).trim(1e-34) == gdual_mp(0.));
    BOOST_CHECK((log(exp(x)) - x).trim(1e-34) == gdual_mp(0.));
}
*/
BOOST_AUTO_TEST_CASE(exponentiation)
{
//    // overload arith^gdual
//    {
//        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
//        gdual_mp x{onethird, "x", 3};
//        auto res = pow(mppp::real128(2.), x);
//        auto log_2 = audi::log(mppp::real128(2.));
//        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(2., onethird)) < 1e-33);
//        BOOST_CHECK(abs(res.find_cf({1}) - audi::pow(2., onethird) * log_2) < 1e-33);
//        BOOST_CHECK(abs(res.find_cf({2}) - audi::pow(2., onethird) * log_2 * log_2 / 2.) < 1e-33);
//        BOOST_CHECK(abs(res.find_cf({3}) - audi::pow(2., onethird) * log_2 * log_2 * log_2 / 6.) < 1e-33);
//
//        x = gdual_mp(4); // exponent is an integer constant
//        res = pow(mppp::real128(2.), x);
//        BOOST_CHECK(abs(res.constant_cf() - 2*2*2*2) < 1e-33);
//    }
//    // overload gdual^arith
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, onethird);
        //BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(onethird, onethird)) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({1}) - onethird * audi::pow(onethird, onethird - 1.)) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({2}) - onethird * (onethird - 1.) * audi::pow(onethird, onethird - 2.) / 2.) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({3}) - onethird * (onethird - 1.) * (onethird - 2.) * audi::pow(onethird, onethird - 3.) / 6.) < 1e-33);

print(res, "\n");

        // exponent is an integer
        //res = pow(x, 3);
        //BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(onethird, 3)) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({1}) - 3 * audi::pow(3, 3 - 1)) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({2}) - 3 * (3 - 1.) * audi::pow(3, 3 - 2)) < 1e-33);
        //BOOST_CHECK(abs(res.find_cf({3}) - 3 * (3 - 1.) * (3 - 2.) * audi::pow(3, 3 - 3)) < 1e-33);
    }

}
/*
BOOST_AUTO_TEST_CASE(trigonometric)
{
    mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
    gdual_mp x{onethird, "x", 3};
    auto res = sin(x);
    BOOST_CHECK(abs(res.find_cf({0}) - audi::sin(onethird)) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({1}) - audi::cos(onethird)) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({2}) + audi::sin(onethird) / 2.) < 1e-33);
    BOOST_CHECK(abs(res.find_cf({3}) + audi::cos(onethird) / 6.) < 1e-33);
    //BOOST_CHECK((sin(asin(x)) - x).trim(1e-34) == gdual_mp(0.));
}
*/