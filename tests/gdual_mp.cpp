#define BOOST_TEST_MODULE audi_gdualld_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>

using namespace audi;

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

BOOST_AUTO_TEST_CASE(exponentiation)
{
    // overload arith^gdual
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(mppp::real128(2.), x);
        auto log_2 = audi::log(mppp::real128(2.));
        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(2., onethird)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - audi::pow(2., onethird) * log_2) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - audi::pow(2., onethird) * log_2 * log_2 / 2.) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - audi::pow(2., onethird) * log_2 * log_2 * log_2 / 6.) < 1e-33);

        x = gdual_mp(4); // exponent is an integer constant
        res = pow(mppp::real128(2.), x);
        BOOST_CHECK(abs(res.constant_cf() - 2 * 2 * 2 * 2) < 1e-33);
    }
    // overload gdual^arith
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, onethird);
        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(onethird, onethird)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - onethird * audi::pow(onethird, onethird - 1.)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - onethird * (onethird - 1.) * audi::pow(onethird, onethird - 2.) / 2.)
                    < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3})
                        - onethird * (onethird - 1.) * (onethird - 2.) * audi::pow(onethird, onethird - 3.) / 6.)
                    < 1e-33);

        // exponent is a double representable as integer
        x = gdual_mp(3, "x", 3);
        res = pow(x, 3.); // note that the constant 3 is double
        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(3, 3)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - 3 * audi::pow(3, 3 - 1)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - 3 * (3 - 1.) * audi::pow(3, 3 - 2) / 2.) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - 3 * (3 - 1.) * (3 - 2.) * audi::pow(3, 3 - 3) / 6.) < 1e-33);
    }
    // overload gdual^int
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, 3);
        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(onethird, 3)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - 3 * audi::pow(onethird, 3 - 1)) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - 3 * (3 - 1.) * audi::pow(onethird, 3 - 2) / 2.) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - 3 * (3 - 1.) * (3 - 2.) * audi::pow(onethird, 3 - 3) / 6.) < 1e-33);
        res = pow(x, -3);
        BOOST_CHECK(abs(res.find_cf({0}) - audi::pow(onethird, -3)) < 1e-32); // we check for absolute error, hence the precision loss
        BOOST_CHECK(abs(res.find_cf({1}) + 3 * audi::pow(onethird, -3 - 1)) < 1e-31); // we check for absolute error, hence the precision loss
        BOOST_CHECK(abs(res.find_cf({2}) + 3 * (-3 - 1) * audi::pow(onethird, -3 - 2) / 2) < 1e-30); // we check for absolute error, hence the precision loss
        BOOST_CHECK(abs(res.find_cf({3}) + 3 * (-3 - 1) * (-3 - 2) * audi::pow(onethird, -3 - 3) / 6) < 1e-29); // we check for absolute error, hence the precision loss

        res = pow(x, 0);
        BOOST_CHECK(abs(res.find_cf({0}) - 1 == 0));
        BOOST_CHECK(abs(res.find_cf({1}) == 0));
        BOOST_CHECK(abs(res.find_cf({2}) == 0));
        BOOST_CHECK(abs(res.find_cf({3}) == 0));
    }
    // overload gdual^gdual
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, x);
        BOOST_CHECK(abs(res.find_cf({0}) - mppp::real128("0.693361274350634704843352274785961795445935113457754036565863693400")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - mppp::real128("-0.0683739421375531895286241688670376535021100832441426057588451227360656")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - mppp::real128("1.043413166985674569326084555526740070128227383299088804220207359748035")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - mppp::real128("-1.142713640471137873959384027190746725008973893684368645673060217249986")) < 1e-33);
    }
}
BOOST_AUTO_TEST_CASE(sqrt_fun)
{
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sqrt(x);
        BOOST_CHECK(abs(res.find_cf({0}) - mppp::real128("0.577350269189625764509148780501957455647601751270126876018602326484")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - mppp::real128("0.8660254037844386467637231707529361834714026269051903140279034897263995")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - mppp::real128("-0.6495190528383289850727923780647021376035519701788927355209276172954492")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - mppp::real128("0.974278579257493477609188567097053206405327955268339103281391425944148")) < 1e-33);
    }
}
BOOST_AUTO_TEST_CASE(cbrt_fun)
{
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = cbrt(x);
        BOOST_CHECK(abs(res.find_cf({0}) - mppp::real128("0.693361274350634704843352274785961795445935113457754036565863693400")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({1}) - mppp::real128("0.6933612743506347048433522747859617954459351134577540365658636934004977")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({2}) - mppp::real128("-0.693361274350634704843352274785961795445935113457754036565863693401191")) < 1e-33);
        BOOST_CHECK(abs(res.find_cf({3}) - mppp::real128("1.155602123917724508072253791309936325743225189096256727609772822336474")) < 1e-33);
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