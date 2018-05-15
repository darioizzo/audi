#define BOOST_TEST_MODULE audi_gdualmp_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <audi/audi.hpp>
#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
    gdual_mp x{onethird, "x", 10};
    BOOST_CHECK(audi::abs(onethird - x.constant_cf()) == 0.);
}

BOOST_AUTO_TEST_CASE(exp_and_log)
{
    mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
    gdual_mp x{onethird, "x", 3};
    auto res = exp(x);
    BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::exp(onethird)) < 1e-33);
    BOOST_CHECK(audi::abs(res.find_cf({1}) - audi::exp(onethird)) < 1e-33);
    BOOST_CHECK(audi::abs(res.find_cf({2}) - audi::exp(onethird) / 2.) < 1e-33);
    BOOST_CHECK(audi::abs(res.find_cf({3}) - audi::exp(onethird) / 6.) < 1e-33);
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
        BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::pow(2., onethird)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - audi::pow(2., onethird) * log_2) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - audi::pow(2., onethird) * log_2 * log_2 / 2.) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - audi::pow(2., onethird) * log_2 * log_2 * log_2 / 6.) < 1e-33);

        x = gdual_mp(4); // exponent is an integer constant
        res = pow(mppp::real128(2.), x);
        BOOST_CHECK(audi::abs(res.constant_cf() - 2 * 2 * 2 * 2) < 1e-33);
    }
    // overload gdual^arith
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, onethird);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::pow(onethird, onethird)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - onethird * audi::pow(onethird, onethird - 1.)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - onethird * (onethird - 1.) * audi::pow(onethird, onethird - 2.) / 2.)
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3})
                              - onethird * (onethird - 1.) * (onethird - 2.) * audi::pow(onethird, onethird - 3.) / 6.)
                    < 1e-33);

        // exponent is a double representable as integer
        x = gdual_mp(3, "x", 3);
        res = pow(x, 3.); // note that the constant 3 is double
        BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::pow(3, 3)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - 3 * audi::pow(3, 3 - 1)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - 3 * (3 - 1.) * audi::pow(3, 3 - 2) / 2.) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - 3 * (3 - 1.) * (3 - 2.) * audi::pow(3, 3 - 3) / 6.) < 1e-33);
    }
    // overload gdual^int
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, 3);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::pow(onethird, 3)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - 3 * audi::pow(onethird, 3 - 1)) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - 3 * (3 - 1.) * audi::pow(onethird, 3 - 2) / 2.) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - 3 * (3 - 1.) * (3 - 2.) * audi::pow(onethird, 3 - 3) / 6.) < 1e-33);
        res = pow(x, -3);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - audi::pow(onethird, -3))
                    < 1e-32); // we check for absolute error, hence the precision loss
        BOOST_CHECK(audi::abs(res.find_cf({1}) + 3 * audi::pow(onethird, -3 - 1))
                    < 1e-31); // we check for absolute error, hence the precision loss
        BOOST_CHECK(audi::abs(res.find_cf({2}) + 3 * (-3 - 1) * audi::pow(onethird, -3 - 2) / 2)
                    < 1e-30); // we check for absolute error, hence the precision loss
        BOOST_CHECK(audi::abs(res.find_cf({3}) + 3 * (-3 - 1) * (-3 - 2) * audi::pow(onethird, -3 - 3) / 6)
                    < 1e-29); // we check for absolute error, hence the precision loss

        res = pow(x, 0);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - 1 == 0));
        BOOST_CHECK(audi::abs(res.find_cf({1}) == 0));
        BOOST_CHECK(audi::abs(res.find_cf({2}) == 0));
        BOOST_CHECK(audi::abs(res.find_cf({3}) == 0));
    }
    // overload gdual^gdual
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = pow(x, x);
        BOOST_CHECK(audi::abs(res.find_cf({0})
                              - mppp::real128("0.693361274350634704843352274785961795445935113457754036565863693400"))
                    < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({1})
                      - mppp::real128("-0.0683739421375531895286241688670376535021100832441426057588451227360656"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({2})
                      - mppp::real128("1.043413166985674569326084555526740070128227383299088804220207359748035"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({3})
                      - mppp::real128("-1.142713640471137873959384027190746725008973893684368645673060217249986"))
            < 1e-33);
    }
}
BOOST_AUTO_TEST_CASE(sqrt_fun)
{
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sqrt(x);
        BOOST_CHECK(audi::abs(res.find_cf({0})
                              - mppp::real128("0.577350269189625764509148780501957455647601751270126876018602326484"))
                    < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({1})
                      - mppp::real128("0.8660254037844386467637231707529361834714026269051903140279034897263995"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({2})
                      - mppp::real128("-0.6495190528383289850727923780647021376035519701788927355209276172954492"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({3})
                      - mppp::real128("0.974278579257493477609188567097053206405327955268339103281391425944148"))
            < 1e-33);
    }
}
BOOST_AUTO_TEST_CASE(cbrt_fun)
{
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = cbrt(x);
        BOOST_CHECK(audi::abs(res.find_cf({0})
                              - mppp::real128("0.693361274350634704843352274785961795445935113457754036565863693400"))
                    < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({1})
                      - mppp::real128("0.6933612743506347048433522747859617954459351134577540365658636934004977"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({2})
                      - mppp::real128("-0.693361274350634704843352274785961795445935113457754036565863693401191"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({3})
                      - mppp::real128("1.155602123917724508072253791309936325743225189096256727609772822336474"))
            < 1e-33);
    }
}
BOOST_AUTO_TEST_CASE(trigonometric)
{
    // sine
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sin(x);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - mppp::real128("0.3271946967961522441733440852676206060643014"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - mppp::real128("0.94495694631473766438828400767588060784585271047"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - mppp::real128("-0.16359734839807612208667204263381030303215068770"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - mppp::real128("-0.157492824385789610731380667945980101307642118412"))
                    < 1e-33);
    }
    // cosine
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = cos(x);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - mppp::real128("0.9449569463147376643882840076758806078458527"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - mppp::real128("-0.32719469679615224417334408526762060606430137540"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - mppp::real128("-0.472478473157368832194142003837940303922926355236"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - mppp::real128("0.054532449466025374028890680877936767677383562566"))
                    < 1e-33);
    }
    // sincos
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sin_and_cos(x);
        BOOST_CHECK(res[0] == sin(x));
        BOOST_CHECK(res[1] == cos(x));
    }
    // tan
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = tan(x);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - mppp::real128("0.3462535495105754910385435656097407745957039"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - mppp::real128("1.11989152054867255286971414281622406522119684361"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - mppp::real128("0.38776641405677346168406834694334831426784259534"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - mppp::real128("0.50756267076436951015232120845633638256350224970"))
                    < 1e-33);
    }
    // inverse functions
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        BOOST_CHECK((sin(asin(x)) - x).trim(1e-34) == gdual_mp(0.));
        BOOST_CHECK((asin(sin(x)) - x).trim(1e-34) == gdual_mp(0.));
        BOOST_CHECK((cos(acos(x)) - x).trim(1e-32)
                    == gdual_mp(0.)); // in this computation we have an accumulating precision loss at high orders
        BOOST_CHECK((acos(cos(x)) - x).trim(1e-32)
                    == gdual_mp(0.)); // in this computation we have an accumulating precision loss at high orders
        BOOST_CHECK((tan(atan(x)) - x).trim(1e-34) == gdual_mp(0.));
        BOOST_CHECK((atan(tan(x)) - x).trim(1e-33) == gdual_mp(0.));
    }
}
BOOST_AUTO_TEST_CASE(hyperbolic)
{
    // hyperbolic sine
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sinh(x);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - mppp::real128("0.339540557256150139101260611338603585072")) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - mppp::real128("1.0560718678299393895268647082639832525253961"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - mppp::real128("0.1697702786280750695506303056693017925360226"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - mppp::real128("0.17601197797165656492114411804399720875423268"))
                    < 1e-33);
    }
    // hyperbolic  cosine
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = cosh(x);
        BOOST_CHECK(audi::abs(res.find_cf({0}) - mppp::real128("1.056071867829939389526864708263983252525")) < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1}) - mppp::real128("0.3395405572561501391012606113386035850720452"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2}) - mppp::real128("0.52803593391496969476343235413199162626269805"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3}) - mppp::real128("0.05659009287602502318354343522310059751200753"))
                    < 1e-33);
    }
    // hyperbolic sine and cosine
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = sinh_and_cosh(x);
        BOOST_CHECK(res[0] == sinh(x));
        BOOST_CHECK(res[1] == cosh(x));
    }
    // hyperbolic tan
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = tanh(x);
        BOOST_CHECK(audi::abs(res.find_cf({0})
                              - mppp::real128("0.3215127375316343447194062224252064660052920025020835116798050"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({1})
                              - mppp::real128("0.89662955960491440420948934034124973587445145348001949661874818114"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({2})
                              - mppp::real128("-0.28827782426035973712472569203365039113005161345723090151632152611"))
                    < 1e-33);
        BOOST_CHECK(audi::abs(res.find_cf({3})
                              - mppp::real128("-0.20619152742069314951779955121670628474663040144796586567125202883"))
                    < 1e-33);
    }
    // hyperbolic inverse functions
    {
        mppp::real128 onethird{"1.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        BOOST_CHECK((sinh(asinh(x)) - x).trim(1e-33) == gdual_mp(0.));
        BOOST_CHECK((asinh(sinh(x)) - x).trim(1e-33) == gdual_mp(0.));
        BOOST_CHECK((cosh(acosh(x)) - x).trim(1e-33) == gdual_mp(0.));
        BOOST_CHECK((acosh(cosh(x)) - x).trim(1e-33) == gdual_mp(0.));
        BOOST_CHECK((tanh(atanh(x - 1)) - (x - 1)).trim(1e-33) == gdual_mp(0.));
        BOOST_CHECK((atanh(tanh(x)) - x).trim(1e-33) == gdual_mp(0.));
    }
}

BOOST_AUTO_TEST_CASE(erf_fun)
{
    {
        mppp::real128 onethird{"0.33333333333333333333333333333333333333333333333"};
        gdual_mp x{onethird, "x", 3};
        auto res = erf(x);
        BOOST_CHECK(audi::abs(res.find_cf({0})
                              - mppp::real128("0.362648111766062933408178640147865879692141590372537239240"))
                    < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({1})
                      - mppp::real128("1.00971804299131606624460336666642604033678715687395758256805086446851098"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({2})
                      - mppp::real128("-0.33657268099710535541486778888880868011226238562465252752268362148947000"))
            < 1e-33);
        BOOST_CHECK(
            audi::abs(res.find_cf({3})
                      - mppp::real128("-0.26177875188663749865600828024685119564287074437472974362875392782518447"))
            < 1e-33);
    }
}
