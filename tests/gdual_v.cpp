#define BOOST_TEST_MODULE audi_gdualv_test

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include <audi/config.hpp>

#if defined(AUDI_WITH_MPPP)
#include <audi/real128.hpp>
#endif

#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/vectorized.hpp>

using namespace audi;

using gdual_v = audi::gdual<audi::vectorized<double>>;


BOOST_AUTO_TEST_CASE(construction)
{
    // Constructing a vectorized constant
    { // 1 - from initializer_list
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{}), std::invalid_argument);
        BOOST_CHECK_EQUAL(gdual_v{}, gdual_v(0.));
        gdual_v x{0.};
        gdual_v y(std::vector<double>{1.2, 2.2});
        BOOST_CHECK_EQUAL(x.get_order(), 0);
        BOOST_CHECK_EQUAL(y.get_order(), 0);
        auto x0 = x.constant_cf();
        auto y0 = y.constant_cf();
        BOOST_CHECK_EQUAL(x0[0], 0.);
        BOOST_CHECK_EQUAL(y0[0], 1.2);
        BOOST_CHECK_EQUAL(y0[1], 2.2);
        BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 0u);
        BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 0u);
    }
    { // 2 - from an std::vector
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{}), std::invalid_argument);
        gdual_v x(std::vector<double>(1, 0.));
        gdual_v y(std::vector<double>{1.2, 2.2});
        BOOST_CHECK_EQUAL(x.get_order(), 0);
        BOOST_CHECK_EQUAL(y.get_order(), 0);
        auto x0 = x.constant_cf();
        auto y0 = y.constant_cf();
        BOOST_CHECK_EQUAL(x0[0], 0.);
        BOOST_CHECK_EQUAL(y0[0], 1.2);
        BOOST_CHECK_EQUAL(y0[1], 2.2);
        BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 0u);
        BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 0u);
    }
    // Constructing a "full" vectorized gdual_v
    { // 1 - from initializer_list
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{}, "x", 3), std::invalid_argument);
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{2, 3}, "dx", 3), std::invalid_argument);
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{2, 3}, "x", std::numeric_limits<unsigned int>::max() - 8u),
                          std::invalid_argument);
        gdual_v x(std::vector<double>{0.}, "x", 3);
        gdual_v y(std::vector<double>{1.2, 2.2}, "y", 3);
        BOOST_CHECK_EQUAL(x.get_order(), 3);
        BOOST_CHECK_EQUAL(y.get_order(), 3);
        auto x0 = x.constant_cf();
        auto y0 = y.constant_cf();
        BOOST_CHECK_EQUAL(x0[0], 0.);
        BOOST_CHECK_EQUAL(y0[0], 1.2);
        BOOST_CHECK_EQUAL(y0[1], 2.2);
        BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 1u);
        BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 1u);
        auto x1 = x.get_derivative({1});
        auto y1 = y.get_derivative({1});
        BOOST_CHECK_EQUAL(x1[0], 1.);
        BOOST_CHECK_EQUAL(y1[0], 1.);
        BOOST_CHECK_EQUAL(y1[1], 1.);
    }
    { // 2 - from an std::vector
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{}, "x", 3), std::invalid_argument);
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{2, 3}, "dx", 3), std::invalid_argument);
        BOOST_CHECK_THROW(gdual_v(std::vector<double>{2, 3}, "x", std::numeric_limits<unsigned int>::max() - 8u),
                          std::invalid_argument);
        gdual_v x(std::vector<double>{0.}, "x", 3);
        gdual_v y(std::vector<double>{1.2, 2.2}, "y", 3);
        BOOST_CHECK_EQUAL(x.get_order(), 3);
        BOOST_CHECK_EQUAL(y.get_order(), 3);
        auto x0 = x.constant_cf();
        auto y0 = y.constant_cf();
        BOOST_CHECK_EQUAL(x0[0], 0.);
        BOOST_CHECK_EQUAL(y0[0], 1.2);
        BOOST_CHECK_EQUAL(y0[1], 2.2);
        BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 1u);
        BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 1u);
        auto x1 = x.get_derivative({1});
        auto y1 = y.get_derivative({1});
        BOOST_CHECK_EQUAL(x1[0], 1.);
        BOOST_CHECK_EQUAL(y1[0], 1.);
        BOOST_CHECK_EQUAL(y1[1], 1.);
    }
}
BOOST_AUTO_TEST_CASE(arithmetic_plus)
{
    { // nominal case
        gdual_v x(std::vector<double>{1., 1.}, "x", 3);
        gdual_v y(std::vector<double>{-1., 1}, "y", 3);
        gdual_v z(std::vector<double>{1.2, 1.3, 1.4}, "z", 3);
        BOOST_CHECK_THROW(x + z, std::invalid_argument);
        auto sum = x + y;
        BOOST_CHECK_EQUAL(sum.get_symbol_set_size(), 2u);
        auto sum0 = sum.constant_cf();
        auto sumdx = sum.get_derivative({1, 0});
        auto sumdy = sum.get_derivative({0, 1});
        BOOST_CHECK_EQUAL(sum0[0], 0.);
        BOOST_CHECK_EQUAL(sum0[1], 2.);
        BOOST_CHECK_EQUAL(sumdx[0], 1.);
        BOOST_CHECK_EQUAL(sumdx[1], 1.);
        BOOST_CHECK_EQUAL(sumdy[0], 1.);
        BOOST_CHECK_EQUAL(sumdy[1], 1.);
    }
    { // scalar case
        auto N = 10u;
        gdual_v x(std::vector<double>(N, 123.), "x", 3);
        auto sum = x + 5.;
        auto sum2 = x + gdual_v(std::vector<double>{5.});
        auto sum3 = 5. + x;
        BOOST_CHECK_EQUAL(sum, sum2);
        BOOST_CHECK_EQUAL(sum2, sum3);
        for (auto i = 0u; i < N; ++i) {
            BOOST_CHECK_EQUAL(sum.constant_cf()[i], 128.);
            BOOST_CHECK_EQUAL(sum.get_derivative({1})[i], 1.);
        }
    }
}
BOOST_AUTO_TEST_CASE(arithmetic_minus)
{
    {
        gdual_v x(std::vector<double>{1., 1.}, "x", 3);
        gdual_v y(std::vector<double>{-1., 1}, "y", 3);
        gdual_v z(std::vector<double>{1.2, 1.3, 1.4}, "z", 3);
        BOOST_CHECK_THROW(x - z, std::invalid_argument);
        auto diff = x - y;
        BOOST_CHECK_EQUAL(diff.get_symbol_set_size(), 2u);
        auto diff0 = diff.constant_cf();
        auto diffdx = diff.get_derivative({1, 0});
        auto diffdy = diff.get_derivative({0, 1});
        BOOST_CHECK_EQUAL(diff0[0], 2.);
        BOOST_CHECK_EQUAL(diff0[1], 0.);
        BOOST_CHECK_EQUAL(diffdx[0], 1.);
        BOOST_CHECK_EQUAL(diffdx[1], 1.);
        BOOST_CHECK_EQUAL(diffdy[0], -1.);
        BOOST_CHECK_EQUAL(diffdy[1], -1.);
    }
    { // scalar case
        auto N = 10u;
        gdual_v x(std::vector<double>(N, 123.), "x", 3);
        auto diff = x - 5.;
        auto diff2 = x - gdual_v(std::vector<double>{5.});
        auto diff3 = 5. - x;
        BOOST_CHECK_EQUAL(diff, diff2);
        BOOST_CHECK_EQUAL(diff2, -diff3);
        for (auto i = 0u; i < N; ++i) {
            BOOST_CHECK_EQUAL(diff.constant_cf()[i], 118.);
            BOOST_CHECK_EQUAL(diff.get_derivative({1})[i], 1.);
        }
    }
}
BOOST_AUTO_TEST_CASE(arithmetic_mul)
{
    {
        gdual_v x(std::vector<double>{1., 1.}, "x", 3);
        gdual_v y(std::vector<double>{-1., 1}, "y", 3);
        gdual_v z(std::vector<double>{1.2, 1.3, 1.4}, "z", 3);
        BOOST_CHECK_THROW(x * z, std::invalid_argument);
        auto mul = x * y;
        BOOST_CHECK_EQUAL(mul.get_symbol_set_size(), 2u);
        auto mul0 = mul.constant_cf();
        auto muldx = mul.get_derivative({1, 0});
        auto muldy = mul.get_derivative({0, 1});
        BOOST_CHECK_EQUAL(mul0[0], -1.);
        BOOST_CHECK_EQUAL(mul0[1], 1.);
        BOOST_CHECK_EQUAL(muldx[0], -1.);
        BOOST_CHECK_EQUAL(muldx[1], 1.);
        BOOST_CHECK_EQUAL(muldy[0], 1.);
        BOOST_CHECK_EQUAL(muldy[1], 1.);
    }
    { // scalar case
        auto N = 10u;
        gdual_v x(std::vector<double>(N, 123.), "x", 3);
        auto mul = x * 5.;
        auto mul2 = x * gdual_v(std::vector<double>{5.});
        auto mul3 = 5. * x;
        BOOST_CHECK_EQUAL(mul, mul2);
        BOOST_CHECK_EQUAL(mul3, mul2);
        for (auto i = 0u; i < N; ++i) {
            BOOST_CHECK_EQUAL(mul.constant_cf()[i], 123. * 5.);
            BOOST_CHECK_EQUAL(mul.get_derivative({1})[i], 5.);
        }
    }
}
BOOST_AUTO_TEST_CASE(arithmetic_div)
{
    {
        gdual_v x(std::vector<double>{1., 1.}, "x", 3);
        gdual_v y(std::vector<double>{-1., 1}, "y", 3);
        gdual_v z(std::vector<double>{1.2, 1.3, 1.4}, "z", 3);
        BOOST_CHECK_THROW(x / z, std::invalid_argument);
        auto div = x / y;
        BOOST_CHECK_EQUAL(div.get_symbol_set_size(), 2u);
        auto div0 = div.constant_cf();
        auto divdx = div.get_derivative({1, 0});
        auto divdy = div.get_derivative({0, 1});
        BOOST_CHECK_EQUAL(div0[0], -1.);
        BOOST_CHECK_EQUAL(div0[1], 1.);
        BOOST_CHECK_EQUAL(divdx[0], -1.);
        BOOST_CHECK_EQUAL(divdx[1], 1.);
        BOOST_CHECK_EQUAL(divdy[0], -1.);
        BOOST_CHECK_EQUAL(divdy[1], -1.);
    }
    { // scalar case
        auto N = 10u;
        gdual_v x(std::vector<double>(N, 50.), "x", 3);
        auto div = x / 5.;
        auto div2 = x / gdual_v(std::vector<double>{5.});
        auto div3 = 5. / x;
        BOOST_CHECK_EQUAL(div, div2);
        BOOST_CHECK_EQUAL(div2, 1. / div3); // this will not work if floats are not representable
        for (auto i = 0u; i < N; ++i) {
            BOOST_CHECK_EQUAL(div.constant_cf()[i], 50. / 5.);
            BOOST_CHECK_EQUAL(div.get_derivative({1})[i], 1 / 5.);
        }
    }
}

BOOST_AUTO_TEST_CASE(comparisons)
{
    gdual_v x(std::vector<double>{1.1, 1.3}, "x", 4);
    gdual_v y(std::vector<double>{3.1, 3.2}, "y", 4);
    BOOST_CHECK(x < y);
    BOOST_CHECK(x * x < y * y);
    BOOST_CHECK(audi::sin(x) < y * y);
    BOOST_CHECK(audi::sin(x) < exp(y));
    BOOST_CHECK(y > x);
    BOOST_CHECK(y * y > x * x);
    BOOST_CHECK(y * y > audi::sin(x));
    BOOST_CHECK(exp(y) > audi::sin(x));
    BOOST_CHECK(!(y>y));
    BOOST_CHECK(!(y<y));
}
