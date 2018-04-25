#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(substitution)
{
    gdual_v x(std::vector<double>{1., 1.}, "x", 1);
    gdual_v y(std::vector<double>{-1., -1}, "y", 1);
    auto res = x * y * x / (x - y); // [-0.75, -0.75]*dx+[-0.5, -0.5]+[0.25, 0.25]*dy
    auto res2 = res.subs("dx", std::vector<double>{1.});
    auto res3 = res.subs("dy", std::vector<double>{1.});
    BOOST_CHECK_EQUAL(res2.constant_cf()[0], -1.25);
    BOOST_CHECK_EQUAL(res2.constant_cf()[1], -1.25);
    BOOST_CHECK_EQUAL(res2.get_derivative({0, 1})[0], 0.25);
    BOOST_CHECK_EQUAL(res2.get_derivative({0, 1})[1], 0.25);
    BOOST_CHECK_EQUAL(res2.get_derivative({1, 0})[0], 0.);
    BOOST_CHECK_EQUAL(res2.get_derivative({1, 0}).size(), 1);
    BOOST_CHECK_EQUAL(res3.constant_cf()[0], -0.25);
    BOOST_CHECK_EQUAL(res3.constant_cf()[1], -0.25);
    BOOST_CHECK_EQUAL(res3.get_derivative({1, 0})[0], -0.75);
    BOOST_CHECK_EQUAL(res3.get_derivative({1, 0})[1], -0.75);
    BOOST_CHECK_EQUAL(res3.get_derivative({0, 1})[0], 0.);
    BOOST_CHECK_EQUAL(res3.get_derivative({0, 1}).size(), 1);
}

BOOST_AUTO_TEST_CASE(is_zero)
{
    // We test some trivial cases where truncation order does not influence the results
    {
        gdual_v x(std::vector<double>{{1, 2, 3, 4, 0.123, -21.211}}, "x", 4);
        gdual_v y(std::vector<double>{{0.123, 1.2, 4.3, 2.4, 0.23, -1.211}}, "y", 4);
        gdual_v z(std::vector<double>{{-0.2, -2.01, 0.123, -0.132, 1.123, -0.211}}, "z", 4);
        gdual_v f = x * x * x + x * y * z + z * x * y;

        BOOST_CHECK((f - f).is_zero(1e-12));
        BOOST_CHECK((f - 1 / (1 / f)).is_zero(1e-12));
        BOOST_CHECK(((f * f) / (f)-f).is_zero(1e-12));
        BOOST_CHECK(!f.is_zero(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(extract_order)
{
    unsigned int order = 8u;
    gdual_v x({0.123, 0.222}, "x", order);
    gdual_v y({0.456, -0.12}, "y", order);
    auto f = audi::sin(x * y);
    // We test that the extracted gduals have the requested order
    for (auto i = 0u; i <= order; ++i) {
        auto fi = f.extract_terms(i);
        BOOST_CHECK_EQUAL(fi.degree(), fi.get_order());
        BOOST_CHECK_EQUAL(fi.get_order(), i);
    }
    // We test that f = f0+f1+f2+...+fn
    std::vector<gdual_v> terms;
    for (auto i = 0u; i <= order; ++i) {
        auto f2 = f.extract_terms(i);
        terms.push_back(f2);
    }
    auto sum = std::accumulate(terms.begin(), terms.end(), gdual_v({0.,0.}));
    BOOST_CHECK((sum - f).is_zero(0.));
    // And we test the throw
    BOOST_CHECK_THROW(f.extract_terms(order + 1), std::invalid_argument);
    BOOST_CHECK_NO_THROW(f.extract_terms(order));
}

BOOST_AUTO_TEST_CASE(trim)
{
    gdual_v x({10., -1e-4, 10., -10.}, "x", 4);
    gdual_v y({1e-3, -1e-4, 1e-9, 1e-5}, "x", 4);
    BOOST_CHECK(x.trim(1e-1) == x);
    BOOST_CHECK(y.trim(1e-2) == gdual_v({0.,0.,0.,0.}, "x", 4));
    BOOST_CHECK_THROW(x.trim(-1e-3), std::invalid_argument);
}