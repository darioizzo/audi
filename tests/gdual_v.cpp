#include "../src/gdual_v.hpp"
#include "helpers.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Constructing a vectorized constant
    { // 1 - from initializer_list
    gdual_v x({0.});
    gdual_v y({1.2, 2.2});
    BOOST_CHECK_EQUAL(x.get_order(), 0);
    BOOST_CHECK_EQUAL(y.get_order(), 0);
    BOOST_CHECK_THROW(gdual_v({}), std::invalid_argument);
    auto x0 = x.constant_cf();
    auto y0 = y.constant_cf();
    BOOST_CHECK_EQUAL(x0[0], 0.);
    BOOST_CHECK_EQUAL(y0[0], 1.2);
    BOOST_CHECK_EQUAL(y0[1], 2.2);
    BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 0u);
    BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 0u);
    }
    { // 2 - from an std::vector
    gdual_v x(std::vector<double>(1,0.));
    gdual_v y(std::vector<double>{1.2, 2.2});
    BOOST_CHECK_EQUAL(x.get_order(), 0);
    BOOST_CHECK_EQUAL(y.get_order(), 0);
    BOOST_CHECK_THROW(gdual_v({}), std::invalid_argument);
    auto x0 = x.constant_cf();
    auto y0 = y.constant_cf();
    BOOST_CHECK_EQUAL(x0[0], 0.);
    BOOST_CHECK_EQUAL(y0[0], 1.2);
    BOOST_CHECK_EQUAL(y0[1], 2.2);
    BOOST_CHECK_EQUAL(x.get_symbol_set_size(), 0u);
    BOOST_CHECK_EQUAL(y.get_symbol_set_size(), 0u);
    }
    // Constructing a "full" vectorized gdual_v
    {
    gdual_v x({0.}, "x", 3);
    gdual_v y({1.2, 2.2}, "y", 3);
    }

}
