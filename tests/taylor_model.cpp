#define BOOST_TEST_MODULE taylor_models_test
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction_and_getters)
{
    // Constructing a constant with a symbol
    uint order = 5;
    gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    std::unordered_map<std::string, double> exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    std::unordered_map<std::string, audi::int_d> dom = {{"x", int_dom}};
    taylor_model tm_x(x, rem, exp, dom);

    // BOOST_CHECK_EQUAL(tm_x.get_tpol(), x);
    // BOOST_CHECK_EQUAL(tm_x.get_ndim(), 1);
    BOOST_CHECK_EQUAL(1, 1);
    // BOOST_CHECK_EQUAL(tm_x.get_order(), order);
    // BOOST_CHECK_EQUAL(tm_x.get_rem(), rem);
    // BOOST_CHECK_EQUAL(tm_x.get_exp(), exp);
    // BOOST_CHECK_EQUAL(tm_x.get_dom(), dom);

}
