#define BOOST_TEST_MODULE taylor_models_test

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction_and_getters)
{
    // Constructing a Taylor model
    uint order = 5;
    gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    std::unordered_map<std::string, double> exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    std::unordered_map<std::string, audi::int_d> dom = {{"x", int_dom}};
    taylor_model tm_x(x, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_x.get_tpol(), x);
    BOOST_CHECK_EQUAL(tm_x.get_ndim(), 1);
    BOOST_CHECK_EQUAL(tm_x.get_order(), order);
    // BOOST_CHECK_EQUAL(tm_x.get_rem(), rem);
    BOOST_CHECK_EQUAL(boost::numeric::lower(tm_x.get_rem()), boost::numeric::lower(rem));
    BOOST_CHECK_EQUAL(boost::numeric::upper(tm_x.get_rem()), boost::numeric::upper(rem));

    // Test the expansion points unordered map
    BOOST_CHECK_EQUAL(tm_x.get_exp().size(), exp.size());

    // Loop and check each key-value pair
    for (const auto &[key, expected_value] : exp) {
        auto it = tm_x.get_exp().find(key);
        BOOST_CHECK(it != tm_x.get_exp().end()); // Ensure key exists
        if (it != tm_x.get_exp().end()) {
            BOOST_CHECK_EQUAL(it->second, expected_value);
        }
    }

    // Test the domain intervals unordered map
    BOOST_CHECK_EQUAL(tm_x.get_dom().size(), dom.size());

    // Loop and check each key-value pair
    for (const auto &[key, expected_value] : dom) {
        auto it = tm_x.get_dom().find(key);
        BOOST_CHECK(it != tm_x.get_dom().end()); // Ensure key exists
        if (it != tm_x.get_dom().end()) {
            BOOST_CHECK_EQUAL(boost::numeric::lower(it->second), boost::numeric::lower(expected_value));
        }
    }
};

BOOST_AUTO_TEST_CASE(comparison)
{
    uint order = 5;
    gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    std::unordered_map<std::string, double> exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    std::unordered_map<std::string, audi::int_d> dom = {{"x", int_dom}};
    taylor_model tm_fx(x, rem, exp, dom);
    int gx = 2;

    auto prod_ans = tm_fx + gx;

    auto exp_x = x + 2;
    taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans.get_tpol(), exp_ans.get_tpol());
    bool exp_point_ans = prod_ans.map_interval_equal<double>(prod_ans.get_exp(), exp_ans.get_exp());
    BOOST_CHECK_EQUAL(exp_point_ans, true);
    bool dom_ans = prod_ans.map_interval_equal<int_d>(prod_ans.get_dom(), exp_ans.get_dom());
    BOOST_CHECK_EQUAL(dom_ans, true);
    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(addition)
{
    uint order = 5;
    gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    std::unordered_map<std::string, double> exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    std::unordered_map<std::string, audi::int_d> dom = {{"x", int_dom}};
    taylor_model tm_fx(x, rem, exp, dom);

    int gx = 2;

    auto prod_ans = tm_fx + gx;

    auto exp_x = x + 2;
    taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
    uint order = 5;
    gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    std::unordered_map<std::string, double> exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    std::unordered_map<std::string, audi::int_d> dom = {{"x", int_dom}};
    taylor_model tm_fx(x, rem, exp, dom);

    int gx = -5;

    auto prod_ans = tm_fx + gx;

    auto exp_x = x - 5;
    taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(test_multiplication)
{
    // ---- Setup expansion point and domain ----
    var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    var_map_i dom = {{"x", int_dom}};

    // Initialize gdual (centered at 0.0, variable "x", order 5)
    gdual_d x(0.0, "x", 5);
    gdual_d f_x = x;

    // Initialize TaylorModel over domain [0,2], expansion point 0
    int_d tpol_range(0.0, 2.0);
    taylor_model tm_fx(f_x, tpol_range, exp, dom);

    // ---- Test multiplication by int ----
    int g_x = 2;
    taylor_model prod_ans = tm_fx * g_x;
    taylor_model exp_ans(f_x * 2.0, tpol_range * 2.0, exp, dom);
    BOOST_CHECK(prod_ans == exp_ans);

    // Right-hand side multiplication
    prod_ans = g_x * tm_fx;
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by double ----
    double h_x = 2.0;
    prod_ans = tm_fx * h_x;
    exp_ans = taylor_model(f_x * 2.0, tpol_range * 2.0, exp, dom);
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by interval ----
    int_d i_x(0.0, 2.0);
    var_map_i dom2 = {{"x", i_x}};
    prod_ans = tm_fx * i_x;
    exp_ans = taylor_model(f_x, tpol_range * 4.0, exp, dom); // [0,8] interval
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by another TaylorModel ----
    prod_ans = tm_fx * tm_fx;
    exp_ans = taylor_model(f_x * f_x, tpol_range * tpol_range.upper(), exp, dom); // [0,8]
    BOOST_CHECK(prod_ans == exp_ans);
}
