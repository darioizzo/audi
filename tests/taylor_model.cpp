#define BOOST_TEST_MODULE taylor_models_test

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model.hpp>

BOOST_AUTO_TEST_CASE(comparison_methods)
{
    audi::var_map_d exp1 = {{"x", 0.0}, {"y", 1.0}};
    audi::var_map_d exp2 = {{"x", 0.0}};

    BOOST_CHECK(!audi::taylor_model::map_equal(exp1, exp2));

    audi::var_map_d exp_2_1 = {{"x", 0.0}, {"y", 1.0}};
    audi::var_map_d exp_2_2 = {{"x", 0.0}, {"y", 1.0}};

    BOOST_CHECK_EQUAL(audi::taylor_model::map_equal(exp_2_1, exp_2_2), true);

    audi::var_map_i dom1 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};
    audi::var_map_i dom2 = {{"x", int_d(2.0, 3.0)}};

    BOOST_CHECK_EQUAL(audi::taylor_model::map_interval_equal(dom1, dom2), false);

    audi::var_map_i dom_2_1 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};
    audi::var_map_i dom_2_2 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};

    BOOST_CHECK_EQUAL(audi::taylor_model::map_interval_equal(dom_2_1, dom_2_2), true);
}

BOOST_AUTO_TEST_CASE(validity_check)
{
    // Constructing a Taylor model
    uint order = 5;
    audi::gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    audi::var_map_i dom = {{"y", int_dom}};
    BOOST_CHECK_THROW(audi::taylor_model tm_x(x, rem, exp, dom), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(construction_and_getters_univariate)
{
    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}};
    audi::var_map_i dom = {{"x", int_d(0.0, 1.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::taylor_model tm_x(x, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_x.get_tpol(), x);
    BOOST_CHECK_EQUAL(tm_x.get_ndim(), 1);
    BOOST_CHECK_EQUAL(tm_x.get_order(), order);
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(tm_x.get_rem(), rem), true);

    // Test the expansion points unordered map
    BOOST_CHECK_EQUAL(tm_x.get_exp().size(), exp.size());
    BOOST_CHECK_EQUAL(audi::taylor_model::map_equal(tm_x.get_exp(), exp), true);
    BOOST_CHECK_EQUAL(audi::taylor_model::map_interval_equal(tm_x.get_dom(), dom), true);
}

BOOST_AUTO_TEST_CASE(construction_and_getters_multivariate)
{
    // Constructing a multivariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}, {"y", 1.0}};
    audi::var_map_i dom = {{"x", int_d(0.0, 1.0)}, {"y", int_d(1.0, 2.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d y(exp.find("y")->second, "y", order);

    audi::gdual_d f_xy = x * x + y + 20;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_xy.get_tpol(), f_xy);
    BOOST_CHECK_EQUAL(tm_xy.get_ndim(), 2);
    BOOST_CHECK_EQUAL(tm_xy.get_order(), order);
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(tm_xy.get_rem(), rem), true);

    // Test the expansion points unordered map
    BOOST_CHECK_EQUAL(tm_xy.get_exp().size(), exp.size());
    BOOST_CHECK_EQUAL(audi::taylor_model::map_equal(tm_xy.get_exp(), exp), true);
    BOOST_CHECK_EQUAL(audi::taylor_model::map_interval_equal(tm_xy.get_dom(), dom), true);
}

BOOST_AUTO_TEST_CASE(construction_and_getters_identity)
{
    audi::taylor_model tm_id = audi::taylor_model::identity();
    BOOST_CHECK_EQUAL(tm_id.get_tpol(), audi::gdual<double>(1.0));
    BOOST_CHECK_EQUAL(tm_id.get_ndim(), 0);
    BOOST_CHECK_EQUAL(tm_id.get_order(), 0);
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(tm_id.get_rem(), audi::int_d(0.0)), true);
}

BOOST_AUTO_TEST_CASE(comparison_of_construction_order)
{
    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}};
    audi::var_map_i dom = {{"x", int_d(0.0, 1.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::taylor_model tm_x(x, rem, exp, dom);
    int gx = 2;
    auto prod_ans = tm_x + gx;

    auto exp_x = x + 2;
    audi::taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans.get_tpol(), exp_ans.get_tpol());
    BOOST_CHECK_EQUAL(audi::taylor_model::map_equal(prod_ans.get_exp(), exp_ans.get_exp()), true);
    BOOST_CHECK_EQUAL(audi::taylor_model::map_interval_equal(prod_ans.get_dom(), exp_ans.get_dom()), true);
    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(increasing_order)
{
    uint order = 5;
    audi::gdual_d x(0.0, "x", order);
    audi::gdual_d y = audi::taylor_model::get_increased_order(x, 10);
    BOOST_CHECK_EQUAL(y.get_order(), 10);
}

BOOST_AUTO_TEST_CASE(flatten)
{
    std::vector<std::vector<double>> vec = {{1, 2, 3}, {4, 5, 6}};
    std::vector<double> exp_vec = {1, 2, 3, 4, 5, 6};

    std::vector<double> prod_vec = audi::taylor_model::flatten(vec);
    for (int i = 0; i < static_cast<int>(exp_vec.size()); ++i) {
        BOOST_CHECK_EQUAL(prod_vec[i], exp_vec[i]);
    }
}

BOOST_AUTO_TEST_CASE(get_bounds)
{
    uint order = 5;
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}, {"y", 1.0}};
    audi::var_map_i dom = {{"x", int_d(0.0, 1.0)}, {"y", int_d(1.0, 2.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d y(exp.find("y")->second, "y", order);

    audi::gdual_d f_xy = x * x + y + 20;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    // Case verified with Python implementation, in turn verified by Makino (1998)
    BOOST_CHECK_EQUAL(audi::taylor_model::interval_equal(tm_xy.get_bounds(), int_d(21.0, 23.0)), true);
}

BOOST_AUTO_TEST_CASE(addition)
{
    uint order = 5;
    audi::gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    audi::var_map_i dom = {{"x", int_dom}};
    audi::taylor_model tm_fx(x, rem, exp, dom);

    int gx = 2;

    auto prod_ans = tm_fx + gx;

    auto exp_x = x + 2;
    audi::taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
    uint order = 5;
    audi::gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    audi::var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    audi::var_map_i dom = {{"x", int_dom}};
    audi::taylor_model tm_fx(x, rem, exp, dom);

    int gx = -5;

    auto prod_ans = tm_fx + gx;

    auto exp_x = x - 5;
    audi::taylor_model exp_ans(exp_x, rem, exp, dom);

    BOOST_CHECK_EQUAL(prod_ans, exp_ans);
}

BOOST_AUTO_TEST_CASE(test_multiplication)
{
    // ---- Setup expansion point and domain ----
    var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    var_map_i dom = {{"x", int_dom}};

    // Initialize gdual (centered at 0.0, variable "x", order 5)
    audi::gdual_d x(0.0, "x", 5);
    audi::gdual_d f_x = x;

    // Initialize TaylorModel over domain [0,2], expansion point 0
    int_d tpol_range(0.0, 2.0);
    audi::taylor_model tm_fx(f_x, tpol_range, exp, dom);

    // ---- Test multiplication by int ----
    int g_x = 2;
    audi::taylor_model prod_ans = tm_fx * g_x;
    audi::taylor_model exp_ans(f_x * 2.0, tpol_range * 2.0, exp, dom);
    BOOST_CHECK(prod_ans == exp_ans);

    // Right-hand side multiplication
    prod_ans = g_x * tm_fx;
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by double ----
    double h_x = 2.0;
    prod_ans = tm_fx * h_x;
    exp_ans = audi::taylor_model(f_x * 2.0, tpol_range * 2.0, exp, dom);
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by interval ----
    int_d i_x(0.0, 2.0);
    prod_ans = tm_fx * i_x;
    exp_ans = audi::taylor_model(f_x, tpol_range * 4.0, exp, dom); // [0,8] interval
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test multiplication by another TaylorModel ----
    prod_ans = tm_fx * tm_fx;
    exp_ans = audi::taylor_model(f_x * f_x, int_d(0.0, 8.0), exp, dom); // [0,8]
    BOOST_CHECK(prod_ans == exp_ans);
}

BOOST_AUTO_TEST_CASE(test_power)
{

    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 0.0);
    audi::var_map_d exp = {{"x", 0.0}};
    audi::var_map_i dom = {{"x", int_d(0.0, 1.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d fx = audi::pow(x, 3);
    audi::taylor_model tm_fx_exp(fx, rem, exp, dom);
    audi::taylor_model tm_x(x, rem, exp, dom);

    audi::taylor_model tm_fx = tm_x.pow(3);
    BOOST_CHECK(tm_fx == tm_fx_exp);
}
