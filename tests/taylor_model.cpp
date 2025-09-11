#include <boost/test/tools/old/interface.hpp>
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

// using int_d = boost::numeric::interval<double>;
using int_d = boost::numeric::interval<
    double, boost::numeric::interval_lib::policies<
                boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double>>,
                boost::numeric::interval_lib::checking_no_nan<double>>>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

BOOST_AUTO_TEST_CASE(comparison_methods)
{
    var_map_d exp1 = {{"x", 0.0}, {"y", 1.0}};
    var_map_d exp2 = {{"x", 0.0}};

    BOOST_CHECK(!audi::taylor_model::map_equal(exp1, exp2));

    var_map_d exp_2_1 = {{"x", 0.0}, {"y", 1.0}};
    var_map_d exp_2_2 = {{"x", 0.0}, {"y", 1.0}};

    BOOST_CHECK(audi::taylor_model::map_equal(exp_2_1, exp_2_2));

    var_map_i dom1 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};
    var_map_i dom2 = {{"x", int_d(2.0, 3.0)}};

    BOOST_CHECK(!audi::taylor_model::map_interval_equal(dom1, dom2));

    var_map_i dom_2_1 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};
    var_map_i dom_2_2 = {{"x", int_d(2.0, 3.0)}, {"y", int_d(2.0, 3.0)}};

    BOOST_CHECK(audi::taylor_model::map_interval_equal(dom_2_1, dom_2_2));
}

BOOST_AUTO_TEST_CASE(validity_check)
{
    // Constructing a Taylor model
    uint order = 5;
    audi::gdual_d x(0.0, "x", order);
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    var_map_i dom = {{"y", int_dom}};
    BOOST_CHECK_THROW(audi::taylor_model tm_x(x, rem, exp, dom), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(construction_and_getters_univariate)
{
    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::taylor_model tm_x(x, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_x.get_tpol(), x);
    BOOST_CHECK_EQUAL(tm_x.get_ndim(), 1);
    BOOST_CHECK_EQUAL(tm_x.get_order(), order);
    BOOST_CHECK(audi::taylor_model::interval_equal(tm_x.get_rem(), rem));

    // Test the expansion points unordered map
    BOOST_CHECK_EQUAL(tm_x.get_exp().size(), exp.size());
    BOOST_CHECK(audi::taylor_model::map_equal(tm_x.get_exp(), exp));
    BOOST_CHECK(audi::taylor_model::map_interval_equal(tm_x.get_dom(), dom));
}

BOOST_AUTO_TEST_CASE(construction_and_getters_multivariate)
{
    // Constructing a multivariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}, {"y", 1.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}, {"y", int_d(1.0, 2.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d y(exp.find("y")->second, "y", order);

    audi::gdual_d f_xy = x * x + y + 20;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);

    BOOST_CHECK_EQUAL(tm_xy.get_tpol(), f_xy);
    BOOST_CHECK_EQUAL(tm_xy.get_ndim(), 2);
    BOOST_CHECK_EQUAL(tm_xy.get_order(), order);
    BOOST_CHECK(audi::taylor_model::interval_equal(tm_xy.get_rem(), rem));

    // Test the expansion points unordered map
    BOOST_CHECK_EQUAL(tm_xy.get_exp().size(), exp.size());
    BOOST_CHECK(audi::taylor_model::map_equal(tm_xy.get_exp(), exp));
    BOOST_CHECK(audi::taylor_model::map_interval_equal(tm_xy.get_dom(), dom));
}

BOOST_AUTO_TEST_CASE(construction_and_getters_identity)
{
    audi::taylor_model tm_id = audi::taylor_model::identity();
    BOOST_CHECK_EQUAL(tm_id.get_tpol(), audi::gdual<double>(1.0));
    BOOST_CHECK_EQUAL(tm_id.get_ndim(), 0);
    BOOST_CHECK_EQUAL(tm_id.get_order(), 0);
    BOOST_CHECK(audi::taylor_model::interval_equal(tm_id.get_rem(), int_d(0.0)));

    // Constructing a multivariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}, {"y", 1.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}, {"y", int_d(1.0, 2.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d y(exp.find("y")->second, "y", order);

    audi::gdual_d f_xy = x * x + y + 20;
    audi::taylor_model tm_xy(f_xy, rem, exp, dom);
    audi::taylor_model tm_id_2 = audi::taylor_model::identity(rem, exp, dom);
    BOOST_CHECK_EQUAL(tm_id_2.get_tpol(), audi::gdual<double>(1.0));
    BOOST_CHECK_EQUAL(tm_id_2.get_ndim(), 0);
    BOOST_CHECK_EQUAL(tm_id_2.get_order(), 0);
    BOOST_CHECK(audi::taylor_model::interval_equal(tm_id_2.get_rem(), rem));
    BOOST_CHECK(audi::taylor_model::map_interval_equal(tm_id_2.get_dom(), dom));
    BOOST_CHECK(audi::taylor_model::map_equal(tm_id_2.get_exp(), exp));
}

BOOST_AUTO_TEST_CASE(comparison_of_construction_order)
{
    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}};
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

BOOST_AUTO_TEST_CASE(get_bounds)
{
    uint order = 5;
    int_d rem(0.0, 2.0);
    var_map_d exp = {{"x", 0.0}, {"y", 1.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}, {"y", int_d(1.0, 2.0)}};
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
    var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    var_map_i dom = {{"x", int_dom}};
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
    var_map_d exp = {{"x", 0.0}};
    int_d int_dom(0.0, 1.0);
    var_map_i dom = {{"x", int_dom}};
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

BOOST_AUTO_TEST_CASE(test_division)
{
    // ---- Setup expansion point and domain ----
    var_map_d exp = {{"x", 1.0}};
    var_map_i dom = {{"x", int_d(0.5, 1.5)}};

    // Initialize gdual (centered at 0.0, variable "x", order 5)
    audi::gdual_d x(exp.find("x")->second, "x", 5);

    // Initialize TaylorModel over domain [0,2], expansion point 0
    int_d rem = int_d(0.0, 0.0);
    audi::taylor_model tm_fx(x, rem, exp, dom);

    // ---- Test division by another TaylorModel ----
    audi::taylor_model prod_ans = 1 / tm_fx;
    audi::taylor_model exp_ans = audi::taylor_model(1 / x, int_d(0.0, 2.0), {{"x", 1.0}}, {{"x", int_d(0.5, 1.5)}});

    BOOST_CHECK(prod_ans == exp_ans);
}

BOOST_AUTO_TEST_CASE(test_division_2)
{
    // ---- Setup expansion point and domain ----
    var_map_d exp = {{"x", 0.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}};

    // Initialize gdual (centered at 0.0, variable "x", order 5)
    audi::gdual_d x(exp.find("x")->second, "x", 5);

    // Initialize TaylorModel over domain [0,2], expansion point 0
    int_d rem = int_d(0.0, 2.0);
    audi::taylor_model tm_fx(x, rem, exp, dom);

    // ---- Test int division ----
    audi::taylor_model prod_ans = tm_fx / 2;
    audi::taylor_model exp_ans = audi::taylor_model(x / 2, int_d(0.0, 1.0), {{"x", 0.0}}, {{"x", int_d(0.0, 1.0)}});
    BOOST_CHECK(prod_ans == exp_ans);

    // ---- Test int division ----
    audi::taylor_model prod_ans_2 = tm_fx / 2.0;
    audi::taylor_model exp_ans_2 = audi::taylor_model(x / 2.0, int_d(0.0, 1.0), {{"x", 0.0}}, {{"x", int_d(0.0, 1.0)}});
    BOOST_CHECK(prod_ans_2 == exp_ans_2);
}

BOOST_AUTO_TEST_CASE(test_power)
{

    // Constructing a univariate Taylor model
    uint order = 5;
    int_d rem(0.0, 0.0);
    var_map_d exp = {{"x", 0.0}};
    var_map_i dom = {{"x", int_d(0.0, 1.0)}};
    audi::gdual_d x(exp.find("x")->second, "x", order);
    audi::gdual_d fx = audi::pow(x, 3);
    audi::taylor_model tm_fx_exp(fx, rem, exp, dom);
    audi::taylor_model tm_x(x, rem, exp, dom);

    audi::taylor_model tm_fx = audi::pow(tm_x, 3);
    BOOST_CHECK(tm_fx == tm_fx_exp);
}

BOOST_AUTO_TEST_CASE(test_makino1998_simpleexample)
{
    // Setup
    int_d rem_bound(0.0, 0.0);
    var_map_d exp_points = {{"x", 2.0}};
    var_map_i domain = {{"x", int_d(1.9, 2.1)}};

    // Expected answers (Table 4.2 p.95)
    std::vector<std::pair<double, double>> exp_ans
        = {{0.0, 1.4579384e-3},  {-7.6733603e-5, 7.6733603e-5},   {0.0, 4.0386107e-6},  {-2.1255845e-7, 2.1255845e-7},
           {0.0, 1.1187287e-8},  {-5.8880459e-10, 5.8880459e-10}, {0.0, 3.0989715e-11}, {-1.6310376e-12, 1.6310376e-12},
           {0.0, 8.5844087e-14}, {-4.5181098e-15, 4.5181098e-15}, {0.0, 2.3779525e-16}, {-1.2515539e-17, 1.2515539e-17},
           {0.0, 6.5871262e-19}, {-3.4669086e-20, 3.4669085e-20}, {0.0, 1.8246887e-21}};

    // Loop over orders
    for (int order = 1; order <= 15; ++order) {
        audi::gdual<double> x(exp_points.find("x")->second, "x", order);
        audi::taylor_model const_fx(x, rem_bound, exp_points, domain);
        audi::taylor_model const_fx_2 = 1 / const_fx;
        audi::taylor_model const_fx_3 = const_fx_2 + const_fx;

        BOOST_CHECK_CLOSE_FRACTION(const_fx_3.get_rem().lower(), exp_ans[order - 1].first, 1e-7);
        BOOST_CHECK_CLOSE_FRACTION(const_fx_3.get_rem().upper(), exp_ans[order - 1].second, 1e-7);
    }

    // Construction test without build-up
    audi::gdual<double> x(exp_points.find("x")->second, "x", 3);
    audi::taylor_model tm_fx(x, rem_bound, exp_points, domain);
    audi::taylor_model tm_fx_2 = 1 / tm_fx + tm_fx;

    var_map_d exp_points_2 = {{"y", 2.0}};
    var_map_i domain_2 = {{"y", int_d(1.9, 2.1)}};
    audi::gdual<double> y(exp_points_2.find("y")->second, "y", 3);
    audi::taylor_model const_gy(y, rem_bound, exp_points_2, domain_2);
    audi::taylor_model const_gy_2 = 1 / const_gy;
    audi::taylor_model const_gy_3 = const_gy_2 + const_gy;

    BOOST_TEST(static_cast<double>(tm_fx_2.get_rem().lower()) == static_cast<double>(const_gy_3.get_rem().lower()));
    BOOST_TEST(static_cast<double>(tm_fx_2.get_rem().upper()) == static_cast<double>(const_gy_3.get_rem().upper()));

    // Verify bound values
    int_d exp_bounds(2.42631, 2.57618);

    BOOST_CHECK_CLOSE_FRACTION(tm_fx_2.get_bounds().lower(), exp_bounds.lower(), 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(tm_fx_2.get_bounds().upper(), exp_bounds.upper(), 1e-5);

    // Verify bound interval (Taylor model order)
    std::vector<double> exp_taymodorder_ans = {0.15145793, 0.15015346, 0.14987903, 0.14987542, 0.14987469, 0.14987468};

    for (int order = 5; order <= 6; ++order) {
        int it = order - 1;
        audi::gdual<double> x(exp_points.find("x")->second, "x", order);
        audi::taylor_model T_x(x, rem_bound, exp_points, domain);
        audi::taylor_model T_fx = 1 / T_x + T_x;

        double interval_size
            = static_cast<double>(T_fx.get_bounds().upper()) - static_cast<double>(T_fx.get_bounds().lower());

        BOOST_CHECK_CLOSE_FRACTION(interval_size, exp_taymodorder_ans[it], 1e-7);
    }
}
