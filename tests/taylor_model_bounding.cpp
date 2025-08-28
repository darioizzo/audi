#define BOOST_TEST_MODULE taylor_models_bounding_test

#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <vector>

#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model_bounding.hpp>

BOOST_AUTO_TEST_SUITE(test1)
BOOST_AUTO_TEST_CASE(test_generate_combinations_small)
{
    std::vector<int> limits = {2, 1};
    auto combos = generate_combinations(limits);

    BOOST_CHECK_EQUAL(combos.size(), 5);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 1}));
    BOOST_CHECK(combos[2] == std::vector<int>({1, 0}));
    BOOST_CHECK(combos[3] == std::vector<int>({1, 1}));
    BOOST_CHECK(combos[4] == std::vector<int>({2, 0}));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test2)
BOOST_AUTO_TEST_CASE(test_generate_combinations_empty_limits)
{
    std::vector<int> limits;
    auto combos = generate_combinations(limits);

    BOOST_CHECK_EQUAL(combos.size(), 0);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test3)
BOOST_AUTO_TEST_CASE(test_generate_combinations_multidim)
{
    std::vector<int> limits = {1, 1, 1}; // 2*2*2 = 8 combos
    auto combos = generate_combinations(limits);

    BOOST_CHECK_EQUAL(combos.size(), 4);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 0, 1}));
    BOOST_CHECK(combos[2] == std::vector<int>({0, 1, 0}));
    BOOST_CHECK(combos[3] == std::vector<int>({1, 0, 0}));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test4)
BOOST_AUTO_TEST_CASE(test_get_poly_univariate)
{
    audi::gdual<double> poly(1.0, "x", 3);

    auto [coeffs, exps] = get_poly(poly);

    auto it = std::find(exps.begin(), exps.end(), std::vector<int>({0}));
    BOOST_CHECK(it != exps.end());
    auto idx = std::distance(exps.begin(), it);
    BOOST_CHECK_EQUAL(coeffs[static_cast<size_t>(idx)], 1.0);

    BOOST_CHECK(!coeffs.empty());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test6)
BOOST_AUTO_TEST_CASE(test_get_poly_multidim)
{
    using namespace audi;
    audi::gdual<double> x(1.0, "x", 2);
    audi::gdual<double> y(1.0, "y", 2);
    audi::gdual<double> poly = x + 2.0 * y + 3.0 * x * y;

    auto [coeffs, exps] = get_poly(poly);

    // Should contain (6, [0,0]), (5, [0,1]), (4, [1,0], (3, [1,1])
    bool corresponding_coeff = true;
    std::vector<double> exp_coeffs = {6.0, 5.0, 4.0, 3.0};
    std::vector<std::vector<int>> exp_exps = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (exps[i] != exp_exps[i] || coeffs[i] != exp_coeffs[i]) corresponding_coeff=false;
    }

    BOOST_CHECK_EQUAL(corresponding_coeff, true);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test7)
BOOST_AUTO_TEST_CASE(test_get_ndim_valid)
{
    std::vector<double> coeffs = {1.0, 2.0, 3.0};
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};

    int ndim = get_ndim(coeffs, exps);
    BOOST_CHECK_EQUAL(ndim, 3);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test8)
BOOST_AUTO_TEST_CASE(test_get_ndim_inconsistent_coeffs_exps)
{
    std::vector<double> coeffs = {1.0, 2.0};
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};

    BOOST_CHECK_THROW(get_ndim(coeffs, exps), std::invalid_argument);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test9)
BOOST_AUTO_TEST_CASE(test_get_max_degrees_ndim_greater_than_1)
{
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};
    int ndim = 3;
    std::vector<int> max_degrees = get_max_degrees(exps, ndim);

    std::vector<int> expected = {4, 3, 5};
    BOOST_CHECK_EQUAL_COLLECTIONS(max_degrees.begin(), max_degrees.end(), expected.begin(), expected.end());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test10)
BOOST_AUTO_TEST_CASE(test_get_max_degrees_ndim_1)
{
    std::vector<std::vector<int>> exps = {{2, 3, 1, 4}};
    int ndim = 1;
    std::vector<int> max_degrees = get_max_degrees(exps, ndim);

    BOOST_CHECK_EQUAL(max_degrees.size(), 1);
    BOOST_CHECK_EQUAL(max_degrees[0], 4);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test11)
BOOST_AUTO_TEST_CASE(test_get_max_degrees_empty)
{
    std::vector<std::vector<int>> exps;
    int ndim = 2;
    BOOST_CHECK_THROW(get_max_degrees(exps, ndim), std::invalid_argument);
}
BOOST_AUTO_TEST_SUITE_END()
