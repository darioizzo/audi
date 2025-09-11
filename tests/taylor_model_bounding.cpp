#define BOOST_TEST_MODULE taylor_models_bounding_test

#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <vector>

#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model_bounding.hpp>

using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

BOOST_AUTO_TEST_CASE(test_generate_combinations_small)
{
    std::vector<int> limits = {2, 1};
    auto combos = audi::generate_combinations(limits, true);

    BOOST_CHECK_EQUAL(combos.size(), 5);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 1}));
    BOOST_CHECK(combos[2] == std::vector<int>({1, 0}));
    BOOST_CHECK(combos[3] == std::vector<int>({1, 1}));
    BOOST_CHECK(combos[4] == std::vector<int>({2, 0}));
}

BOOST_AUTO_TEST_CASE(test_generate_combinations_small_2)
{
    std::vector<int> limits = {2, 1};
    auto combos = audi::generate_combinations(limits);
    BOOST_CHECK_EQUAL(combos.size(), 6);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 1}));
    BOOST_CHECK(combos[5] == std::vector<int>({2, 1}));
}

BOOST_AUTO_TEST_CASE(test_generate_combinations_empty_limits)
{
    std::vector<int> limits;
    auto combos = audi::generate_combinations(limits);

    BOOST_CHECK_EQUAL(combos.size(), 0);
}

BOOST_AUTO_TEST_CASE(test_generate_combinations_multidim)
{
    std::vector<int> limits = {1, 1, 1}; // 2*2*2 = 8 combos
    auto combos = audi::generate_combinations(limits, true);

    BOOST_CHECK_EQUAL(combos.size(), 4);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 0, 1}));
    BOOST_CHECK(combos[2] == std::vector<int>({0, 1, 0}));
    BOOST_CHECK(combos[3] == std::vector<int>({1, 0, 0}));
}

BOOST_AUTO_TEST_CASE(test_generate_combinations_multidim_2)
{
    std::vector<int> limits = {1, 1, 1}; // 2*2*2 = 8 combos
    auto combos = audi::generate_combinations(limits);

    BOOST_CHECK_EQUAL(combos.size(), 8);
    BOOST_CHECK(combos[0] == std::vector<int>({0, 0, 0}));
    BOOST_CHECK(combos[1] == std::vector<int>({0, 0, 1}));
    BOOST_CHECK(combos[2] == std::vector<int>({0, 1, 0}));
    BOOST_CHECK(combos[3] == std::vector<int>({0, 1, 1}));
    BOOST_CHECK(combos.back() == std::vector<int>({1, 1, 1}));
}

BOOST_AUTO_TEST_CASE(test_get_poly_univariate)
{
    audi::gdual<double> poly(1.0, "x", 3);

    auto [coeffs, exps] = audi::get_poly(poly);

    auto it = std::find(exps.begin(), exps.end(), std::vector<int>({0}));
    BOOST_CHECK(it != exps.end());
    auto idx = std::distance(exps.begin(), it);
    BOOST_CHECK_EQUAL(coeffs[static_cast<size_t>(idx)], 1.0);

    BOOST_CHECK(!coeffs.empty());
}

BOOST_AUTO_TEST_CASE(test_get_poly_multidim)
{
    using namespace audi;
    audi::gdual<double> x(1.0, "x", 2);
    audi::gdual<double> y(1.0, "y", 2);
    audi::gdual<double> poly = x + 2.0 * y + 3.0 * x * y;

    auto [coeffs, exps] = audi::get_poly(poly);

    // Should contain (6, [0,0]), (5, [0,1]), (4, [1,0], (3, [1,1])
    bool corresponding_coeff = true;
    std::vector<double> exp_coeffs = {6.0, 5.0, 4.0, 3.0};
    std::vector<std::vector<int>> exp_exps = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (exps[i] != exp_exps[i] || coeffs[i] != exp_coeffs[i]) corresponding_coeff = false;
    }

    BOOST_CHECK_EQUAL(corresponding_coeff, true);
}

BOOST_AUTO_TEST_CASE(test_get_ndim_valid)
{
    std::vector<double> coeffs = {1.0, 2.0, 3.0};
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};

    uint ndim = audi::get_ndim(coeffs, exps);
    BOOST_CHECK_EQUAL(ndim, 3);
}

BOOST_AUTO_TEST_CASE(test_get_ndim_inconsistent_coeffs_exps)
{
    std::vector<double> coeffs = {1.0, 2.0};
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};

    BOOST_CHECK_THROW(audi::get_ndim(coeffs, exps), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_get_max_degrees_ndim_greater_than_1)
{
    std::vector<std::vector<int>> exps = {{2, 3, 1}, {4, 1, 5}, {3, 2, 2}};
    uint ndim = 3;
    std::vector<int> max_degrees = audi::get_max_degrees(exps, ndim);

    std::vector<int> expected = {4, 3, 5};
    BOOST_CHECK_EQUAL_COLLECTIONS(max_degrees.begin(), max_degrees.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(test_get_max_degrees_ndim_1)
{
    std::vector<std::vector<int>> exps = {{2, 3, 1, 4}};
    uint ndim = 1;
    std::vector<int> max_degrees = audi::get_max_degrees(exps, ndim);

    BOOST_CHECK_EQUAL(max_degrees.size(), 1);
    BOOST_CHECK_EQUAL(max_degrees[0], 4);
}

BOOST_AUTO_TEST_CASE(test_get_coefficient_found)
{
    std::vector<double> coeffs = {22.0, 16.0};
    std::vector<std::vector<int>> exps = {{1, 2}, {4, 3}};

    double c1 = audi::get_coefficient<double>({1, 2}, coeffs, exps);
    double c2 = audi::get_coefficient<double>({4, 3}, coeffs, exps);

    BOOST_CHECK_EQUAL(c1, 22.0);
    BOOST_CHECK_EQUAL(c2, 16.0);
}

BOOST_AUTO_TEST_CASE(test_get_coefficient_not_found)
{
    std::vector<double> coeffs = {22.0, 16.0};
    std::vector<std::vector<int>> exps = {{1, 2}, {4, 3}};

    double c = audi::get_coefficient<double>({0, 0}, coeffs, exps);
    BOOST_CHECK_EQUAL(c, 0.0); // default return
}

BOOST_AUTO_TEST_CASE(test_a_matrix_2d)
{
    std::vector<double> coeffs = {22.0, 16.0};
    std::vector<std::vector<int>> exps = {{1, 2}, {4, 3}};

    auto vec_A = audi::get_a_matrix_vec<double>(coeffs, exps);

    std::vector<std::vector<double>> expected = {
        {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 22.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 16.0}};

    BOOST_CHECK_EQUAL(vec_A.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        BOOST_CHECK_EQUAL(vec_A[i].size(), expected[i].size());
        for (size_t j = 0; j < expected[i].size(); ++j) {
            BOOST_CHECK_CLOSE(vec_A[i][j], expected[i][j], 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_binomial)
{
    BOOST_CHECK_EQUAL(audi::binomial(5, 0), 1ULL);
    BOOST_CHECK_EQUAL(audi::binomial(5, 5), 1ULL);
    BOOST_CHECK_EQUAL(audi::binomial(5, 2), 10ULL);
    BOOST_CHECK_EQUAL(audi::binomial(6, 3), 20ULL);
    BOOST_CHECK_EQUAL(audi::binomial(10, 3), 120ULL);
    BOOST_CHECK_EQUAL(audi::binomial(0, 0), 1ULL);
    BOOST_CHECK_EQUAL(audi::binomial(5, -1), 0ULL); // invalid
    BOOST_CHECK_EQUAL(audi::binomial(5, 6), 0ULL);  // invalid
}

BOOST_AUTO_TEST_CASE(test_get_d_acc_matrix)
{
    std::vector<int> max_degrees = {3};
    auto Dacc = audi::get_d_acc_matrix(max_degrees, 0);

    BOOST_REQUIRE_EQUAL(Dacc.size(), 4);
    BOOST_CHECK_CLOSE(Dacc[0][0], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(Dacc[1][1], 1.0 / 3.0, 1e-12);
    BOOST_CHECK_CLOSE(Dacc[2][2], 1.0 / 3.0, 1e-12);
    BOOST_CHECK_CLOSE(Dacc[3][3], 1.0, 1e-12);

    // off-diagonal must be zero
    BOOST_CHECK_EQUAL(Dacc[0][1], 0.0);
}

BOOST_AUTO_TEST_CASE(test_get_d_matrix)
{
    std::vector<int> max_degrees = {3};
    double t = 2.0;
    auto D = audi::get_d_matrix(max_degrees, 0, t);

    BOOST_REQUIRE_EQUAL(D.size(), 4);
    BOOST_CHECK_CLOSE(D[0][0], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(D[1][1], 2.0, 1e-12);
    BOOST_CHECK_CLOSE(D[2][2], 4.0, 1e-12);
    BOOST_CHECK_CLOSE(D[3][3], 8.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_get_lower_pascal_matrix)
{
    std::vector<int> max_degrees = {4};
    auto P = audi::get_lower_pascal_matrix(max_degrees, 0);

    BOOST_REQUIRE_EQUAL(P.size(), 5);

    BOOST_CHECK_EQUAL(P[0][0], 1.0);
    BOOST_CHECK_EQUAL(P[4][2], 6.0); // binomial(4,2)
    BOOST_CHECK_EQUAL(P[3][1], 3.0); // binomial(3,1)

    // upper part must be zero
    BOOST_CHECK_EQUAL(P[1][2], 0.0);
}

BOOST_AUTO_TEST_CASE(test_build_dim_to_var_map)
{
    var_map_i domain;
    domain["x"] = int_d(0.0, 1.0);
    domain["y"] = int_d(2.0, 3.0);

    auto dim_map = audi::build_dim_to_var_map(domain);

    BOOST_REQUIRE_EQUAL(dim_map.size(), 2);

    bool found0 = (dim_map.at(0).first == "x" || dim_map.at(0).first == "y");
    bool found1 = (dim_map.at(1).first == "x" || dim_map.at(1).first == "y");
    BOOST_CHECK(found0);
    BOOST_CHECK(found1);
}

BOOST_AUTO_TEST_CASE(test_get_q_matrix)
{
    var_map_i domain;
    domain["x"] = int_d(0.0, 2.0);
    domain["y"] = int_d(1.0, 3.0);
    std::vector<int> max_degrees = {2, 2};

    // Case 1: a == 0.0 → should return D(b)
    auto Q0 = audi::get_q_matrix(max_degrees, domain, 0);
    BOOST_CHECK_CLOSE(Q0[1][1], 2.0, 1e-12); // t^1 = 2
    BOOST_CHECK_CLOSE(Q0[2][2], 4.0, 1e-12); // t^2 = 4

    // Case 2: a != 0.0
    auto Q1 = audi::get_q_matrix(max_degrees, domain, 1);
    BOOST_REQUIRE_EQUAL(Q1.size(), 3);
    BOOST_REQUIRE_EQUAL(Q1[0].size(), 3);

    // We don’t hardcode exact values (since it involves matmul),
    // but at least check it's invertible-looking and diagonal entries > 0
    BOOST_CHECK(Q1[0][0] > 0.0);
    BOOST_CHECK(Q1[1][1] > 0.0);
}

BOOST_AUTO_TEST_CASE(test_shift_indices_empty)
{
    std::vector<int> input{};
    auto result = audi::shift_indices(input);
    BOOST_TEST(result.empty());
}

BOOST_AUTO_TEST_CASE(test_shift_indices_basic)
{
    std::vector<int> input{1, 2, 3, 4};
    auto result = audi::shift_indices(input);
    std::vector<int> expected{2, 3, 4, 1};
    BOOST_TEST(result == expected, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(test_a_matrix_vec_and_cycle)
{
    // coefficients and exponents
    std::vector<double> coeffs{1.0, 1.0, -1.1, 1.0};
    std::vector<std::vector<int>> exps{{0, 0, 0}, {1, 2, 0}, {1, 0, 0}, {1, 0, 2}};

    // compute ndim and max_degrees
    unsigned int ndim = audi::get_ndim(coeffs, exps);
    std::vector<int> max_degrees = audi::get_max_degrees(exps, ndim);

    // get A-matrix
    auto A = audi::get_a_matrix_vec(coeffs, exps);

    // check A-matrix values (row-major)
    BOOST_CHECK_CLOSE(A[0][0], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(A[1][0], -1.1, 1e-12);
    BOOST_CHECK_CLOSE(A[1][2], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(A[1][6], 1.0, 1e-12);

    // expected shape: (2 x 9)
    BOOST_CHECK_EQUAL(A.size(), 2);
    BOOST_CHECK_EQUAL(A[0].size(), 9);

    // now test cycling along dimension 0
    auto cycled = audi::get_cycled_2d_array_vec(A, max_degrees, 0);

    // check a few key entries according to your example
    BOOST_CHECK_CLOSE(cycled[0][0], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(cycled[0][3], -1.1, 1e-12);
    BOOST_CHECK_CLOSE(cycled[0][5], 1.0, 1e-12);
    BOOST_CHECK_CLOSE(cycled[2][3], 1.0, 1e-12);

    // shape check: 3 rows, 6 columns (as per example)
    BOOST_CHECK_EQUAL(cycled.size(), 3);
    BOOST_CHECK_EQUAL(cycled[0].size(), 6);
}

// from Titi (2018) "Matrix methods for tensorial Bernstein form"
BOOST_AUTO_TEST_CASE(test_lambda_generalbox_known_case)
{
    // coefficients
    std::vector<double> coeffs{22.0, 16.0};

    // exponents
    std::vector<std::vector<int>> exps{{1, 2}, {4, 3}};

    // domain [0,1] for x and y
    var_map_i domain{{"x", int_d(0.0, 1.0)}, {"y", int_d(0.0, 1.0)}};

    auto Lambda = audi::get_titi_base_lambda_generalbox(coeffs, exps, domain);

    // expected result (verified from Python implementation)
    std::vector<std::vector<double>> expected = {{0.0, 0.0, 0.0, 0.0},
                                                 {0.0, 0.0, 1.83333333, 0.0},
                                                 {0.0, 0.0, 0.0, 0.0},
                                                 {0.0, 0.0, 0.0, 0.0},
                                                 {0.0, 0.0, 0.0, 16.0}};

    BOOST_REQUIRE_EQUAL(Lambda.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        BOOST_REQUIRE_EQUAL(Lambda[i].size(), expected[i].size());
        for (std::size_t j = 0; j < expected[i].size(); ++j) {
            BOOST_CHECK_CLOSE(Lambda[i][j], expected[i][j], 1e-6); // tolerance
        }
    }
}

// Mag 6 sub-problem from Titi (2018) "Matrix methods for tensorial Bernstein form"
BOOST_AUTO_TEST_CASE(test_lambda_generalbox_selective)
{
    std::vector<double> coeffs{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -1.0};

    std::vector<std::vector<int>> exps{{0, 0, 0, 2, 0, 0}, {0, 0, 2, 0, 0, 0}, {0, 0, 0, 0, 2, 0}, {2, 0, 0, 0, 0, 0},
                                       {0, 2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 0, 1}};

    var_map_i domain{{"x1", int_d(-5.0, 5.0)}, {"x2", int_d(-5.0, 5.0)}, {"x3", int_d(-5.0, 5.0)},
                     {"x4", int_d(-5.0, 5.0)}, {"x5", int_d(-5.0, 5.0)}, {"x6", int_d(-5.0, 5.0)}};

    auto Lambda = audi::get_titi_base_lambda_generalbox(coeffs, exps, domain);

    // selectively check a few non-zero entries
    BOOST_CHECK_CLOSE(Lambda[0][0], 305.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][1], -100.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][2], 200.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[1][0], -100.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[2][0], 200.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][27], -100.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][54], 200.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][81], -105.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[0][162], 200.0, 1e-6);
    BOOST_CHECK_CLOSE(Lambda[1][12], 0.0, 1e-6); // one zero check
    BOOST_CHECK_CLOSE(Lambda[2][36], 0.0, 1e-6); // another zero check
}

// Worked out example from Titi (2018) "Matrix methods for tensorial Bernstein form"
BOOST_AUTO_TEST_CASE(test_titi_bernstein_patch_ndim_full_matrix)
{
    // Coefficients and exponents
    std::vector<double> coeffs = {170.0, 1.0, -22.0, 2.0, 1.0, -21.0, 2.0, -13.0, -14.0};
    std::vector<std::vector<int>> exps = {{0, 0}, {4, 0}, {0, 1}, {2, 1}, {0, 4}, {2, 0}, {1, 2}, {0, 2}, {1, 0}};

    // Domain
    var_map_i domain = {{"x", int_d(-5.0, 5.0)}, {"y", int_d(-5.0, 5.0)}};

    // Compute Lambda patch
    std::vector<std::vector<double>> lambda_patch
        = audi::get_titi_bernstein_patch_ndim_generalbox(coeffs, exps, domain);

    // Expected Bernstein patch
    std::vector<std::vector<double>> expected_patch
        = {{250.0, -355.0, 1156.66666667, -215.0, 530.0},
           {-135.0, -990.0, 355.0, -1100.0, -355.0},
           {1463.33333333, 441.66666667, 1703.33333333, 248.33333333, 1076.66666667},
           {45.0, -1060.0, 201.66666667, -1170.0, -175.0},
           {610.0, -495.0, 850.0, -355.0, 890.0}};

    // Tolerance for floating-point comparison
    const double tol = 1e-6;

    // Compare element-wise
    BOOST_REQUIRE_EQUAL(lambda_patch.size(), expected_patch.size());
    for (size_t i = 0; i < lambda_patch.size(); ++i) {
        BOOST_REQUIRE_EQUAL(lambda_patch[i].size(), expected_patch[i].size());
        for (size_t j = 0; j < lambda_patch[i].size(); ++j) {
            BOOST_TEST(lambda_patch[i][j] == expected_patch[i][j], boost::test_tools::tolerance(tol));
        }
    }
}

// L.V. 3 problem from Titi (2018) "Matrix methods for tensorial Bernstein form"
BOOST_AUTO_TEST_CASE(test_titi_bernstein_patch_ndim_3d_example)
{
    // Coefficients and exponents
    std::vector<double> coeffs = {1.0, 1.0, -1.1, 1.0};
    std::vector<std::vector<int>> exps = {{0, 0, 0}, {1, 2, 0}, {1, 0, 0}, {1, 0, 2}};

    // Domain
    var_map_i domain = {{"x", int_d(-1.5, 2.0)}, {"y", int_d(-1.5, 2.0)}, {"z", int_d(-1.5, 2.0)}};

    // Compute Lambda patch
    std::vector<std::vector<double>> lambda_patch
        = audi::get_titi_bernstein_patch_ndim_generalbox(coeffs, exps, domain);

    // Expected Bernstein patch
    std::vector<std::vector<double>> expected_patch = {{-4.1, 3.775, -6.725, 3.775, 11.65, 1.15, -6.725, 1.15, -9.35},
                                                       {7.8, -2.7, 11.3, -2.7, -13.2, 0.8, 11.3, 0.8, 14.8}};

    // Tolerance for floating-point comparison
    const double tol = 1e-6;

    // Compare element-wise
    BOOST_REQUIRE_EQUAL(lambda_patch.size(), expected_patch.size());
    for (size_t i = 0; i < lambda_patch.size(); ++i) {
        BOOST_REQUIRE_EQUAL(lambda_patch[i].size(), expected_patch[i].size());
        for (size_t j = 0; j < lambda_patch[i].size(); ++j) {
            BOOST_TEST(lambda_patch[i][j] == expected_patch[i][j], boost::test_tools::tolerance(tol));
        }
    }
}

// Mag 6 problem from Titi (2018) "Matrix methods for tensorial Bernstein form"
BOOST_AUTO_TEST_CASE(test_titi_bernstein_patch_ndim_6d_example_partial)
{

    // Coefficients and exponents
    std::vector<double> coeffs = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -1.0};
    std::vector<std::vector<int>> exps
        = {{0, 0, 0, 2, 0, 0}, {0, 0, 2, 0, 0, 0}, {0, 0, 0, 0, 2, 0}, {2, 0, 0, 0, 0, 0},
           {0, 2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 0, 1}};

    // Domain
    var_map_i domain = {{"x1", int_d(-5.0, 5.0)}, {"x2", int_d(-5.0, 5.0)}, {"x3", int_d(-5.0, 5.0)},
                        {"x4", int_d(-5.0, 5.0)}, {"x5", int_d(-5.0, 5.0)}, {"x6", int_d(-5.0, 5.0)}};

    // Compute Lambda patch
    std::vector<std::vector<double>> lambda_patch
        = audi::get_titi_bernstein_patch_ndim_generalbox(coeffs, exps, domain);

    // Select a few elements to test (row, col, expected value)
    struct ElemCheck {
        size_t row;
        size_t col;
        double value;
    };
    std::vector<ElemCheck> checks = {
        {0, 0, 305.0},  {0, 9, 205.0},  {0, 18, 305.0}, {0, 27, 205.0}, {0, 36, 105.0}, {0, 45, 205.0}, {0, 54, 305.0},
        {0, 63, 205.0}, {0, 72, 305.0}, {0, 81, 200.0}, {0, 4, 105.0},  {0, 13, 5.0},   {0, 22, 105.0}, {0, 31, 5.0},
        {0, 40, -95.0}, {0, 49, 5.0},   {0, 58, 105.0}, {0, 67, 5.0},   {0, 76, 105.0}, {0, 85, 0.0},
    };

    // Tolerance for floating-point comparison
    const double tol = 1e-6;

    for (const auto &c : checks) {
        BOOST_TEST(lambda_patch[c.row][c.col] == c.value, boost::test_tools::tolerance(tol));
    }
}
