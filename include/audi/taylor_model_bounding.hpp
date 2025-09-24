#ifndef AUDI_TM_BOUNDING_HPP
#define AUDI_TM_BOUNDING_HPP

#include <algorithm>
#include <boost/numeric/interval.hpp>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

#include <audi/gdual.hpp>
#include <audi/taylor_model_utilities.hpp>

// using int_d = boost::numeric::interval<double>;
using int_d = boost::numeric::interval<
    double, boost::numeric::interval_lib::policies<
                boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double>>,
                boost::numeric::interval_lib::checking_no_nan<double>>>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

namespace audi
{

/// Function for the combinations of given maximum indices
/**
 * Function gives either all possible combinations given a list of maximum limits, or (if
 * cap_sum_indices = true), the total sum of the indices of a given combination are limited by the
 * maximum index given by the input maximum indices.
 *
 * The cap_sum_indices defaults to off, and is generally only set to true when retrieving the
 * polynomial from a gdual (to prevent a monomial from exceeding the degree of the gdual).
 *
 * @param limits std::vector<int> the maximum degree limits per dimension
 * @param cap_sum_indices bool whether or not to limit the combinations by the maximum of the
 * maximum degrees
 *
 * @return a nested vector of combinations MxN where M is the number of combinations, and N is the
 * number of dimensions that then have a given degree.
 */
std::vector<std::vector<int>> generate_combinations(const std::vector<int> &limits, bool cap_sum_indices = false)
{
    if (limits.empty()) return {};

    int cap = *std::max_element(limits.begin(), limits.end());
    std::vector<std::vector<int>> result = {{}};

    for (int limit : limits) {
        std::vector<std::vector<int>> tmp;
        for (const auto &comb : result) {
            for (int i = 0; i <= limit; ++i) {
                auto copy = comb;
                copy.push_back(i);

                // Check sum of indices
                int s = std::accumulate(copy.begin(), copy.end(), 0);
                if ((s <= cap && cap_sum_indices == true) || (cap_sum_indices == false)) {
                    tmp.push_back(std::move(copy));
                }
            }
        }
        result = std::move(tmp);
    }

    return result;
}

/// Extract the polynomial representation of a gdual
/**
 * Retrieves the coefficients and corresponding exponents of the nonzero monomials
 * contained in a gdual.
 *
 * For univariate polynomials, the exponents are generated from 0 up to the degree
 * of the gdual. For multivariate polynomials, all valid exponent combinations are
 * generated such that the sum of indices does not exceed the degree of the `gdual`.
 *
 * @tparam T numeric type of the coefficients (e.g., `double`, `long double`)
 *
 * @param tpol the input truncated polynomial (`audi::gdual<T>`) from which the
 * polynomial representation is extracted
 *
 * @return a pair `(coeffs, exponents)` where:
 *   - `coeffs` is a vector of coefficients of type `T`
 *   - `exponents` is a vector of integer vectors, each containing the exponents
 *      of a corresponding monomial term in `coeffs`
 */
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<int>>> get_poly(const audi::gdual<T> &tpol)
{

    auto ndim = tpol.get_symbol_set_size();

    std::vector<std::vector<int>> combs;
    if (ndim == 1) {
        int deg = static_cast<int>(tpol.degree());
        for (int i = 0; i <= deg; ++i) {
            combs.push_back({i});
        }
    } else if (ndim > 1) {
        std::vector<int> limits(static_cast<size_t>(ndim), static_cast<const int>(tpol.degree()));
        combs = generate_combinations(limits, true);
    }

    std::vector<T> coeffs;
    std::vector<std::vector<int>> exponents;

    for (const auto &exps : combs) {
        T cf = tpol.find_cf(exps); // coefficient for monomial x^exps
        if (cf != T(0)) {
            coeffs.push_back(cf);
            exponents.push_back(exps);
        }
    }

    return {coeffs, exponents};
}

/// Get the number of dimensions in a polynomial representation
/**
 * Determines the dimensionality of a polynomial given its coefficients and
 * exponent matrix.
 *
 * @param coeffs vector of coefficients
 * @param exps vector of exponent rows, each row representing the exponents
 *             of a monomial term
 *
 * @return the number of dimensions (variables) of the polynomial
 *
 * @throws std::invalid_argument if the number of coefficients does not match
 *         the number of exponent rows, or if exponent rows are inconsistent
 *         in length
 */
unsigned int get_ndim(const std::vector<double> &coeffs, const std::vector<std::vector<int>> &exps)
{
    if (coeffs.empty() || exps.empty()) {
        return 0;
    }

    // number of terms must match
    if (coeffs.size() != exps.size() && coeffs.size() != 0 && exps.size() != 0) {
        throw std::invalid_argument("Number of coefficients and number of exponent rows must match.");
    }

    size_t ndim = exps[0].size();
    for (const auto &row : exps) {
        if (row.size() != ndim) {
            throw std::invalid_argument("All exponent rows must have the same length.");
        }
    }

    return static_cast<unsigned int>(ndim);
}

/// Get the maximum degree in each dimension
/**
 * Computes the maximum degree reached in each variable of a polynomial,
 * given its exponent representation. The result is a vector of per-dimension
 * maximum exponents.
 *
 * - For multivariate polynomials (`ndim > 1`), returns one maximum per variable.
 * - For univariate polynomials (`ndim == 1`), returns a single-element vector
 *   containing the maximum degree.
 * - For zero-dimensional input (`ndim == 0`), returns `{0}`.
 *
 * @param exponents vector of exponent rows, each row representing the exponents
 *                  of a monomial term
 * @param ndim number of dimensions (variables) of the polynomial
 *
 * @return vector of maximum degrees, one per dimension
 *
 * @throws std::invalid_argument if exponent rows are inconsistent in length
 *         or if `ndim < 0`
 */
std::vector<int> get_max_degrees(const std::vector<std::vector<int>> &exponents, unsigned int ndim)
{
    if (ndim > 1) {
        size_t cols = exponents[0].size();
        std::vector<int> max_degrees(cols, 0);

        for (const auto &row : exponents) {
            if (row.size() != cols) {
                throw std::invalid_argument("All inner vectors must have the same size.");
            }
            for (size_t i = 0; i < cols; ++i) {
                max_degrees[i] = std::max(max_degrees[i], row[i]);
            }
        }
        return max_degrees;
    } else if (ndim == 1) {
        int max_degree = 0;
        for (const auto &row : exponents) {
            for (int val : row) {
                max_degree = std::max(max_degree, val);
            }
        }
        return {max_degree}; // single element vector
    } else if (ndim == 0) {
        return {0};
    } else {
        throw std::invalid_argument("ndim must be >= 1");
    }
}

/// Retrieve a coefficient by exponent combination
/**
 * Finds the coefficient of a monomial term in a polynomial representation
 * given its exponent combination. If the exponent vector is not found in
 * the exponent list, the coefficient is assumed to be zero.
 *
 * @tparam T numeric type of the coefficients (e.g., `double`, `long double`)
 *
 * @param comb exponent combination to look up (multi-index)
 * @param coeffs vector of coefficients
 * @param exps vector of exponent rows corresponding to the coefficients
 *
 * @return the coefficient corresponding to the given exponent combination,
 *         or `0` if not found
 */
template <typename T>
T get_coefficient(const std::vector<int> &comb, const std::vector<T> &coeffs, const std::vector<std::vector<int>> &exps)
{
    for (size_t k = 0; k < exps.size(); ++k) {
        if (exps[k] == comb) {
            return coeffs[k];
        }
    }
    return 0.0; // if not found, coefficient is zero
}

/// Construct the coefficient matrix representation of a polynomial
/**
 * Builds a 2D matrix `A` from a polynomial given in coefficient–exponent form.
 *
 * - The **rows** of `A` correspond to the degree of the first variable.
 * - The **columns** of `A` represent the flattened index of all remaining variables.
 *
 * More precisely, the matrix has size:
 *   - `(max_degrees[0] + 1)` rows
 *   - `∏_{i=1..ndim-1} (max_degrees[i] + 1)` columns
 *
 * Each entry `A[i][j]` contains the coefficient of the monomial with multi-index
 * `(i, *rest)` mapped into column index `j`.
 *
 * @tparam T numeric type of the coefficients (e.g., `double`, `long double`)
 *
 * @param coeffs vector of coefficients
 * @param exps vector of exponent rows corresponding to the coefficients
 *
 * @return a 2D vector `A` holding the coefficients arranged as a matrix
 */
template <typename T>
std::vector<std::vector<T>> get_a_matrix_vec(const std::vector<T> &coeffs, const std::vector<std::vector<int>> &exps)
{
    // number of dimensions
    unsigned int ndim = get_ndim(coeffs, exps);

    // max degrees for each dimension
    std::vector<int> max_degrees = get_max_degrees(exps, ndim);

    // compute l_star = ∏_{i=1..ndim-1} (max_degrees[i] + 1)
    int l_star = 1;
    for (size_t i = 1; i < static_cast<size_t>(ndim); ++i) {
        l_star *= (max_degrees[i] + 1);
    }

    // initialize A with zeros
    std::vector<std::vector<T>> A(static_cast<size_t>(max_degrees[0]) + 1,
                                  std::vector<T>(static_cast<size_t>(l_star), 0.0));

    // generate all index combinations (multi-indices up to max_degrees)
    std::vector<std::vector<int>> combs = generate_combinations(max_degrees);

    for (std::vector<int> &comb : combs) {
        // row index (first component)
        size_t i_val = static_cast<size_t>(comb[0]);

        // column index (flatten the remaining components into one index)
        size_t j_val = 0;
        int stride = 1;
        for (size_t q = 1; q < static_cast<size_t>(ndim); ++q) {
            j_val += static_cast<size_t>(comb[q] * stride);
            stride *= (max_degrees[q] + 1);
        }

        // coefficient lookup
        T a_val = get_coefficient<T>(comb, coeffs, exps);

        // assign
        A[i_val][j_val] = a_val;
    }

    return A;
}

/// Compute a binomial coefficient
/**
 * Computes the number of ways to choose `k` elements from `n` without repetition:
 *
 * \f[
 *   \binom{n}{k} = \frac{n!}{k!(n-k)!}
 * \f]
 *
 * Uses an iterative approach that avoids intermediate factorial computation
 * and exploits the symmetry relation \f$\binom{n}{k} = \binom{n}{n-k}\f$.
 *
 * @param n nonnegative integer
 * @param k integer satisfying \f$0 \leq k \leq n\f$
 *
 * @return the binomial coefficient as an unsigned long long
 */
unsigned long long binomial(int n, int k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k; // symmetry

    unsigned long long res = 1;
    for (int i = 1; i <= k; ++i) {
        res *= static_cast<unsigned long long>(n - (k - i));
        res /= static_cast<unsigned long long>(i);
    }
    return res;
}

/// Compute the product of binomial coefficients
/**
 * Given two vectors `n` and `k` of the same length, computes
 * the product
 *
 * \f[
 *   \prod_{i=1}^m \binom{n_i}{k_i}
 * \f]
 *
 * where \f$m = n.size() = k.size()\f$.
 *
 * @param n vector of nonnegative integers
 * @param k vector of nonnegative integers of the same size as `n`
 *
 * @return the product of binomial coefficients as an unsigned long long
 *
 * @throws std::runtime_error if `n` and `k` do not have the same size
 */
unsigned long long binom_product(const std::vector<int> &n, const std::vector<int> &k)
{

    if (n.size() != k.size()) throw std::runtime_error("Binomial lists must be of equal size.");
    unsigned long long result = 1;
    for (std::size_t i = 0; i < n.size(); ++i) {
        result *= binomial(n[i], k[i]);
    }
    return result;
}

/// Construct the diagonal accumulation matrix for a given dimension
/**
 * Builds a diagonal matrix of size `(n+1) x (n+1)` where `n` is the maximum
 * degree in the specified dimension. The diagonal entries are given by
 *
 * \f[
 *   (D_{\text{acc}})_{kk} = \frac{1}{\binom{n}{k}}, \quad k = 0, \dots, n.
 * \f]
 *
 * @param max_degrees vector of maximum degrees per dimension
 * @param dim the dimension index to construct the matrix for
 *
 * @return a `(n+1) x (n+1)` diagonal matrix of scaling factors
 *
 * @throws std::out_of_range if `dim` is negative or not less than
 *         `max_degrees.size()`
 */
std::vector<std::vector<double>> get_d_acc_matrix(const std::vector<int> &max_degrees, int dim)
{
    if (dim < 0 || dim >= static_cast<int>(max_degrees.size())) throw std::out_of_range("Dimension out of range.");

    int n = max_degrees[static_cast<size_t>(dim)];

    // initialize diagonal matrix
    std::vector<std::vector<double>> d_acc(static_cast<size_t>(n + 1),
                                           std::vector<double>(static_cast<size_t>(n + 1), 0.0));

    for (int k = 0; k <= n; ++k) {
        double val = 1.0 / static_cast<double>(binomial(n, k));
        d_acc[static_cast<size_t>(k)][static_cast<size_t>(k)] = val;
    }

    return d_acc;
}

/// Construct the diagonal evaluation matrix for a given dimension
/**
 * Builds a diagonal matrix of size `(n+1) x (n+1)` where `n` is the maximum
 * degree in the specified dimension. The diagonal entries are
 *
 * \f[
 *   (D)_{kk} =
 *     \begin{cases}
 *       1 & \text{if } k = 0, \\
 *       t^k & \text{if } k \geq 1.
 *     \end{cases}
 * \f]
 *
 * This corresponds to evaluating powers of `t` up to the maximum degree.
 *
 * @param max_degrees vector of maximum degrees per dimension
 * @param dim the dimension index to construct the matrix for
 * @param t the evaluation point
 *
 * @return a `(n+1) x (n+1)` diagonal matrix of powers of `t`
 */
std::vector<std::vector<double>> get_d_matrix(const std::vector<int> &max_degrees, int dim, double t)
{
    int n = max_degrees[static_cast<size_t>(dim)];
    std::vector<std::vector<double>> D(static_cast<size_t>(n + 1),
                                       std::vector<double>(static_cast<size_t>(n + 1), 0.0));

    D[0][0] = 1.0;
    for (int k = 1; k <= n; ++k) {
        D[static_cast<size_t>(k)][static_cast<size_t>(k)] = std::pow(t, k);
    }

    return D;
}

/// Construct the lower Pascal matrix for a given dimension
/**
 * Builds a lower-triangular Pascal matrix of size `(n+1) x (n+1)` where `n`
 * is the maximum degree in the specified dimension. The entries are defined by
 *
 * \f[
 *   P_{ij} =
 *     \begin{cases}
 *       \binom{i}{j} & \text{if } j \leq i, \\
 *       0            & \text{if } j > i.
 *     \end{cases}
 * \f]
 *
 * This matrix encodes binomial coefficients in triangular form.
 *
 * @param max_degrees vector of maximum degrees per dimension
 * @param dim the dimension index to construct the matrix for
 *
 * @return a `(n+1) x (n+1)` lower-triangular Pascal matrix
 */
std::vector<std::vector<double>> get_lower_pascal_matrix(const std::vector<int> &max_degrees, int dim)
{
    int n = max_degrees[static_cast<size_t>(dim)];
    std::vector<std::vector<double>> P(static_cast<size_t>(n + 1),
                                       std::vector<double>(static_cast<size_t>(n + 1), 0.0));

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= i; ++j) {
            P[static_cast<size_t>(i)][static_cast<size_t>(j)] = static_cast<double>(binomial(i, j));
        }
    }

    return P;
}

/// Build a dimension-to-variable mapping
/**
 * Constructs a mapping from dimension index to `(variable name, interval)`
 * pairs based on the given domain.
 *
 * @param domain a mapping from variable name to its interval
 *
 * @return unordered map from dimension index (0-based) to `(var_name, interval)`
 */
std::unordered_map<int, std::pair<std::string, int_d>> build_dim_to_var_map(const var_map_i &domain)
{
    std::unordered_map<int, std::pair<std::string, int_d>> dim_map;
    int dim = 0;
    for (const auto &[var_name, interval] : domain) {
        dim_map[dim] = {var_name, interval};
        ++dim;
    }
    return dim_map;
}

/// Construct the Q-matrix for a given dimension
/**
 * Builds the Q-matrix (from Titi (2018)) associated with a given variable domain and maximum degree. This function
 * basically maps the domain to [0, 1] for the subsequent matrix compuations.
 *
 * - If the lower bound `a` of the variable interval is zero, the matrix reduces
 *   to a diagonal matrix \f$D(b)\f$ where diagonal entries are \f$b^k\f$.
 * - Otherwise, the matrix is constructed as
 *
 * \f[
 *   Q = D\!\left(\frac{b-a}{a}\right) \cdot P^T \cdot D(a),
 * \f]
 *
 * where:
 *   - `P` is the lower Pascal matrix of size `(n+1) x (n+1)`
 *   - `D(x)` is the diagonal evaluation matrix with entries \f$x^k\f$
 *   - `n = max_degrees[dim]`
 *
 * @param max_degrees vector of maximum degrees per dimension
 * @param domain variable domain mapping (variable name → interval)
 * @param dim the dimension index
 *
 * @return a `(n+1) x (n+1)` Q-matrix for the specified dimension
 *
 * @throws std::out_of_range if `dim` is invalid
 * @throws std::runtime_error if `dim` is not found in the domain map
 */
std::vector<std::vector<double>> get_q_matrix(const std::vector<int> &max_degrees, const var_map_i &domain, int dim)
{
    if (dim < 0 || dim >= static_cast<int>(max_degrees.size())) throw std::out_of_range("Dimension out of range");

    // build dim -> {var, interval} map
    std::unordered_map<int, std::pair<std::string, int_d>> dim_map = build_dim_to_var_map(domain);

    auto it = dim_map.find(dim);
    if (it == dim_map.end()) throw std::runtime_error("Dimension not found in domain map");

    double a = it->second.second.lower();
    double b = it->second.second.upper();

    if (a == 0.0) {
        return get_d_matrix(max_degrees, dim, b);
    } else {
        double t1 = (b - a) / a;
        auto D2 = get_d_matrix(max_degrees, dim, t1);
        auto P = get_lower_pascal_matrix(max_degrees, dim);
        auto D1 = get_d_matrix(max_degrees, dim, a);

        auto PT = transpose(P);
        auto temp = matmul(D2, PT);
        auto Q = matmul(temp, D1);
        return Q;
    }
}

/// Shift indices cyclically to the left
/**
 * Rotates a vector of indices one position to the left.
 *
 * Example:
 * - Input: `{2, 3, 4}`
 * - Output: `{3, 4, 2}`
 *
 * @param indices vector of indices
 *
 * @return shifted vector with elements rotated left by one
 */
std::vector<int> shift_indices(const std::vector<int> &indices)
{
    if (indices.empty()) return {};
    std::vector<int> shifted = indices;
    std::rotate(shifted.begin(), shifted.begin() + static_cast<std::ptrdiff_t>(1), shifted.end());
    return shifted;
}

/// Cycle a 2D array representation along dimensions
/**
 * Rearranges a 2D array representation of polynomial coefficients by cycling
 * the role of dimensions. This operation effectively rotates the "leading"
 * dimension used in the row index.
 *
 * If the number of dimensions is one, the input array is returned unchanged.
 *
 * @param arr 2D array of shape `(rows x cols)` representing coefficients
 * @param max_degrees vector of maximum degrees per dimension
 * @param dim the current dimension index to cycle
 *
 * @return a 2D array with cycled dimension ordering; new shape is
 *         `(max_degrees[next_dim]+1, ∏_{m≠next_dim} (max_degrees[m]+1))`
 */
std::vector<std::vector<double>> get_cycled_2d_array_vec(const std::vector<std::vector<double>> &arr,
                                                         const std::vector<int> &max_degrees, int dim)
{
    int num_rows = static_cast<int>(arr.size());
    int num_cols = static_cast<int>(arr[0].size());
    int ndim = static_cast<int>(max_degrees.size());
    if (ndim == 1) {
        return arr;
    }

    int next_dim = (dim + 1) % ndim;

    // new shape
    int new_rows = max_degrees[static_cast<size_t>(next_dim)] + 1;
    int new_cols = 1;
    for (int m = 0; m < ndim; ++m) {
        if (m != next_dim) {
            new_cols *= (max_degrees[static_cast<size_t>(m)] + 1);
        }
    }

    std::vector<std::vector<double>> cycled_arr(static_cast<size_t>(new_rows),
                                                std::vector<double>(static_cast<size_t>(new_cols), 0.0));

    // Compute index mapping
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            int i_new = j % (max_degrees[static_cast<size_t>(next_dim)] + 1);

            int stride = 1;
            for (int m = 0; m < ndim; ++m) {
                if (m != (dim % ndim) && m != next_dim) {
                    stride *= (max_degrees[static_cast<size_t>(m)] + 1);
                }
            }

            int j_new = j / (max_degrees[static_cast<size_t>(next_dim)] + 1) + i * stride;

            cycled_arr[static_cast<size_t>(i_new)][static_cast<size_t>(j_new)]
                = arr[static_cast<size_t>(i)][static_cast<size_t>(j)];
        }
    }

    return cycled_arr;
}

/// Construct the base Lambda matrix for the general box domain
/**
 * Computes the base \f$\Lambda\f$ matrix representation for a polynomial
 * expressed in the Bernstein basis over a general (non-unit) domain.
 *
 * Iteratively constructs \f$\Lambda\f$ by applying, for each dimension:
 *  - diagonal accumulation matrix \f$D_{\text{acc}}\f$
 *  - Q-matrix (interval scaling/translation)
 *  - cycling of dimensions
 *
 * Intermediate Lambda matrices are stored in a dictionary by dimension.
 *
 * @param coeffs vector of polynomial coefficients
 * @param exps exponent matrix (each row a monomial’s exponents)
 * @param domain variable domain mapping (variable name → interval)
 *
 * @return the final \f$\Lambda\f$ matrix of shape
 *         `∏ (deg_i+1) x ∏ (deg_i+1)` for all dimensions
 */
std::vector<std::vector<double>> get_titi_base_lambda_generalbox(const std::vector<double> &coeffs,
                                                                 const std::vector<std::vector<int>> &exps,
                                                                 const var_map_i &domain)
{
    unsigned int ndim = get_ndim(coeffs, exps);

    // L_dict[0] = A-matrix
    std::unordered_map<int, std::vector<std::vector<double>>> L_dict;
    L_dict[0] = get_a_matrix_vec(coeffs, exps);

    // max_degrees for later
    std::vector<int> max_degrees = get_max_degrees(exps, ndim);

    // shape_nd = (deg+1, deg+1, ...)
    std::vector<int> shape_nd;
    shape_nd.reserve(max_degrees.size());
    for (int deg : max_degrees) {
        shape_nd.push_back(deg + 1);
    }

    std::vector<std::vector<double>> D_acc;
    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> temp;
    std::vector<std::vector<double>> lam_2d_next_uncycled;
    std::vector<std::vector<double>> cycled_lam;
    for (int dim = 0; dim < static_cast<int>(ndim); ++dim) {
        // lam_2d_next_uncycled = D_acc(dim) * Q(dim) * L_dict[dim]
        D_acc = get_d_acc_matrix(max_degrees, dim);
        Q = get_q_matrix(max_degrees, domain, dim);
        temp = matmul(D_acc, Q);
        lam_2d_next_uncycled = matmul(temp, L_dict[dim]);

        // Cycle order
        shape_nd = shift_indices(shape_nd);

        cycled_lam = get_cycled_2d_array_vec(lam_2d_next_uncycled, max_degrees, dim);

        L_dict[dim + 1] = cycled_lam;
    }

    return L_dict[static_cast<int>(ndim)];
}

/// Construct the Bernstein patch for an n-dimensional polynomial (with a general box domain)
/**
 * Computes the Bernstein patch (control points) of a multivariate polynomial
 * defined on a general hyper-rectangular domain. Uses audi::get_titi_base_lambda_generalbox(const std::vector<double>
 * &coeffs, const std::vector<std::vector<int>> &exps, const var_map_i &domain) followed by repeated multiplication with
 * lower Pascal matrices and cycling of dimensions.
 *
 * @tparam T numeric type of coefficients
 *
 * @param coeffs vector of polynomial coefficients
 * @param exps exponent matrix (each row a monomial’s exponents)
 * @param domain variable domain mapping (variable name → interval)
 *
 * @return a 2D array of control points in the Bernstein basis,
 *         reshaped according to `(deg_1+1, deg_2+1, ..., deg_n+1)`
 */
template <typename T>
std::vector<std::vector<T>> get_titi_bernstein_patch_ndim_generalbox(const std::vector<T> &coeffs,
                                                                     const std::vector<std::vector<int>> &exps,
                                                                     const var_map_i &domain)
{
    // Determine number of dimensions
    unsigned int ndim = get_ndim(coeffs, exps);

    // Compute base Lambda using the generalbox function
    std::vector<std::vector<T>> L_base = get_titi_base_lambda_generalbox(coeffs, exps, domain);

    // Prepare L_dict: map of iteration -> Lambda
    std::unordered_map<int, std::vector<std::vector<T>>> L_dict;
    L_dict[0] = L_base;

    // max_degrees for later
    std::vector<int> max_degrees = get_max_degrees(exps, ndim);

    // shape_nd = (deg+1, deg+1, ...)
    std::vector<int> shape_nd;
    shape_nd.reserve(max_degrees.size());
    for (int deg : max_degrees) {
        shape_nd.push_back(deg + 1);
    }

    for (int dim = 0; dim < static_cast<int>(ndim); ++dim) {
        // Compute lam_2d_next_uncycled = P_s @ L_dict[dim]
        auto P = get_lower_pascal_matrix(max_degrees, dim);
        auto lam_2d_next_uncycled = matmul(P, L_dict[dim]);

        // Cycle order
        auto old_shape_nd = shape_nd;
        shape_nd = shift_indices(shape_nd);
        auto cycled_lam = get_cycled_2d_array_vec(lam_2d_next_uncycled, max_degrees, dim);

        L_dict[dim + 1] = cycled_lam;
    }

    return L_dict[static_cast<int>(ndim)];
}

} // namespace audi

#endif // !AUDI_TM_BOUNDING_HPP
