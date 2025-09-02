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

using int_d = boost::numeric::interval<double>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

namespace audi
{

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

// Return type: pair of (coefficients, exponents)
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

uint get_ndim(const std::vector<double> &coeffs, const std::vector<std::vector<int>> &exps)
{
    if (coeffs.empty() || exps.empty()) {
        return 0;
    }
    //     throw std::invalid_argument("Coefficients and exponents cannot be empty.");
    // } else if (!coeffs.empty() && exps.empty()) {
    // }

    // number of terms must match
    if (coeffs.size() != exps.size() && coeffs.size() != 0 && exps.size() != 0) {
        throw std::invalid_argument("Number of coefficients and number of exponent rows must match.");
    }

    size_t ndim = exps[0].size();

    // // check all rows in exps have the same length
    // if (ndim == 0) {
    //     throw std::invalid_argument("Exponent rows cannot be empty.");
    // }

    for (const auto &row : exps) {
        if (row.size() != ndim) {
            throw std::invalid_argument("All exponent rows must have the same length.");
        }
    }

    return static_cast<uint>(ndim);
}

std::vector<int> get_max_degrees(const std::vector<std::vector<int>> &exponents, uint ndim)
{
    // if (exponents.empty()) {
    //     throw std::invalid_argument("Exponents array cannot be empty.");
    // }

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

// coefficients lookup by matching multi-index in exps
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

template <typename T>
std::vector<std::vector<T>> get_a_matrix_vec(const std::vector<T> &coeffs, const std::vector<std::vector<int>> &exps)
{
    // number of dimensions
    uint ndim = get_ndim(coeffs, exps);

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

unsigned long long binom_product(const std::vector<int> &n, const std::vector<int> &k)
{

    if (n.size() != k.size()) throw std::runtime_error("Binomial lists must be of equal size.");
    unsigned long long result = 1;
    for (std::size_t i = 0; i < n.size(); ++i) {
        result *= binomial(n[i], k[i]);
    }
    return result;
}

std::vector<std::vector<double>> get_d_acc_matrix(const std::vector<int> &max_degrees, int dim)
{
    if (dim < 0 || dim >= static_cast<int>(max_degrees.size())) throw std::out_of_range("Dimension out of range");

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

// matrix multiplication helper
std::vector<std::vector<double>> matmul(const std::vector<std::vector<double>> &A,
                                        const std::vector<std::vector<double>> &B)
{
    size_t m = A.size();
    size_t n = B[0].size();
    size_t p = A[0].size(); // must equal B.size()

    std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < p; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// transpose helper
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &M)
{
    size_t rows = M.size();
    size_t cols = M[0].size();
    std::vector<std::vector<double>> T(cols, std::vector<double>(rows, 0.0));
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            T[j][i] = M[i][j];
    return T;
}

// dim_to_var map helper
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

std::vector<int> shift_indices(const std::vector<int> &indices)
{
    if (indices.empty()) return {};
    std::vector<int> shifted = indices;
    std::rotate(shifted.begin(), shifted.begin() + static_cast<std::ptrdiff_t>(1), shifted.end());
    return shifted;
}

std::vector<std::vector<double>> get_cycled_2d_array_vec(const std::vector<std::vector<double>> &arr,
                                                         const std::vector<int> &max_degrees, int dim)
{
    int num_rows = static_cast<int>(arr.size());
    int num_cols = static_cast<int>(arr[0].size());
    int ndim = static_cast<int>(max_degrees.size());

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

std::vector<std::vector<double>> get_titi_base_lambda_generalbox(const std::vector<double> &coeffs,
                                                                 const std::vector<std::vector<int>> &exps,
                                                                 const var_map_i &domain)
{
    uint ndim = get_ndim(coeffs, exps);

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
    std::vector<int> old_shape_nd;
    std::vector<std::vector<double>> cycled_lam;
    for (int dim = 0; dim < static_cast<int>(ndim); ++dim) {
        // lam_2d_next_uncycled = D_acc(dim) * Q(dim) * L_dict[dim]
        D_acc = get_d_acc_matrix(max_degrees, dim);
        Q = get_q_matrix(max_degrees, domain, dim);
        temp = matmul(D_acc, Q);
        lam_2d_next_uncycled = matmul(temp, L_dict[dim]);

        // Cycle order
        old_shape_nd = shape_nd;
        shape_nd = shift_indices(shape_nd);

        cycled_lam = get_cycled_2d_array_vec(lam_2d_next_uncycled, max_degrees, dim);

        L_dict[dim + 1] = cycled_lam;
    }

    return L_dict[static_cast<int>(ndim)];
}

template <typename T>
std::vector<std::vector<T>> get_titi_bernstein_patch_ndim_generalbox(const std::vector<T> &coeffs,
                                                                     const std::vector<std::vector<int>> &exps,
                                                                     const var_map_i &domain)
{
    // Determine number of dimensions
    uint ndim = get_ndim(coeffs, exps);

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

// For univariate functions, calculate the bernstein coefficients exhaustively (no guarantee for the
// best performance but same method as matrix one.

template <typename T>
std::vector<T> get_bernstein_coefficients(const std::vector<T> coeffs, const std::vector<std::vector<int>> exponents,
                                          uint ndim)
{

    std::vector<int> max_degrees = get_max_degrees(exponents, ndim);
    std::vector<std::vector<int>> bernstein_combs = generate_combinations(max_degrees);

    std::vector<T> b_coeffs(bernstein_combs.size(), 0.0);
    for (int i = 0; i < static_cast<int>(bernstein_combs.size()); ++i) {
        std::vector<int> b_comb = bernstein_combs[i];
        std::vector<std::vector<int>> current_coeff_combs = generate_combinations(b_comb);

        double b_coeff = 0;
        for (std::vector<int> &comb : current_coeff_combs) {
            b_coeff += get_coefficient(comb, coeffs, exponents) * static_cast<double>(binom_product(b_comb, comb))
                       / static_cast<double>(binom_product(max_degrees, comb));
        }
        b_coeffs[i] = b_coeff;
    }
    return b_coeffs;
}

} // namespace audi

#endif // !AUDI_TM_BOUNDING_HPP
