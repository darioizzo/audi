#ifndef AUDI_TM_BOUNDING_HPP
#define AUDI_TM_BOUNDING_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <audi/gdual.hpp>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

std::vector<std::vector<int>> generate_combinations(const std::vector<int> &limits)
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
                if (s <= cap) {
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
        combs = generate_combinations(limits);
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
        throw std::invalid_argument("Coefficients and exponents cannot be empty.");
    }

    // number of terms must match
    if (coeffs.size() != exps.size()) {
        throw std::invalid_argument("Number of coefficients and number of exponent rows must match.");
    }

    // check all rows in exps have the same length
    size_t ndim = exps[0].size();
    if (ndim == 0) {
        throw std::invalid_argument("Exponent rows cannot be empty.");
    }

    for (const auto &row : exps) {
        if (row.size() != ndim) {
            throw std::invalid_argument("All exponent rows must have the same length.");
        }
    }

    return static_cast<uint>(ndim);
}

std::vector<int> get_max_degrees(const std::vector<std::vector<int>> &exponents, uint ndim)
{
    if (exponents.empty()) {
        throw std::invalid_argument("Exponents array cannot be empty.");
    }

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
std::vector<std::vector<T>> get_a_matrix_vec(const std::vector<T> &coeffs,
                                                  const std::vector<std::vector<int>> &exps)
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

#endif // !AUDI_TM_BOUNDING_HPP
