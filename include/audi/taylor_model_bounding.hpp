#ifndef AUDI_TM_BOUNDING_HPP
#define AUDI_TM_BOUNDING_HPP

#include <Eigen/Dense>
#include <unordered_map>
#include <audi/gdual.hpp>

std::vector<std::vector<int>> generate_combinations(const std::vector<int>& limits) {
    if (limits.empty()) return {};

    int cap = *std::max_element(limits.begin(), limits.end());
    std::vector<std::vector<int>> result = {{}};

    for (int limit : limits) {
        std::vector<std::vector<int>> tmp;
        for (const auto& comb : result) {
            for (int i = 0; i <= limit; ++i) {
                auto copy = comb;
                copy.push_back(i);

                // Check sum of indices
                int s = std::accumulate(copy.begin(), copy.end(), 0);
                if (s <= cap) {
                    tmp.push_back(std::move(copy));
                }    }
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
        auto deg = tpol.degree();
        for (int i = 0; i <= static_cast<int>(deg); ++i) {
            combs.push_back({i});
        }
    } else if (ndim > 1) {
        std::vector<int> limits(ndim, static_cast<const int>(tpol.degree()));
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

int get_ndim(const std::vector<double>& coeffs, const std::vector<std::vector<int>>& exps) {
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

    for (const auto& row : exps) {
        if (row.size() != ndim) {
            throw std::invalid_argument("All exponent rows must have the same length.");
        }
    }

    return static_cast<int>(ndim);
}

std::vector<int> get_max_degrees(const std::vector<std::vector<int>>& exponents, int ndim) {
    if (exponents.empty()) {
        throw std::invalid_argument("Exponents array cannot be empty.");
    }

    if (ndim > 1) {
        size_t cols = exponents[0].size();
        std::vector<int> max_degrees(cols, 0);

        for (const auto& row : exponents) {
            if (row.size() != cols) {
                throw std::invalid_argument("All inner vectors must have the same size.");
            }
            for (size_t i = 0; i < cols; ++i) {
                max_degrees[i] = std::max(max_degrees[i], row[i]);
            }
        }
        return max_degrees;
    } 
    else if (ndim == 1) {
        int max_degree = 0;
        for (const auto& row : exponents) {
            for (int val : row) {
                max_degree = std::max(max_degree, val);
            }
        }
        return {max_degree};  // single element vector
    } 
    else {
        throw std::invalid_argument("ndim must be >= 1");
    }
}


// Eigen::MatrixXd get_a_matrix_vec()
// {
//     /// Generates the A matrix for a n-dimensional polynomial in Bernstein form.
//     /// This is a vectorized version of the original get_a_matrix function.
// }
//
// std::unordered_map<uint, Eigen::MatrixXd> taylor_model_bounding(const audi::taylor_model &tm)
// {
//     std::unordered_map<uint, Eigen::MatrixXd> L_dict;
//     L_dict[0] = get_a_matrix_vec()
// }

#endif // !AUDI_TM_BOUNDING_HPP
