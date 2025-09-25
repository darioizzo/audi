// Copyright © 2018–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com),
// Sean Cowan (lambertarc@icloud.com)
//
// This file is part of the audi library.
//
// The audi library is free software: you can redistribute it and/or modify
// it under the terms of either:
//   - the GNU General Public License as published by the Free Software
//     Foundation, either version 3 of the License, or (at your option)
//     any later version, or
//   - the GNU Lesser General Public License as published by the Free
//     Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
// The audi library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License and the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser General Public License along with the audi library.
// If not, see <https://www.gnu.org/licenses/>.

#ifndef AUDI_TM_UTILS_HPP
#define AUDI_TM_UTILS_HPP

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <unordered_map>

using int_d = boost::numeric::interval<
    double, boost::numeric::interval_lib::policies<
                boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double>>,
                boost::numeric::interval_lib::checking_no_nan<double>>>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

/// Matrix multiplication
/**
 * Computes the matrix product \f$ C = A \cdot B \f$ where
 * - `A` is of shape `(m x p)`
 * - `B` is of shape `(p x n)`
 *
 * The result `C` has shape `(m x n)`.
 *
 * @param A left-hand matrix of shape `(m x p)`
 * @param B right-hand matrix of shape `(p x n)`
 *
 * @return matrix product `C` of shape `(m x n)`
 *
 * @throws std::out_of_range if `A` and `B` have incompatible dimensions
 */
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

/// Matrix transpose
/**
 * Computes the transpose of a given matrix \f$M\f$.
 *
 * - If `M` has shape `(rows x cols)`, the result has shape `(cols x rows)`.
 *
 * @param M input matrix
 *
 * @return transpose of `M`
 */
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

/// Overloaded stream insertion operator for var_map_d
/**
 * Outputs the contents of a var_map_d (mapping variable names to doubles)
 * to the given output stream in a human-readable format.
 *
 * Example output: {"x": 1.23, "y": 4.56}
 * @param os the output stream
 * @param m the var_map_d mapping variable names to double values
 *
 * @return a reference to the output stream with the formatted contents inserted
 */
std::ostream &operator<<(std::ostream &os, const var_map_d &m)
{
    os << "{";
    bool first = true;
    for (const auto &[key, value] : m) {
        if (!first) os << ", ";
        os << key << ": " << value;
        first = false;
    }
    os << "}";
    return os;
}

/// Overloaded stream insertion operator for var_map_i
/**
 * Outputs the contents of a var_map_i (mapping variable names to intervals)
 * to the given output stream in a human-readable format.
 *
 * Each interval is displayed as [lower, upper].
 * Example output: {"x":[0, 1], "y":[2, 5]}
 *
 * @param os the output stream
 * @param m the var_map_i mapping variable names to int_d intervals
 *
 * @return a reference to the output stream with the formatted contents inserted
 */
std::ostream &operator<<(std::ostream &os, const var_map_i &m)
{
    os << "{";
    bool first = true;
    for (const auto &[key, value] : m) {
        if (!first) os << ", ";
        os << key << ":" << "[" << value.lower() << ", " << value.upper() << "]";

        first = false;
    }
    os << "}";
    return os;
}

#endif
