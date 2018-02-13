#ifndef AUDI_INVERT_MAP_HPP
#define AUDI_INVERT_MAP_HPP

#include <Eigen/Dense>
#include <stdexcept>
#include <string>
#include <vector>

#include <audi/audi.hpp>
#include <audi/exceptions.hpp>
#include <audi/io.hpp>

namespace audi
{
using taylor_map = std::vector<gdual_d>;

namespace detail
{
// This is the composition operator.
taylor_map operator&(const taylor_map &A, const taylor_map &B)
{
    taylor_map retval;
    auto symbols = B[0].get_symbol_set();
    for (decltype(B.size()) i = 0u; i < B.size(); ++i) {
        auto tmp = B[i];
        for (decltype(B.size()) j = 0u; j < B.size(); ++j) {
            tmp = tmp.subs("d" + symbols[j], A[j]);
        }
        retval.push_back(tmp);
    }
    return retval;
}

// This is the sum operator.
taylor_map operator+(const taylor_map &A, const taylor_map &B)
{
    taylor_map retval(B.size());
    for (decltype(B.size()) i = 0u; i < B.size(); ++i) {
        retval[i] = A[i] + B[i];
    }
    return retval;
}

// This is the multiplication for a scalar operator.
taylor_map operator*(double c, const taylor_map &A)
{
    taylor_map retval(A.size());
    for (decltype(A.size()) i = 0u; i < A.size(); ++i) {
        retval[i] = A[i] * c;
    }
    return retval;
}
} // ends detail

// An Introduction to Beam Physics pag. 138-139
taylor_map invert_map(const taylor_map &map_in)
{
    using namespace detail; // To find the overloaded operators
    // TODO CHECKS: at least order 1, same symbols, same order
    // Some preliminary definitions
    auto map_size = map_in.size();
    auto map_order = map_in[0].get_order();

    // We write the map, at order n, as M = M1 + Mn. n_order_map contains M1 and then Mn for n=2...map_order
    std::vector<taylor_map> n_order_map(map_order);
    // We write the inverse, at order n as M* = M1* + Mn*. n_order_map_inv contains M1* and then Mn* for n=2...map_order
    taylor_map dummy1(map_size, gdual_d(0.));
    std::vector<taylor_map> n_order_map_inv(map_order, dummy1);

    // Decompose map_in into n_order_map
    for (decltype(map_order) j = 0u; j < map_order; ++j) {
        for (decltype(map_size) i = 0u; i < map_size; ++i) {
            n_order_map[j].push_back(map_in[i].extract_terms(j + 1));
        }
        if (j > 1) {
            n_order_map[j] = n_order_map[j] + n_order_map[j - 1];
        }
    }

    // Invert the linear part
    using namespace Eigen;
    MatrixXd mat(map_size, map_size);
    for (unsigned i = 0; i < mat.rows(); ++i) {
        for (unsigned j = 0; j < mat.cols(); ++j) {
            std::vector<unsigned> monomial(map_size, 0u);
            monomial[j] = 1u;
            mat(i, j) = n_order_map[0][i].find_cf(monomial);
        }
    }
    auto det = mat.determinant();
    if (std::abs(det) < 1e-8) {
        audi_throw(std::invalid_argument, "The map you are trying to invert has a non ivertable linear part");
    }
    auto invm = mat.inverse();

    // Populate n_order_map_inv[0]
    for (decltype(invm.rows()) i = 0; i < invm.rows(); ++i) {
        for (decltype(invm.cols()) j = 0; j < invm.cols(); ++j) {
            // note that the order must be already the final one as during iterations terms would otherwise disappear
            gdual_d dummy(0., "p" + std::to_string(j), map_order);
            dummy *= invm(i, j);
            n_order_map_inv[0][i] += dummy;
        }
    }

    taylor_map A, B, BC, ABC;
    // Perform the main iterations to build M* at the final order
    for (auto n = 2u; n <= map_order; ++n) {
        unsigned m = n / 2;
        taylor_map C = n_order_map_inv[0];
        if (m > 1) {
            C = C + n_order_map_inv[m - 1];
        }
        A = n_order_map_inv[0];
        B = n_order_map[n - 1];
        BC = B & C;
        ABC = -1 * (A & BC);
        n_order_map_inv[n - 1] = ABC;
    }
    // TEST
    taylor_map inv2 = A + ABC;
    taylor_map dir2 = n_order_map[0] + n_order_map[1];
    print(dir2 & inv2, "\n");
    return n_order_map_inv[map_order - 1];
}

} // end of namespace audi

#endif
