#ifndef AUDI_INVERT_MAP_HPP
#define AUDI_INVERT_MAP_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <audi/exceptions.hpp>
#include <audi/gdual.hpp>
#include <audi/io.hpp>

namespace audi
{
// The basic type of a Taylor map, that is a "square" collection of gduals
template <typename Cf, typename M>
using taylor_map = std::vector<gdual<Cf, M>>;

namespace detail
{
// This is the composition operator. It substitutes all symbol in A
// with the gduals in B. Careful that A and B do not share symbols otherwise
// unexpexted effect may arise
template <typename Cf, typename M>
inline taylor_map<Cf, M> operator&(const taylor_map<Cf, M> &A, const taylor_map<Cf, M> &B)
{
    taylor_map<Cf, M> retval;
    auto ss = A[0].get_symbol_set();
    for (decltype(A.size()) i = 0u; i < A.size(); ++i) {
        auto tmp = A[i];
        for (decltype(B.size()) j = 0u; j < B.size(); ++j) {
            tmp = tmp.subs("d" + ss[j], B[j]);
        }
        retval.push_back(tmp);
    }
    return retval;
}

// This is the sum operator.
template <typename Cf, typename M>
inline taylor_map<Cf, M> operator+(const taylor_map<Cf, M> &A, const taylor_map<Cf, M> &B)
{
    taylor_map<Cf, M> retval(B.size());
    for (decltype(B.size()) i = 0u; i < B.size(); ++i) {
        retval[i] = A[i] + B[i];
    }
    return retval;
}

// This is the diff operator.
template <typename Cf, typename M>
inline taylor_map<Cf, M> operator-(const taylor_map<Cf, M> &A, const taylor_map<Cf, M> &B)
{
    taylor_map<Cf, M> retval(B.size());
    for (decltype(B.size()) i = 0u; i < B.size(); ++i) {
        retval[i] = A[i] - B[i];
    }
    return retval;
}

// This is the "multiplication" operator.
template <typename T, typename Cf, typename M>
inline taylor_map<Cf, M> operator*(const T &c, const taylor_map<Cf, M> &A)
{
    taylor_map<Cf, M> retval(A.size());
    for (decltype(A.size()) i = 0u; i < A.size(); ++i) {
        retval[i] = A[i] * c;
    }
    return retval;
}

template <typename Cf, typename M>
inline taylor_map<Cf, M> trim(const taylor_map<Cf, M> &A, double epsilon)
{
    taylor_map<Cf, M> retval(A.size());
    for (decltype(A.size()) i = 0u; i < A.size(); ++i) {
        retval[i] = A[i].trim(epsilon);
    }
    return retval;
}

// Eigen stores indexes and sizes as signed types, while audi
// uses STL containers thus sizes and indexes are unsigned. To
// make the conversion as painless as possible this template is provided
// allowing, for example, syntax of the type M[_(i)] to adress an element
// of a vector when i is an Eigen::DenseIndex
template <typename I>
static unsigned _(I n)
{
    return static_cast<unsigned>(n);
}

} // namespace detail

/// Inversion of Taylor maps
/**
 * Consider a system of non linear equations in \f$n\f$ variables:
 * \f[
 * \left\{
 * \begin{array}{l}
 * f_1(x_1, x_2, .., x_n) = 0 \\
 * ...                        \\
 * f_n(x_1, x_2, .., x_n) = 0 \\
 * \end{array}
 * \right.
 * \f]
 * and represent each of it via a Taylor polynomial truncated at order \f$m\f$ (i.e. \f$\in \mathcal P_m\f$) around
 * some point \f$\overline {\mathbf x}\f$:
 * \f[
 * \left\{
 * \begin{array}{l}
 * T_{f_1}(dx_1, dx_2, .., dx_n) = 0 \\
 * ...                        \\
 * T_{f_n}(dx_1, dx_2, .., dx_n) = 0 \\
 * \end{array}
 * \right.
 * \f]
 * The above system of, now, polynomial equations can be written as \f$C + M + N = C + \mathcal M = 0\f$, where we have
 * explicitly indicated the constant term \f$C\f$, the linear term \f$M\f$ and the rest \f$N\f$. The symbol
 * \f$\mathcal M\f$ indicates what we call a Taylor map and its inverse is found by this function by performing
 * some Picard iterations (fixed point) over the following simple formal equalities:
 *
 * \f[
 * \mathcal M \circ \mathcal M^{-1} = \mathcal I \rightarrow  (M + N) \circ \mathcal M^{-1}  = \mathcal I \rightarrow M
 * \circ \mathcal M^{-1} + N \circ \mathcal M^{-1} = \mathcal I \f] and hence: \f[ \mathcal M^{-1} = M^{-1} \circ
 * \left(\mathcal I - N \circ \mathcal M^{-1}\right) \f] which can be used to make a fixed point iterative scheme once
 * the inverse of the linear part \f$M^{-1}\f$ is found using any linear algebraic tool. Note that the inverse is
 * guaranteed to exist if the linear part is invertible.
 *
 * @param map_in The input map represented as an <tt> std::vector </tt> of <tt> gduals </tt>. They can have a constant
 * coefficient which will be neglected
 * @param verbose when true some output is shown during the Picard iterations
 *
 * @return The inverse map \f$\mathcal M^{-1}\f$ represented as an <tt> std::vector </tt> of <tt> gduals </tt> with
 * symbol set [p0, p1, .. pn]
 *
 * @throws std::invalid_argument if \p map_in is empty, if its linear part is not invertable, if the symbol set of the
 * map components are not all equal, if the map is not square (i.e. it is of size n with n symbols) or if the order of
 * the gduals are not all the same.
 *
 *  * Example usage:
 *
 * @code{.cpp}
 * gdual<double> x(0., "x1", 9);
 * auto sinx = sin(x);
 * auto asinx = asin(x);
 * auto inv_sinx = invert_map({sinx});
 * std::cout << "sin x:\t\t" << sinx << std::endl;
 * std::cout << "asin x:\t\t" << asinx << std::endl;
 * std::cout << "inv sinx:\t" << inv_sinx[0] << std::endl;
 * @endcode
 *
 * The code above will produce the output:
 *
 * @code{.txt}
 * sin x:        dx1+2.75573e-06*dx1**9-0.166667*dx1**3+0.00833333*dx1**5-0.000198413*dx1**7
 * asin x:       dx1+0.0303819*dx1**9+0.166667*dx1**3+0.075*dx1**5+0.0446429*dx1**7
 * inv sinx:     dp0+0.0303819*dp0**9+0.166667*dp0**3+0.075*dp0**5+0.0446429*dp0**7
 * @endcode
 */
template <typename Cf, typename Monomial>
inline taylor_map<Cf, Monomial> invert_map(const taylor_map<Cf, Monomial> &map_in, bool verbose = false)
{
    // To find the overloaded operators
    using namespace detail;
    // Some preliminary definitions and consistency checks
    // The map cannot be empty
    auto map_size = map_in.size();
    if (map_size == 0) {
        audi_throw(std::invalid_argument, "Impossible to invert the Tayor Map as it is empty!");
    }
    // The map has to contain equal order gduals. (We could allow for different order implementing some promotion)
    auto map_order = map_in[0].get_order();
    if (!std::all_of(map_in.cbegin(), map_in.cend(),
                     [&map_order](const gdual<Cf, Monomial> &el) { return (el.get_order() == map_order); })) {
        audi_throw(std::invalid_argument,
                   "Impossible to invert the Tayor Map as the order of its components are mismatching.");
    }
    // The map needs to contain gdual having the same symbols. (We could allow different symbols and take the union)
    auto ss = map_in[0].get_symbol_set();
    if (!std::all_of(map_in.cbegin(), map_in.cend(),
                     [&ss](const gdual<Cf, Monomial> &el) { return (el.get_symbol_set() == ss); })) {
        audi_throw(std::invalid_argument,
                   "Impossible to invert the Tayor Map as the symbol sets of its components are mismatching.");
    }
    // The map must be "square"
    if (ss.size() != map_size) {
        audi_throw(std::invalid_argument,
                   "Impossible to invert the Tayor Map as the number of symbols is not the same as the map size.");
    }
    // We add the differential symbol to the variable symbols
    for (auto &s : ss) {
        s = "d" + s;
    }
    // We create the symbol set of the inverted map
    std::vector<std::string> ss_inv;
    for (decltype(map_size) i = 0u; i < map_size; ++i) {
        ss_inv.push_back("dp" + std::to_string(i));
    }
    /// ---------------------- Algorithm Start ---------------------------------------///
    // We decompose the map as map = M + N (linear plus non linear part). The inverse of the linear part will be Minv
    taylor_map<Cf, Monomial> M(map_size, gdual<Cf, Monomial>(0)), Minv(map_size, gdual<Cf, Monomial>(0)),
        N(map_size, gdual<Cf, Monomial>(0)), I(map_size, gdual<Cf, Monomial>(0)),
        Minv2(map_size, gdual<Cf, Monomial>(0));
    taylor_map<Cf, Monomial> retval;

    // We construct the identity map [dp0, dp1, dp2, ... ]
    for (decltype(map_size) i = 0u; i < map_size; ++i) {
        I[i] = gdual<Cf, Monomial>(0., "p" + std::to_string(i), map_order);
        I[i].extend_symbol_set(ss_inv);
    }

    // Populate M, the linear part of the input map
    for (decltype(map_size) i = 0u; i < map_size; ++i) {
        M[i] = map_in[i].extract_terms(1);
        M[i].extend_symbol_set(ss);
    }

    // Populate N, the non linear part of the input map
    for (decltype(map_size) i = 0u; i < map_size; ++i) {
        for (decltype(map_order) j = 2u; j <= map_order; ++j) {
            N[i] = N[i] + map_in[i].extract_terms(j);
        }
        N[i].extend_symbol_set(ss);
    }

    // Invert M (using eigen)
    using namespace Eigen;
    MatrixXd mat(map_size, map_size);
    for (decltype(mat.rows()) i = 0u; i < mat.rows(); ++i) {
        for (decltype(mat.rows()) j = 0u; j < mat.cols(); ++j) {
            std::vector<unsigned> monomial(map_size, 0u);
            monomial[_(j)] = 1u;
            mat(i, j) = M[_(i)].find_cf(monomial);
        }
    }
    auto det = mat.determinant();
    if (std::abs(det) < 1e-10) {
        audi_throw(std::invalid_argument,
                   "The map you are trying to invert has no inverse, the determinant of the linear part is: "
                       + std::to_string(det));
    }
    auto invm = mat.inverse();

    if (verbose) {
        print("\nLinear part inverted successfully, determinant was: " + std::to_string(det), "\n");
    }

    // Populate Minv and Minv2 the inverse of the linear part (they are the same map, only with different
    // symbols to avoid mess when subs symbols later)
    for (decltype(invm.rows()) i = 0; i < invm.rows(); ++i) {
        for (decltype(invm.cols()) j = 0; j < invm.cols(); ++j) {
            gdual<Cf, Monomial> dummy(0., "p" + std::to_string(j), map_order);
            gdual<Cf, Monomial> dummy2(0., "pdummy" + std::to_string(j), map_order);
            dummy *= invm(i, j);
            dummy2 *= invm(i, j);
            // this must have the same symbols of the identity and the final symbols of the onverse map returned
            Minv[_(i)] += dummy;
            // this must have different symbols as to avoid mess when calling subs (in the & operator)
            Minv2[_(i)] += dummy2;
        }
    }

    // Here is where the magic happens!!!
    // Picard Iterations (each iteration will increase the order by one).
    // See http://bt.pa.msu.edu/pub/papers/AIEP108book/AIEP108book.pdf page 102 for the formal derivation
    retval = Minv;
    if (verbose) {
        print("\nAt Order 1: ", retval, "\n");
    }
    for (decltype(map_order) i = 1u; i < map_order; ++i) {
        retval = Minv2 & (I - (N & retval));
        if (verbose) {
            print("\nAt Order " + std::to_string(i + 1) + ": ", retval, "\n");
        }
    }

    return retval;
}

} // end of namespace audi

#endif
