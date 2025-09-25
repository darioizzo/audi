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

#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits> // For std::enable_if, std::is_same, etc.
#include <utility>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <obake/key/key_degree.hpp>
#include <obake/math/degree.hpp>
#include <obake/math/diff.hpp>
#include <obake/math/evaluate.hpp>
#include <obake/math/integrate.hpp>
#include <obake/math/subs.hpp>
#include <obake/math/trim.hpp>
#include <obake/math/truncate_degree.hpp>
#include <obake/polynomials/d_packed_monomial.hpp>
#include <obake/polynomials/polynomial.hpp>
#include <obake/series.hpp>
#include <obake/symbols.hpp>

#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>            //for audi::abs
#include <audi/detail/overloads_vectorized.hpp> //for audi::abs

/// Root namespace for AuDi symbols
namespace audi
{

/// Generalized dual number class.
/**
 * This class represents a generalized dual number, in a nutshell, a truncated multivariate Taylor polynomial.
 * Using the multi-index notation, a generalized dual number may be written as:
 * \f[
 * T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)
 * \f]
 *
 * and thus depends on the order \f$m\f$ as well as on the expansion point \f$\mathbf a\f$. All arithmetic
 * operators +,*,/,- are overloaded so that the Taylor expansion of arithmetic computations is obtained.
 * A basic use case, where \f$m = 2\f$, \f$\mathbf a = [1.2, -0.1]\f$ and \f$f = \frac{x1+x2}{x1-x2}\f$ is thus:
 *
 * @code
 * gdual<double> x1(1.2, "x1", 2);
 * gdual<double> x2(-0.1, "x2", 2);
 * std::cout << (x1+x2) / (x1-x2) << "\n";
 * @endcode
 *
 * resulting in the output:
 *
 * @code
 * 0.118343*dx1+0.846154+1.42012*dx2+1.0924*dx2**2-0.0910332*dx1**2-1.00137*dx1*dx2
 * @endcode
 *
 *
 * Integration and differentiation are also implemented so that the generalized dual computations are
 * formally made in a differential algebra.
 *
 * @note The class can be instantiated with any type that is suitable to be a coefficient in a obake polynomial.
 *
 */
template <typename Cf, typename Monomial = obake::d_packed_monomial<std::uint64_t, 8>>
class gdual
{
    // Static checks.
    static_assert(
        obake::is_cf_v<Cf> && obake::is_differentiable_v<Cf>,
        "A gdual must be constructed from a coefficient satisfying obake's conditions is_cf and is_differentiable.");

public:
    using cf_type = Cf;
    using key_type = Monomial;
    using p_type = obake::polynomial<key_type, cf_type>;

private:
    // We enable the overloads of the +,-,*,/ operators only in the following cases:
    // - at least one operand is a gdual,
    // - the other operand, if not gdual<Cf>, must be double, int or unsigned int or Cf
    template <typename T, typename U>
    using gdual_if_enabled = typename std::enable_if<
        (std::is_same<T, gdual<Cf, Monomial>>::value && std::is_same<U, gdual<Cf, Monomial>>::value)
            || (std::is_same<T, gdual<Cf, Monomial>>::value && std::is_same<U, Cf>::value)
            || (std::is_same<T, gdual<Cf, Monomial>>::value && std::is_same<U, double>::value)
            || (std::is_same<T, gdual<Cf, Monomial>>::value && std::is_same<U, int>::value)
            || (std::is_same<T, gdual<Cf, Monomial>>::value && std::is_same<U, unsigned int>::value)
            || (std::is_same<U, gdual<Cf, Monomial>>::value && std::is_same<T, Cf>::value)
            || (std::is_same<U, gdual<Cf, Monomial>>::value && std::is_same<T, double>::value)
            || (std::is_same<U, gdual<Cf, Monomial>>::value && std::is_same<T, int>::value)
            || (std::is_same<U, gdual<Cf, Monomial>>::value && std::is_same<T, unsigned int>::value),
        gdual>::type;

    // Enable the generic ctor only if T is not a gdual<Cf,Monomial> (after removing
    // const/reference qualifiers).
    template <typename T>
    using generic_ctor_enabler = enable_if_t<
        !std::is_same<gdual<Cf, Monomial>, decay_t<T>>::value && std::is_constructible<p_type, decay_t<T>>::value, int>;

    void check_order() const
    {
        if (m_order >= std::numeric_limits<unsigned int>::max() - 10u) {
            throw std::invalid_argument("polynomial truncation order is too large");
        }
    }

    void check_var_name(const std::string &name) const
    {
        if (name.at(0) == 'd') {
            throw std::invalid_argument("symbol names cannot start with the letter d");
        }
    }

    void check_var_name_has_d(const std::string &name) const
    {
        if (name.at(0) != 'd') {
            throw std::invalid_argument("symbol variations must start with the letter d");
        }
    }

    // A private constructor to move-initialise a gdual from a polynomial. Used
    // in the implementation of the operators.
    explicit gdual(p_type &&p, unsigned int order) : m_p(std::move(p)), m_order(order) {}
    // A private constructor used in the implementation of the operators.
    template <typename T>
    explicit gdual(T value, unsigned int order) : m_p(Cf(value)), m_order(order)
    {
    }

    // Basic overloads for the addition
    static gdual add(const gdual &d1, const gdual &d2)
    {
        return gdual(d1.m_p + d2.m_p, std::max(d1.get_order(), d2.get_order()));
    }

    template <typename T>
    static gdual add(const T &d1, const gdual &d2)
    {
        return gdual(d1 + d2.m_p, d2.get_order());
    }

    template <typename T>
    static gdual add(const gdual &d1, const T &d2)
    {
        return gdual(d2 + d1.m_p, d1.get_order());
    }

    // Basic overloads for the subtraction
    static gdual sub(const gdual &d1, const gdual &d2)
    {
        return gdual(d1.m_p - d2.m_p, std::max(d1.get_order(), d2.get_order()));
    }

    template <typename T>
    static gdual sub(const T &d1, const gdual &d2)
    {
        return gdual(d1 - d2.m_p, d2.get_order());
    }

    template <typename T>
    static gdual sub(const gdual &d1, const T &d2)
    {
        return gdual(d1.m_p - d2, d1.get_order());
    }

    // Dual * dual.
    static gdual mul(const gdual &d1, const gdual &d2)
    {
        const unsigned int order = std::max(d1.get_order(), d2.get_order());
        return gdual(truncated_mul(d1.m_p, d2.m_p, order), order);
    }

    template <typename T>
    static gdual mul(const T &d1, const gdual &d2)
    {
        return gdual(d1, d2.get_order()) * d2;
    }

    template <typename T>
    static gdual mul(const gdual &d1, const T &d2)
    {
        return d1 * gdual(d2, d1.get_order());
    }

    // Basic overloads for the division
    static gdual div(const gdual &d1, const gdual &d2)
    {
        gdual retval(1);
        Cf fatt(-1);
        auto p0 = d2.constant_cf();
        auto phat = (d2 - p0);
        phat = phat / p0;
        gdual tmp(phat);

        retval = retval - phat;
        for (auto i = 2u; i <= d2.m_order; ++i) {
            fatt = -1. * fatt;
            phat *= tmp;
            retval = retval + fatt * phat;
        }

        return (d1 * retval) / p0;
    }

    template <typename T>
    static gdual div(const T &d1, const gdual &d2)
    {
        gdual retval(1.);
        Cf fatt(-1.);
        auto p0 = d2.constant_cf();

        auto phat = (d2 - p0);
        phat = phat / p0;
        gdual tmp(phat);

        retval = retval - phat;
        for (auto i = 2u; i <= d2.m_order; ++i) {
            fatt = -1. * fatt;
            phat *= tmp;
            retval = retval + fatt * phat;
        }
        return (d1 * retval) / p0;
    }

    template <typename T>
    static gdual div(const gdual &d1, const T &d2)
    {
        // Note: we create the reciprocal of d2 constructing a Cf otherwise we could loose precision, e.g. when T is
        // double and Cf is real128
        return d1 * (1. / Cf(d2));
    }

    p_type m_p;
    unsigned int m_order;

public:
    /// Defaulted copy constructor.
    gdual(const gdual &) = default;
    /// Defaulted move constructor.
    gdual(gdual &&) = default;
    /// Default constuctor
    gdual() : m_p(0.), m_order(0u) {}
    /// Destructor (contains a sanity check)
    ~gdual()
    {
        assert(obake::degree(m_p) >= 0);
        assert(static_cast<unsigned>(obake::degree(m_p)) <= m_order);
    }

    template <typename T, generic_ctor_enabler<T> = 0>
    explicit gdual(const T &value) : m_p(), m_order(0u)
    {
        m_p += Cf(value);
    }

    template <typename T>
    explicit gdual(const std::initializer_list<T> &value) : m_p(value), m_order(0u)
    {
    }

    template <typename T, generic_ctor_enabler<T> = 0>
    explicit gdual(const T &value, const std::string &symbol, unsigned order) : m_p(), m_order(order)
    {
        check_order();
        check_var_name(symbol);
        if (order == 0) {
            extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
        } else {
            m_p = obake::make_polynomials<p_type>(std::string("d") + symbol)[0];
        }
        m_p += Cf(value);
    }

    template <typename T>
    explicit gdual(const std::initializer_list<T> &value, const std::string &symbol, unsigned order)
        : m_p(), m_order(order)
    {
        check_order();
        check_var_name(symbol);
        if (order == 0) {
            extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
        } else {
            m_p = obake::make_polynomials<p_type>(std::string("d") + symbol)[0];
        }
        m_p += Cf(value);
    }

    /// Defaulted assignment operator
    gdual &operator=(const gdual &) = default;
    /// Defaulted assignment operator
    gdual &operator=(gdual &&) = default;

    /// Gets the symbol set size
    /**
     * Returns the size of the symbol set.
     *
     * @return the size of the symbol set.
     */
    auto get_symbol_set_size() const
    {
        return m_p.get_symbol_set().size();
    }

    /// Returns the symbol set
    /**
     * Returns the symbol set (not the differentials)
     *
     * @return the symbol set
     */
    std::vector<std::string> get_symbol_set() const
    {
        std::vector<std::string> retval;
        for (auto s : m_p.get_symbol_set()) {
            std::string tmp(s);
            retval.push_back(tmp.erase(0, 1));
        }
        return retval;
    }

    /// Extends the symbol set
    /**
     * Adds some symbolic variables to the current polynomial
     * This is useful in situations where some variable \f$ dx\f$ does not
     * appear in the final polynomial but we still want to
     * treat the Taylor expansion as a function of \f$ dx\f$ too (for example
     * when extracting the relative coefficient)
     *
     * @param sym_vars list of symbolic names. It must contain all symbolic names of
     * the current polynomial. It may contain more. All symbols must start with the letter "d".
     *
     * @throws std::invalid_argument:
     * - if any symbol in \p sym_vars does not start with the letter "d"
     */
    void extend_symbol_set(const std::vector<std::string> &sym_vars)
    {
        for (auto sym_var : sym_vars) {
            check_var_name_has_d(sym_var);
        }
        m_p = obake::add_symbols(m_p, obake::symbol_set(sym_vars.begin(), sym_vars.end()));
    }

    /// Integration
    /**
     * Performs the integration of the gdual with respect to the symbol.
     *
     * \note If the \p var_name differential is not in the symbol set, then it is added.
     *
     * \note Information may be lost as the truncation order is preserved.
     *
     * @param var_name Symbol name (cannot start with "d").
     *
     * @throws std::invalid_argument:
     * - if \p var_name starts with the letter "d" (this avoid creating confusing names for symbol's differentials)
     *
     * @throws unspecified any exception thrown by:
     * - obake::truncate_degree,
     * - obake::integrate
     */
    template <enable_if_t<obake::is_differentiable<Cf>::value, int> = 0>
    gdual integrate(const std::string &var_name)
    {
        check_var_name(var_name);
        auto new_p = obake::integrate(m_p, "d" + var_name);
        obake::truncate_degree(new_p, static_cast<int>(m_order));
        return gdual(std::move(new_p), m_order);
    }

    /// Partial derivative
    /**
     * Performs the partial derivative of the gdual with respect to the symbol
     *
     * \note If the \p var_name differential is not in the symbol set, then it is added.
     *
     * @param var_name Symbol name (cannot start with "d").
     *
     * @throws std::invalid_argument:
     * - if \p symbol starts with the letter "d" (this avoid creating confusing names for symbol's differentials)
     *
     * @throws unspecified any exception thrown by:
     * - obake::diff,
     */
    template <enable_if_t<obake::is_differentiable<Cf>::value, int> = 0>
    gdual partial(const std::string &var_name)
    {
        check_var_name(var_name);
        auto new_p = obake::diff(m_p, "d" + var_name);
        return gdual(std::move(new_p), m_order);
    }

    /// Substitute symbol with value
    /**
     * Substitute the symbol \p sym with the value \p val
     *
     * @throws unspecified any exception thrown by:
     * - obake::subs,
     */
    template <typename T>
    gdual subs(const std::string &sym, const T &val) const
    {
        auto new_p = obake::subs(m_p, obake::symbol_map<Cf>{{sym, Cf(val)}});
        return gdual(std::move(new_p), m_order);
    }

    /// Substitute a symbol with a gdual
    /**
     * Substitute the symbol \p sym with the gdual \p val
     * The returned gdual has the same order, its symbol set will be
     * expanded with the symbol set of \p val and then trimmed of the substituded symbol and
     * truncated to the order of the gdual.
     *
     * @param sym Symbol to be substituded (i.e. "dx")
     * @param val The gdual to substitute \p sym with.
     *
     * @throws unspecified any exception thrown by:
     * - obake::subs,
     * - obake::trim
     * - obake::truncate_degree
     */
    gdual subs(const std::string &sym, const gdual &val) const
    {
        auto new_p = obake::subs(m_p, obake::symbol_map<decltype(val.m_p)>{{sym, val.m_p}});
        obake::truncate_degree(new_p, static_cast<decltype(obake::degree(m_p))>(m_order));
        new_p = obake::trim(new_p);
        return gdual(std::move(new_p), m_order);
    }

    /// Trims small coefficients
    /**
     * All coefficients smaller than a set magnitude \p epsilon are set to zero
     *
     * @param epsilon Tolerance used to trim coefficients.
     *
     * @return The new trimmed gdual.
     */
    gdual trim(double epsilon) const
    {
        if (epsilon < 0) {
            throw std::invalid_argument(
                "When trimming a gdual the trim tolerance must be positive, you seem to have used a negative value: "
                + std::to_string(epsilon));
        }
        auto new_p = m_p;
        obake::filter(new_p, [epsilon](const auto &p) { return !(audi::abs(p.second) < epsilon); });
        return gdual(std::move(new_p), m_order);
    }

    /// Extracts all terms of some order
    /**
     * Extracts all the terms with the given \p order into a new gdual
     *
     * @param order Order of the terms to be extracted
     *
     * @return A new gdual containing the terms extracted, but preserving the order of the original gdual.
     *
     * @throws std::invalid_argument
     * - if the \p order is higher than the gdual order
     */
    gdual extract_terms(unsigned order) const
    {
        if (order > m_order) {
            throw std::invalid_argument("requested order is beyond the truncation order.");
        }
        auto res = m_p;
        obake::filter(res, [order, &ss = m_p.get_symbol_set()](const auto &p) {
            return static_cast<unsigned>(obake::key_degree(p.first, ss)) == order;
        });
        return gdual(std::move(res), order);
    }

    /// Evaluates the Taylor polynomial
    /**
     * Evaluates the Taylor polynomial using the values in \p dict for all the
     * differentials (variables variations)
     *
     * @throws unspecified any exception thrown by:
     * - obake::evaluate,
     */
    auto evaluate(const obake::symbol_map<double> &dict) const
    {
        return obake::evaluate(m_p, dict);
    }

    /// Current degree
    /**
     * Returns the current degree of the polynomial represented as an audi::gdual.
     * This may be different from the truncation order \f$m\f$ and,
     * in particular will be smaller or equal.
     *
     * @return the current degree of the polynomial
     */
    auto degree() const
    {
        return obake::degree(m_p);
    }

    /// Getter for the truncation order
    /**
     * Returns the truncation order of the underlying \f$\mathcal P_{n,m}\f$ algebra.
     *
     * @return the truncation order
     */
    unsigned int get_order() const
    {
        return m_order;
    }

    /// Finds the coefficient of a particular monomial
    /**
     * Returns the coefficient of the monomial specified in the container \p c
     *
     * \note The container contains the exponents of the requested monomial. In a three
     * variable polynomial with "x", "y" and "z" as symbols, [1, 3, 2] would denote
     * the term x y^3 z^2.
     *
     * \note Alphabetical order is used to order the symbol set and thus specify
     * the coefficients.
     *
     * This method will first construct a term with zero coefficient
     * and key initialised from the begin/end iterators of c and the
     * symbol set of this, and it will then try to locate the term inside
     * this. If the term is found, its coefficient will be returned.
     * Otherwise, a coefficient initialised from 0 will be returned.
     *
     * @return the coefficient
     *
     * @throws std::invalid_argument:
     * - if the requested coefficient is beyond the truncation order
     */
    template <typename T>
    Cf find_cf(const T &c) const
    {
        if (std::accumulate(c.begin(), c.end(), 0u) > m_order) {
            throw std::invalid_argument("requested coefficient is beyond the truncation order.");
        }
        if (c.size() != this->get_symbol_set_size()) {
            throw std::invalid_argument(
                "requested monomial does not exist, check the length of the input with respect to the symbol set size");
        }
        const auto it = m_p.find(key_type(c.begin(), c.end()));
        if (it == m_p.end()) {
            return Cf(0);
        } else {
            return it->second;
        }
    }

    /// Finds the coefficient of a particular monomial
    /**
     * Returns the coefficient of the monomial specified in the
     * initializer list \p l
     *
     * \note This method is identical to the other overload with the same name, and it is provided for convenience.
     *
     * @return the coefficient
     *
     * @throws std::invalid_argument:
     * - if the requested coefficient is beyond the truncation order
     */
    template <typename T>
    Cf find_cf(std::initializer_list<T> l) const
    {
        if (std::accumulate(l.begin(), l.end(), 0u) > m_order) {
            throw std::invalid_argument("requested coefficient is beyond the truncation order.");
        }
        if (l.size() != this->get_symbol_set_size()) {
            throw std::invalid_argument(
                "requested monomial does not exist, check the length of the input with respect to the symbol set size");
        }
        const auto it = m_p.find(key_type(l.begin(), l.end()));
        if (it == m_p.end()) {
            return Cf(0);
        } else {
            return it->second;
        }
    }

    /// Gets the derivative value
    /**
     * Returns the (mixed) derivative value of order specified
     * by the container \p c
     *
     * \note The container contains the order requested. In a three
     * variable polynomial with "x", "y" and "z" as symbols, [1, 3, 2] would denote
     * the sixth order derivative \f$ \frac{d^6}{dxdy^3dz^2}\f$.
     *
     * \note No computations are made at this points as all derivatives are already
     * contained in the Taylor expansion
     *
     * @return the value of the derivative
     *
     * @throws std::invalid_argument:
     * - if the requested coefficient is beyond the truncation order
     */
    template <typename T>
    Cf get_derivative(const T &c) const
    {
        double cumfact = 1.;
        for (auto i = c.begin(); i < c.end(); ++i) {
            cumfact *= boost::math::factorial<double>(static_cast<unsigned>(*i));
        }
        return this->find_cf(c) * cumfact;
    }

    /// Gets the derivative value
    /**
     * Returns the (mixed) derivative value of order specified
     * by the initializer list \p l
     *
     * \note This method is identical to the other overload with the same name,
     * and it is provided for convenience.
     *
     * \note No computations are made at this points as all derivatives are already
     * contained in the Taylor expansion
     *
     * @return the value of the derivative
     *
     * @throws std::invalid_argument:
     * - if the requested coefficient is beyond the truncation order
     */
    template <typename T>
    Cf get_derivative(std::initializer_list<T> l) const
    {
        double cumfact = 1.;
        for (auto i = l.begin(); i < l.end(); ++i) {
            cumfact *= boost::math::factorial<double>((unsigned)(*i));
        }
        return this->find_cf(l) * cumfact;
    }

    /// Gets the derivative value
    /**
     * Returns the (mixed) derivative value of order specified
     * by the unordered_map \p dict
     *
     * \note To get the following derivative: \f$ \frac{d^6}{dxdy^3dz^2}\f$
     * the input should be {{"dx", 1u},{"dy",3u},{"dz",2u}}
     *
     * \note The current implementation call internally the other templated
     * implementations. When obake will implement the sparse monomial
     * this will change and be more efficient.
     *
     * @return the value of the derivative
     *
     * @throws unspecified all exceptions thrown by the templated version call.
     * @throws std::invalid_argument: if one of the symbols is not found in the expression
     */
    auto get_derivative(const std::unordered_map<std::string, unsigned> &dict) const
    {
        const auto &ss = m_p.get_symbol_set();
        std::vector<unsigned> coeff(ss.size(), 0);
        for (const auto &entry : dict) {
            auto idx = ss.index_of(ss.find(entry.first));
            if (idx == ss.size()) {
                return cf_type(0.);
            }
            coeff[idx] = entry.second;
        }
        return get_derivative(coeff);
    }

    /// Finds the constant coefficient
    /**
     * Returns the coefficient of the of the constant part of the polynomial
     * so that if \f$T_{f} = f_0 + \hat f\f$, \f$f_0\f$ is returned.
     *
     * @return the coefficient
     */
    Cf constant_cf() const
    {
        using v_size_type = std::vector<int>::size_type;
        return find_cf(std::vector<int>(boost::numeric_cast<v_size_type>(get_symbol_set_size()), 0));
    }

    /// Determines if a gdual is zero within tolerance
    /**
     *
     * @return true if all coefficients of the gdual are zero within a tolerance \p tol
     */
    bool is_zero(double tol) const
    {
        for (const auto &t : m_p) {
            if (audi::abs(t.second) > tol) // call to audi abs has precedence
            {
                return false;
            }
        }
        return true;
    }

    /// Info on the gdual
    /**
     * Will assemble a string containing information on the gdual.
     * It uses the obake streaming operator for the type obake::series. Refer to that
     * documentation for further details 
     */
    std::string info() const
    {
        std::ostringstream ss;
        ss << m_p << "\n";
        ss << "Order: " << get_order() << "\n";
        ss << "Degree: " << degree() << "\n";
        return ss.str();
    }

    /// Overloaded stream operator
    /**
     * Will direct to stream a human-readable representation of the generalized dual number.
     *
     * \note The print order of the terms will be undefined.
     *
     * @param[in,out] os target stream.
     * @param d audi::gdual argument.
     *
     * @return reference to \p os.
     *
     */
    friend std::ostream &operator<<(std::ostream &os, const gdual &d)
    {
        obake::detail::series_stream_terms_impl<false>(os, d.m_p);
        return os;
    }

    /** @name Comparison operators
     *
     */
    //@{

    /// Overloaded Equality operator.
    /**
     * Compares the single polynomial coefficients of
     * two audi::gdual objects and returns true if equal.
     *
     * /note The truncation order of \p d1 and \p d2 may be different
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return The result of the comparison
     */
    friend bool operator==(const gdual &d1, const gdual &d2)
    {
        return (d1.m_p == d2.m_p);
    }

    /// Overloaded Non equality operator.
    /**
     * Compares the single polynomial coefficients of
     * two audi::gdual objects and returns false if equal.
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return The result of the comparison
     */
    friend bool operator!=(const gdual &d1, const gdual &d2)
    {
        return !(d1 == d2);
    }

    /// Overloaded less than operator
    /**
     * Compares the constant coefficients of two gduals
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return The result of the comparison
     */
    friend bool operator<(const gdual &d1, const gdual &d2)
    {
        return d1.constant_cf() < d2.constant_cf();
    }

    /// Overloaded more than operator
    /**
     * Compares the constant coefficients of two gduals
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return The result of the comparison
     */
    friend bool operator>(const gdual &d1, const gdual &d2)
    {
        return d1.constant_cf() > d2.constant_cf();
    }
    //@}

    /** @name Algebraic operators
     *
     */
    //@{

    /// Add and assignment operator
    template <typename T>
    auto operator+=(const T &d1) -> decltype(*this = *this + d1)
    {
        return *this = *this + d1;
    }

    /// Subtract and assignment operator
    template <typename T>
    auto operator-=(const T &d1) -> decltype(*this = *this - d1)
    {
        return *this = *this - d1;
    }

    /// Multiply and assignment operator
    template <typename T>
    auto operator*=(const T &d1) -> decltype(*this = *this * d1)
    {
        return *this = *this * d1;
    }

    /// Divide and assignment operator
    template <typename T>
    auto operator/=(const T &d1) -> decltype(*this = *this / d1)
    {
        return *this = *this / d1;
    }

    /// Negate operator
    gdual operator-() const
    {
        return gdual(-m_p, m_order);
    }

    /// Identity operator
    gdual operator+() const
    {
        return *this;
    }

    /// Overloaded addition operator
    /**
     * Implements the sum operation between truncated Taylor polynomials.
     * \note In order for this overload to be active (SFINAE rules), at least one
     * of the arguments must be an audi::gdual, while the second argument
     * may only be a double or int.
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return the sum between d1 and d2
     */
    template <typename T, typename U>
    friend gdual_if_enabled<T, U> operator+(const T &d1, const U &d2)
    {
        return add(d1, d2);
    }

    /// Overloaded difference operator
    /**
     * Implements the difference operation between truncated Taylor polynomials.
     * \note In order for this overload to be active (SFINAE rules), at least one
     * of the arguments must be an audi::gdual, while the second argument
     * may only be a double or int.
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return the difference between d1 and d2
     */
    template <typename T, typename U>
    friend gdual_if_enabled<T, U> operator-(const T &d1, const U &d2)
    {
        return sub(d1, d2);
    }

    /// Overloaded multiplication operator
    /**
     * Implements the multiplication operation between truncated Taylor polynomials.
     * \note In order for this overload to be active (SFINAE rules), at least one
     * of the arguments must be an audi::gdual, while the second argument
     * may only be a double or int.
     *
     * \note The truncated polynomial multiplication operator is at the very heart of AuDi
     * and its details / performances are those of the obake multiplication
     * algorithm which is, essentially, used.
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return the (truncated) multiplication between d1 and d2
     */
    template <typename T, typename U>
    friend gdual_if_enabled<T, U> operator*(const T &d1, const U &d2)
    {
        return mul(d1, d2);
    }

    /// Overloaded division operator
    /**
     * Implements the division operation between truncated Taylor polynomials.
     * Essentially, defined (in case of two audi::gdual) by a multiplication and
     * the reciprocal rule:
     *
     * \f[
     * T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
     * \f]
     *
     * where \f$T_f = f_0 + \hat f\f$.
     *
     * \note In order for this overload to be active (SFINAE rules), at least one
     * of the arguments must be an audi::gdual, while the second argument
     * may only be a double or int.
     *
     * @param d1 first audi::gdual argument
     * @param d2 second audi::gdual argument
     *
     * @return an audi:gdual containing the (truncated) multiplication between d1 and d2
     */
    template <typename T, typename U>
    friend gdual_if_enabled<T, U> operator/(const T &d1, const U &d2)
    {
        return div(d1, d2);
    }
    //@}

    /** @name Low-level interface
     *
     */
    //@{
    /// Get a mutable reference to the container of terms.
    /**
     * @return a reference to the internal container of terms.
     */
    auto _container() -> decltype(m_p._get_s_table())
    {
        return m_p._get_s_table();
    }

    /// Get a const reference to the container of terms.
    /**
     * @return a const reference to the internal container of terms.
     */
    auto _container() const -> decltype(m_p._get_s_table())
    {
        return m_p._get_s_table();
    }

    /// Get a const reference to the polynomial.
    /**
     * @return a const reference to the internal polynomial.
     */
    const p_type &_poly() const
    {
        return m_p;
    }

    /// Get a mutable reference to the polynomial.
    /**
     * @return a reference to the internal polynomial.
     */
    p_type &_poly()
    {
        return m_p;
    }
    //@}
private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_p;
        ar &m_order;
    }
};

} // end of namespace audi

#endif
