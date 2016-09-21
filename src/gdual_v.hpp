#ifndef AUDI_GDUAL_V_HPP
#define AUDI_GDUAL_V_HPP

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <piranha/monomial.hpp>
#include <piranha/math.hpp>
#include <piranha/polynomial.hpp>
#include <piranha/safe_cast.hpp>
#include <piranha/series_multiplier.hpp>
#include <piranha/symbol.hpp>
#include <stdexcept>
#include <string>
#include <type_traits> // For std::enable_if, std::is_same, etc.
#include <utility>
#include <vector>
#include <cassert>

#include "detail/coefficient_v.hpp"

/// Root namespace for AuDi symbols
namespace audi
{

/// Vectorized Generalized dual number class.
/**
 * This class represents a vectorized audi::gdual, or more formally, an element
 * of the truncated polynomial algebra \f$\mathcal P_{n,m}\f$ where the polynomial
 * coefficients are vectors and all operations operate element-wise
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class gdual_v
{
        using p_type = piranha::polynomial<detail::coefficient_v,piranha::monomial<char>>;

        // We enable the overloads of the +,-,*,/ operators only in the following cases:
        // - at least one operand is a gdual_v,
        // - the other operand, if not gdual_v, must be double, int or unsigned int or detail::coefficient_v
        template <typename T, typename U>
        using gdual_v_if_enabled = typename std::enable_if<
        (std::is_same<T,gdual_v>::value && std::is_same<U,gdual_v>::value) ||
        (std::is_same<T,gdual_v>::value && std::is_same<U,detail::coefficient_v>::value) ||
        (std::is_same<T,gdual_v>::value && std::is_same<U,double>::value) ||
        (std::is_same<T,gdual_v>::value && std::is_same<U,int>::value) ||
        (std::is_same<T,gdual_v>::value && std::is_same<U,unsigned int>::value) ||
        (std::is_same<U,gdual_v>::value && std::is_same<T,detail::coefficient_v>::value) ||
        (std::is_same<U,gdual_v>::value && std::is_same<T,double>::value) ||
        (std::is_same<U,gdual_v>::value && std::is_same<T,int>::value) ||
        (std::is_same<U,gdual_v>::value && std::is_same<T,unsigned int>::value),
        gdual_v>::type;

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

        // A private constructor to move-initialise a gdual_v from a polynomial. Used
        // in the implementation of the operators.
        explicit gdual_v(p_type &&p, unsigned int order):m_p(std::move(p)),m_order(order) {}
        // A private constructor used in the implementation of the operators (is it necessary?)
        explicit gdual_v(detail::coefficient_v value, unsigned int order):m_p(value),m_order(order) {}

        // Basic overloads for the addition
        static gdual_v add(const gdual_v &d1, const gdual_v &d2)
        {
            return gdual_v(d1.m_p + d2.m_p, std::max(d1.get_order(), d2.get_order()));
        }

        template <typename T>
        static gdual_v add(const T &d1, const gdual_v &d2)
        {
            return gdual_v(d1 + d2.m_p, d2.get_order());
        }

        template <typename T>
        static gdual_v add(const gdual_v &d1, const T &d2)
        {
            return gdual_v(d2 + d1.m_p, d1.get_order());
        }

        // Basic overloads for the subtraction
        static gdual_v sub(const gdual_v &d1, const gdual_v &d2)
        {
            return gdual_v(d1.m_p - d2.m_p, std::max(d1.get_order(), d2.get_order()));
        }

        template <typename T>
        static gdual_v sub(const T &d1, const gdual_v &d2)
        {
            return gdual_v(d1 - d2.m_p, d2.get_order());
        }

        template <typename T>
        static gdual_v sub(const gdual_v &d1, const T &d2)
        {
            return gdual_v(d1.m_p - d2, d1.get_order());
        }

        // Dual * dual.
        static gdual_v mul(const gdual_v &d1, const gdual_v &d2)
        {
            const unsigned int order = std::max(d1.get_order(), d2.get_order());
            return gdual_v(p_type::truncated_multiplication(d1.m_p,d2.m_p,order),order);
        }

        template <typename T>
        static gdual_v mul(const T &d1, const gdual_v &d2)
        {
            return gdual_v(d1, d2.get_order()) * d2;
        }

        template <typename T>
        static gdual_v mul(const gdual_v &d1, const T &d2)
        {
            return d1 * gdual_v(d2, d1.get_order());
        }

        // Basic overloads for the division
        static gdual_v div(const gdual_v &d1, const gdual_v &d2)
        {
            gdual_v retval({1.});
            double fatt = -1.;
            auto p0 = d2.constant_cf();
            auto phat = (d2 - p0);
            phat = phat / p0;
            gdual_v tmp(phat);

            retval = retval - phat;
            for (auto i = 2u; i <= d2.m_order; ++i) {
                fatt *= -1.;
                phat*=tmp;
                retval =  retval + fatt * phat;
            }

            return (d1 * retval) / p0;
        }

        template <typename T>
        static gdual_v div(const T &d1, const gdual_v &d2)
        {
            gdual_v retval({1.});
            double fatt = -1.;
            auto p0 = d2.constant_cf();
            if (p0 == 0) {
                throw std::domain_error("gdual_v: divide by zero");
            }
            auto phat = (d2 - p0);
            phat = phat / p0;
            gdual_v tmp(phat);

            retval = retval - phat;
            for (auto i = 2u; i <= d2.m_order; ++i) {
                fatt *= -1.;
                phat*=tmp;
                retval =  retval + fatt * phat;
            }
            return retval / (d1 * p0);
        }

        template <typename T>
        static gdual_v div(const gdual_v &d1, const T &d2)
        {
            return d1 * (1. / d2);
        }

        p_type	m_p;
        unsigned int m_order;

    public:
        /// Defaulted copy constructor.
        gdual_v(const gdual_v &) = default;
        /// Defaulted move constructor.
        gdual_v(gdual_v &&) = default;
        /// Default constuctor
        explicit gdual_v() : m_order(0u) {}
        /// Destructor (contains a sanity check)
        ~gdual_v()
         {
            assert(m_p.degree() >= 0);
            assert(static_cast<unsigned>(m_p.degree()) <= m_order);
         }

        /// Constructs a constant gdual_v
        /**
         * Will construct a generalized dual number of order 0 containing a constant coefficient
         *
         * @param[in] value initializer_list containing the various vaues of the vectorized coefficient
         */
        explicit gdual_v(std::initializer_list<double> value):m_p(value), m_order(0u) {}
        /// Constructs a constant gdual_v
        /**
         * Will construct a generalized dual number of order 0 containing a constant coefficient
         *
         * @param[in] value std::vector<double> containing the vectorized coefficient
         */
        explicit gdual_v(std::vector<double> value):m_p(value), m_order(0u) {}

        /// Constructs a gdual_v
        /**
         * Will construct a vectorized generalized dual number representing the expansion around \p value
         * of the symbolic variable \p symbol. The truncation order is also set to \p order.
         *
         * @note If the \p order is requested to be zero, this will instead construct a constant, while
         * keeping in the symbol set the requested symbol name. If, later on,
         * any derivative will be requested with respect to that symbol, it will be zero.
         *
         * The type of \p symbol must be a string type (either C or C++) and its variation will be indicated prepending the letter "d"
         * so that "x" -> "dx".
         *
         * @param[in] value std::vector<double> containing all the values of the variable
         * @param[in] symbol symbolic name
         * @param[in] order truncation order
         *
         * @throws std::invalid_argument:
         * - if \p order is not in [0, std::numeric_limits<int>::max() - 10u]
         * - if \p symbol already starts with the letter "d" (this avoids to create confusing variation symbols of the form "ddname")
         */
        explicit gdual_v(std::vector<double> value, const std::string &symbol, unsigned int order):m_p(),m_order(order)
        {
            check_order();
            check_var_name(symbol);
            if (order == 0) {
                extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
            } else {
                m_p = p_type(std::string("d") + symbol) * std::vector<double>(value.size(),1.);
            }
            m_p+=detail::coefficient_v(value);
        }
        /// Constructs a gdual_v
        /**
         * Will construct a vectorized generalized dual number representing the expansion around \p value
         * of the symbolic variable \p symbol. The truncation order is also set to \p order.
         *
         * @note If the \p order is requested to be zero, this will instead construct a constant, while
         * keeping in the symbol set the requested symbol name. If, later on,
         * any derivative will be requested with respect to that symbol, it will be zero.
         *
         * The type of \p symbol must be a string type (either C or C++) and its variation will be indicated prepending the letter "d"
         * so that "x" -> "dx".
         *
         * @param[in] value std::initializer_list<double> containing all the values of the variable
         * @param[in] symbol symbolic name
         * @param[in] order truncation order
         *
         * @throws std::invalid_argument:
         * - if \p order is not in [0, std::numeric_limits<int>::max() - 10u]
         * - if \p symbol already starts with the letter "d" (this avoids to create confusing variation symbols of the form "ddname")
         */
        explicit gdual_v(std::initializer_list<double> value, const std::string &symbol, unsigned int order):m_p(),m_order(order)
        {
            check_order();
            check_var_name(symbol);
            if (order == 0) {
                extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
            } else {
                m_p = p_type(std::string("d") + symbol) * std::vector<double>(value.size(),1.);
            }
            m_p+=detail::coefficient_v(value);
        }

        /// Defaulted assignment operator
        gdual_v &operator=(const gdual_v &) = default;
        /// Defaulted assignment operator
        gdual_v &operator=(gdual_v &&) = default;

        /// Gets the symbol set size
        /**
         * Returns the size of the symbol set.
         *
         * @return the size of the symbol set.
         *
         */
        auto get_symbol_set_size() const -> decltype(m_p.get_symbol_set().size())
        {
            return m_p.get_symbol_set().size();
        }

        /// Returns the symbol set size
        /**
         * Returns the symbol set (not the differentials)
         *
         * @return the symbol set
         */
        std::vector<std::string> get_symbol_set() const
        {
            std::vector<std::string> retval;
            for (auto s : m_p.get_symbol_set()) {
                std::string tmp(s.get_name());
                retval.push_back(tmp.erase(0,1));
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
         * @param[in] sym_vars list of symbolic names. It must contain all symbolic names of
         * the current polynomial. It may contain more. All symbols must start with the letter "d".
         *
         * @throws std::invalid_argument:
         * - if any symbol in \p sym_vars does not start with the letter "d"
         *
         * @throws unspecified any exception thrown by:
         * - piranha::series::extend_symbol_set,
         */
        void extend_symbol_set(const std::vector<std::string>& sym_vars)
        {
            piranha::symbol_set ss;
            for (auto sym_var : sym_vars) {
                check_var_name_has_d(sym_var);
                ss.add(sym_var);
            }
            m_p = m_p.extend_symbol_set(ss);
        }

        /// Integration
        /**
         * Performs the integration of the gdual_v with respect to the symbol.
         *
         * \note If the \p var_name differential is not in the symbol set, then it is added.
         *
         * \note Information may be lost as the truncation order is preserved.
         *
         * @param[in] var_name Symbol name (cannot start with "d").
         *
         * @throws std::invalid_argument:
         * - if \p var_name starts with the letter "d" (this avoid creating confusing names for symbol's differentials)
         *
         * @throws unspecified any exception thrown by:
         * - piranha::series::truncate_degree,
         * - piranha::series::integrate
         */
        gdual_v integrate(const std::string& var_name)
        {
            check_var_name(var_name);
            auto new_p = m_p.integrate("d" + var_name);
            new_p = new_p.truncate_degree(static_cast<int>(m_order));
            return gdual_v(std::move(new_p), m_order);
        }

        /// Partial derivative
        /**
         * Performs the partial derivative of the gdual_v with respect to the symbol
         *
         * \note If the \p var_name differential is not in the symbol set, then it is added.
         *
         * @param[in] var_name Symbol name (cannot start with "d").
         *
         * @throws std::invalid_argument:
         * - if \p symbol starts with the letter "d" (this avoid creating confusing names for symbol's differentials)
         *
         * @throws unspecified any exception thrown by:
         * - piranha::series::partial,
         */
        gdual_v partial(const std::string& var_name)
        {
            check_var_name(var_name);
            auto new_p = m_p.partial("d" + var_name);
            return gdual_v(std::move(new_p), m_order);
        }

        /// Substitute symbol with value
        /**
         * Substitute the symbol \p sym with the value \p val
         *
         * @throws unspecified any exception thrown by:
         * - piranha::series::subs,
         */
        gdual_v subs( const std::string sym, detail::coefficient_v val)
        {
            auto new_p = m_p.subs(sym, val);
            return gdual_v(std::move(new_p), m_order);
        }

        /// Evaluates the Taylor polynomial
        /**
         * Evaluates the Taylor polynomial using the values in \p dict for all the
         * differentials (variables variations)
         *
         * @throws unspecified any exception thrown by:
         * - piranha::math::evaluate,
         */
        auto evaluate(const std::unordered_map< std::string, double> &dict) const -> decltype(piranha::math::evaluate(m_p, dict))
        {
            auto retval = piranha::math::evaluate(m_p, dict);
            return retval;
        }
        /// Current degree
        /**
         * Returns the current degree of the polynomial represented as an audi::gdual_v.
         * This may be different from the truncation order \f$m\f$ and,
         * in particular will be smaller or equal.
         *
         * @return the current degree of the polynomial
         */
        auto degree() const -> decltype(m_p.degree())
        {
            return m_p.degree();
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
        auto find_cf(const T &c) const -> decltype(m_p.find_cf(c))
        {
            if (std::accumulate(c.begin(),c.end(),0u) > m_order) {
                throw std::invalid_argument("requested coefficient is beyond the truncation order.");
            }
            return m_p.find_cf(c);
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
        auto find_cf(std::initializer_list<T> l) const -> decltype(m_p.find_cf(l))
        {
            if (std::accumulate(l.begin(),l.end(),0u) > m_order) {
                throw std::invalid_argument("requested coefficient is beyond the truncation order.");
            }
            return m_p.find_cf(l);
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
        auto get_derivative(const T &c) const -> decltype(m_p.find_cf(c))
        {
            double cumfact = 1;
            for (auto i = c.begin(); i < c.end(); ++i)
            {
                cumfact*=boost::math::factorial<double>(static_cast<unsigned int>(*i));
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
        auto get_derivative(std::initializer_list<T> l) const -> decltype(m_p.find_cf(l))
        {
            double cumfact = 1;
            for (auto i = l.begin(); i < l.end(); ++i)
            {
                cumfact*=boost::math::factorial<double>((unsigned int)(*i));
            }
            return this->find_cf(l) * cumfact;
        }

        /// Gets the derivative value
        /**
         * Returns the (mixed) derivative value of order specified
         * by the container \p c
         *
         * \note To get the following derivative: \f$ \frac{d^6}{dxdy^3dz^2}\f$
         * the input should be {{"x", 1u},{"y",3u},{"z",2u}}
         *
         * \note The current implementation call internally the other templated
         * implementations. WHen piranha will implement the sparse monomial
         * this will change and be more efficient.
         *
         * @return the value of the derivative
         *
         * @throws unspecified all exceptions thrown by the templated version call.
         * @throws std::invalid_argument: if one of the symbols is not found in the expression
         */
        auto get_derivative(const std::unordered_map<std::string, unsigned int> &dict) const -> decltype(get_derivative(std::vector<double>{}))
        {
            const auto &ss = m_p.get_symbol_set();
            std::vector<double> coeff(ss.size(), 0);
            for (const auto &entry : dict) {
                auto idx = ss.index_of(piranha::symbol{entry.first});
                if (idx == ss.size()) {
                    throw std::invalid_argument("Symbol not found in the symbol set, cannot return a derivative");
                }
                coeff[idx] = entry.second;
            }
            return get_derivative(coeff);
        }

        /// Finds the constant coefficient
        /**
         * Returns the coefficient of the of the constant part of the polynomial
         * so that if \f$T_{f} = f_0 + \hat f\f$,m \f$f_0\f$ is returned
         *
         * \note This method is identical to the other overload with the same name, and it is provided for convenience.
         * @return the coefficient
         */
        detail::coefficient_v constant_cf() const
        {
            using v_size_type = std::vector<int>::size_type;
            return find_cf(std::vector<int>(boost::numeric_cast<v_size_type>(get_symbol_set_size()),0));
        }


        /// Overloaded stream operator
        /**
         * Will direct to stream a human-readable representation of the generalized dual number.
         * It uses the piranha overload for the type piranha::series. Refer to that
         * documentation for further details
         *
         * \note The print order of the terms will be undefined.
         * \note At most piranha::settings::get_max_term_output() terms
         * are printed, and terms in excess are represented with ellipsis "..."
         * at the end of the output; if piranha::settings::get_max_term_output()
         * is zero, all the terms will be printed. piranha::settings::set_max_term_output()
         * is used to set this parameter.
         *
         * @param[in,out] os target stream.
         * @param[in] d audi::gdual_v argument.
         *
         * @return reference to \p os.
         *
        */
        friend std::ostream& operator<<(std::ostream& os, const gdual_v& d)
        {
            os << d.m_p;
            return os;
        }

        /// Overloaded Equality operator.
        /**
         * Compares the single polynomial coefficients of
         * two audi::gdual_v objects and returns true if equal.
         *
         * /note The truncation order of \p d1 and \p d2 may be different
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return The result of the cmparison
        */
        friend bool operator==(const gdual_v &d1, const gdual_v &d2)
        {
            return (d1.m_p == d2.m_p);
        }

        /// Overloaded Non equality operator.
        /**
         * Compares the truncation order and the single polynomial coefficients of
         * two audi::gdual_v objects and returns false if equal.
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return The result of the cmparison
        */
        friend bool operator!=(const gdual_v &d1, const gdual_v &d2)
        {
            return !(d1 == d2);
        }


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
        gdual_v operator-() const
        {
            return gdual_v(-m_p, m_order);
        }

        /// Identity operator
        gdual_v operator+() const
        {
            return *this;
        }

        /// Overloaded addition operator
        /**
         * Implements the sum operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual_v, while the second argument
         * may only be a double or int.
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return the sum between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_v_if_enabled<T, U> operator+(const T &d1, const U &d2)
        {
            return add(d1,d2);
        }

        /// Overloaded difference operator
        /**
         * Implements the difference operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual_v, while the second argument
         * may only be a double or int.
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return the difference between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_v_if_enabled<T, U> operator-(const T &d1, const U &d2)
        {
            return sub(d1,d2);
        }

        /// Overloaded multiplication operator
        /**
         * Implements the multiplication operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual_v, while the second argument
         * may only be a double or int.
         *
         * \note The truncated polynomial multiplication operator is at the very heart of AuDi
         * and its details / performances are those of the piranha multiplication
         * algorithm which is, essentially, used.
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return the (truncated) multiplication between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_v_if_enabled<T, U> operator*(const T &d1, const U &d2)
        {
            return mul(d1,d2);
        }

        /// Overloaded division operator
        /**
         * Implements the division operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials. Essentially defined (in case of two audi:gdual_v) by a multiplication and
         * the reciprocal rule:
         *
         * \f[
         * T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
         * \f]
         *
         * where \f$T_f = f_0 + \hat f\f$.
         *
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual_v, while the second argument
         * may only be a double or int.
         *
         * @param[in] d1 first audi::gdual_v argument
         * @param[in] d2 second audi::gdual_v argument
         *
         * @return an audi:gdual_v containing the (truncated) multiplication between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_v_if_enabled<T, U> operator/(const T &d1, const U &d2)
        {
            return div(d1,d2);
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
        auto _container() -> decltype(m_p._container())
        {
            return m_p._container();
        }

        /// Get a const reference to the container of terms.
        /**
         * @return a const reference to the internal container of terms.
         */
        auto _container() const -> decltype(m_p._container())
        {
            return m_p._container();
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

        /// Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int)
        {
            ar & m_p;
            ar & m_order;
        }
};



} // end of namespace audi

#endif
