#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

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
#include <stdexcept>
#include <string>
#include <type_traits> // For std::enable_if, std::is_same, etc.
#include <utility>
#include <vector>
#include <cassert>

/// Root namespace for AuDi symbols
namespace audi
{

/// Generalized dual number class.
/**
 * This class represents a generalized dual number, or more formally, an element
 * of the truncated polynomial algebra \f$\mathcal P_{n,m}\f$. 
 *
 * The basic operations defined in the algebra \f$\mathcal P_{n,m}\f$ are
 * implemented as operators overloads, thus the new audi::gdual type can be used 
 * in substitution to the simple double type to also compute derivatives.
 *
 * The order of truncation \f$m\f$ is determined upon construction and cannot be later
 * modified. The number of variables \f$n\f$ will instead by determined dynamically
 * when operations are performed on the audi::gdual type.
 *
 * The actual truncated polynomial is contained in audi::gdual as a data member 
 * of type piranha::polynomial<double,piranha::monomial<char> >, 
 * allowing to support a high number of monomials.
 * 
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class gdual
{
        using p_type = piranha::polynomial<double,piranha::monomial<char>>;

        // We enable the overloads of the +,-,*,/ operators only in the following cases:
        // - at least one operand is a dual,
        // - the other operand, if not dual, must be double, int or unsigned int
        template <typename T, typename U>
        using gdual_if_enabled = typename std::enable_if<
        (std::is_same<T,gdual>::value && std::is_same<U,gdual>::value) ||
        (std::is_same<T,gdual>::value && std::is_same<U,double>::value) ||
        (std::is_same<T,gdual>::value && std::is_same<U,int>::value) ||
        (std::is_same<T,gdual>::value && std::is_same<U,unsigned int>::value) ||
        (std::is_same<U,gdual>::value && std::is_same<T,double>::value) ||
        (std::is_same<U,gdual>::value && std::is_same<T,int>::value) ||
        (std::is_same<U,gdual>::value && std::is_same<T,unsigned int>::value),
        gdual>::type;

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
        explicit gdual(p_type &&p, unsigned int order):m_p(std::move(p)),m_order(order) {}

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

        // Basic overloads for the multiplication
        // This is the low-level multiplication between two duals when they have
        // identical symbol set.
        static gdual mul_impl(const p_type &p1, const p_type &p2, unsigned int order)
        {
            using degree_type = decltype(p1.degree());
            piranha::series_multiplier<p_type> m(p1,p2);
            return gdual(m.truncated_multiplication(piranha::safe_cast<degree_type>(order)),order);
        }

        // Dual * dual.
        static gdual mul(const gdual &d1, const gdual &d2)
        {
            const unsigned int order = std::max(d1.get_order(), d2.get_order());
            const auto &ss1 = d1.m_p.get_symbol_set(), &ss2 = d2.m_p.get_symbol_set();
            if (ss1 == ss2) {
                return mul_impl(d1.m_p,d2.m_p,order);
            }
            // If the symbol sets are not the same, we need to merge them and make
            // copies of the original operands as needed.
            auto merge = ss1.merge(ss2);
            const bool need_copy_1 = (merge != ss1), need_copy_2 = (merge != ss2);
            if (need_copy_1) {
                p_type copy_1(d1.m_p.extend_symbol_set(merge));
                if (need_copy_2) {
                    p_type copy_2(d2.m_p.extend_symbol_set(merge));
                    return mul_impl(copy_1,copy_2,order);
                }
                return mul_impl(copy_1,d2.m_p,order);
            } else {
                p_type copy_2(d2.m_p.extend_symbol_set(merge));
                return mul_impl(d1.m_p, copy_2,order);
            }
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
            gdual retval(1.);
            double fatt = -1.;
            auto p0 = d2.constant_cf();
            auto phat = (d2 - p0);
            phat = phat / p0;
            gdual tmp(phat);

            retval = retval - phat;
            for (auto i = 2u; i <= d2.m_order; ++i) {
                fatt *= -1.;
                phat*=tmp;
                retval =  retval + fatt * phat;
            }

            return (d1 * retval) / p0;
        }

        template <typename T>
        static gdual div(const T &d1, const gdual &d2)
        {
            gdual retval(1.);
            double fatt = -1.;
            auto p0 = d2.constant_cf();
            if (p0 == 0) {
                throw std::domain_error("gdual: divide by zero");
            }
            auto phat = (d2 - p0);
            phat = phat / p0;
            gdual tmp(phat);

            retval = retval - phat;
            for (auto i = 2u; i <= d2.m_order; ++i) {
                fatt *= -1.;
                phat*=tmp;
                retval =  retval + fatt * phat;
            }
            return retval / (d1 * p0);
        }

        template <typename T>
        static gdual div(const gdual &d1, const T &d2)
        {
            return d1 * (1. / d2);
        }

        p_type	m_p;
        unsigned int m_order;

    public:
        /// Defaulted copy constructor.
        gdual(const gdual &) = default;
        /// Defaulted move constructor.
        gdual(gdual &&) = default;
        /// Default constuctor
        explicit gdual() : m_order(0u) {}
        /// Destructor (contains a sanity check)
        ~gdual() {assert(m_p.degree() <= (int)m_order);}

        /// Constructor from symbol and truncation order
        /**
         *
         * Will construct a generalized dual number made of a constant and a single term with unitary coefficient and exponent,
         * representing the expansion around zero of the symbolic variable \p symbol. The truncation order
         * is also set to \p order. 
         *
         * @note If the \p order is requested to be zero, this will instead construct a constant, while
         * keeping in the symbol set the requested symbol name. If, later on,
         * any derivative will be requested with respect to that symbol, it will be zero.
         * 
         * The type of \p symbol must be a string type (either C or C++) and its variation will be indicated prepending the letter "d"
         * so that "x" -> "dx". 
         * 
         * @param[in] symbol symbolic name
         * @param[in] order truncation order
         * 
         * @throws std::invalid_argument:
         * - if \p order is not in [0, std::numeric_limits<int>::max() - 10u]
         * - if \p symbol already starts with the letter "d" (this avoids to create confusing variation symbols of the form "ddname")
         */
        explicit gdual(const std::string &symbol, unsigned int order):m_p(),m_order(order)
        {
            check_var_name(symbol);
            if (order == 0) {
                extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
            } else {
                m_p = p_type(std::string("d") + symbol);
            }
        }

        /// Constructor from value and truncation order
        /**
         *
         * Will construct a generalized dual number representing a constant number
         * 
         * @param[in] value value of the constant
         * @param[in] order truncation order of the underlying algebra
         * 
         * @throws std::invalid_argument:
         * - if \p order is not in [0, std::numeric_limits<int>::max() - 10u]
         */
        explicit gdual(double value, unsigned int order):m_p(value),m_order(order)
        {
            check_order();
        }

        /// Constructor from value 
        /**
         *
         * Will construct a generalized dual number of order 0 representing
         * a constant number
         * 
         * @param[in] value value of the constant
         * 
         */
        explicit gdual(double value):m_p(value), m_order(0u) {}

        /// Constructor from value, symbol and truncation order
        /**
         *
         * Will construct a generalized dual number representing the expansion around \p value 
         * of the symbolic variable \p symbol. The truncation order is also set to \p order. 
         *
         * @note If the \p order is requested to be zero, this will instead construct a constant, while
         * keeping in the symbol set the requested symbol name. If, later on,
         * any derivative will be requested with respect to that symbol, it will be zero.
         * 
         * The type of \p symbol must be a string type (either C or C++) and its variation will be indicated prepending the letter "d"
         * so that "x" -> "dx". 
         * 
         * @param[in] value value of the variable
         * @param[in] symbol symbolic name
         * @param[in] order truncation order
         * 
         * @throws std::invalid_argument:
         * - if \p order is not in [0, std::numeric_limits<int>::max() - 10u]
         * - if \p symbol already starts with the letter "d" (this avoids to create confusing variation symbols of the form "ddname")
         */
        explicit gdual(double value, const std::string &symbol, unsigned int order):m_p(),m_order(order)
        {
            check_var_name(symbol);
            if (order == 0) {
                extend_symbol_set(std::vector<std::string>{std::string("d") + symbol});
            } else {
                m_p = p_type(std::string("d") + symbol);
            }
            m_p+=value;
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
         * Performs the integration of the gdual with respect to the symbol.
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
        gdual integrate(const std::string& var_name)
        {
            check_var_name(var_name);
            auto new_p = m_p.integrate("d" + var_name);
            new_p = new_p.truncate_degree(m_order);
            return gdual(std::move(new_p), m_order);
        }

        /// Partial derivative
        /**
         * Performs the partial derivative of the gdual with respect to the symbol
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
        gdual partial(const std::string& var_name)
        {
            check_var_name(var_name);
            auto new_p = m_p.partial("d" + var_name);
            return gdual(std::move(new_p), m_order);
        }

        /// Substitute symbol with value
        /**
         * Substitute the symbol \p sym with the value \p val
         *
         * @throws unspecified any exception thrown by:
         * - piranha::series::subs,
         */
        gdual subs( const std::string sym, double val)
        {
            auto new_p = m_p.subs(sym, val);
            return gdual(std::move(new_p), m_order);
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
         * Returns the current degree of the polynomial represented as an audi::gdual.
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
                cumfact*=boost::math::factorial<double>(*i);
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

        /// Finds the constant coefficient
        /**
         * Returns the coefficient of the of the constant part of the polynomial
         * so that if \f$T_{f} = f_0 + \hat f\f$,m \f$f_0\f$ is returned
         * 
         * \note This method is identical to the other overload with the same name, and it is provided for convenience.
         * @return the coefficient
         */
        double constant_cf() const
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
         * @param[in] d audi::gdual argument.
         * 
         * @return reference to \p os.
         * 
        */
        friend std::ostream& operator<<(std::ostream& os, const gdual& d)
        {
            os << d.m_p;
            return os;
        }

        /// Overloaded Equality operator.
        /**
         * Compares the single polynomial coefficients of
         * two audi::gdual objects and returns true if equal. 
         *
         * /note The truncation order of \p d1 and \p d2 may be different
         *
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         * 
         * @return The result of the cmparison
        */
        friend bool operator==(const gdual &d1, const gdual &d2)
        {
            return (d1.m_p == d2.m_p);
        }

        /// Overloaded Non equality operator.
        /**
         * Compares the truncation order and the single polynomial coefficients of
         * two audi::gdual objects and returns false if equal.
         *
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         * 
         * @return The result of the cmparison
        */
        friend bool operator!=(const gdual &d1, const gdual &d2)
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
         * Implements the sum operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual, while the second argument
         * may only be a double or int. 
         *
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         *
         * @return the sum between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator+(const T &d1, const U &d2)
        {
            return add(d1,d2);
        }

        /// Overloaded difference operator
        /**
         * Implements the difference operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual, while the second argument
         * may only be a double or int. 
         *
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         *
         * @return the difference between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator-(const T &d1, const U &d2)
        {
            return sub(d1,d2);
        }

        /// Overloaded multiplication operator
        /**
         * Implements the multiplication operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials.
         * \note In order for this overload to be active (SFINAE rules), at least one
         * of the arguments must be an audi::gdual, while the second argument
         * may only be a double or int. 
         *
         * \note The truncated polynomial multiplication operator is at the very heart of AuDi
         * and its details / performances are those of the piranha multiplication
         * algorithm which is, essentially, used. 
         *
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         *
         * @return the (truncated) multiplication between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator*(const T &d1, const U &d2)
        {
            return mul(d1,d2);
        }

        /// Overloaded division operator
        /**
         * Implements the division operation in the algebra \f$\mathcal P_{n,m}\f$
         * of truncated polynomials. Essentially defined (in case of two audi:gdual) by a multiplication and
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
         * @param[in] d1 first audi::gdual argument
         * @param[in] d2 second audi::gdual argument
         *
         * @return an audi:gdual containing the (truncated) multiplication between d1 and d2
        */
        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator/(const T &d1, const U &d2)
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
