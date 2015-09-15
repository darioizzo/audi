#ifndef AUDI_GDUAL_HPP
#define AUDI_GDUAL_HPP

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <piranha/mp_integer.hpp>
#include <piranha/polynomial.hpp>
#include <piranha/series_multiplier.hpp>
#include <stdexcept>
#include <string>
#include <type_traits> // For std::enable_if, std::is_same, etc.
#include <utility>
#include <vector>

/// Root namespace for AuDi symbols
namespace audi
{

namespace detail
{

// Main definition of the Piranha polynomial type.
using p_type = piranha::polynomial<double,piranha::monomial<char> >;

// Multiplier functor for the the gdual class, based on
// Piranha's polynomial multiplier.
class gdual_multiplier: piranha::series_multiplier<p_type>
{
        using base = piranha::series_multiplier<p_type>;
    public:
        explicit gdual_multiplier(const p_type &p1, const p_type &p2, int max_degree):base(p1,p2),m_max_degree(max_degree) {}
        p_type operator()() const
        {
            using term_type = p_type::term_type;
            using degree_type = piranha::integer;
            using size_type = base::size_type;
            // First let's order the terms in the second series according to the degree.
            std::stable_sort(this->m_v2.begin(),this->m_v2.end(),[this](term_type const *p1, term_type const *p2) {
                return p1->m_key.degree(this->m_ss) < p2->m_key.degree(this->m_ss);
            });
            // Next we create two vectors with the degrees of the terms in the two series. In the second series,
            // we negate and add the max degree in order to avoid adding in the skipping functor.
            std::vector<degree_type> v_d1, v_d2;
            std::transform(this->m_v1.begin(),this->m_v1.end(),std::back_inserter(v_d1),[this](term_type const *p) {
                return p->m_key.degree(this->m_ss);
            });
            std::transform(this->m_v2.begin(),this->m_v2.end(),std::back_inserter(v_d2),[this](term_type const *p) {
                return this->m_max_degree - p->m_key.degree(this->m_ss);
            });
            // The skipping functor.
            auto sf = [&v_d1,&v_d2](const size_type &i, const size_type &j) -> bool {
                using d_size_type = decltype(v_d1.size());
                return v_d1[static_cast<d_size_type>(i)] > v_d2[static_cast<d_size_type>(j)];
            };
            // The filter functor: will return 1 if the degree of the term resulting from the multiplication of i and j
            // is greater than the max degree, zero otherwise.
            auto ff = [&sf](const size_type &i, const size_type &j) {
                return static_cast<unsigned>(sf(i,j));
            };
            return this->plain_multiplication(sf,ff);
        }
    private:
        const int m_max_degree;
};

} // end of namespace detail


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
        // Lift the polynomial type from the detail namespace.
        using p_type = detail::p_type;

        // We enable the overloads of the +,-,*,/ operators only in the following cases:
        // - at least one operand is a dual,
        // - the other operand, if not dual, must be double or int.
        template <typename T, typename U>
        using gdual_if_enabled = typename std::enable_if<
        (std::is_same<T,gdual>::value && std::is_same<U,gdual>::value) ||
        (std::is_same<T,gdual>::value && std::is_same<U,double>::value) ||
        (std::is_same<T,gdual>::value && std::is_same<U,int>::value) ||
        (std::is_same<U,gdual>::value && std::is_same<T,double>::value) ||
        (std::is_same<U,gdual>::value && std::is_same<T,int>::value),
        gdual>::type;

        void check_order() const
        {
            if (m_order < 1) {
                throw std::invalid_argument("polynomial truncation order must be >= 1");
            }
            if (m_order == std::numeric_limits<int>::max()) {
                throw std::invalid_argument("polynomial truncation order is too large");
            }
        }

        void check_var_name(const std::string &name) const
        {
            if (name.at(0) == 'd') {
                throw std::invalid_argument("variable names cannot start with the letter d");
            }
        }


    private:
        // A private constructor to move-initialise a gdual from a polynomial. Used
        // in the implementation of the operators.
        explicit gdual(p_type &&p, int order):m_p(std::move(p)),m_order(order) {}

        // Basic overloads for the addition
        static gdual add(const gdual &d1, const gdual &d2)
        {
            if (d1.get_order() != d2.get_order()) {
                throw std::invalid_argument("different truncation limit");
            }
            return gdual(d1.m_p + d2.m_p,d1.get_order());
        }

        template <typename T>
        static gdual add(const T &d1, const gdual &d2)
        {
            return gdual(d1 + d2.m_p,d2.get_order());
        }

        template <typename T>
        static gdual add(const gdual &d1, const T &d2)
        {
            return add(d2, d1);
        }

        // Basic overloads for the subtraction
        static gdual sub(const gdual &d1, const gdual &d2)
        {
            if (d1.get_order() != d2.get_order()) {
                throw std::invalid_argument("different truncation limit");
            }
            return gdual(d1.m_p - d2.m_p,d1.get_order());
        }

        template <typename T>
        static gdual sub(const T &d1, const gdual &d2)
        {
            return gdual(d1 - d2.m_p,d2.get_order());
        }

        template <typename T>
        static gdual sub(const gdual &d1, const T &d2)
        {
            return gdual(d1.m_p - d2,d1.get_order());
        }

        // Basic overloads for the multiplication
        // This is the low-level multiplication between two duals when they have
        // identical symbol set.
        static gdual mul_impl(const p_type &p1, const p_type &p2, int order)
        {
            detail::gdual_multiplier gdm(p1,p2,order);
            return gdual(gdm(),order);
        }

        // Dual * dual.
        static gdual mul(const gdual &d1, const gdual &d2)
        {
            const int order = d1.get_order();
            if (order != d2.get_order()) {
                throw std::invalid_argument("different truncation limit");
            }
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
            if (d1.get_order() != d2.get_order()) {
                throw std::invalid_argument("different truncation limit");
            }

            gdual retval(1, d2.get_order());
            double fatt = -1;
            auto p0 = d2.constant_cf();
            if (p0 == 0) {
                throw std::domain_error("divide by zero");
            }
            auto phat = (d2 - p0);
            phat = phat/p0;
            gdual tmp(phat);

            retval = retval - phat;
            for (auto i = 2; i <= d2.m_order; ++i) {
                fatt *= -1;
                phat*=tmp;
                retval =  retval + fatt * phat;
            }

            return (d1 * retval) / p0;
        }

        template <typename T>
        static gdual div(const T &d1, const gdual &d2)
        {
            gdual retval(1, d2.get_order());
            double fatt = -1;
            auto p0 = d2.constant_cf();
            if (p0 == 0) {
                throw std::domain_error("divide by zero");
            }
            auto phat = (d2 - p0);
            phat = phat/p0;
            gdual tmp(phat);

            retval = retval - phat;
            for (auto i = 2; i <= d2.m_order; ++i) {
                fatt *= -1;
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
        int	m_order;

    public:
        /// Defaulted copy constructor.
        gdual(const gdual &) = default;
        /// Defaulted move constructor.
        gdual(gdual &&) = default;

        /// Constructor from symbol name and truncation order
        /**
         *
         * Will construct a generalized dual number made of a single term with unitary coefficient and exponent,
         * representing the expansion around zero of the symbolic variable \p name. The truncation order
         * is also set to \p order. 
         * The type of \p name must be a string type (either C or C++) and its variation will be indicated prepending the letter "d"
         * so that "x" -> "dx". 
         * 
         * @param[in] name name of the symbolic variable
         * 
         * @throws std::invalid_argument:
         * - if \p order is not in [1, std::numeric_limits<int>::max()]
         * - if \p name already starts with the letter "d" (this avoids to create confusing variation symbols of the form "ddname")
         */
        explicit gdual(const std::string &name, int order):m_p(std::string("d") + name),m_order(order)
        {
            check_order();
            check_var_name(name);
        }

        explicit gdual(double x, int order):m_p(x),m_order(order)
        {
            check_order();
        }

        explicit gdual(double value, const std::string &name, int order):m_p(std::string("d") + name),m_order(order)
        {
            check_order();
            check_var_name(name);
            m_p+=value;
        }

        // Defaulted assignment operators.
        gdual &operator=(const gdual &) = default;
        gdual &operator=(gdual &&) = default;

        std::vector<std::string> get_symbols() const
        {
            std::vector<std::string> retval;
            for (const auto &symbol : m_p.get_symbol_set()) {
                retval.push_back(symbol.get_name());
            }
            return retval;
        }

        auto _container() -> decltype(m_p._container())
        {
            return m_p._container();
        }

        auto degree() const -> decltype(m_p.degree())
        {
            return m_p.degree();
        }

        auto get_n_variables() const -> decltype(m_p.get_symbol_set().size())
        {
            return m_p.get_symbol_set().size();
        }

        int get_order() const
        {
            return m_order;
        }

        template <typename T>
        auto find_cf(const T &c) const -> decltype(m_p.find_cf(c))
        {
            return m_p.find_cf(c);
        }

        template <typename T>
        auto find_cf(std::initializer_list<T> l) const -> decltype(m_p.find_cf(l))
        {
            return m_p.find_cf(l);
        }

        double constant_cf() const
        {
            using v_size_type = std::vector<int>::size_type;
            return find_cf(std::vector<int>(boost::numeric_cast<v_size_type>(get_n_variables()),0));
        }

        friend std::ostream& operator<<(std::ostream& os, const gdual& d)
        {
            os << d.m_p;
            return os;
        }

        friend bool operator==(const gdual &d1, const gdual &d2)
        {
            return (d1.m_p == d2.m_p) && (d1.m_order == d2.m_order);
        }

        friend bool operator!=(const gdual &d1, const gdual &d2)
        {
            return !(d1 == d2);
        }

        template <typename T>
        auto operator+=(const T &d1) -> decltype(*this = *this + d1)
        {
            return *this = *this + d1;
        }

        template <typename T>
        auto operator-=(const T &d1) -> decltype(*this = *this - d1)
        {
            return *this = *this - d1;
        }

        template <typename T>
        auto operator*=(const T &d1) -> decltype(*this = *this * d1)
        {
            return *this = *this * d1;
        }

        template <typename T>
        auto operator/=(const T &d1) -> decltype(*this = *this / d1)
        {
            return *this = *this / d1;
        }

        gdual operator-() const
        {
            return gdual(-m_p,m_order);
        }

        gdual operator+() const
        {
            return *this;
        }

        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator+(const T &d1, const U &d2)
        {
            return add(d1,d2);
        }

        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator-(const T &d1, const U &d2)
        {
            return sub(d1,d2);
        }

        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator*(const T &d1, const U &d2)
        {
            return mul(d1,d2);
        }

        template <typename T, typename U>
        friend gdual_if_enabled<T, U> operator/(const T &d1, const U &d2)
        {
            return div(d1,d2);
        }

    
};



} // end of namespace audi 

#endif
