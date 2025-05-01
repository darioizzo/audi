#ifndef AUDI_TM_HPP
#define AUDI_TM_HPP

#include <boost/numeric/interval.hpp>
#include <cassert>
#include <iostream>
#include <ranges>
#include <string>
#include <unordered_map>
#include <unordered_set>

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

#include <audi/audi.hpp>
#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>            //for audi::abs
#include <audi/detail/overloads_vectorized.hpp> //for audi::abs
#include <audi/gdual.hpp>

/// Root namespace for AuDi symbols
namespace audi
{

using int_d = boost::numeric::interval<double>;

/// Taylor model class.
/**
 *
 */
class taylor_model
{

public:
    using var_map_d = std::unordered_map<std::string, double>;
    using var_map_i = std::unordered_map<std::string, int_d>;

private:
    // The Taylor polynomial
    audi::gdual<double> m_tpol;
    // The order of the Taylor polynomial
    uint m_order;
    // The dimension
    uint m_ndim;
    // The remainder bound
    int_d m_rem;
    // The expansion point(s)
    var_map_d m_exp;
    // The domain
    var_map_i m_domain;

    void check_input_validity()
    {
        // For domain
        auto domain_keys = m_domain | std::views::keys;
        std::unordered_set<std::string> domain_set(domain_keys.begin(), domain_keys.end());

        // For exp_point
        auto exp_point_keys = m_exp | std::views::keys;
        std::unordered_set<std::string> exp_point_set(exp_point_keys.begin(), exp_point_keys.end());

        // For symbol set
        std::unordered_set<std::string> symbol_set(m_tpol.get_symbol_set().begin(), m_tpol.get_symbol_set().end());

        if (domain_set != exp_point_set || exp_point_set != symbol_set) {
            throw std::invalid_argument("Maps have different number of items.");
        }
    }

public:
    /// Defaulted copy constructor.
    taylor_model(const taylor_model &) = default;
    /// Defaulted move constructor.
    taylor_model(taylor_model &&) = default;
    /// Default constuctor
    taylor_model() : m_tpol(0.), m_order(0u), m_rem(0.) {}
    ~taylor_model() {}

    explicit taylor_model(const audi::gdual<double> &tpol, const int_d &rem_bound, const var_map_d &exp_point,
                          const var_map_i &domain)
    {

        m_tpol = tpol;
        m_order = tpol.get_order();
        m_ndim = static_cast<uint>(tpol.get_symbol_set_size());
        m_rem = rem_bound;
        m_exp = exp_point;
        m_domain = domain;

        check_input_validity();
    }

    ///////////////
    /// Getters ///
    ///////////////

    const audi::gdual<double>& get_tpol() const
    {
        return m_tpol;
    }

    const int_d& get_rem() const
    {
        return m_rem;
    }

    const var_map_d& get_exp() const
    {
        return m_exp;
    }

    const var_map_i& get_dom() const
    {
        return m_domain;
    }

    const uint& get_order() const
    {
        return m_order;
    }

    const uint& get_ndim() const
    {
        return m_ndim;
    }

    ///////////////
    /// Setters ///
    ///////////////

    void set_tpol(const audi::gdual<double> &tpol)
    {
        m_tpol = tpol;
        m_order = tpol.get_order();
        m_ndim = static_cast<uint>(tpol.get_symbol_set_size());
        check_input_validity();
    }

    void set_rem(const int_d &rem_bound)
    {
        m_rem = rem_bound;
    }

    void set_exp(const var_map_d &exp_point)
    {
        m_exp = exp_point;
        check_input_validity();
    }

    void set_domain(const var_map_i &domain)
    {
        m_domain = domain;
        check_input_validity();
    }

private:
    template <typename K, typename V>
    static std::unordered_map<K, V> get_common_map(const std::unordered_map<K, V> &a, const std::unordered_map<K, V> &b)
    {
        std::unordered_set<K> all_keys;

        // Insert keys from both maps using views::keys
        all_keys.insert((a | std::views::keys).begin(), (a | std::views::keys).end());
        all_keys.insert((b | std::views::keys).begin(), (b | std::views::keys).end());

        std::unordered_map<K, V> result;

        for (const K &key : all_keys) {
            auto it_a = a.find(key);
            auto it_b = b.find(key);

            if (it_a != a.end() && it_b != b.end()) {
                if (it_a->second != it_b->second) {
                    throw std::runtime_error("Conflicting values for key: " + key);
                }
                result[key] = it_a->second;
            } else if (it_a != a.end()) {
                result[key] = it_a->second;
            } else {
                result[key] = it_b->second;
            }
        }

        return result;
    }

    //////////////////////////////////////////////
    /// Static member functions for arithmetic ///
    //////////////////////////////////////////////

    template <typename T>
    static taylor_model add(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 + d2.m_tpol, d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T>
    static taylor_model add(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol + d2, d1.m_rem, d1.m_exp, d1.m_domain);
    }

    template <typename T>
    static taylor_model add(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d new_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i new_domain = get_common_map(d1.m_domain, d2.m_domain);
        return taylor_model(d1.m_tpol + d2.m_tpol, d1.m_rem + d2.m_rem, new_exp, new_domain);
    }

public:
    ////////////////////////////////////
    /// Overloaded addition operator ///
    ////////////////////////////////////

    template <typename T, typename U>
    // TODO: This needs to become taylor_model_if_enabled to specify  exactly what types the + operator can accept
    friend taylor_model operator+(const T &d1, const U &d2)
    {
        return add(d1, d2);
    }

    friend std::ostream &operator<<(std::ostream &os, const taylor_model &tm)
    {
        os << "Taylor Polynomial:\n" << tm.get_tpol(); // << "\nRemainder Bound: " << tm.get_rem().lower() << tm.get_rem().upper();
           // << "\nExpansion Point: " << tm.get_exp() << "\nDomain: " << tm.get_dom();
        return os;
    }


}; // end of taylor_model class


} // end of namespace audi

#endif
