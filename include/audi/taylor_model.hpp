#ifndef AUDI_TM_HPP
#define AUDI_TM_HPP

#include <boost/numeric/interval.hpp>
#include <cassert>
#include <ranges>
#include <unordered_map>
#include <unordered_set>

#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>            //for audi::abs
#include <audi/detail/overloads_vectorized.hpp> //for audi::abs
#include <audi/gdual.hpp>

#include <audi/taylor_model_bounding.hpp>

/// Root namespace for AuDi symbols
namespace audi
{

using int_d = boost::numeric::interval<double>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

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

/// Taylor model class.
/**
 *
 */
class taylor_model
{

public:
    void check_input_validity()
    {
        // For domain
        auto domain_keys = m_domain | std::views::keys;
        std::unordered_set<std::string> domain_set(domain_keys.begin(), domain_keys.end());

        // For exp_point
        auto exp_point_keys = m_exp | std::views::keys;
        std::unordered_set<std::string> exp_point_set(exp_point_keys.begin(), exp_point_keys.end());

        if (domain_set != exp_point_set) {
            throw std::invalid_argument("Maps have different number of items.");
        }
    }

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

public:
    /// Defaulted copy constructor & assignment
    taylor_model(const taylor_model &) = default;
    taylor_model &operator=(const taylor_model &) = default;

    /// Defaulted move constructor & assignment
    taylor_model(taylor_model &&) = default;
    taylor_model &operator=(taylor_model &&) = default;

    /// Default constructor
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

    taylor_model(double constant)
    {
        m_tpol = audi::gdual<double>(constant);
        m_order = m_tpol.get_order();
        m_ndim = static_cast<uint>(m_tpol.get_symbol_set_size());
        m_rem = int_d(0.0);
        m_exp = {};
        m_domain = {};
    }

    ///////////////
    /// Getters ///
    ///////////////

    const audi::gdual<double> &get_tpol() const
    {
        return m_tpol;
    }

    const int_d &get_rem() const
    {
        return m_rem;
    }

    const var_map_d &get_exp() const
    {
        return m_exp;
    }

    const var_map_i &get_dom() const
    {
        return m_domain;
    }

    const uint &get_order() const
    {
        return m_order;
    }

    const uint &get_ndim() const
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
                if (!interval_equal(it_a->second, it_b->second)) {
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

public:
    template <typename T>
    static audi::gdual<T> get_increased_order(const audi::gdual<T> &poly, uint new_order)
    {
        // Create a dummy polynomial with a higher order
        audi::gdual<T> temp(0.0, "temp", new_order);

        // Substitute poly into temp
        return temp.subs("dtemp", poly);
    }

    static taylor_model identity()
    {
        taylor_model tm;
        tm.m_tpol = audi::gdual<double>(1.0);
        tm.m_order = 0u;
        tm.m_ndim = 0u;
        tm.m_rem = 0.0;
        tm.m_exp = {};
        tm.m_domain = {};
        return tm;
    }

    template <typename T>
    static std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
    {
        std::vector<T> ret;
        for (const auto &v : orig)
            ret.insert(ret.end(), v.begin(), v.end());
        return ret;
    }

    int_d get_bounds() const
    {
        return get_bounds(m_tpol, m_exp, m_domain);
    }

private:
    template <typename T>
    static int_d get_bounds(const audi::gdual<T> &tpol, const var_map_d &exp_points, const var_map_i &domain)
    {
        // Build shifted domain relative to expansion points
        var_map_i domain_shifted;

        auto [coeffs, exps] = audi::get_poly(tpol);
        auto ndim = audi::get_ndim(coeffs, exps); // necessary because static function

        std::vector<T> flat;
        if (ndim > 1) {

            for (size_t i = 0; i < tpol.get_symbol_set_size(); ++i) {
                std::string var_name = tpol.get_symbol_set()[i];

                T a_shifted = domain.at(var_name).lower() - exp_points.at(var_name);
                T b_shifted = domain.at(var_name).upper() - exp_points.at(var_name);
                domain_shifted[var_name] = int_d(a_shifted, b_shifted);
            }

            // Compute Bernstein patch using generalbox (takes var_map_i directly)
            std::vector<std::vector<T>> bern_patch
                = audi::get_titi_bernstein_patch_ndim_generalbox(coeffs, exps, domain_shifted);

            flat = flatten(bern_patch); // flatten the nested vector
        } else if (ndim == 1) {
            flat = get_bernstein_coefficients(coeffs, exps, ndim);
        } else {
            throw std::runtime_error("The dimension cannot be negative.");
        }

        auto [min_val, max_val] = std::minmax_element(flat.begin(), flat.end());

        return int_d(*min_val, *max_val);
    }

    //////////////////////////////////////////////
    /// Static member functions for arithmetic ///
    //////////////////////////////////////////////

    // Addition
    static taylor_model add(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d new_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i new_domain = get_common_map(d1.m_domain, d2.m_domain);
        return taylor_model(d1.m_tpol + d2.m_tpol, d1.m_rem + d2.m_rem, new_exp, new_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 + d2.m_tpol, d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol + d2, d1.m_rem, d1.m_exp, d1.m_domain);
    }

    template <typename T>
    static taylor_model add(const boost::numeric::interval<T> &d1, const taylor_model &d2)
    {
        return taylor_model(d2.m_tpol, d1 + d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T>
    static taylor_model add(const taylor_model &d1, const boost::numeric::interval<T> &d2)
    {
        return taylor_model(d1.m_tpol, d1.m_rem + d2, d1.m_exp, d1.m_domain);
    }

    // Subtraction
    static taylor_model sub(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d new_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i new_domain = get_common_map(d1.m_domain, d2.m_domain);
        return taylor_model(d1.m_tpol - d2.m_tpol, d1.m_rem - d2.m_rem, new_exp, new_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 - d2.m_tpol, d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol - d2, d1.m_rem, d1.m_exp, d1.m_domain);
    }

    template <typename T>
    static taylor_model sub(const boost::numeric::interval<T> &d1, const taylor_model &d2)
    {
        return taylor_model(d2.m_tpol, d1 - d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T>
    static taylor_model sub(const taylor_model &d1, const boost::numeric::interval<T> &d2)
    {
        return taylor_model(d1.m_tpol, d1.m_rem - d2, d1.m_exp, d1.m_domain);
    }

    // // Multiplication

    // TODO: Should be templated as well I guess, but raises error with ambiguity
    static taylor_model multiply(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d new_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i new_domain = get_common_map(d1.m_domain, d2.m_domain);

        // Get untruncated polynomial product
        uint new_order = d1.m_tpol.get_order() + d2.m_tpol.get_order();

        audi::gdual<double> d1_ho = get_increased_order(d1.m_tpol, new_order);
        audi::gdual<double> d2_ho = get_increased_order(d2.m_tpol, new_order);

        audi::gdual<double> polynomial_product = d1_ho * d2_ho;

        // Get agreeable polynomial
        const uint agreeable_degree = std::max(d1.m_tpol.get_order(), d2.m_tpol.get_order());

        audi::gdual<double> agreeable_pol(0.0);
        for (uint i = 0; i <= agreeable_degree; ++i) {
            agreeable_pol += polynomial_product.extract_terms(i);
        }

        audi::gdual<double> pol_e = polynomial_product - agreeable_pol;

        int_d pol_e_bounds;
        if (pol_e.is_zero(1e-16)) {
            pol_e_bounds = boost::numeric::interval<double>(0.0, 0.0);
        } else {
            // extract symbols used in pol_e
            std::vector<std::string> symbs_e = pol_e.get_symbol_set();

            // extract relevant domain and expansion points
            var_map_i domain_e;
            var_map_d exp_point_e;
            for (const auto &sym : symbs_e) {
                domain_e[sym] = new_domain[sym];
                exp_point_e[sym] = new_exp[sym];
            }

            pol_e_bounds = get_bounds(pol_e, exp_point_e, domain_e);
        }

        int_d interval_of_product
            = pol_e_bounds + d1.get_bounds() * d2.m_rem + d1.m_rem * d2.get_bounds() + d1.m_rem * d2.m_rem;

        return taylor_model(agreeable_pol, interval_of_product, new_exp, new_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 * d2.m_tpol, int_d(d1) * d2.m_rem, d2.m_exp, d2.m_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol * d2, d1.m_rem * int_d(d2), d1.m_exp, d1.m_domain);
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const boost::numeric::interval<T> &d1, const taylor_model &d2)
    {
        audi::gdual<T> gd_one(1.0);
        const taylor_model tm_one(gd_one, d1, d2.m_exp, d2.m_domain);
        return tm_one * d2;
    }

    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const taylor_model &d1, const boost::numeric::interval<T> &d2)
    {
        audi::gdual<T> gd_one(1.0, "x", 0);
        const taylor_model tm_one(gd_one, d2, d1.m_exp, d1.m_domain);
        return d1 * tm_one;
    }

public:
    //////////////////////////////////////
    /// Overloaded arithmetic operator ///
    //////////////////////////////////////

    // TODO: These need to become taylor_model_if_enabled (analogous to gdual_if_enabled) to specify exactly what
    // types the operators can accept

    template <typename T, typename U>
    friend taylor_model operator+(const T &d1, const U &d2)
    {
        return add(d1, d2);
    }

    template <typename T, typename U>
    friend taylor_model operator-(const T &d1, const U &d2)
    {
        return sub(d1, d2);
    }

    template <typename T, typename U>
    friend taylor_model operator*(const T &d1, const U &d2)
    {
        return multiply(d1, d2);
    }

    /////////////////////////////////
    /// Overloaded other operator ///
    /////////////////////////////////

    friend std::ostream &operator<<(std::ostream &os, const taylor_model &tm)
    {
        os << "Taylor Polynomial:\n"
           << tm.get_tpol() << "\nRemainder Bound: [" << tm.get_rem().lower() << ", " << tm.get_rem().upper() << "]"
           << "\nExpansion Point: " << tm.get_exp() << "\nDomain: " << tm.get_dom();
        return os;
    }

    // Compare scalars
    template <typename T>
    static bool interval_equal(const T &a, const T &b)
    {
        if constexpr (std::is_floating_point_v<T>) {
            auto eps = std::numeric_limits<T>::epsilon() * 10;
            return std::fabs(a - b) < eps;
        } else {
            return a == b;
        }
    }

    // Compare intervals
    template <typename T>
    static bool interval_equal(const boost::numeric::interval<T> &a, const boost::numeric::interval<T> &b)
    {
        if constexpr (std::is_floating_point_v<T>) {
            auto eps = std::numeric_limits<T>::epsilon() * 10;
            return std::fabs(a.lower() - b.lower()) < eps && std::fabs(a.upper() - b.upper()) < eps;
        } else {
            return a.lower() == b.lower() && a.upper() == b.upper();
        }
    }

    // Generic map comparison
    template <typename T>
    static bool map_interval_equal(const std::unordered_map<std::string, T> &a,
                                   const std::unordered_map<std::string, T> &b)
    {
        if (a.size() != b.size()) return false;
        for (const auto &[key, val_a] : a) {
            auto it = b.find(key);
            if (it == b.end() || !taylor_model::interval_equal(val_a, it->second)) return false;
        }
        return true;
    }

    // Generic map comparison
    template <typename T>
    static bool map_equal(const std::unordered_map<std::string, T> &a, const std::unordered_map<std::string, T> &b)
    {
        if (a.size() != b.size()) return false;
        for (const auto &[key, val_a] : a) {
            auto it = b.find(key);
            if (it == b.end() || !(val_a == it->second)) return false;
        }
        return true;
    }

    bool operator==(const taylor_model &other) const
    {
        return m_tpol == other.m_tpol && interval_equal(m_rem, other.m_rem) && map_equal(m_exp, other.m_exp)
               && map_interval_equal(m_domain, other.m_domain);
    }

    ///////////////////////////////////
    /// Custom arithmetic functions ///
    ///////////////////////////////////

    // pow method
    taylor_model pow(int exponent) const
    {
        if (exponent < 0) {
            throw std::invalid_argument("Negative integer exponents are not supported yet.");
        }

        if (exponent == 0) {
            return taylor_model::identity();
        }

        taylor_model product = 1.0 * (*this); // like "1"
        if (exponent > 1) {
            for (int i = 1; i < exponent; ++i) {
                product = product * (*this);
            }
        } else if (exponent < 0) {
            throw std::logic_error("Not yet implemented.");
            // for (int i = 0; i < exponent + 1; ++i) {
            //     product = product * 1 / (*this);
            // }
        }
        return product;
    }

}; // end of taylor_model class

} // end of namespace audi

#endif
