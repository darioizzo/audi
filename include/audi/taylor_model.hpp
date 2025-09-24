#ifndef AUDI_TM_HPP
#define AUDI_TM_HPP

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <cassert>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <unordered_set>

#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>            //for audi::abs
#include <audi/detail/overloads_vectorized.hpp> //for audi::abs
#include <audi/gdual.hpp>

#include <audi/taylor_model_bounding.hpp>
#include <audi/taylor_model_utilities.hpp>

/// Root namespace for AuDi symbols
namespace audi
{

/// Taylor model class.
/**
 * This class represents a Taylor model, containing a generalized dual number in combination with an
 * interval. This implementation is derived from Makino (1998) and is built on top of a gdual object
 * representing the Taylor polynomial (generalized dual number) part of the Taylor model.
 *
 * A key concept here is Taylor's theorem (formulation used from Makino (1998) p.80. Taylors theorem allows a
 * quantitative estimate of the error that is to be expected when approximating a function by its Taylor polynomial
 * Furthermore it even offers a way to obtain bounds for the error in practice based on bounding the \f$(n+1)\f$st
 * derivative a method that has sometimes been employed in interval calculations.
 *
 * As a result, you get \f$ \forall \vec{x} \in [\vec{a}, \vec{b}] \text{, a given order } n \text{, and an expansion
 * point} \vec{x_o} \f$:
 *
 * \f$ f(\vec{x}) \in P_{\alpha, f}(\vec{x} - \vec{x_0}) + I_{\alpha, f} \f$
 *
 * where f is the function you're creating a Taylor model for, P is the Taylor polynomial, and I is
 * the interval remainder.
 *
 * TODO: multiply() needs to be templated as well, but raises error with ambiguity.
 * TODO: operator+ etc. need to become taylor_model_if_enabled (analogous to gdual_if_enabled) to specify exactly what
 * types the operators can accept
 */
class taylor_model
{

private:
    // The Taylor polynomial
    audi::gdual<double> m_tpol;
    // The order of the Taylor polynomial
    unsigned int m_order;
    // The dimension
    unsigned int m_ndim;
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

        if (domain_set != exp_point_set) {
            throw std::invalid_argument("Domain and expansion point maps have different items.");
        }

        for (const std::string &key : domain_set) {
            auto exp_it = m_exp.find(key);
            auto dom_it = m_domain.find(key);

            if (exp_it == m_exp.end() || dom_it == m_domain.end()) {
                throw std::invalid_argument("Key not found in expansion or domain map: " + key);
            }

            if (!boost::numeric::in(exp_it->second, dom_it->second)) {
                throw std::invalid_argument("The expansion point falls outside of the valid domain for key: " + key);
            }
        }
        // For ndim
        if (domain_set.size() != m_ndim) {
            throw std::invalid_argument("The number of variables in the domain is not equal to the dimension.");
        }
    }

public:
    ////////////////
    /// Constructors
    ////////////////

    /// Defaulted copy constructor & assignment
    taylor_model(const taylor_model &) = default;
    taylor_model &operator=(const taylor_model &) = default;

    /// Defaulted move constructor & assignment
    taylor_model(taylor_model &&) = default;
    taylor_model &operator=(taylor_model &&) = default;

    /// Default constructor
    taylor_model() : m_tpol(0.), m_order(0u), m_rem(0.) {}
    ~taylor_model() {}

    // Standard complete construction
    explicit taylor_model(const audi::gdual<double> &tpol, const int_d &rem_bound, const var_map_d &exp_point,
                          const var_map_i &domain)
    {

        m_tpol = tpol;
        m_order = tpol.get_order();
        m_ndim = static_cast<unsigned int>(tpol.get_symbol_set_size());
        m_rem = rem_bound;
        m_exp = exp_point;
        m_domain = domain;

        check_input_validity();
    }

    // Construction of Taylor model representing a constant
    taylor_model(double constant)
    {
        m_tpol = audi::gdual<double>(constant);
        m_order = 0;
        m_ndim = 0;
        m_rem = int_d(0.0, 0.0);
        m_exp = {};
        m_domain = {};
    }

    // Specific member function for a taylor model representing the constant 1.0
    // Identical to taylor_model(1.0)
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

    // Specific member function for a taylor model representing the constant 1.0 with custom remainder bound, and
    // variables. Identical to taylor_model(1.0).
    static taylor_model identity(int_d rem, var_map_d exp, var_map_i domain)
    {
        taylor_model tm;
        tm.m_tpol = audi::gdual<double>(1.0);
        tm.m_order = 0u;
        tm.m_ndim = 0u;
        tm.m_rem = rem;
        tm.m_exp = exp;
        tm.m_domain = domain;
        return tm;
    }

    ///////////
    /// Getters
    ///////////

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

    const unsigned int &get_order() const
    {
        return m_order;
    }

    const unsigned int &get_ndim() const
    {
        return m_ndim;
    }

    ///////////
    /// Setters
    ///////////

    void set_tpol(const audi::gdual<double> &tpol)
    {
        m_tpol = tpol;
        m_order = tpol.get_order();
        m_ndim = static_cast<unsigned int>(tpol.get_symbol_set_size());
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

    void set_dom(const var_map_i &domain)
    {
        m_domain = domain;
        check_input_validity();
    }

    void extend_symbol_set(const std::vector<std::string> &sym_vars, const var_map_d &exp_point,
                           const var_map_i &domain)
    {
        std::vector<std::string> d_sym_vars;
        d_sym_vars.reserve(sym_vars.size());

        for (const auto &var : sym_vars) {
            d_sym_vars.push_back("d" + var);

            auto it = exp_point.find(var);
            if (it != exp_point.end()) {
                m_exp.insert(*it);
            } else {
                throw std::logic_error("A symbol is passed that is not present in the exp_point argument.");
            }

            auto it2 = domain.find(var);
            if (it2 != domain.end()) {
                m_domain.insert(*it2);
            } else {
                throw std::logic_error("A symbol is passed that is not present in the domain argument.");
            }
        }

        m_tpol.extend_symbol_set(d_sym_vars);
        m_ndim = static_cast<unsigned int>(m_tpol.get_symbol_set_size());
        check_input_validity();
    }

    /////////////////
    /// Miscellaneous
    /////////////////

private:
    /// Increase the order of a gdual polynomial
    /**
     * Constructs a new gdual with the specified higher order and substitutes
     * the given polynomial into it. This effectively increases the truncation
     * order of the input gdual.
     *
     * @tparam T the coefficient type of the gdual
     * @param poly the input gdual polynomial
     * @param new_order the new truncation order to assign
     *
     * @return a gdual polynomial with the specified higher order
     */
    template <typename T>
    static audi::gdual<T> get_increased_order(const audi::gdual<T> &poly, unsigned int new_order)
    {
        // Create a dummy polynomial with a higher order
        audi::gdual<T> temp(0.0, "temp", new_order);

        // Substitute poly into temp
        return temp.subs("dtemp", poly);
    }

    /// Flatten a nested vector
    /**
     * Converts a vector of vectors into a single vector by concatenating
     * all inner vectors in order.
     *
     * @tparam T the element type of the vector
     * @param orig the original nested vector
     *
     * @return a flat vector containing all elements of the inner vectors
     */
    template <typename T>
    static std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
    {
        std::vector<T> ret;
        for (const auto &v : orig)
            ret.insert(ret.end(), v.begin(), v.end());
        return ret;
    }

public:
    /// Compute the bounds of the polynomial over the stored domain
    /**
     * Evaluates the current gdual polynomial over its associated domain
     * and returns the interval containing its minimum and maximum values.
     *
     * @return an int_d interval representing the lower and upper bounds
     */
    int_d get_bounds() const
    {
        return get_bounds(m_tpol, m_exp, m_domain);
    }

    /// Compute the bounds of a gdual polynomial over a given domain
    /**
     * Computes the Bernstein enclosure of a gdual polynomial on a specified domain.
     * The domain is shifted relative to the expansion points before evaluation.
     * The result is the interval containing its minimum and maximum values.
     *
     * @tparam T the coefficient type of the gdual
     * @param tpol the input gdual polynomial
     * @param exp_points the expansion points (variable → double)
     * @param domain the domain of variables as intervals (variable → int_d)
     *
     * @return an int_d interval representing the minimum and maximum polynomial values
     */
    template <typename T>
    static int_d get_bounds(const audi::gdual<T> &tpol, const var_map_d &exp_points, const var_map_i &domain)
    {
        // Build shifted domain relative to expansion points
        var_map_i domain_shifted;

        auto [coeffs, exps] = audi::get_poly(tpol);
        auto ndim = audi::get_ndim(coeffs, exps); // necessary because static function

        std::vector<T> flat;
        if (ndim > 1 || ndim == 1) {

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
        } else if (ndim == 0) {
            return int_d(tpol.constant_cf(), tpol.constant_cf());
        } else {
            throw std::runtime_error("The dimension cannot be negative.");
        }

        auto [min_val, max_val] = std::minmax_element(flat.begin(), flat.end());

        return int_d(*min_val, *max_val);
    }

    //////////////////////////////////////////
    /// Static member functions for arithmetic
    //////////////////////////////////////////

    /// Merge two maps with consistency checking
    /**
     * Combines two unordered maps into one by taking the union of their key–value pairs.
     * If a key is present in both maps, the values must be equal (checked using
     * `interval_equal`), otherwise a runtime error is thrown.
     *
     * If the key exists in only one map, that key–value pair is added directly.
     *
     * @tparam K the key type of the maps
     * @tparam V the value type of the maps
     * @param a the first unordered map
     * @param b the second unordered map
     *
     * @throws std::runtime_error if a key is present in both maps with conflicting values
     *
     * @return a new unordered_map containing the merged key–value pairs
     */
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

    /// Filter a map by a given set of keys with subset checking
    /**
     * Constructs a new unordered map containing only the entries from the input
     * map `a` whose keys are listed in `symbol_set`.
     *
     * Each key in `symbol_set` must be present in `a`. If any key is missing,
     * the function throws a logic error, since the symbol set is expected to be
     * a strict subset of the map's keys.
     *
     * @tparam K the key type of the map
     * @tparam V the value type of the map
     * @param symbol_set the set of keys to extract (must all exist in `a`)
     * @param a the input unordered map
     *
     * @throws std::logic_error if a key in `symbol_set` does not exist in `a`
     *
     * @return a new unordered_map containing only the key–value pairs from `a`
     *         whose keys are listed in `symbol_set`
     */
    template <typename K, typename V>
    static std::unordered_map<K, V> trim_map(const std::vector<K> &symbol_set, const std::unordered_map<K, V> &a)
    {
        std::unordered_map<K, V> result;
        result.reserve(symbol_set.size()); // optimization

        for (const K &symb : symbol_set) {
            auto it = a.find(symb);
            if (it != a.end()) {
                result.emplace(symb, it->second);
            } else {
                throw std::logic_error("The symbol set is supposed to be a subset of the map provided, "
                                       "but it is not because a symbol was found in the symbol set "
                                       "that doesn't exist in the map.");
            }
        }
        return result;
    }

    /// Add two objects, producing a new taylor_model
    /**
     * Adds a taylor_model to another object and returns the result as a new taylor_model.
     * Supported combinations include:
     *
     * - taylor_model + taylor_model
     * - scalar (arithmetic type) + taylor_model
     * - taylor_model + scalar (arithmetic type)
     * - interval + taylor_model
     * - taylor_model + interval
     *
     * In all cases, the domains and expansion points are merged consistently.
     *
     * @param d1 the first operand
     * @param d2 the second operand
     *
     * @return a new taylor_model representing the sum
     */

    // Add two Taylor models
private:
    static taylor_model add(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d common_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i common_domain = get_common_map(d1.m_domain, d2.m_domain);
        audi::gdual<double> new_tpol = d1.m_tpol + d2.m_tpol;
        var_map_d new_exp = trim_map(new_tpol.get_symbol_set(), common_exp);
        var_map_i new_domain = trim_map(new_tpol.get_symbol_set(), common_domain);
        return taylor_model(new_tpol, d1.m_rem + d2.m_rem, new_exp, new_domain);
    }

    // Add value (double, int, etc.) to Taylor model
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 + d2.m_tpol, d2.m_rem, d2.m_exp, d2.m_domain);
    }

    // Add Taylor model to value (double, int, etc.)
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol + d2, d1.m_rem, d1.m_exp, d1.m_domain);
    }

    // Add interval to Taylor model
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const boost::numeric::interval<T, Policies> &d1, const taylor_model &d2)
    {
        return taylor_model(d2.m_tpol, d1 + d2.m_rem, d2.m_exp, d2.m_domain);
    }

    // Add Taylor model to interval
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model add(const taylor_model &d1, const boost::numeric::interval<T, Policies> &d2)
    {
        return taylor_model(d1.m_tpol, d1.m_rem + d2, d1.m_exp, d1.m_domain);
    }

    /// Subtract two objects, producing a new taylor_model
    /**
     * Subtracts a taylor_model from another object and returns the result as a new taylor_model.
     * Supported combinations include:
     *
     * - taylor_model - taylor_model
     * - scalar (arithmetic type) - taylor_model
     * - taylor_model - scalar (arithmetic type)
     * - interval - taylor_model
     * - taylor_model - interval
     *
     * In all cases, the domains and expansion points are merged consistently.
     *
     * @param d1 the left operand
     * @param d2 the right operand
     *
     * @return a new taylor_model representing the difference
     */

    // Subtract Taylor model from another Taylor model
    static taylor_model sub(const taylor_model &d1, const taylor_model &d2)
    {

        // Get union of domains
        var_map_d common_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i common_domain = get_common_map(d1.m_domain, d2.m_domain);
        audi::gdual<double> new_tpol = d1.m_tpol - d2.m_tpol;
        var_map_d new_exp = trim_map(new_tpol.get_symbol_set(), common_exp);
        var_map_i new_domain = trim_map(new_tpol.get_symbol_set(), common_domain);
        return taylor_model(new_tpol, d1.m_rem - d2.m_rem, new_exp, new_domain);
    }

    // Subtract Taylor model from value (double, int, etc.)
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 - d2.m_tpol, d2.m_rem, d2.m_exp, d2.m_domain);
    }

    // Subtract value (double, int, etc.) from Taylor model
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol - d2, d1.m_rem, d1.m_exp, d1.m_domain);
    }

    // Subtract Taylor model from interval
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const boost::numeric::interval<T, Policies> &d1, const taylor_model &d2)
    {
        return taylor_model(d2.m_tpol, d1 - d2.m_rem, d2.m_exp, d2.m_domain);
    }

    // Subtract interval from Taylor model
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model sub(const taylor_model &d1, const boost::numeric::interval<T, Policies> &d2)
    {
        return taylor_model(d1.m_tpol, d1.m_rem - d2, d1.m_exp, d1.m_domain);
    }

    /// Multiply two objects, producing a new taylor_model
    /**
     * Multiplies a taylor_model with another object and returns the result
     * as a new taylor_model.
     *
     * Supported combinations include:
     *
     * - taylor_model × taylor_model
     * - scalar (arithmetic type) × taylor_model
     * - taylor_model × scalar (arithmetic type)
     * - interval × taylor_model
     * - taylor_model × interval
     *
     * For taylor_model × taylor_model:
     * - A higher-order untruncated polynomial product is formed.
     * - The result is truncated to the maximum degree of the operands.
     * - The truncation error is enclosed in an interval remainder using bounds.
     *
     * In all cases, domains and expansion points are merged consistently.
     *
     * @param d1 the left operand
     * @param d2 the right operand
     *
     * @return a new taylor_model representing the product
     */

    // Multiply two Taylor models
    static taylor_model multiply(const taylor_model &d1, const taylor_model &d2)
    {

        if (d1.get_tpol() == audi::gdual<double>(0.0, "irrelevant", 0)
            || d2.get_tpol() == audi::gdual<double>(0.0, "irrelevant", 0)) {
            return taylor_model(0.0);
        }

        // Get union of domains
        var_map_d common_exp = get_common_map(d1.m_exp, d2.m_exp);
        var_map_i common_domain = get_common_map(d1.m_domain, d2.m_domain);

        // Get untruncated polynomial product
        unsigned int new_order = d1.m_tpol.get_order() + d2.m_tpol.get_order();

        audi::gdual<double> d1_ho = get_increased_order(d1.m_tpol, new_order);
        audi::gdual<double> d2_ho = get_increased_order(d2.m_tpol, new_order);

        audi::gdual<double> polynomial_product = d1_ho * d2_ho;

        // Get agreeable polynomial
        const unsigned int agreeable_degree = std::max(d1.m_tpol.get_order(), d2.m_tpol.get_order());

        audi::gdual<double> agreeable_pol(0.0);
        for (unsigned int i = 0; i <= agreeable_degree; ++i) {
            agreeable_pol += polynomial_product.extract_terms(i);
        }
        var_map_d new_exp = trim_map(agreeable_pol.get_symbol_set(), common_exp);
        var_map_i new_domain = trim_map(agreeable_pol.get_symbol_set(), common_domain);

        audi::gdual<double> pol_e = polynomial_product - agreeable_pol;

        int_d pol_e_bounds;
        if (pol_e.is_zero(1e-16)) {
            pol_e_bounds = int_d(0.0, 0.0);
        } else {
            // extract symbols used in pol_e
            std::vector<std::string> symbs_e = pol_e.get_symbol_set();

            // extract relevant domain and expansion points
            var_map_i domain_e;
            var_map_d exp_point_e;
            for (const auto &sym : symbs_e) {
                domain_e[sym] = common_domain[sym];
                exp_point_e[sym] = common_exp[sym];
            }

            pol_e_bounds = get_bounds(pol_e, exp_point_e, domain_e);
        }

        int_d interval_of_product
            = pol_e_bounds + d1.get_bounds() * d2.m_rem + d1.m_rem * d2.get_bounds() + d1.m_rem * d2.m_rem;

        return taylor_model(agreeable_pol, interval_of_product, new_exp, new_domain);
    }

    // Multiply value (double, int, etc.) with Taylor model
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const T &d1, const taylor_model &d2)
    {
        return taylor_model(d1 * d2.m_tpol, int_d(d1, d1) * d2.m_rem, d2.m_exp, d2.m_domain);
    }

    // Multiply Taylor model with value (double, int, etc.)
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol * d2, d1.m_rem * int_d(d2, d2), d1.m_exp, d1.m_domain);
    }

    // Multiply interval with Taylor model
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const boost::numeric::interval<T, Policies> &d1, const taylor_model &d2)
    {
        audi::gdual<T> gd_one(1.0);
        const taylor_model tm_one(gd_one, d1, d2.m_exp, d2.m_domain);
        return tm_one * d2;
    }

    // Multiply Taylor model with interval
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model multiply(const taylor_model &d1, const boost::numeric::interval<T, Policies> &d2)
    {
        audi::gdual<T> gd_one(1.0, "x", 0);
        const taylor_model tm_one(gd_one, d2, d1.m_exp, d1.m_domain);
        return d1 * tm_one;
    }

    /// Divide two objects, producing a new taylor_model
    /**
     * Divides a taylor_model by another object and returns the result as a new taylor_model.
     *
     * Supported combinations include:
     *
     * - scalar (arithmetic type) ÷ taylor_model
     * - taylor_model ÷ taylor_model
     * - taylor_model ÷ scalar (arithmetic type)
     * - interval ÷ taylor_model
     * - taylor_model ÷ interval
     *
     * For scalar ÷ taylor_model:
     * - The denominator’s polynomial and remainder are analyzed to ensure that
     *   zero is not contained in the range (division by zero is forbidden).
     * - Bounds are computed for the remainder, and a validated enclosure of the
     *   division is returned.
     *
     * In all cases, domains and expansion points are merged consistently.
     *
     * @param d1 the numerator
     * @param d2 the denominator
     *
     * @throws std::runtime_error if the denominator’s range includes zero
     *
     * @return a new taylor_model representing the quotient
     */

    // Division between two Taylor models
    static taylor_model div(const taylor_model &d1, const taylor_model &d2)
    {
        return d1 * (1 / d2);
    }

    // Division of value (double, int, etc.) by Taylor model
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model div(const T &d1, const taylor_model &d2)
    {
        if (boost::numeric::zero_in(d2.get_bounds() + d2.get_rem())) {
            throw std::runtime_error("The range of the function includes 0, and "
                                     "therefore is ill-defined when performing a division (cannot divide by 0).");
        }
        double const_term = d2.get_tpol().constant_cf();
        audi::gdual<double> f_bar = d2.get_tpol() - const_term;
        int_d f_bar_bounds = d2.get_bounds(f_bar, d2.get_exp(), d2.get_dom());
        int_d f_bar_remainder = d2.get_rem();
        int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
        unsigned int k = d2.get_tpol().get_order();
        int_d total_rem_bound
            = std::pow(-1, k + 1) * boost::numeric::pow(f_bar_remainder_bounds, static_cast<int>(k + 1))
              / std::pow(const_term, static_cast<int>(k + 2))
              * (int_d(1.0, 1.0)
                 / boost::numeric::pow(int_d(1.0, 1.0) + f_bar_remainder_bounds / const_term * int_d(0.0, 1.0),
                                       static_cast<int>(k + 2)));
        return taylor_model(d1 / d2.get_tpol(), total_rem_bound, d2.get_exp(), d2.get_dom());
    }

    // Division of Taylor model by value (double, int, etc.)
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model div(const taylor_model &d1, const T &d2)
    {
        return taylor_model(d1.m_tpol / d2, d1.m_rem / int_d(d2, d2), d1.m_exp, d1.m_domain);
    }

    // Division of interval by Taylor model
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model div(const boost::numeric::interval<T, Policies> &d1, const taylor_model &d2)
    {
        audi::gdual<T> gd_one(1.0);
        const taylor_model tm_one(gd_one, d1, d2.m_exp, d2.m_domain);
        return tm_one * d2;
    }

    // Division of Taylor model by interval
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static taylor_model div(const taylor_model &d1, const boost::numeric::interval<T, Policies> &d2)
    {
        audi::gdual<T> gd_one(1.0, "x", 0);
        const taylor_model tm_one(gd_one, d2, d1.m_exp, d1.m_domain);
        return d1 * tm_one;
    }

public:
    //////////////////////////////////
    /// Overloaded arithmetic operator
    //////////////////////////////////

    /// Overloaded addition operator
    /**
     * Adds two objects and returns the result as a taylor_model.
     * This is equivalent to calling `add(d1, d2)`.
     *
     * @tparam T the type of the left operand
     * @tparam U the type of the right operand
     * @param d1 the left operand
     * @param d2 the right operand
     *
     * @return a new taylor_model representing the sum
     */
    template <typename T, typename U>
    friend taylor_model operator+(const T &d1, const U &d2)
    {
        return add(d1, d2);
    }

    /// Add and assign operator
    /**
     * Adds an object to this taylor_model in place.
     * This is equivalent to `*this = *this + d1`.
     *
     * @tparam T the type of the right operand
     * @param d1 the right operand
     *
     * @return reference to this taylor_model after modification
     */
    /// Add and assignment operator
    template <typename T>
    auto operator+=(const T &d1) -> decltype(*this = *this + d1)
    {
        return *this = *this + d1;
    }

    /// Overloaded subtraction operator
    /**
     * Subtracts one object from another and returns the result as a taylor_model.
     * This is equivalent to calling `sub(d1, d2)`.
     *
     * @tparam T the type of the left operand
     * @tparam U the type of the right operand
     * @param d1 the left operand
     * @param d2 the right operand
     *
     * @return a new taylor_model representing the difference
     */
    template <typename T, typename U>
    friend taylor_model operator-(const T &d1, const U &d2)
    {
        return sub(d1, d2);
    }

    /// Subtract and assign operator
    /**
     * Subtracts an object from this taylor_model in place.
     * This is equivalent to `*this = *this - d1`.
     *
     * @tparam T the type of the right operand
     * @param d1 the right operand
     *
     * @return reference to this taylor_model after modification
     */
    template <typename T>
    auto operator-=(const T &d1) -> decltype(*this = *this - d1)
    {
        return *this = *this - d1;
    }

    /// Unary negation operator
    /**
     * Returns the negation of this taylor_model.
     * This is equivalent to `sub(0.0, *this)`.
     *
     * @return a new taylor_model representing the negated value
     */
    taylor_model operator-() const
    {
        return sub(0.0, *this);
    }

    /// Overloaded multiplication operator
    /**
     * Multiplies two objects and returns the result as a taylor_model.
     * This is equivalent to calling `multiply(d1, d2)`.
     *
     * @tparam T the type of the left operand
     * @tparam U the type of the right operand
     * @param d1 the left operand
     * @param d2 the right operand
     *
     * @return a new taylor_model representing the product
     */
    template <typename T, typename U>
    friend taylor_model operator*(const T &d1, const U &d2)
    {
        return multiply(d1, d2);
    }

    /// Multiply and assign operator
    /**
     * Multiplies this taylor_model by an object in place.
     * This is equivalent to `*this = *this * d1`.
     *
     * @tparam T the type of the right operand
     * @param d1 the right operand
     *
     * @return reference to this taylor_model after modification
     */
    template <typename T>
    auto operator*=(const T &d1) -> decltype(*this = *this * d1)
    {
        return *this = *this * d1;
    }

    /// Overloaded division operator
    /**
     * Divides one object by another and returns the result as a taylor_model.
     * This is equivalent to calling `div(d1, d2)`.
     *
     * @tparam T the type of the numerator
     * @tparam U the type of the denominator
     * @param d1 the numerator
     * @param d2 the denominator
     *
     * @return a new taylor_model representing the quotient
     */
    template <typename T, typename U>
    friend taylor_model operator/(const T &d1, const U &d2)
    {
        return div(d1, d2);
    }

    /// Divide and assign operator
    /**
     * Divides this taylor_model by an object in place.
     * This is equivalent to `*this = *this / d1`.
     *
     * @tparam T the type of the right operand
     * @param d1 the right operand (denominator)
     *
     * @return reference to this taylor_model after modification
     */
    template <typename T>
    auto operator/=(const T &d1) -> decltype(*this = *this / d1)
    {
        return *this = *this / d1;
    }

    /////////////////////////////////
    /// Quantity comparison utilities
    /////////////////////////////////

    /// Compare two scalar values for equality
    /**
     * Checks whether two arithmetic values are equal within a specified tolerance.
     * For floating-point types, if no tolerance is provided, a default of
     * 10 × machine epsilon is used.
     * For integral types, strict equality is used.
     *
     * @tparam T the scalar type (must be arithmetic)
     * @param a the first value
     * @param b the second value
     * @param tol optional tolerance for floating-point comparison
     *
     * @return true if the values are considered equal, false otherwise
     */
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static bool interval_equal(const T &a, const T &b, std::optional<double> tol = std::nullopt)
    {
        if constexpr (std::is_floating_point_v<T>) {
            if (!tol.has_value()) {
                tol = std::numeric_limits<T>::epsilon() * 10;
            }
            return std::fabs(a - b) < tol;
        } else {
            return a == b;
        }
    }

    /// Compare two intervals for equality
    /**
     * Checks whether two boost::numeric::interval values are equal within a
     * specified tolerance.
     * For floating-point types, lower and upper bounds are compared with tolerance.
     * For integral types, strict equality of bounds is used.
     *
     * @tparam T the interval bound type (must be arithmetic)
     * @param a the first interval
     * @param b the second interval
     * @param tol optional tolerance for floating-point comparison
     *
     * @return true if both intervals are considered equal, false otherwise
     */
    template <typename T, typename Policies, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static bool interval_equal(const boost::numeric::interval<T, Policies> &a,
                               const boost::numeric::interval<T, Policies> &b, std::optional<double> tol = std::nullopt)
    {
        if constexpr (std::is_floating_point_v<T>) {
            if (!tol.has_value()) {
                tol = std::numeric_limits<T>::epsilon() * 10;
            }
            return std::fabs(a.lower() - b.lower()) < tol && std::fabs(a.upper() - b.upper()) < tol;
        } else {
            return a.lower() == b.lower() && a.upper() == b.upper();
        }
    }

    /// Compare two maps of intervals for equality
    /**
     * Checks whether two unordered_maps mapping variable names to int_d intervals
     * are equal within a specified tolerance.
     * Both maps must have the same size and the same keys.
     *
     * @param a the first map
     * @param b the second map
     * @param tol optional tolerance for floating-point comparison of interval bounds
     *
     * @return true if the maps are considered equal, false otherwise
     */
    static bool map_interval_equal(const std::unordered_map<std::string, int_d> &a,
                                   const std::unordered_map<std::string, int_d> &b,
                                   std::optional<double> tol = std::nullopt)
    {
        if (a.size() != b.size()) return false;
        for (const auto &[key, val_a] : a) {
            auto it = b.find(key);
            if (it == b.end() || !taylor_model::interval_equal(val_a, it->second, tol)) return false;
        }
        return true;
    }

    /// Compare two maps of scalar values for equality
    /**
     * Checks whether two unordered_maps mapping variable names to arithmetic values
     * are equal within a specified tolerance.
     * Both maps must have the same size and the same keys.
     * If no tolerance is provided, strict equality is used.
     *
     * @tparam T the value type (must be arithmetic)
     * @param a the first map
     * @param b the second map
     * @param tol optional tolerance for floating-point comparison
     *
     * @return true if the maps are considered equal, false otherwise
     */
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    static bool map_equal(const std::unordered_map<std::string, T> &a, const std::unordered_map<std::string, T> &b,
                          std::optional<double> tol = std::nullopt)
    {
        if (a.size() != b.size()) return false;
        for (const auto &[key, val_a] : a) {
            auto it = b.find(key);

            if (!tol.has_value()) {
                if (it == b.end() || !(val_a == it->second)) return false;
            } else {
                if (it == b.end() || !(std::abs(val_a - it->second) < tol)) return false;
            }
        }
        return true;
    }

    /////////////////////////////
    /// Overloaded other operator
    /////////////////////////////

    /// Stream insertion operator for taylor_model
    /**
     * Outputs the contents of a taylor_model to the given output stream in a
     * human-readable format.
     *
     * @param os the output stream
     * @param tm the taylor_model to output
     *
     * @return a reference to the output stream with the formatted contents inserted
     */
    friend std::ostream &operator<<(std::ostream &os, const taylor_model &tm)
    {
        os << "Taylor Polynomial:\n"
           << tm.get_tpol() << "\nRemainder Bound: [" << tm.get_rem().lower() << ", " << tm.get_rem().upper() << "]"
           << "\nExpansion Point: " << tm.get_exp() << "\nDomain: " << tm.get_dom();
        return os;
    }
    /// Equality operator for taylor_model
    /**
     * Compares two taylor_model instances for equality.
     * Two taylor_models are considered equal if:
     * - Their polynomials are equal (`m_tpol`)
     * - Their remainder intervals are equal (`interval_equal`)
     * - Their expansion points are equal (`map_equal`)
     * - Their domains are equal (`map_interval_equal`)
     *
     * @param other the taylor_model to compare with
     *
     * @return true if all components are equal, false otherwise
     */
    bool operator==(const taylor_model &other) const
    {
        return m_tpol == other.m_tpol && interval_equal(m_rem, other.m_rem) && map_equal(m_exp, other.m_exp)
               && map_interval_equal(m_domain, other.m_domain);
    }

    /// Inequality operator for taylor_model
    /**
     * Compares two taylor_model instances for inequality.
     * Returns the negation of the equality operator.
     *
     * @param other the taylor_model to compare with
     *
     * @return true if any component differs, false otherwise
     */
    bool operator!=(const taylor_model &other) const
    {
        return m_tpol != other.m_tpol || !interval_equal(m_rem, other.m_rem) || !map_equal(m_exp, other.m_exp)
               || !map_interval_equal(m_domain, other.m_domain);
    }

}; // end of taylor_model class

} // end of namespace audi

#endif
