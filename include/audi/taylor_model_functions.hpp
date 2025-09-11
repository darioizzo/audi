#ifndef AUDI_TM_FUNCTIONS_HPP
#define AUDI_TM_FUNCTIONS_HPP

#include <audi/config.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <cmath>
#include <stdexcept>

#include <audi/back_compatibility.hpp>
#include <audi/detail/overloads.hpp>
#include <audi/detail/overloads_vectorized.hpp>
#include <audi/gdual.hpp>
#include <audi/taylor_model.hpp>
#include <audi/type_traits.hpp>
#include <audi/vectorized.hpp>

// using int_d = boost::numeric::interval<double>;
using int_d = boost::numeric::interval<
    double, boost::numeric::interval_lib::policies<
                boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double>>,
                boost::numeric::interval_lib::checking_no_nan<double>>>;
using var_map_d = std::unordered_map<std::string, double>;
using var_map_i = std::unordered_map<std::string, int_d>;

namespace audi
{

/// Overload for the double factorial function
/**
 * Computes the double factorial of an integer \p n, denoted as \f$n!!\f$.
 * The double factorial is defined as:
 *
 * \f[
 * n!! =
 * \begin{cases}
 *   1 & \text{if } n = 0 \text{ or } n = -1, \\
 *   n \cdot (n-2) \cdot (n-4) \cdot \ldots \cdot 2 & \text{if } n > 0 \text{ and even}, \\
 *   n \cdot (n-2) \cdot (n-4) \cdot \ldots \cdot 1 & \text{if } n > 0 \text{ and odd}.
 * \end{cases}
 * \f]
 *
 * @param n the integer argument (must satisfy \f$ n \geq -1 \f$)
 *
 * @throws std::invalid_argument if \p n < -1
 *
 * @return the double factorial of \p n
 */
constexpr long long factorial2(int n)
{
    if (n < -1) {
        throw std::invalid_argument("factorial2 is undefined for n < -1");
    }
    if (n == 0 || n == -1) {
        return 1;
    }
    long long result = 1;
    for (int k = n; k > 0; k -= 2) {
        result *= k;
    }
    return result;
}

/// Overload for the power of a taylor_model with integer exponent
/**
 * Computes the power of a taylor_model with an integer exponent. The Taylor polynomial is
 * computed as \f$(T_f)^n\f$ using repeated multiplication, while the remainder is propagated
 * according to the interval arithmetic rules of the taylor_model class.
 *
 * @param tm the taylor_model base
 * @param n integer exponent
 *
 * @return a new taylor_model representing \f$tm^n\f$
 */
inline audi::taylor_model pow(const audi::taylor_model &tm, int n)
{

    if (n == 0) {
        return audi::taylor_model::identity(tm.get_rem(), tm.get_exp(), tm.get_dom());
    }

    audi::taylor_model product;
    if (n > 1) {
        product = 1.0 * tm; // like "1"
        for (int i = 1; i < n; ++i) {
            product = product * tm;
        }
    } else if (n < 0) {
        product = 1.0 / tm; // like "1"
        for (int i = 0; i < -n - 1; ++i) {
            product = product * 1.0 / tm;
        }
    }
    return product;
}

/// \overload
/**
 * Overload for the sqrt and inverse sqrt of a taylor_model with a real-valued arithmetic exponent. A different, generic
 * real-valued exponent is not implemented.
 *
 * The remainder for the square root is:
 *
 * \f[
 * R = (-1)^k \sqrt(f_0) \cdot \frac{(2k - 1)!}{(k + 1)!2^{k+1}} \frac{\hat{f}^{k+1}}{f_0^{k+1}}
 * \frac{1}{(1 + \theta \cdot \frac{\hat{f}}{f_0})^{k + \frac{1}{2}}}
 * \f]
 *
 * And for the inverse square root:
 *
 * \f[
 * R = (-1)^{k+1} \frac{1}{\sqrt(f_0)} \cdot \frac{(2k + 1)!!}{(k + 1)!2^{k+1}} \frac{\hat{f}^{k+1}}{f_0^{k+1}}
 * \frac{1}{(1 + \theta \cdot \frac{\hat{f}}{f_0})^{k + \frac{3}{2}}}
 * \f]
 *
 * where \f$ \forall \vec{x} \in [\vec{a}, \vec{b}] \text{, } B(\hat{f}) + I_f \subset (0, \infty) \f$
 *
 * @tparam T arithmetic type of the exponent (e.g. double, float)
 * @param tm the taylor_model base
 * @param n real-valued exponent
 *
 * @return a new taylor_model representing \f$tm^n\f$
 */
template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
inline audi::taylor_model pow(const audi::taylor_model &tm, T n)
{

    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d total_rem_bound;
    if (n == 1.0 / 2.0) {
        if (!boost::numeric::subset(tm.get_bounds() + tm.get_rem(),
                                    int_d(0, std::numeric_limits<double>::infinity()))) {
            throw std::runtime_error("The range of the function is not a subset of the positive Real numbers, and "
                                     "therefore is ill-defined when performing a sqrt operation.");
        }

        int_d k_int_denominator
            = boost::numeric::pow(int_d(1.0) + (f_bar_remainder_bounds / const_term) * int_d(0, 1), 2 * k + 1);
        int_d k_half_denominator = boost::numeric::sqrt(k_int_denominator);
        total_rem_bound = std::pow(-1, k) * audi::pow(const_term, 1.0 / 2.0)
                          * static_cast<double>(factorial2(2 * k - 1))
                          / std::pow(boost::math::factorial<double>(k + 1) * 2, k + 1)
                          * boost::numeric::pow(f_bar_remainder_bounds, k + 1) / std::pow(const_term, k + 1)
                          * int_d(1.0) / k_half_denominator;
    } else if (n == -1.0 / 2.0) {
        if (!boost::numeric::subset(tm.get_bounds() + tm.get_rem(),
                                    int_d(0, std::numeric_limits<double>::infinity()))) {
            throw std::runtime_error("The range of the function is not a subset of the positive Real numbers, and "
                                     "therefore is ill-defined when performing an inverse sqrt operation.");
        }
        int_d k_int_denominator
            = boost::numeric::pow(int_d(1.0) + (f_bar_remainder_bounds / const_term) * int_d(0, 1), 2 * k + 3);
        int_d k_half_denominator = boost::numeric::sqrt(k_int_denominator);
        total_rem_bound = std::pow(-1, k + 1) * audi::pow(const_term, 1.0 / 2.0)
                          * static_cast<double>(factorial2(2 * k + 1))
                          / std::pow(boost::math::factorial<double>(k + 1) * 2, k + 1)
                          * boost::numeric::pow(f_bar_remainder_bounds, k + 1) / std::pow(const_term, k + 1)
                          * int_d(1.0) / k_half_denominator;
    } else if (n - std::round(n) < std::numeric_limits<T>::epsilon() * 10) {
        return audi::pow(tm, static_cast<int>(n));
    } else {
        throw std::runtime_error("This exponent is not implemented for Taylor model arithmetic.");
    }
    return audi::taylor_model(audi::pow(tm.get_tpol(), n), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Overload for the exponential of a taylor_model
/**
 * Computes the exponential of a taylor_model. For this, the Taylor polynomial is calculated using
 * audi::exp(const gdual<T, M> &d). The remainder term is represented as follows:
 *
 * \f[
 * R = \exp(f_0) \frac{\hat{f}^{k+1}}{(k+1)!} \exp(\theta \cdot \hat{f})
 * \f]
 * for any \f$ x \in [\vec{a}, \vec{b}] \text{, } \hat{f} \in B(\hat{f}) \text{, and } 0 < \theta <
 * 1 \f$.
 *
 * @param tm the taylor_model whose exponential is to be computed
 *
 * @return a new taylor_model representing \f$ exp(tm) \f$
 */
inline audi::taylor_model exp(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d total_rem_bound = std::pow(std::exp(1.0), const_term) * (1.0 / boost::math::factorial<double>(k + 1))
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1)
                            * boost::numeric::exp(f_bar_remainder_bounds * int_d(0, 1));

    return audi::taylor_model(audi::exp(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Overload for the log of a taylor_model.
/**
 * Function calculating the log of a taylor_model.
 *
 * The remainder for the square root is:
 *
 * \f[
 * R = (-1)^{k+2} \frac{1}{k+1} \cdot \frac{\hat{f}^{k+1}}{f_0^{k+1}}
 * \frac{1}{(1 + \theta \cdot \frac{\hat{f}}{f_0})^{k + 1}}
 * \f]
 *
 * where \f$ \forall \vec{x} \in [\vec{a}, \vec{b}] \text{, } B(\hat{f}) + I_f \subset (0, \infty) \f$
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$log(tm)\f$
 */
inline audi::taylor_model log(const audi::taylor_model &tm)
{
    if (!boost::numeric::subset(tm.get_bounds() + tm.get_rem(), int_d(0, std::numeric_limits<double>::infinity()))) {
        throw std::runtime_error("The range of the function is not a subset of the positive Real numbers, and "
                                 "therefore is ill-defined when performing a log operation.");
    }
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d total_rem_bound
        = int_d(-1 * (k + 2) * (1.0 / (k + 1.0))) * boost::numeric::pow(f_bar_remainder_bounds, k + 1)
          / std::pow(const_term, k + 1)
          * boost::numeric::pow(int_d(1.0) / (int_d(1.0) + (f_bar_remainder_bounds / const_term) * int_d(0, 1)), k + 1);

    return audi::taylor_model(audi::log(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/**
 * Overload for the sqrt of a taylor_model. This function runs `audi::pow(tm, 1/2)`.
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$sqrt(tm)\f$
 */
inline audi::taylor_model sqrt(const audi::taylor_model &tm)
{
    return audi::pow(tm, 1.0 / 2.0);
}

/// Overload for the sine of a taylor_model.
/**
 * Function that calculates the sine of a taylor_model.
 *
 * The remainder for the sine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} \bigl(\bar f(\vec{x})\bigr)^{k+1} \cdot J ,
 * \f]
 *
 * where
 *
 * \f[
 * J =
 * \begin{cases}
 *   -J_0 & \text{if } \operatorname{mod}(k, 4) = 1, 2, \\
 *   J_0  & \text{else},
 * \end{cases}
 * \f]
 *
 * and
 *
 * \f[
 * J_0 =
 * \begin{cases}
 *   \cos\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{if } k \text{ is even}, \\
 *   \sin\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{else}.
 * \end{cases}
 * \f]
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$sin(tm)\f$
 */
inline audi::taylor_model sin(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func;
    if (k % 2 == 0) {
        conditional_func = boost::numeric::cos(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func = boost::numeric::sin(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d conditional_func_2;
    if (k % 4 == 1 || k % 4 == 2) {
        conditional_func_2 = -conditional_func;
    } else {
        conditional_func_2 = conditional_func;
    }

    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_2
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1);

    return audi::taylor_model(audi::sin(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Overload for the cos of a taylor_model.
/**
 * Function that calculates the cos of a taylor_model.
 *
 * The remainder for the cosine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} \bigl(\bar f(\vec{x})\bigr)^{k+1} \cdot J ,
 * \f]
 *
 * where
 *
 * \f[
 * J =
 * \begin{cases}
 *   -J_0 & \text{if } \operatorname{mod}(k, 4) = 0, 1, \\
 *   J_0  & \text{else},
 * \end{cases}
 * \f]
 *
 * and
 *
 * \f[
 * J_0 =
 * \begin{cases}
 *   \sin\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{if } k \text{ is even}, \\
 *   \cos\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{else}.
 * \end{cases}
 * \f]
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$cos(tm)\f$
 */
inline audi::taylor_model cos(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func;
    if (k % 2 == 0) {
        conditional_func = boost::numeric::sin(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func = boost::numeric::cos(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d conditional_func_2;
    if (k % 4 == 0 || k % 4 == 1) {
        conditional_func_2 = -conditional_func;
    } else {
        conditional_func_2 = conditional_func;
    }

    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_2
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1);

    return audi::taylor_model(audi::cos(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Function for both sine and cosine of a taylor_model
/**
 * Function that calculates both the sine and cosine of a taylor_model. Most computations are
 * overlapping, so this saves time. This structure is inspired by sin_and_cos(const gdual<T, M> &d).
 *
 * @param tm the taylor_model base
 *
 * @return an array of two new taylor_models representing \f${sin(tm), cos(tm)}\f$
 */
inline std::array<audi::taylor_model, 2> sin_and_cos(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func_sin;
    int_d conditional_func_cos;
    if (k % 2 == 0) {
        conditional_func_sin = boost::numeric::cos(const_term + f_bar_remainder_bounds * int_d(0, 1));
        conditional_func_cos = boost::numeric::sin(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func_sin = boost::numeric::sin(const_term + f_bar_remainder_bounds * int_d(0, 1));
        conditional_func_cos = boost::numeric::cos(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d conditional_func_2_sin;
    if (k % 4 == 1 || k % 4 == 2) {
        conditional_func_2_sin = -conditional_func_sin;
    } else {
        conditional_func_2_sin = conditional_func_sin;
    }

    int_d conditional_func_2_cos;
    if (k % 4 == 0 || k % 4 == 1) {
        conditional_func_2_cos = -conditional_func_cos;
    } else {
        conditional_func_2_cos = conditional_func_cos;
    }

    int_d total_rem_bound_sin = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_2_sin
                                * boost::numeric::pow(f_bar_remainder_bounds, k + 1);
    int_d total_rem_bound_cos = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_2_cos
                                * boost::numeric::pow(f_bar_remainder_bounds, k + 1);
    audi::taylor_model tm_sin
        = audi::taylor_model(audi::sin(tm.get_tpol()), total_rem_bound_sin, tm.get_exp(), tm.get_dom());
    audi::taylor_model tm_cos
        = audi::taylor_model(audi::cos(tm.get_tpol()), total_rem_bound_cos, tm.get_exp(), tm.get_dom());

    return std::array<audi::taylor_model, 2>{{std::move(tm_sin), std::move(tm_cos)}};
}

/// Overload for tangent of taylor_model
/**
 * Function that calculates the tangent of a taylor_model. This just uses the usual identity of `tan(x) = sin(x) /
 * cos(x)`
 *
 * @param tm the taylor_model base
 *
 * @return an array of two new taylor_models representing \f$tan(tm)\f$
 */
inline audi::taylor_model tan(const audi::taylor_model &tm)
{
    return audi::sin(tm) / audi::cos(tm);
}

/// Overload for the sinh of a taylor_model.
/**
 * Function that calculates the sinh of a taylor_model.
 *
 * The remainder for the hyperbolic sine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} \bigl(\bar f(\vec{x})\bigr)^{k+1} \cdot J ,
 * \f]
 *
 * where
 *
 * \f[
 * J =
 * \begin{cases}
 *   \cosh\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{if } k \text{ is even}, \\
 *   \sinh\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{else}.
 * \end{cases}
 * \f]
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$sinh(tm)\f$
 */
inline audi::taylor_model sinh(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func;
    if (k % 2 == 0) {
        conditional_func = boost::numeric::cosh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func = boost::numeric::sinh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1);

    return audi::taylor_model(audi::sinh(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Overload for the cosh of a taylor_model.
/**
 * Function that calculates the cosh of a taylor_model.
 *
 * The remainder for the hyperbolic cosine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} \bigl(\bar f(\vec{x})\bigr)^{k+1} \cdot J ,
 * \f]
 *
 * where
 *
 * \f[
 * J =
 * \begin{cases}
 *   \sinh\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{if } k \text{ is even}, \\
 *   \cosh\!\bigl(f_0 + \theta \cdot \bar f(\vec{x})\bigr) & \text{else}.
 * \end{cases}
 * \f]
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$cosh(tm)\f$
 */
inline audi::taylor_model cosh(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func;
    if (k % 2 == 0) {
        conditional_func = boost::numeric::sinh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func = boost::numeric::cosh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1);

    return audi::taylor_model(audi::cosh(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

/// Function for both sine and cosine of a taylor_model
/**
 * Function that calculates both the hyperbolic sine and cosine of a taylor_model. Most computations are
 * overlapping, so this saves time. This structure is inspired by sinh_and_cosh(const gdual<T, M> &d) and
 * analogous to sinh_and_cosh(const audi::taylor_model &tm).
 *
 * @param tm the taylor_model base
 *
 * @return an array of two new taylor_models representing \f${sinh(tm), cosh(tm)}\f$
 */
inline std::array<audi::taylor_model, 2> sinh_and_cosh(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d conditional_func_sinh;
    int_d conditional_func_cosh;
    if (k % 2 == 0) {
        conditional_func_sinh = boost::numeric::cosh(const_term + f_bar_remainder_bounds * int_d(0, 1));
        conditional_func_cosh = boost::numeric::sinh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    } else {
        conditional_func_sinh = boost::numeric::sinh(const_term + f_bar_remainder_bounds * int_d(0, 1));
        conditional_func_cosh = boost::numeric::cosh(const_term + f_bar_remainder_bounds * int_d(0, 1));
    }

    int_d total_rem_bound_sinh = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_sinh
                                 * boost::numeric::pow(f_bar_remainder_bounds, k + 1);
    int_d total_rem_bound_cosh = int_d(1.0 / boost::math::factorial<double>(k + 1)) * conditional_func_cosh
                                 * boost::numeric::pow(f_bar_remainder_bounds, k + 1);
    audi::taylor_model tm_sinh
        = audi::taylor_model(audi::sinh(tm.get_tpol()), total_rem_bound_sinh, tm.get_exp(), tm.get_dom());
    audi::taylor_model tm_cosh
        = audi::taylor_model(audi::cosh(tm.get_tpol()), total_rem_bound_cosh, tm.get_exp(), tm.get_dom());

    return std::array<audi::taylor_model, 2>{{std::move(tm_sinh), std::move(tm_cosh)}};
}

/// Overload of the hyperbolic tangent of a taylor_model
/**
 * Function that calculates the hyperbolic tangent of a taylor_model. This just uses the usual identity of `tanh(x) =
 * sinh(x) / cosh(x)`
 *
 * @param tm the taylor_model base
 *
 * @return an array of two new taylor_models representing \f$tanh(tm)\f$
 */
inline audi::taylor_model tanh(const audi::taylor_model &tm)
{
    return audi::sinh(tm) / audi::cosh(tm);
}

/// Function for higher-order derivatives of arcsin
/**
 * Computes the \f$n\f$-th derivative of the arcsine function \f$\arcsin(a)\f$
 * using a recurrence relation. This function is used internally in the
 * implementation of \c asin() for taylor_model.
 *
 * The recurrence is based on:
 *
 * \f[
 * R = \frac{d}{dx} \arcsin(x) = \frac{1}{\sqrt{1 - x^2}}
 * \f]
 *
 * and higher-order derivatives are obtained via:
 *
 * \f[
 * f^{(k+2)}(x) =
 * \frac{(2k+1)x f^{(k+1)}(x) + k^2 f^{(k)}(x)}{1 - x^2},
 * \quad k \geq 1.
 * \f]
 *
 * Base cases:
 * \f[
 * f^{(1)}(x) = \frac{1}{\sqrt{1-x^2}}, \qquad
 * f^{(2)}(x) = \frac{x}{(1-x^2)^{3/2}}.
 * \f]
 *
 * @param a the interval argument \f$a\f$
 * @param order the derivative order (must be \f$\geq 1\f$)
 *
 * @return the \p order -th derivative of \f$\arcsin(a)\f$
 *
 * @throws std::out_of_range if \p order < 1
 */
int_d asin_derivative(int_d a, int order)
{
    if (order == 0) {
        return boost::numeric::asin(a);
    }
    std::vector<int_d> derivs(order + 1);

    // Base cases
    derivs[1] = int_d(1.0) / boost::numeric::sqrt(int_d(1.0) - a * a); // arcsin'(a)
    if (order >= 2)
        derivs[2] = a * boost::numeric::pow(int_d(1.0) - a * a, 2)
                    / boost::numeric::pow(int_d(1.0) - a * a, 3); // arcsin''(a)

    // Recursive computation for higher orders
    for (int k = 1; k <= order - 2; ++k) {
        derivs[k + 2] = 1.0 / (int_d(1.0) - a * a) * (int_d(2 * k + 1) * a * derivs[k + 1] + int_d(k * k) * derivs[k]);
    }

    return derivs[order];
}

/// Overload for the asin of a taylor_model.
/**
 * Function that calculates the asin of a taylor_model.
 *
 * Here, some custom arithmetic is used from Makino (1998):
 *
 * \f$ \arcsin(f) = \arcsin(f_0) + \arcsin(g) \f$
 *
 * where \f$ g = f \cdot \sqrt{1 - f_0^2} - f_0 \cdot \sqrt{1 - f^2} \f$
 *
 * and where \f$ \forall \vec{x} \in [\vec{a}, \vec{b}] \text{, } B(f) + I_f \subset (-1, 1) \f$ applies.
 *
 * The remainder for the arcsine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} g^{k+1} \cdot \arcsin^{(k+1)}(\theta \cdot g)
 * \f]
 *
 * Here, \f$ \arcsin^{(k+1)}(\theta \cdot g) \f$ is the higher order derivative for an arcsine: see
 * Makino (1998).
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$asin(tm)\f$
 */
inline audi::taylor_model asin(const audi::taylor_model &tm)
{
    if (!boost::numeric::subset(tm.get_bounds() + tm.get_rem(), int_d(-1, 1))) {
        throw std::runtime_error("The range of the function is not a subset of [-1, 1], and "
                                 "therefore is ill-defined when performing an asin operation.");
    }

    double const_term = tm.get_tpol().constant_cf();
    audi::taylor_model tm_adapted
        = tm * std::sqrt(1 - std::pow(const_term, 2)) - const_term * audi::sqrt(1 - audi::pow(tm, 2));
    int_d tm_adapted_bounds
        = audi::taylor_model::get_bounds(tm_adapted.get_tpol(), tm_adapted.get_exp(), tm_adapted.get_dom());
    int_d tm_adapted_remainder = tm_adapted.get_rem();
    int_d tm_adapted_remainder_bounds = tm_adapted_bounds + tm_adapted_remainder;
    int k = tm.get_tpol().get_order();

    int_d asin_deriv = asin_derivative(int_d(0, 1) * tm_adapted_remainder_bounds, k + 1);
    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1))
                            * boost::numeric::pow(tm_adapted_remainder_bounds, k + 1) * asin_deriv;

    return audi::taylor_model(audi::asin(tm.get_tpol()), total_rem_bound, tm.get_exp(),
                              tm.get_dom());
}

/// Overload for the acos of a taylor_model.
/**
 * Function that calculates the acos of a taylor_model.
 *
 * Here, we use the following identity utilizing the arcsine. So see asin(const audi::taylor_model &tm) for the
 * implementation.
 *
 * \f$ \arccos(f) = \frac{\pi}{2} - \arcsin(f) \f$
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$acos(tm)\f$
 */
inline audi::taylor_model acos(const audi::taylor_model &tm)
{
    return std::numbers::pi / 2 - audi::asin(tm);
}

/// Overload for the atan of a taylor_model.
/**
 * Function that calculates the atan of a taylor_model.
 *
 * Here, some custom arithmetic is used from Makino (1998):
 *
 * \f$ \arctan(f) = \arctan(f_0) + \arctan(g) \f$
 *
 * where \f$ g = \frac{\hat{f}}{1 + f_0 \cdot f} \f$
 *
 * The remainder for the arcsine is:
 *
 * \f[
 * R = \frac{1}{(k+1)!} g^{k+1} \cdot \arccos^{(k+1)}(\arctan(\theta \cdot g)) \cdot \sin((k + 1) \cdot
 * (\arctan(\theta \cdot g) + \frac{\pi}{2}))
 * \f]
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$atan(tm)\f$
 */
inline audi::taylor_model atan(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::taylor_model tm_bar = tm - const_term;
    audi::taylor_model tm_adapted = tm_bar / (1 + const_term * tm);
    int_d tm_adapted_bounds
        = audi::taylor_model::get_bounds(tm_adapted.get_tpol(), tm_adapted.get_exp(), tm_adapted.get_dom());
    int_d tm_adapted_remainder = tm_adapted.get_rem();
    int_d tm_adapted_remainder_bounds = tm_adapted_bounds + tm_adapted_remainder;
    int k = tm.get_tpol().get_order();

    int_d total_rem_bound
        = int_d(1.0 / (k + 1.0)) * boost::numeric::pow(tm_adapted_remainder_bounds, k + 1)
          * boost::numeric::pow(boost::numeric::cos(boost::numeric::atan(int_d(0, 1) * tm_adapted_remainder_bounds)),
                                k + 1)
          * boost::numeric::sin(
              int_d(k + 1)
              * (boost::numeric::atan(int_d(0, 1) * tm_adapted_remainder_bounds) + std::numbers::pi / 2.0));

    return audi::taylor_model(audi::atan(tm.get_tpol()), total_rem_bound, tm.get_exp(),
                              tm.get_dom());
}

/// Overload for the asinh of a taylor_model.
/**
 * Function that calculates the asinh of a taylor_model.
 *
 * Here, the following identity is used:
 *
 * \f$ arcsinh(f) = \log(f + \sqrt(f^2 + 1)) \f$
 *
 * For the remainder, see log(const audi::taylor_model &tm), sqrt(const audi::taylor_model &tm), and pow(const
 * audi::taylor_model &tm).
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$arcsinh(tm)\f$
 */
inline audi::taylor_model asinh(const audi::taylor_model &tm)
{
    return audi::log(tm + audi::sqrt(audi::pow(tm, 2) + 1));
}

/// Overload for the acosh of a taylor_model.
/**
 * Function that calculates the acosh of a taylor_model.
 *
 * Here, the following identity is used:
 *
 * \f$ arccosh(f) = \log(f + \sqrt(f^2 - 1)) \f$
 *
 * where \f$ B(f) + I_f \subset (1, \infty) \f$ to prevent complex values from the sqrt.
 *
 * For the remainder, see log(const audi::taylor_model &tm), sqrt(const audi::taylor_model &tm), and pow(const
 * audi::taylor_model &tm).
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$arccosh(tm)\f$
 */
inline audi::taylor_model acosh(const audi::taylor_model &tm)
{
    audi::taylor_model tm_square = audi::pow(tm, 2);
    if (!boost::numeric::subset(tm_square.get_bounds() + tm_square.get_rem(),
                                int_d(1, std::numeric_limits<double>::infinity()))) {
        throw std::runtime_error("The range of the function is not a subset of the [1, infty], and "
                                 "therefore is ill-defined when performing an acosh operation.");
    }

    return audi::log(tm + audi::sqrt(tm_square - 1));
}

/// Overload for the atanh of a taylor_model.
/**
 * Function that calculates the atanh of a taylor_model.
 *
 * Here, the following identity is used:
 *
 * \f$ arctanh(f) = \frac{1}{2} \cdot \log(\frac{1 + f}{1 - f}) \f$
 *
 * For the remainder, see log(const audi::taylor_model &tm).
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$arctanh(tm)\f$
 */
inline audi::taylor_model atanh(const audi::taylor_model &tm)
{
    return (1.0 / 2.0) * audi::log((1 + tm) / (1 - tm));
}

/// Overload for the abs of a taylor_model.
/**
 * UNVERIFIED. Function that calculates the absolute value of a taylor_model.
 *
 * For the remainder, it is assumed that one just takes the absolute value of the remainder bound.
 *
 * @param tm the taylor_model base
 *
 * @return a new taylor_model representing \f$abs(tm)\f$
 */
inline audi::taylor_model abs(const audi::taylor_model &tm)
{
    return audi::taylor_model(audi::abs(tm.get_tpol()), boost::numeric::abs(tm.get_rem()), tm.get_exp(), tm.get_dom());
}

} // end of namespace audi

#endif
