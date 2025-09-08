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

// pow method
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

// pow method 1/2 and -1/2
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
                                     "therefore is ill-defined when performing a log operation.");
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
                                     "therefore is ill-defined when performing a log operation.");
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

inline audi::taylor_model sqrt(const audi::taylor_model &tm)
{
    return audi::pow(tm, 1.0 / 2.0);
}

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

inline audi::taylor_model tan(const audi::taylor_model &tm)
{
    return audi::sin(tm) / audi::cos(tm);
}

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

inline audi::taylor_model tanh(const audi::taylor_model &tm)
{
    return audi::sinh(tm) / audi::cosh(tm);
}

// Compute higher-order derivatives of arcsin(a) using the recurrence relation
// Only used for asin() function below
int_d asin_derivative(int_d a, int order)
{
    std::vector<int_d> derivs(order + 1);

    // Base cases
    derivs[1] = int_d(1.0) / boost::numeric::sqrt(int_d(1.0) - a * a); // arcsin'(a)
    if (order >= 2)
        derivs[2] = a * boost::numeric::pow(int_d(1.0) - a * a, 2)
                    / boost::numeric::pow(int_d(1.0) - a * a, 3); // arcsin''(a)

    // Recursive computation for higher orders
    for (int k = 1; k <= order - 2; ++k) {
        derivs[k + 2] = (int_d(2 * k + 1) * a * derivs[k + 1] + int_d(k * k) * derivs[k]) / int_d(1.0) - a * a;
    }

    return derivs[order];
}

inline audi::taylor_model asin(const audi::taylor_model &tm)
{
    if (!boost::numeric::subset(tm.get_bounds() + tm.get_rem(), int_d(-1, 1))) {
        throw std::runtime_error("The range of the function is not a subset of the positive Real numbers, and "
                                 "therefore is ill-defined when performing a log operation.");
    }
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d asin_deriv = asin_derivative(int_d(0, 1) * f_bar_remainder_bounds, k + 1);
    int_d total_rem_bound = int_d(1.0 / boost::math::factorial<double>(k + 1))
                            * boost::numeric::pow(f_bar_remainder_bounds, k + 1) * asin_deriv;

    return audi::taylor_model(audi::asin(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

inline audi::taylor_model acos(const audi::taylor_model &tm)
{
    return std::numbers::pi / 2 - audi::asin(tm);
}

inline audi::taylor_model atan(const audi::taylor_model &tm)
{
    double const_term = tm.get_tpol().constant_cf();
    audi::gdual<double> f_bar = tm.get_tpol() - const_term;
    int_d f_bar_bounds = audi::taylor_model::get_bounds(f_bar, tm.get_exp(), tm.get_dom());
    int_d f_bar_remainder = tm.get_rem();
    int_d f_bar_remainder_bounds = f_bar_bounds + f_bar_remainder;
    int k = tm.get_tpol().get_order();

    int_d total_rem_bound
        = int_d(1.0 / (k + 1.0)) * boost::numeric::pow(f_bar_remainder_bounds, k + 1)
          * boost::numeric::pow(boost::numeric::cos(boost::numeric::atan(int_d(0, 1) * f_bar_remainder_bounds)), k + 1)
          * boost::numeric::sin(int_d(k + 1)
                                * (boost::numeric::atan(int_d(0, 1) * f_bar_remainder_bounds) + std::numbers::pi / 2));

    return audi::taylor_model(audi::atan(tm.get_tpol()), total_rem_bound, tm.get_exp(), tm.get_dom());
}

inline audi::taylor_model asinh(const audi::taylor_model &tm)
{
    return audi::log(tm + audi::sqrt(audi::pow(tm, 2) + 1));
}

inline audi::taylor_model acosh(const audi::taylor_model &tm)
{
    // TODO: Add check on validity of acosh (1 <= x <= \infty)
    return audi::log(tm + audi::sqrt(audi::pow(tm, 2) - 1));
}

inline audi::taylor_model atanh(const audi::taylor_model &tm)
{
    return (1.0 / 2.0) * audi::log((1 + tm) / (1 - tm));
}

// This is not given in Makino (1998), but follows directly from the definition. However, this
// function is NOT derived mathematically.
inline audi::taylor_model abs(const audi::taylor_model &tm)
{
    return audi::taylor_model(audi::abs(tm.get_tpol()), tm.get_rem(), tm.get_exp(), tm.get_dom());
}

} // end of namespace audi

#endif
