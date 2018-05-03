#ifndef AUDI_REAL128_HPP
#define AUDI_REAL128_HPP

#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <exception>
#include <vector>

#include <mp++/real128.hpp>

// This header adds to the mppp::real128 class  the necessary methods that allow it to be
// considered as a type in gdual (piranha::is_Cf, piranha::is_differentiable)

namespace piranha
{
namespace math
{
template <>
struct partial_impl<mppp::real128> {
    /// Call operator.
    /**
     * @return an instance of piranha::real128 constructed from zero.
     */
    mppp::real128 operator()(const mppp::real128 &, const std::string &) const
    {
        return mppp::real128{};
    }
};
}} // end of piranha namespace

// Specialization for numeric limits. Implementing the becessary ones for bernoulli numbers to work in boost
namespace std {
template<> 
class numeric_limits<mppp::real128> {
public:
constexpr static int radix = 2u;
constexpr static int digits = 113u;
constexpr static int max_exponent = 16384u;
constexpr static bool is_specialized = true;
constexpr static bool has_infinity = true;
constexpr static bool has_quiet_NaN = true;
constexpr static bool has_signaling_NaN = false;
static mppp::real128 min() noexcept { return mppp::real128("3.36210314311209350626267781732175260e-4932"); }
static mppp::real128 max() noexcept { return mppp::real128("1.18973149535723176508575932662800702e4932"); }
static mppp::real128 lowest() {return mppp::real128_min();};
static mppp::real128 quiet_NaN() {return mppp::real128(0.)/mppp::real128(0.);}; //todo: use nanq
static mppp::real128 infinity() {return mppp::real128(1.)/mppp::real128(0.);}; //todo: can this be made better?
};
}

#endif
