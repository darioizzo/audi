#ifndef AUDI_REAL128_HPP
#define AUDI_REAL128_HPP

#include <audi/config.hpp>

#if defined(AUDI_WITH_MPPP)

#include <algorithm>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <exception>
#include <mp++/real128.hpp>
#include <piranha/math.hpp>
#include <piranha/pow.hpp>
#include <string>

#include <vector>


// This header adds to the mppp::real128 class  the necessary methods that allow it to be
// considered as a type in gdual (piranha::is_Cf, piranha::is_differentiable)

namespace boost
{
namespace serialization
{
template <class Archive>
void save(Archive &ar, const mppp::real128 &t, unsigned int version)
{
    std::string s(t.to_string());
    ar << s;
}
template <class Archive>
void load(Archive &ar, mppp::real128 &t, unsigned int version)
{
    std::string s;
    ar >> s;
    t = mppp::real128(s);
}
} // namespace serialization
} // namespace boost
BOOST_SERIALIZATION_SPLIT_FREE(mppp::real128)

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

template <typename U>
struct pow_impl<mppp::real128, U> {
    mppp::real128 operator()(const mppp::real128 &c, const U &exp) const
    {
        return mppp::pow(c, exp);
    };
};
} // namespace math

template <typename Archive>
struct boost_save_impl<Archive, mppp::real128> : boost_save_via_boost_api<Archive, mppp::real128> {
};

template <typename Archive>
struct boost_load_impl<Archive, mppp::real128> : boost_load_via_boost_api<Archive, mppp::real128> {
};

template <>
struct zero_is_absorbing<mppp::real128> : std::false_type {
};
 
} // namespace piranha



#else

#error The real128.hpp header was included but audi was not configured with the AUDI_WITH_MPPP option.

#endif

#endif
