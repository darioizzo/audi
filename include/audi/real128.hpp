#ifndef AUDI_REAL128_HPP
#define AUDI_REAL128_HPP

#include <audi/config.hpp>

#if defined(AUDI_WITH_MPPP)

#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <exception>
#include <mp++/real128.hpp>
#include <piranha/math.hpp>
#include <vector>

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
} // namespace math
} // namespace piranha

#else

#error The real128.hpp header was included but audi was not configured with the AUDI_WITH_MPPP option.

#endif

#endif
