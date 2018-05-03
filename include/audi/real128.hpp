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
} // namespace math
} // namespace piranha



#endif
