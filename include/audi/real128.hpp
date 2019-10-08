#ifndef AUDI_REAL128_HPP
#define AUDI_REAL128_HPP

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)

#include <algorithm>
#include <exception>
#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/split_free.hpp>

#include <mp++/real128.hpp>

// This header adds to the mppp::real128 class  the necessary methods that allow it to be
// considered as a type in gdual (obake::is_cf, obake::is_differentiable)

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

#else

#error The real128.hpp header was included but audi was not configured with the AUDI_WITH_QUADMATH option.

#endif

#endif
