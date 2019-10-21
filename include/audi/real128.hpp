#ifndef AUDI_REAL128_HPP
#define AUDI_REAL128_HPP

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)

#include <algorithm>
#include <exception>
#include <string>
#include <vector>

#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/serialization.hpp>

#include <mp++/real128.hpp>

// This header adds to the mppp::real128 class  the necessary methods that allow it to be
// considered as a type in gdual (obake::is_cf, obake::is_differentiable)

namespace boost::serialization
{

template <class Archive>
void serialize(Archive &ar, mppp::real128 &t, unsigned int)
{
    ar &serialization::make_binary_object(&t, sizeof(t));
}

} // namespace boost::serialization

#else

#error The real128.hpp header was included but audi was not configured with the AUDI_WITH_QUADMATH option.

#endif

#endif
