// Copyright © 2018–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com),
// Sean Cowan (lambertarc@icloud.com)
//
// This file is part of the audi library.
//
// The audi library is free software: you can redistribute it and/or modify
// it under the terms of either:
//   - the GNU General Public License as published by the Free Software
//     Foundation, either version 3 of the License, or (at your option)
//     any later version, or
//   - the GNU Lesser General Public License as published by the Free
//     Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
// The audi library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License and the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser General Public License along with the audi library.
// If not, see <https://www.gnu.org/licenses/>.

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
