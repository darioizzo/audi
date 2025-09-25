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

#ifndef AUDI_HPP
#define AUDI_HPP

#include <complex>

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)
#include <audi/real128.hpp>
#endif

#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>
#include <audi/gdual.hpp>
#include <audi/invert_map.hpp>
#include <audi/vectorized.hpp>
#include <audi/taylor_model_bounding.hpp>
#include <audi/taylor_model_utilities.hpp>
#include <audi/taylor_model.hpp>
#include <audi/taylor_model_functions.hpp>

namespace audi
{
using gdual_d = audi::gdual<double>;
using gdual_v = audi::gdual<audi::vectorized<double>>;
using gdual_c = audi::gdual<std::complex<double>>;
#if defined(AUDI_WITH_QUADMATH)
using gdual_mp = audi::gdual<mppp::real128>;
#endif
} // namespace audi

#endif
