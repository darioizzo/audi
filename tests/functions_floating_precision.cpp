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

#define BOOST_TEST_MODULE audi_functions_floating_precision_test
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>
#include <audi/gdual.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(division_precision)
{
    {
        gdual<double> x(0.1, "x", 4);
        gdual<double> y(0.13, "y", 4);

        auto p1 = (x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x);
        auto p2 = gdual<double>(1.) / p1;
        BOOST_CHECK(EPSILON_COMPARE(p1 * p2, gdual<double>(1), 1e-10) == true);
    }
}
