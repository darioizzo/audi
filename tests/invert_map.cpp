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

#define BOOST_TEST_MODULE audi_invert_map_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include <obake/polynomials/packed_monomial.hpp>

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <audi/invert_map.hpp>

#include "helpers.hpp"

using namespace audi;

BOOST_AUTO_TEST_CASE(computations)
{
    // We need access to operators & and other hidden stuff :)
    using namespace detail;
    // 1 - D cases
    {
        // We invert (around 0) asin(x), hence checking to get the MacLaurin series
        // for sin(x) = dx - dx^3/3! + dx^5/5! + H.O.T.
        gdual<double, obake::packed_monomial<std::uint64_t>> x(0, "x", 5);
        std::vector<gdual<double, obake::packed_monomial<std::uint64_t>>> map, map_inv;
        map.push_back(asin(x));
        map_inv = invert_map(map);
        BOOST_CHECK_CLOSE(map_inv[0].constant_cf(), 0, 1e-12);        // constant
        BOOST_CHECK_CLOSE(map_inv[0].find_cf({1}), 1, 1e-12);         // x
        BOOST_CHECK_CLOSE(map_inv[0].find_cf({2}), 0., 1e-12);        // x^2
        BOOST_CHECK_CLOSE(map_inv[0].find_cf({3}), -1. / 6., 1e-12);  // x^3
        BOOST_CHECK_CLOSE(map_inv[0].find_cf({4}), 0., 1e-12);        // x^4
        BOOST_CHECK_CLOSE(map_inv[0].find_cf({5}), 1. / 120., 1e-12); // x^5
    }
    {
        // We invert (around 0.1) some complex mathematical formula and check that its inversion
        // computes correctly
        gdual<double, obake::packed_monomial<std::uint64_t>> x(0.1, "x", 10);
        std::vector<gdual<double, obake::packed_monomial<std::uint64_t>>> map, map_inv;
        auto p0 = x * 3. * sin(x / 2.) + x + 1. * pow(1 + x, cos(x * 3));
        map.push_back(p0);
        map_inv = invert_map(map);
        auto dp0 = 0.01;
        auto dx = map_inv[0].evaluate({{"dp0", dp0}});
        auto p0_in_dx = p0.evaluate({{"dx", dx}});
        BOOST_CHECK_CLOSE(p0_in_dx, p0.constant_cf() + dp0, 1e-10);
    }
    // 2 - D cases
    {
        // We invert a 2-D map in x,y (polynomial) and check that the map composition
        // is the identity f(g(x)) = x (since no constant_cf is present)
        gdual<double, obake::packed_monomial<std::uint64_t>> x(1., "x", 4);
        gdual<double, obake::packed_monomial<std::uint64_t>> y(1., "y", 4);
        std::vector<gdual<double, obake::packed_monomial<std::uint64_t>>> map, map_inv;
        map.push_back(x + y + x * x * y - 3.);
        map.push_back(x - y + y * y * x - 1.);
        map_inv = invert_map(map);
        auto I = map_inv & map;
        BOOST_CHECK(EPSILON_COMPARE(I[0], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "x", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[1], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "y", 4), 1e-10));
    }
    {
        // We invert a 2-D map in x,y and check that the map composition
        // is the identity f(g(x)) = x (since no constant_cf is present)
        gdual<double, obake::packed_monomial<std::uint64_t>> x(1., "x", 4);
        gdual<double, obake::packed_monomial<std::uint64_t>> y(1., "y", 4);
        std::vector<gdual<double, obake::packed_monomial<std::uint64_t>>> map, map_inv;
        map.push_back(sin(x * y - 1.) + x + y - 2.);
        map.push_back(cos(x * y * y - 1) - y + x - 1.);
        map_inv = invert_map(map);
        auto I = map_inv & map;
        BOOST_CHECK(EPSILON_COMPARE(I[0], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "x", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[1], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "y", 4), 1e-10));
    }
    // 3 - D cases
    {
        // We invert a 3-D map in x,y,z and check that the map composition
        // is the identity f(g(x)) = x (since no constant_cf is present)
        gdual<double, obake::packed_monomial<std::uint64_t>> x(1., "x", 4);
        gdual<double, obake::packed_monomial<std::uint64_t>> y(1., "y", 4);
        gdual<double, obake::packed_monomial<std::uint64_t>> z(1., "z", 4);
        std::vector<gdual<double, obake::packed_monomial<std::uint64_t>>> map, map_inv;
        map.push_back(x * y * z + x * y - 2 * z);
        map.push_back(y * z * x * x + z * y - z - y);
        map.push_back(x - 2 * z + cos(x * y * z - 1.));
        map_inv = invert_map(map);
        auto I = map_inv & map;
        BOOST_CHECK(EPSILON_COMPARE(I[0], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "x", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[1], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "y", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[2], gdual<double, obake::packed_monomial<std::uint64_t>>(0, "z", 4), 1e-10));
    }
}
