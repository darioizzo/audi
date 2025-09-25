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

#define BOOST_TEST_MODULE audi_gdualv_test

#include <boost/test/unit_test.hpp>

#include <fstream>
#include <stdexcept>
#include <vector>

#include <audi/config.hpp>

#if defined(AUDI_WITH_QUADMATH)
#include <audi/real128.hpp>
#endif

#include <audi/functions.hpp>
#include <audi/gdual.hpp>
#include <audi/vectorized.hpp>

using namespace audi;

using gdual_v = audi::gdual<audi::vectorized<double>>;

template <typename T>
T f(const T &x, const T &y)
{
    return (x * y) / (y - x);
}

template <typename T>
T dfdxy(const T &x, const T &y)
{
    return -2 * x * y / (y - x) / (y - x) / (y - x);
}

BOOST_AUTO_TEST_CASE(bug_issue_57)
{
    double a = 0.;
    double b = 1.;

    gdual_v x(std::vector<double>{1, a}, "x", 2);
    gdual_v y(std::vector<double>{2, b}, "y", 2);
    // BOOST_CHECK_EQUAL(f(x, y).find_cf({1, 1})[1], dfdxy(a, b));
    std::cout << "\nWRONG: " << std::endl;
    auto wrong = x * y / (y - x);
    std::cout << "Wrong: " << wrong << std::endl;
    std::cout << "\nRIGHT: " << std::endl;
    auto right = (1. / (y - x)) * (x * y) ;
    std::cout << "Right: " << right << std::endl;

}
