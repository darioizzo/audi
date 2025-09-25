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

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <vector>

#include <boost/optional.hpp>

#include <audi/gdual.hpp>

using namespace audi;

void scalable_div(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> denom(1.);
    gdual<double> num(1.);
    for (int i = 0u; i < n; ++i) {
        num += variables[i] * variables[i];
    }
    for (int i = 0u; i < n; ++i) {
        denom += variables[0] * variables[i];
    }
    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    auto result = num / denom;
}

BOOST_AUTO_TEST_CASE(division_performance)
{
    std::cout << "Testing division of (x1^2 + ... + xn^2) / (x1 * (x1 + x2 + .. + xn) ): " << std::endl;
    for (auto m = 5; m < 10; ++m) {
        for (auto n = 5; n < 12; ++n) {
            scalable_div(m, n);
        }
    }
}
