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

#include <audi/functions.hpp>
#include <audi/gdual.hpp>

using namespace audi;

void scalable_test_sin_cos(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> p1(1);
    for (int i = 0u; i < n; ++i) {
        p1 += variables[i];
    } // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    auto res = sin_and_cos(p1);
}

void scalable_test_sin_and_cos(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> p1(1);
    for (int i = 0u; i < n; ++i) {
        p1 += variables[i];
    } // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    auto res = sin_and_cos(p1);
}

void scalable_test_tan(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> p1(1);
    gdual<double> tangent(p1);
    gdual<double> sine(p1);
    gdual<double> cosine(p1);
    for (int i = 0u; i < n; ++i) {
        p1 += variables[i];
    } // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    tangent = tan(p1);
}

void scalable_test_sin_over_cos(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> p1(1.);
    gdual<double> tangent(p1);
    for (int i = 0u; i < n; ++i) {
        p1 += variables[i];
    } // 1 + x1 + x2 + ...

    auto res = sin_and_cos(p1);
    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    tangent = res[0] / res[1];
}

BOOST_AUTO_TEST_CASE(trigonometry_perf)
{
    unsigned int low = 9, high = 10;

    // sin and cos
    std::cout << "Computing sin(1 + x1 + x2 + ...) and cos(1 + x1 + x2 + ...) separately: " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_cos(m, n);
        }
    }

    // sin and cos together
    std::cout << "\nComputing sin(1 + x1 + x2 + ...) and cos(1 + x1 + x2 + ...) at once: " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_and_cos(m, n);
        }
    }

    // tan
    std::cout << "\nTesting performance of tan(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_tan(m, n);
        }
    }

    // sin / cos
    std::cout << "\nTesting division of sin(1 + x1 + x2 + ...) / cos(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_over_cos(m, n);
        }
    }
}
