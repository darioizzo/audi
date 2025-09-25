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

#include <DACE/DA.h>
#include <algorithm>
#include <boost/timer/timer.hpp>
#include <functional>
#include <iomanip>
#include <random>
#include <vector>

#include <audi/audi.hpp>

// compile with g++ multiplication_perf_DACE.cpp -std=c++14 -I/usr/local/include -ldace -lboost_system
// -lboost_unit_test_framework -lboost_timer -lpthread -lgmp -lmpfr -O3 -DNDEBUG
//
using namespace DACE;

std::vector<double> random_vector_double(unsigned int N, double lb, double ub)
{
    std::vector<double> retval(N);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(lb, ub);
    for (auto &x : retval) {
        x = dist(mt);
    }
    return retval;
}

void measure_speedup(int points, unsigned int order, unsigned int size)
{
    auto lb = -0.5;
    auto ub = 0.5;
    // We generate random vector coefficients
    auto xcoeff = random_vector_double(points, lb, ub);
    auto ycoeff = random_vector_double(points, lb, ub);
    auto zcoeff = random_vector_double(points, lb, ub);
    auto coeff = random_vector_double(points, lb, ub);
    auto coeff_v = audi::gdual_v(coeff);

    boost::timer::cpu_timer t1;
    audi::gdual_v x{xcoeff, "x", order}, y{ycoeff, "y", order}, z{zcoeff, "z", order}, w{zcoeff, "w", order},
        q{zcoeff, "q", order}, r{zcoeff, "r", order}, s{zcoeff, "s", order};
    auto foo = (coeff_v + x + y + z + w + q + r + s) / (coeff_v - x - y - z - w - q - r - s);
    for (int i = 1; i < size; ++i) {
        foo *= (coeff_v + x + y + z + w + q + r + s) / (coeff_v - x - y - z - w - q - r - s);
    }
    foo *= foo;
    auto wall1 = t1.elapsed().wall;

    DA::init(order, 7); // 7 variables xyzwqrs
    boost::timer::cpu_timer t2;
    for (auto i = 0u; i < points; ++i) {
        DA x{1, xcoeff[i]}, y{2, ycoeff[i]}, z{3, zcoeff[i]}, w{4, zcoeff[i]}, q{5, zcoeff[i]}, r{6, zcoeff[i]},
            s{7, zcoeff[i]};
        auto foo = (coeff[i] + x + y + z + w + q + r + s) / (coeff[i] - x - y - z - w - q - r - s);
        for (int j = 1; j < size; ++j) {
            foo *= (coeff[i] + x + y + z + w + q + r + s) / (coeff[i] - x - y - z - w - q - r - s);
        }
        foo *= foo;
    }
    auto wall2 = t2.elapsed().wall;
    std::cout << "Number of points in vector: " << std::setw(10) << points << "\t\tOrder: " << std::setw(3) << order
              << "\tSpeedup: " << std::setw(4) << (1. * wall2) / wall1 << "\n";
}

int main()
{
    for (auto j = 1u; j < 6; ++j) {
        for (auto i = 2u; i < 8; ++i) {
            measure_speedup(pow(2, 2 * i), j, 13);
        }
    }
    for (auto j = 6u; j < 10; ++j) {
        for (auto i = 2u; i < 8; ++i) {
            measure_speedup(pow(2, 2 * i), j, 5);
        }
    }
}
