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

#include <algorithm>
#include <boost/timer/timer.hpp>
#include <functional>
#include <iomanip>
#include <random>
#include <vector>

#include <audi/audi.hpp>

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

void measure_speedup(int points, unsigned int order)
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
    auto foo = coeff_v + x + y + z - w + q + r + s;
    for (int i = 1; i < 20; ++i) {
        foo *= coeff_v + x + y + z - w + q + r + s;
    }
    foo *= foo;
    auto wall1 = t1.elapsed().wall;

    boost::timer::cpu_timer t2;
    for (auto i = 0u; i < points; ++i) {
        audi::gdual_d x{xcoeff[i], "x", order}, y{ycoeff[i], "y", order}, z{zcoeff[i], "z", order},
            w{zcoeff[i], "w", order}, q{zcoeff[i], "q", order}, r{zcoeff[i], "r", order}, s{zcoeff[i], "s", order};
        auto foo = coeff[i] + x + y + z - w + q + r + s;
        for (int j = 1; j < 20; ++j) {
            foo *= coeff[i] + x + y + z - w + q + r + s;
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
        for (auto i = 1u; i < 4; ++i) {
            measure_speedup(pow(10, i), j);
        }
    }
}
