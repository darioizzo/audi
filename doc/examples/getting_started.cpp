// Copyright © 2018–2026 Dario Izzo (dario.izzo@gmail.com),
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

#include <audi/gdual.hpp>
#include <audi/functions.hpp>
#include <iostream>

using namespace audi;

int main()
{
    // We want to compute the Taylor expansion of a function f (and thus all derivatives) at x=2, y=3
    using gdual = gdual<double>;
    // 1 - Define the generalized dual numbers (over doubles, 7 is the truncation order, i.e. the maximum
    // order of derivation we will need)
    gdual x(2, "x", 7);
    gdual y(3, "y", 7);

    // 2 - Compute your function as usual
    auto f = exp(x * x + cbrt(y) / log(x * y));

    // 3 - Inspect the results (this has a constant complexity now as all computations have been made already)
    std::cout << "Taylor polynomial: " << f
              << std::endl; // This is the Taylor expansion of f (truncated at the 7th order)
    std::cout << "Derivative value: " << f.get_derivative({1, 0})
              << std::endl; // This is the value of the derivative (d / dx)
    std::cout << "Derivative value: " << f.get_derivative({4, 3})
              << std::endl; // This is the value of the mixed derivative (d^7 / dx^4dy^3)

    // 4 - Using the dictionary interface (note the presence of the "d" before all variables)
    std::cout << "Derivative value: " << f.get_derivative({{"dx", 1}})
              << std::endl; // This is the value of the derivative (d / dx)
    std::cout << "Derivative value: " << f.get_derivative({{"dx", 4}, {"dy", 3}})
              << std::endl; // This is the value of the mixed derivative (d^7 / dx^4dy^3)
}
