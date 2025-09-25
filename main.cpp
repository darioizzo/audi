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

#include <audi/audi.hpp>

using namespace audi;

int main()
{
    gdual_d x(0.1, "x", 2);
    gdual_d y(0.1, "y", 2);
    gdual_d z(0.1, "z", 2);
    // Order extraction
    auto f = audi::sin(x * y);
    auto f2 = audi::sin(z);
    std::cout << "Full map: " << f << std::endl;
    auto g = f.extract_terms(2);
    std::cout << "Only order 2: " << g << std::endl;
    // Map concatenation
    std::cout << "\n\nFull map (x,y): " << f << std::endl;
    std::cout << "Full map (z): " << f2 << std::endl;
    auto new_f = f.subs("dy", f2);
    std::cout << "Full map (x,z): " << new_f << std::endl;
    std::cout << new_f.info() << std::endl;
}
