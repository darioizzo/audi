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

#include <boost/test/tools/old/interface.hpp>
#define BOOST_TEST_MODULE taylor_models_test

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>
#include <audi/io.hpp>
#include <audi/taylor_model_utilities.hpp>

BOOST_AUTO_TEST_CASE(test_matmul)
{
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}};
    std::vector<std::vector<double>> B = {{5, 6}, {7, 8}};

    auto C = matmul(A, B);

    BOOST_CHECK_EQUAL(C[0][0], 19.0);
    BOOST_CHECK_EQUAL(C[0][1], 22.0);
    BOOST_CHECK_EQUAL(C[1][0], 43.0);
    BOOST_CHECK_EQUAL(C[1][1], 50.0);
}

BOOST_AUTO_TEST_CASE(test_transpose)
{
    std::vector<std::vector<double>> M = {{1, 2, 3}, {4, 5, 6}};
    auto T = transpose(M);

    BOOST_REQUIRE_EQUAL(T.size(), 3);
    BOOST_REQUIRE_EQUAL(T[0].size(), 2);

    BOOST_CHECK_EQUAL(T[0][0], 1.0);
    BOOST_CHECK_EQUAL(T[1][0], 2.0);
    BOOST_CHECK_EQUAL(T[2][1], 6.0);
}


