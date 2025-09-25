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
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/vectorized_double.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Default Constructor
    {
        vectorized_double x;
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 0);
    }
    // Constructor from int.
    {
        vectorized_double x{123};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123);
    }
    // Constructor from double.
    {
        vectorized_double x{123.};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123.);
    }
    // Constructor from an std::vector
    {
        std::vector<double> in{1, 2, 3, 4, 2.2};
        vectorized_double x(in);
        BOOST_CHECK(x.size() == 5);
        BOOST_CHECK(std::equal(in.begin(), in.end(), x.begin()));
        std::vector<double> empty;
        // rvalue
        BOOST_CHECK_THROW(vectorized_double(std::vector<double>{}), std::invalid_argument);
        // lvalue
        BOOST_CHECK_THROW(vectorized_double{empty}, std::invalid_argument);
        // initializer list
        BOOST_CHECK_THROW(vectorized_double({}), std::invalid_argument);
    }
}
BOOST_AUTO_TEST_CASE(math)
{
    vectorized_double x1{1., 1., 2., -3., 4.};
    vectorized_double x2{1., 1., 2., 3., 4.};
    vectorized_double x3{-100., -100., -100., -100., -100.};
    vectorized_double x4{100., 100., 100., 100., 100.};
    BOOST_CHECK(abs(x1) == x2);
    BOOST_CHECK(x3 < x1);
    BOOST_CHECK(x1 > x3);
    BOOST_CHECK(x2 < x4);
    BOOST_CHECK(x4 > x2);
    BOOST_CHECK(x1 < 100.);
    BOOST_CHECK(-100 < x1);
    BOOST_CHECK(100. > x1);
    BOOST_CHECK(x1 > -100);
    BOOST_CHECK(x1 != x2);
    BOOST_CHECK(x1 == x1);
    BOOST_CHECK(x3 == -100);
    BOOST_CHECK(x4 == 100);
    BOOST_CHECK(x3 != 32);
    BOOST_CHECK(x4 != 32);
    BOOST_CHECK(abs(x1) == x2);
    BOOST_CHECK(abs(x2) != x1);
    BOOST_CHECK(x1 - 1. == (vectorized_double{0., 0., 1., -4., 3.}));
}

BOOST_AUTO_TEST_CASE(more_math)
{
    vectorized_double x1{1., 2.};
    vectorized_double x2{3., 3.};
    // Interoperability against arith types
    BOOST_CHECK(x1 * 2. == (vectorized_double{2., 4.}));
    BOOST_CHECK(x1 / 2. == (vectorized_double{0.5, 1.}));
    BOOST_CHECK(x1 + 2. == (vectorized_double{3., 4.}));
    BOOST_CHECK(x1 - 2. == (vectorized_double{-1., 0.}));
    BOOST_CHECK(x1 * 2 == (vectorized_double{2., 4.}));
    BOOST_CHECK(x1 / 2 == (vectorized_double{0.5, 1.}));
    BOOST_CHECK(x1 + 2 == (vectorized_double{3., 4.}));
    BOOST_CHECK(x1 - 2 == (vectorized_double{-1., 0.}));
    BOOST_CHECK(2 * x1 == (vectorized_double{2., 4.}));
    BOOST_CHECK(2 / x1 == (vectorized_double{2., 1.}));
    BOOST_CHECK(2 + x1 == (vectorized_double{3., 4.}));
    BOOST_CHECK(2 - x1 == (vectorized_double{1., 0.}));
    BOOST_CHECK(2 * x1 == (vectorized_double{2., 4.}));
    BOOST_CHECK(2 / x1 == (vectorized_double{2., 1.}));
    BOOST_CHECK(2 + x1 == (vectorized_double{3., 4.}));
    BOOST_CHECK(2 - x1 == (vectorized_double{1., 0.}));
    BOOST_CHECK(x1 != 1.);
    BOOST_CHECK(x2 == 3.);
    BOOST_CHECK(x1 > -1.);
    BOOST_CHECK(x1 < 3.);
    // operators with vectorized-vectorized
    BOOST_CHECK(x1 + x2 == (vectorized_double{4., 5.}));
    BOOST_CHECK(x1 - x2 == (vectorized_double{-2., -1.}));
    BOOST_CHECK(x1 / x2 == (vectorized_double{1./3., 2./3.}));
    BOOST_CHECK(x1 * x2 == (vectorized_double{3., 6.}));
}
