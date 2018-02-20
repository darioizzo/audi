#define BOOST_TEST_MODULE audi_invert_map_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <audi/invert_map.hpp>

#include "helpers.hpp"

using namespace audi;

BOOST_AUTO_TEST_CASE(stress_test)
{
    using namespace detail;    
    {
        // We invert a 3-D map in x,y,z and check that the map composition
        // is the identity f(g(x)) = x (since no constant_cf is present)
        gdual_d x(1., "x", 4);
        gdual_d y(1., "y", 4);
        gdual_d z(1., "z", 4);
        std::vector<gdual_d> map, map_inv;
        map.push_back(x * y * z + x * y - 2 * z);
        map.push_back(y * z * x * x + z * y - z - y);
        map.push_back(x - 2*z + cos(x * y * z - 1.));
        map_inv = invert_map(map, false);
        auto I = map_inv & map;
        BOOST_CHECK(EPSILON_COMPARE(I[0], gdual_d(0, "x", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[1], gdual_d(0, "y", 4), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[2], gdual_d(0, "z", 4), 1e-10));
    }
    {
        // We invert a 3-D map in x,y,z and check that the map composition
        // is the identity f(g(x)) = x (since no constant_cf is present)
        gdual_d x(1., "x", 7);
        gdual_d y(1., "y", 7);
        gdual_d z(1., "z", 7);
        std::vector<gdual_d> map, map_inv;
        map.push_back(x * y * z + x * y - 2 * z);
        map.push_back(y * z * x * x + z * y - z - y);
        map.push_back(x - 2*z + cos(x * y * z - 1.));
        map_inv = invert_map(map, true);
        auto I = map_inv & map;
        BOOST_CHECK(EPSILON_COMPARE(I[0], gdual_d(0, "x", 7), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[1], gdual_d(0, "y", 7), 1e-10));
        BOOST_CHECK(EPSILON_COMPARE(I[2], gdual_d(0, "z", 7), 1e-10));
    }
}
