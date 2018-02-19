#define BOOST_TEST_MODULE audi_invert_map_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include <audi/audi.hpp>
#include <audi/invert_map.hpp>
#include <audi/functions.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(computations)
{
    gdual_d x(0, "x", 3);
    gdual_d y(0, "y", 3);
    std::vector<gdual_d> map;
    map.push_back(x + y*y);
    map.push_back(y + x*x);
    invert_map(map);
}
