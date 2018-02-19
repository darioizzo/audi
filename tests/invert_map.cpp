#define BOOST_TEST_MODULE audi_invert_map_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <audi/invert_map.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(computations)
{
    gdual_d x(0.3, "x", 4);
    gdual_d y(0.3, "y", 4);
    std::vector<gdual_d> map, map_inv;
    map.push_back(sin(x * y) * x * y + x - y);
    map.push_back(x - y + cos(x * y));
    map_inv = invert_map(map);
    map[0] -= map[0].constant_cf();
    map[1] -= map[1].constant_cf();
    using namespace detail;
    auto res = trim(map & map_inv, 1e-8);
    print(res);
}
