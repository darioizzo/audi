#define BOOST_TEST_MODULE audi_functions_floating_precision_test
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>
#include <audi/gdual.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(division_precision)
{
    {
        gdual<double> x(0.1, "x", 4);
        gdual<double> y(0.13, "y", 4);

        auto p1 = (x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x);
        auto p2 = gdual<double>(1.) / p1;
        BOOST_CHECK(EPSILON_COMPARE(p1 * p2, gdual<double>(1), 1e-10) == true);
    }
}
