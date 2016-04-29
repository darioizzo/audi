#include "../src/gdual.hpp"
#include "helpers.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Constructing a constant with a symbol
    gdual x("x", 0);
    gdual y(3., "y", 0);
    gdual z(2., "z", 3);
    auto p = (x+y)*z*(x+y)*z*(x+y)*z*(x+y)*z*(x+y)*z;
    BOOST_CHECK_EQUAL(x.get_order(), 0);
    BOOST_CHECK_EQUAL(y.get_order(), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1,0,0}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({0,1,0}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({0,1,1}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1,1,1}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1,0,2}), 0);

    // TODO: other constructors
}

BOOST_AUTO_TEST_CASE(order_promotion)
{
    gdual c(1);
    gdual x("x", 4);
    gdual y("y", 2);
    BOOST_CHECK_EQUAL(c.get_order(), 0);
    BOOST_CHECK_EQUAL((x + y).get_order(), 4);
    BOOST_CHECK_EQUAL((x - y).get_order(), 4);
    BOOST_CHECK_EQUAL((x * y).get_order(), 4);
    BOOST_CHECK_EQUAL((x / y).get_order(), 4);
    BOOST_CHECK_EQUAL((x + c).get_order(), 4);
    BOOST_CHECK_EQUAL((x - c).get_order(), 4);
    BOOST_CHECK_EQUAL((x * c).get_order(), 4);
    BOOST_CHECK_EQUAL((x / c).get_order(), 4);
    BOOST_CHECK_EQUAL((y + x).get_order(), 4);
    BOOST_CHECK_EQUAL((y - x).get_order(), 4);
    BOOST_CHECK_EQUAL((y * x).get_order(), 4);
    BOOST_CHECK_EQUAL((y / x).get_order(), 4);
    BOOST_CHECK_EQUAL((c + x).get_order(), 4);
    BOOST_CHECK_EQUAL((c - x).get_order(), 4);
    BOOST_CHECK_EQUAL((c * x).get_order(), 4);
    BOOST_CHECK_EQUAL((c / x).get_order(), 4);
}
BOOST_AUTO_TEST_CASE(addition)
{
    gdual p1(1, 4);
    gdual p2("x", 4);
    BOOST_CHECK_EQUAL(p1 + p2, p2 + p1);
    BOOST_CHECK_EQUAL(1 + p2, p2 + 1);
    BOOST_CHECK_EQUAL(1. + p2, p2 + 1.);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
    gdual p1(1, 4);
    gdual p2("x", 4);
    BOOST_CHECK_EQUAL(p1 - p1, gdual(0,4));
    BOOST_CHECK_EQUAL(p2 - p2, gdual(0,4));
    BOOST_CHECK_EQUAL(p1 - p2, - (p2 - p1));
    BOOST_CHECK_EQUAL(1 + p1 - p1, gdual(1,4));
    BOOST_CHECK_EQUAL(1 + p2 - p2, gdual(1,4));
    BOOST_CHECK_EQUAL(1. + p1 - p1, gdual(1,4));
    BOOST_CHECK_EQUAL(1. + p2 - p2, gdual(1,4));
    BOOST_CHECK_EQUAL((1 - p1) + (p1 + p2), 1 + p2);
}

BOOST_AUTO_TEST_CASE(multiplication)
{
    // we test normalized derivatives computations of f = x + 3*x*y + y*y
    gdual x(3, "x",2);
    gdual y(7, "y",2);

    auto p1 = x + 3*x*y + y*y;
    BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 115);
    BOOST_CHECK_EQUAL(p1.find_cf({1,0}), 22);
    BOOST_CHECK_EQUAL(p1.find_cf({0,1}), 23);
    BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 3);
    BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 1);

    // we test the truncation order (5) on a 3 variables polynomial p2 = (1 + x0 + x1 + x2)^10
    std::vector<gdual> vars;
    for (auto i = 0u; i<3; ++i) {
        vars.emplace_back("x" + std::to_string(i), 5);
    }
    gdual p2(1, 5);
    for (auto var : vars) {
        p2+=var;
    }
    p2 = p2*p2*p2*p2*p2*p2;
    BOOST_CHECK_EQUAL(p2.get_order(), 5);
}

BOOST_AUTO_TEST_CASE(division)
{
    // we test normalized derivatives computations of f = 1 / (x + 2*x*y + y*y)
    {
        gdual x(0, "x",2);
        gdual y(1, "y",2);
        auto p1 = 1 / (x + 2*x*y + y*y);
        BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 1);
        BOOST_CHECK_EQUAL(p1.find_cf({1,0}), -3);
        BOOST_CHECK_EQUAL(p1.find_cf({0,1}), -2);
        BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 10);
        BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 9);
        BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 3);
    }

    // the same but calling the gdouble / gdouble overload
    {
        gdual x(0, "x",2);
        gdual y(1, "y",2);

        auto p1 = gdual(1, 2) / (x + 2*x*y + y*y);
        BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 1);
        BOOST_CHECK_EQUAL(p1.find_cf({1,0}), -3);
        BOOST_CHECK_EQUAL(p1.find_cf({0,1}), -2);
        BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 10);
        BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 9);
        BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 3);
    }
}

BOOST_AUTO_TEST_CASE(identities)
{
    // we test some trivial identities
    gdual x(2, "x",3);
    gdual y(3, "y",3);

    auto p1 = x*x+y-x*x*x*x*y-y*x*x;
    auto p2 = y*y-x+y*y*y*y*x-2*x;
    BOOST_CHECK_EQUAL((x + y)*(x + y), x*x + y*y + 2*x*y);
    BOOST_CHECK_EQUAL((p1 + p2)*(p1 + p2), p1*p1 + p2*p2 + 2*p1*p2);
    BOOST_CHECK_EQUAL(x*x*x*x-y*y*y*y, (x-y)*(x+y)*(x*x+y*y));
    BOOST_CHECK_EQUAL(p1*p1*p1*p1-p2*p2*p2*p2, (p1-p2)*(p1+p2)*(p1*p1+p2*p2));

    BOOST_CHECK(EPSILON_COMPARE((p1/p2) * (p2/p1), gdual(1,3), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE((p1/p2) * p2, p1, 1e-12) == true);

    //Uncomment if boost 1.59 macro defined in helpers.hpp is active
    //BOOST_EQUAL_GDUALS_TOL((p1/p2) * (p2/p1), gdual(1,3));
    //BOOST_EQUAL_GDUALS_TOL((p1/p2) * p2, p1);
}


BOOST_AUTO_TEST_CASE(find_cf)
{
    // Exceeding truncation order
    BOOST_CHECK_THROW(gdual("x",4).find_cf({5}),std::invalid_argument);
    BOOST_CHECK_THROW(gdual("x",4).find_cf(std::vector<int>{5}),std::invalid_argument);
    // Mismtach between monomial requested and number of variables
    BOOST_CHECK_THROW(gdual(1,4).find_cf({5}),std::invalid_argument);
    BOOST_CHECK_THROW(gdual(1,4).find_cf(std::vector<int>{5}),std::invalid_argument);

    BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({0}),3);
    BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({1}),1);

    BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({0,1}),1);
    BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({1,0}),1);
}

BOOST_AUTO_TEST_CASE(get_derivative)
{
    // we test some trivial derivatives
    {
    gdual x(1, "x",4);
    gdual y(1, "y",4);
    gdual z(1, "z",4);
    gdual f = x*x*x*x*x + x*y*z*x*x + z*x*y*y*y;
    BOOST_CHECK_EQUAL(f.get_derivative({1,1,1}), 6.);
    BOOST_CHECK_EQUAL(f.get_derivative({2,1,1}), 6.);
    BOOST_CHECK_EQUAL(f.get_derivative({1,2,1}), 6.);
    BOOST_CHECK_EQUAL(f.get_derivative({1,1,2}), 0.);

    BOOST_CHECK_EQUAL(f.get_derivative({4,0,0}),120.);
    BOOST_CHECK_THROW(f.get_derivative({4,1,1}),std::invalid_argument);
    }
    // and a less trivial case
    {
    gdual x(1, "x",8);
    gdual y(1, "y",8);
    gdual f = (x * y + 2 * x * x * y) / (1 + x + y);
    BOOST_CHECK_EQUAL(f.get_derivative({1,1}), 1.);
    BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({2,2}), 0, 1e-14) == true);
    BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({3,3}), 0, 1e-14) == true);
    BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({4,4}), 0, 1e-14) == true);
    }
}

BOOST_AUTO_TEST_CASE(integrate_partial)
{
    // We test some trivial cases where truncation order does not influence the results
    {
    gdual x(1, "x", 4);
    gdual y(1, "y", 4);
    gdual z(1, "z", 4);
    gdual f = x*x*x + x*y*z + z*x*y;
    gdual fx = 3*x*x + y*z + z*y;
    gdual fy = x*z + z*x;
    gdual fz = y*x + x*y;

    BOOST_CHECK_EQUAL(f.partial("x"), fx);
    BOOST_CHECK_EQUAL(f.partial("y"), fy);
    BOOST_CHECK_EQUAL(f.partial("z"), fz);
    }
}

