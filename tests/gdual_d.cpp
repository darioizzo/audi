#define BOOST_TEST_MODULE audi_gduald_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/audi.hpp>
#include <audi/functions.hpp>


using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Constructing a constant with a symbol
    gdual_d x(0, "x", 0);
    gdual_d y(3., "y", 0);
    gdual_d z(2., "z", 3);
    auto p = (x + y) * z * (x + y) * z * (x + y) * z * (x + y) * z * (x + y) * z;
    BOOST_CHECK_EQUAL(x.get_order(), 0);
    BOOST_CHECK_EQUAL(y.get_order(), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1, 0, 0}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({0, 1, 0}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({0, 1, 1}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1, 1, 1}), 0);
    BOOST_CHECK_EQUAL(p.get_derivative({1, 0, 2}), 0);

    // TODO: other constructors
}

BOOST_AUTO_TEST_CASE(order_promotion)
{
    gdual_d c(1);
    gdual_d x(0, "x", 4);
    gdual_d y(0, "y", 2);
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
    gdual_d p1(1);
    gdual_d p2(0, "x", 4);
    BOOST_CHECK_EQUAL(p1 + p2, p2 + p1);
    BOOST_CHECK_EQUAL(1 + p2, p2 + 1);
    BOOST_CHECK_EQUAL(1. + p2, p2 + 1.);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
    gdual_d p1(1);
    gdual_d p2(0, "x", 4);
    BOOST_CHECK_EQUAL(p1 - p1, gdual_d(0));
    BOOST_CHECK_EQUAL(p2 - p2, gdual_d(0));
    BOOST_CHECK_EQUAL(p1 - p2, -(p2 - p1));
    BOOST_CHECK_EQUAL(1 + p1 - p1, gdual_d(1));
    BOOST_CHECK_EQUAL(1 + p2 - p2, gdual_d(1));
    BOOST_CHECK_EQUAL(1. + p1 - p1, gdual_d(1));
    BOOST_CHECK_EQUAL(1. + p2 - p2, gdual_d(1));
    BOOST_CHECK_EQUAL((1 - p1) + (p1 + p2), 1 + p2);
}

BOOST_AUTO_TEST_CASE(multiplication)
{
    // we test normalized derivatives computations of f = x + 3*x*y + y*y
    gdual_d x(3, "x", 2);
    gdual_d y(7, "y", 2);

    auto p1 = x + 3 * x * y + y * y;
    BOOST_CHECK_EQUAL(p1.find_cf({0, 0}), 115);
    BOOST_CHECK_EQUAL(p1.find_cf({1, 0}), 22);
    BOOST_CHECK_EQUAL(p1.find_cf({0, 1}), 23);
    BOOST_CHECK_EQUAL(p1.find_cf({1, 1}), 3);
    BOOST_CHECK_EQUAL(p1.find_cf({2, 0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0, 2}), 1);

    // we test the truncation order (5) on a 3 variables polynomial p2 = (1 + x0 + x1 + x2)^10
    std::vector<gdual_d> vars;
    for (auto i = 0u; i < 3; ++i) {
        vars.emplace_back(0, "x" + std::to_string(i), 5);
    }
    gdual_d p2(1);
    for (auto var : vars) {
        p2 += var;
    }
    p2 = p2 * p2 * p2 * p2 * p2 * p2;
    BOOST_CHECK_EQUAL(p2.get_order(), 5);
}

BOOST_AUTO_TEST_CASE(division)
{
    // we test normalized derivatives computations of f = 1 / (x + 2*x*y + y*y)
    {
        gdual_d x(0, "x", 2);
        gdual_d y(1, "y", 2);
        auto p1 = 1 / (x + 2 * x * y + y * y);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 0}), 1);
        BOOST_CHECK_EQUAL(p1.find_cf({1, 0}), -3);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 1}), -2);
        BOOST_CHECK_EQUAL(p1.find_cf({1, 1}), 10);
        BOOST_CHECK_EQUAL(p1.find_cf({2, 0}), 9);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 2}), 3);
    }

    // the same but calling the gdouble / gdouble overload
    {
        gdual_d x(0, "x", 2);
        gdual_d y(1, "y", 2);

        auto p1 = gdual_d(1) / (x + 2 * x * y + y * y);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 0}), 1);
        BOOST_CHECK_EQUAL(p1.find_cf({1, 0}), -3);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 1}), -2);
        BOOST_CHECK_EQUAL(p1.find_cf({1, 1}), 10);
        BOOST_CHECK_EQUAL(p1.find_cf({2, 0}), 9);
        BOOST_CHECK_EQUAL(p1.find_cf({0, 2}), 3);
    }
    { // added when a bug was fixed in the corresponding called template
        gdual_d x(123., "x", 2);
        BOOST_CHECK_EQUAL(7.1 / x, gdual_d(7.1) / x);
        BOOST_CHECK_EQUAL(7.1 / x, 1. / (x / 7.1));
    }
}

BOOST_AUTO_TEST_CASE(identities)
{
    // we test some trivial identities
    gdual_d x(2, "x", 3);
    gdual_d y(3, "y", 3);

    auto p1 = x * x + y - x * x * x * x * y - y * x * x;
    auto p2 = y * y - x + y * y * y * y * x - 2 * x;
    BOOST_CHECK_EQUAL((x + y) * (x + y), x * x + y * y + 2 * x * y);
    BOOST_CHECK_EQUAL((p1 + p2) * (p1 + p2), p1 * p1 + p2 * p2 + 2 * p1 * p2);
    BOOST_CHECK_EQUAL(x * x * x * x - y * y * y * y, (x - y) * (x + y) * (x * x + y * y));
    BOOST_CHECK_EQUAL(p1 * p1 * p1 * p1 - p2 * p2 * p2 * p2, (p1 - p2) * (p1 + p2) * (p1 * p1 + p2 * p2));

    BOOST_CHECK(EPSILON_COMPARE((p1 / p2) * (p2 / p1), gdual_d(1), 1e-12) == true);
    BOOST_CHECK(EPSILON_COMPARE((p1 / p2) * p2, p1, 1e-12) == true);
}

BOOST_AUTO_TEST_CASE(find_cf)
{
    // Exceeding truncation order
    BOOST_CHECK_THROW(gdual_d(0, "x", 4).find_cf({5}), std::invalid_argument);
    BOOST_CHECK_THROW(gdual_d(0, "x", 4).find_cf(std::vector<int>{5}), std::invalid_argument);
    // Mismtach between monomial requested and number of variables
    BOOST_CHECK_THROW(gdual_d(1).find_cf({5}), std::invalid_argument);
    BOOST_CHECK_THROW(gdual_d(1).find_cf(std::vector<int>{5}), std::invalid_argument);

    BOOST_CHECK_EQUAL((gdual_d(3) + gdual_d(0, "x", 4)).find_cf({0}), 3);
    BOOST_CHECK_EQUAL((gdual_d(3) + gdual_d(0, "x", 4)).find_cf({1}), 1);

    BOOST_CHECK_EQUAL((gdual_d(0, "x", 4) + gdual_d(0, "y", 4)).find_cf({0, 1}), 1);
    BOOST_CHECK_EQUAL((gdual_d(0, "x", 4) + gdual_d(0, "y", 4)).find_cf({1, 0}), 1);
}

BOOST_AUTO_TEST_CASE(get_derivative)
{
    // we test some trivial derivatives
    {
        gdual_d x(1, "x", 4);
        gdual_d y(1, "y", 4);
        gdual_d z(1, "z", 4);
        gdual_d f = x * x * x * x * x + x * y * z * x * x + z * x * y * y * y;
        BOOST_CHECK_EQUAL(f.get_derivative({1, 1, 1}), 6.);
        BOOST_CHECK_EQUAL(f.get_derivative({2, 1, 1}), 6.);
        BOOST_CHECK_EQUAL(f.get_derivative({1, 2, 1}), 6.);
        BOOST_CHECK_EQUAL(f.get_derivative({1, 1, 2}), 0.);

        BOOST_CHECK_EQUAL(f.get_derivative({4, 0, 0}), 120.);
        BOOST_CHECK_THROW(f.get_derivative({4, 1, 1}), std::invalid_argument);
    }
    // and a less trivial case
    {
        gdual_d x(1, "x", 8);
        gdual_d y(1, "y", 8);
        gdual_d f = (x * y + 2 * x * x * y) / (1 + x + y);
        BOOST_CHECK_EQUAL(f.get_derivative({1, 1}), 1.);
        BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({2, 2}), 0, 1e-14) == true);
        BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({3, 3}), 0, 1e-14) == true);
        BOOST_CHECK(EPSILON_COMPARE(f.get_derivative({4, 4}), 0, 1e-14) == true);
    }
    /// we test the dictionary interface
    gdual_d x(1, "x", 4);
    gdual_d y(1, "y", 4);
    gdual_d z(1, "z", 4);
    gdual_d f = x * x * x * x * x + x * y * z * x * x + z * x * y * y * y;
    std::unordered_map<std::string, unsigned int> dict{{"dx", 1u}, {"dy", 1u}, {"dz", 1u}};
    BOOST_CHECK_EQUAL(f.get_derivative(dict), 6.);
    dict = {{"dx", 4u}, {"dy", 1u}, {"dz", 1u}};
    BOOST_CHECK_THROW(f.get_derivative(dict), std::invalid_argument);
    dict = {{"dx", 1u}, {"dr", 1u}, {"dz", 1u}};
    BOOST_CHECK_EQUAL(f.get_derivative(dict), 0.);
}

BOOST_AUTO_TEST_CASE(integrate_partial)
{
    // We test some trivial cases where truncation order does not influence the results
    {
        gdual_d x(1, "x", 4);
        gdual_d y(1, "y", 4);
        gdual_d z(1, "z", 4);
        gdual_d f = x * x * x + x * y * z + z * x * y;
        gdual_d fx = 3 * x * x + y * z + z * y;
        gdual_d fy = x * z + z * x;
        gdual_d fz = y * x + x * y;

        BOOST_CHECK_EQUAL(f.partial("x"), fx);
        BOOST_CHECK_EQUAL(f.partial("y"), fy);
        BOOST_CHECK_EQUAL(f.partial("z"), fz);
    }
}

BOOST_AUTO_TEST_CASE(is_zero)
{
    // We test some trivial cases where truncation order does not influence the results
    {
        gdual_d x(1, "x", 4);
        gdual_d y(1, "y", 4);
        gdual_d z(1, "z", 4);
        gdual_d f = x * x * x + x * y * z + z * x * y;

        BOOST_CHECK((f - f).is_zero(1e-12));
        BOOST_CHECK((f - 1 / (1 / f)).is_zero(1e-12));
        BOOST_CHECK(((f * f) / (f)-f).is_zero(1e-12));
        BOOST_CHECK(!f.is_zero(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(subs)
{
    // Trivial substitution 1 + dx -> 1 + dy
    {
        gdual_d x(1, "x", 3);
        gdual_d y(0., "y", 5);
        auto z = x.subs("dx", y);
        BOOST_CHECK((z - 1. - y).is_zero(1e-14));
        BOOST_CHECK(z.get_order() == 3);
        BOOST_CHECK(z.get_symbol_set_size() == y.get_symbol_set_size());
    }
    // Function composition works on truncated Taylor series if the constant coefficient is zero
    // and we then take advantage of the nihilpotency property.
    {
        gdual_d x(0., "x", 7);
        gdual_d y(0., "y", 7);
        auto sx = log(x + 1) / cos(x);
        auto sy = log(y + 1) / cos(y);
        auto ssy = sx.subs("dx", sy);
        auto ssy2 = log(log(y + 1) / cos(y) + 1) / cos(log(y + 1) / cos(y));
        BOOST_CHECK((ssy - ssy2).is_zero(1e-14));
    }
    // Testing the substitution symbol->value. The Taylor expansion of sin(0+dx) evaluated in dx=1 is
    // compared against the value of sin(1)
    {
        gdual_d x(0, "x", 100);
        auto sx = sin(x);
        auto sx1 = sx.subs("dx", 1.);
        BOOST_CHECK_CLOSE(sx1.constant_cf(), std::sin(1.), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(extract_order)
{
    unsigned int order = 8u;
    gdual_d x(0.123, "x", order);
    gdual_d y(0.456, "y", order);
    auto f = audi::sin(x * y);
    // We test that the extracted gduals have the requested order
    for (auto i = 0u; i <= order; ++i) {
        auto fi = f.extract_terms(i);
        BOOST_CHECK_EQUAL(fi.degree(), fi.get_order());
        BOOST_CHECK_EQUAL(fi.get_order(), i);
    }
    // We test that f = f0+f1+f2+...+fn
    std::vector<gdual_d> terms;
    for (auto i = 0u; i <= order; ++i) {
        auto f2 = f.extract_terms(i);
        terms.push_back(f2);
    }
    auto sum = std::accumulate(terms.begin(), terms.end(), gdual_d(0.));
    BOOST_CHECK((sum - f).is_zero(0.));
    // And we test the throw
    BOOST_CHECK_THROW(f.extract_terms(order + 1), std::invalid_argument);
    BOOST_CHECK_NO_THROW(f.extract_terms(order));
}

BOOST_AUTO_TEST_CASE(serialization_test)
{
    gdual_d x(1, "x", 4);
    gdual_d y(1, "y", 4);
    gdual_d z(1, "z", 4);
    gdual_d f = (x * x * x + x * y * z + z * x * y) * (x * x * x + x * y * z + z * x * y)
                * (x * x * x + x * y * z + z * x * y) * (x * x * x + x * y * z + z * x * y);

    // create and open a character archive for output
    std::ofstream ofs("test.dump");

    // save data to archive
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << f;
        // archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    gdual_d newf;
    {
        // create and open an archive for input
        std::ifstream ifs("test.dump");
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> newf;
        // archive and stream closed when destructors are called
    }
    BOOST_CHECK(newf == f);
    BOOST_CHECK(newf.get_order() == f.get_order());
}

BOOST_AUTO_TEST_CASE(trim)
{
    gdual_d x(1e-5, "x", 4);
    gdual_d y(1e-3, "x", 4);
    BOOST_CHECK(x.trim(1e-6) == x);
    BOOST_CHECK(y.trim(1e-2) == gdual_d(0., "x", 1));
    BOOST_CHECK(y.trim(100) == gdual_d(0., "x", 0));
    BOOST_CHECK_THROW(x.trim(-1e-3), std::invalid_argument);
}