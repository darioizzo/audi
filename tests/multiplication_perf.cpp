#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_mul(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double> > variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x"+std::to_string(i), m);
    }
    gdual<double> p1(1.);
    gdual<double> p2(1.);
    gdual<double> p3(0.);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...
    for (int i = 0u; i < n; ++i) {p2 -= variables[i];} // 1 - x1 - x2 + ...

    p1 = pow(p1,10);
    p2 = pow(p2,10);
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        p3=p1*p2;
    }
}

BOOST_AUTO_TEST_CASE(multiplication_performance)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }
    std::cout << "Testing multiplication of (1 + x1 + .. + xn)^m * (1 - x1 - .. - xn)^m: " << std::endl;
    for (auto m = 10; m < 11; ++m) {
        for (auto n = 10; n < 11; ++n) {
            scalable_mul(m,n);
        }
    }
}
