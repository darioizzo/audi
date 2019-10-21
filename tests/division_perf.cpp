#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <vector>

#include <boost/optional.hpp>

#include <tbb/task_scheduler_init.h>

#include <audi/gdual.hpp>

using namespace audi;

void scalable_div(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double>> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x" + std::to_string(i), m);
    }
    gdual<double> denom(1.);
    gdual<double> num(1.);
    for (int i = 0u; i < n; ++i) {
        num += variables[i] * variables[i];
    }
    for (int i = 0u; i < n; ++i) {
        denom += variables[0] * variables[i];
    }
    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    auto result = num / denom;
}

BOOST_AUTO_TEST_CASE(division_performance)
{
    boost::optional<tbb::task_scheduler_init> tinit;
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        tinit.emplace(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }
    std::cout << "Testing division of (x1^2 + ... + xn^2) / (x1 * (x1 + x2 + .. + xn) ): " << std::endl;
    for (auto m = 5; m < 10; ++m) {
        for (auto n = 5; n < 12; ++n) {
            scalable_div(m, n);
        }
    }
}
