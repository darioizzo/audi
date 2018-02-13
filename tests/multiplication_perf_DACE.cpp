#define BOOST_TEST_MODULE dace_gdual_test
#include <DACE/DA.h>
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <vector>

// compile with g++ multiplication_perf_DACE.cpp -std=c++14 -I/usr/local/include -ldace -lboost_system
// -lboost_unit_test_framework -lboost_timer

using namespace DACE;

void scalable_mul(unsigned int m, unsigned int n, double value)
{
    // m is order, n is the number of variables
    DA::init(m, n);
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<DA> variables;
    for (auto i = 1u; i <= n; ++i) {
        variables.emplace_back(i, value);
    } // [x1, x2,...,xn]

    DA p1(1.);
    DA p2(1.);
    DA res;
    for (int i = 0; i < n; ++i) {
        p1 += variables[i];
    } // 1 + x1 + x2 + ...
    for (int i = 0; i < n; ++i) {
        p2 -= variables[i];
    } // 1 - x1 - x2 + ...
    p1 = pow(p1, 10);
    p2 = pow(p2, 10);
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        res = p1 * p2;
    }
}

int main()
{
    std::cout << "Testing multiplication of (1 + x1 + .. + xn)^m * (1 - x1 - .. - xn)^m: " << std::endl;
    for (auto m = 10u; m < 11u; ++m) {
        for (auto n = 10u; n < 11u; ++n) {
            scalable_mul(m, n, 1.);
        }
    }
    return 0;
}
