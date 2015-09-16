#define BOOST_TEST_MODULE dace_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <DACE/DA.h>

#include <vector>

// compile with g++ multiplication_perf_DACE.cpp -std=c++11 -I/usr/local/include -ldace -lboost_system 
// -lboost_unit_test_framework -lboost_timer

using namespace DACE;

void scalable_mul(int m, int n)
{
    // m is order, n is the number of variables
    DA::init(m,n);
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<DA> variables;
    for (int i = 1; i <= n; ++i) {variables.emplace_back(i,1.);} // 1 + x1 + x2 + ...

    DA p1(1.);
    DA p2(1.);
    for (int i = 0; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...
    for (int i = 0; i < n; ++i) {p2 -= variables[i];} // 1 - x1 - x2 + ...
    auto result = p1 * p2;
    auto factor = result;
    for (auto i = 1; i < m; ++i) {
        result*=factor;
    }
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        result*=factor;
    }
    int cicco;
    std::cin >> cicco;
}

int main()
{
    std::cout << "Testing multiplication of (1 + x1 + .. + xn)^m * (1 - x1 - .. - xn)^m: " << std::endl;
    for (auto m = 5; m < 11; ++m) {
        for (auto n = 5; n < 11; ++n) {
            scalable_mul(m,n);
        }
    }
    return 0;
}
