#define BOOST_TEST_MODULE dace_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <DACE/DA.h>

#include <vector>

using namespace DACE;

void scalable_mul(int m, int n)
{
    DA::init(m,n);
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<DA> variables(n, 1./n);
    DA p1(variables[0]+variables[1]);
    DA p2(2 - variables[0]);
    for (int i = 2u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...
    for (int i = 1u; i < n; ++i) {p2 -= variables[i];} // 1 - x1 - x2 + ...
std::cout<< p1 << std::endl;
    auto result = p1 * p2;
    for (auto i = 1; i < m-1; ++i) {
        result*=result;
    }
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        result*=result;
    }
std::cout << result << std::endl;
}

int main()
{
    std::cout << "Testing multiplication of (1 + x1 + .. + xn)^m * (1 - x1 - .. - xn)^m: " << std::endl;
    for (auto m = 5; m < 6; ++m) {
        for (auto n = 5; n < 6; ++n) {
            scalable_mul(m,n);
        }
    }
    return 0;
}