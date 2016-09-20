#include <piranha/polynomial.hpp>
#include <piranha/type_traits.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <boost/timer/timer.hpp>
#include "src/gdual_v.hpp"
#include "src/gdual.hpp"



std::vector<double> random_vector_double(unsigned int N, double lb, double ub)
{
    std::vector<double> retval(N);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(lb, ub);
    for (auto &x : retval) {
        x = dist(mt);
    }
    return retval;
}

using namespace piranha;

int main() {
    using p_type_vector = piranha::polynomial<coefficient_v,piranha::monomial<char>>;
    using p_type = piranha::polynomial<double,piranha::monomial<char>>;

    auto c1  = 0.5;
    auto c2  = 0.23;
    auto c3  = -0.44;
    auto c4  = 0.1;


    unsigned int N = 100u;
    auto lb = -0.5;
    auto ub = 0.5;
    // We generate random vector coefficients
    auto xcoeff = random_vector_double(N, lb, ub);
    auto ycoeff = random_vector_double(N, lb, ub);
    auto zcoeff = random_vector_double(N, lb, ub);
    auto coeff = random_vector_double(N, lb, ub);
    auto coeff_v = audi::gdual_v(coeff);

    {
        audi::gdual_v x{xcoeff, "x",11}, y{ycoeff, "y",11}, z{zcoeff, "z",11}, w{zcoeff, "w",11}, q{zcoeff, "q",11}, r{zcoeff, "r",11}, s{zcoeff, "s",11} ;
        auto foo = coeff_v + x + y + z - w + q + r + s;
        for (int i = 1; i < 25; ++i) {
            foo *= coeff_v + x + y + z - w + q + r + s;
        }

        boost::timer::auto_cpu_timer t;
        foo *= foo;
    }
    {
        boost::timer::auto_cpu_timer t;
        for (auto i = 0u; i < N; ++i) {
            audi::gdual x{xcoeff[i], "x",11}, y{ycoeff[i], "y",11}, z{zcoeff[i], "z",11}, w{zcoeff[i], "w",11}, q{zcoeff[i], "q",11}, r{zcoeff[i], "r",11}, s{zcoeff[i], "s",11} ;

            auto foo = coeff[i] + x + y + z - w + q + r + s;
            for (int j = 1; j < 25; ++j) {
                foo *= coeff[i] + x + y + z - w + q + r + s;
            }
            foo *= foo;
        }
    }

}
