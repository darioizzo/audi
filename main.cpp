#include <audi/audi.hpp>

using namespace audi;

int main()
{
    gdual_d x(0.1, "x", 2);
    gdual_d y(0.1, "y", 2);
    gdual_d z(0.1, "z", 2);
    // Order extraction
    auto f = audi::sin(x * y);
    auto f2 = audi::sin(z);
    std::cout << "Full map: " << f << std::endl;
    auto g = f.extract_terms(2);
    std::cout << "Only order 2: " << g << std::endl;
    // Map concatenation
    std::cout << "\n\nFull map (x,y): " << f << std::endl;
    std::cout << "Full map (z): " << f2 << std::endl;
    auto new_f = f.subs("dy", f2);
    std::cout << "Full map (x,z): " << new_f << std::endl;
    std::cout << new_f.info() << std::endl;
}
