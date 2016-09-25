#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual_c xc(std::complex<double>(1., -1.), "x", 3);
	gdual_d xd(0.5, "x", 3);
	gdual_v xv(std::vector<double>{0.5,-0.5}, "x", 3);

	auto fc = audi::atanh_d(xc);
	auto fd = audi::atanh(xc);
	auto fv = audi::atanh_d(xv);

	std::cout << fc << std::endl;
	std::cout << fd << std::endl;
	std::cout << fv << std::endl;

	gdual_d x(1.2, "x", 2);
    gdual_d y(-0.1, "y", 2);
    std::cout << (x+y) / (x-y) << "\n";

}
