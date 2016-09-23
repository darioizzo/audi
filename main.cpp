#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual<std::complex<double>> xc(std::complex<double>(1., -1.), "x", 3);
	gdual<double> xd(0.5, "x", 3);
	gdual<vectorized_double> xv(std::vector<double>{0.5,-0.5}, "x", 3);

	//auto fc = audi::atanh_d(xc);
	auto fd = audi::atanh_d(xd);
	auto fv = audi::atanh_d(xv);
	//std::cout << fc << std::endl;
	std::cout << fd << std::endl;
	std::cout << fv << std::endl;
}
