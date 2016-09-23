#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual<std::complex<double>> xc(std::complex<double>(1., -1.), "x", 3);
	gdual<double> xd(-1, "x", 3);
	gdual<vectorized_double> xv(std::vector<double>{1,-1}, "x", 3);

	//auto fc = audi::abs(xc);
	auto fd = audi::abs(xd);
	auto fv = audi::abs(xv);
	//std::cout << fc << std::endl;
	std::cout << fd << std::endl;
	std::cout << fv << std::endl;
}
