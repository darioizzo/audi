#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual<std::complex<double>> x(std::complex<double>(1., 1.), "x", 3);

	auto f = exp(x);
	std::cout << f << std::endl;
}
