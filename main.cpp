#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual x(0.3, "y", 3);
	gdual y(0.2, "x", 3);
	gdual z(0.2, "z", 3);

	auto f = 0.1+x+y-z;
	auto g = atanh(f);
	std::cout << atanh(f) << std::endl;
	std::cout << atanh_d(f) << std::endl;
}