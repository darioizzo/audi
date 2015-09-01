#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual p1(1, 4);
	gdual p2("x", 4);
	std::cout << p1+p2 << std::endl;
	std::cout << 1+p2 << std::endl;
	std::cout << p2+1 << std::endl;
}