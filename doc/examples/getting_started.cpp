#include <audi.hpp>
#include <iostream>

using namespace audi;

int main() {
	// Define the generalized dal numbers (maximum order is 4)
	gdual x("dx", 4);
	gdual y("dy", 4);
	x+=2;
	y+=5;

	// Compute your function as usual
	gdual f = exp((x*x+cbtr(y)/log(x*y)));

	// Inspect the results
	std::cout << f << std::endl;
	std::cout << f.derivative({2,3}) << std::endl;
	std::cout << p2+1 << std::endl;
}