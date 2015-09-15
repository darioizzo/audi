#include <audi/audi.hpp>
#include <iostream>

using namespace audi;

int main() {
	// Define the generalized dal numbers (4 is the truncation order, i.e. the maximum order of derivation we will need)
	gdual x("dx", 7);
	gdual y("dy", 7);
	x+=2;
	y+=5;

	// Compute your function as usual
	gdual f = exp((x*x+cbtr(y)/log(x*y)));

	// Inspect the results
	std::cout << f << std::endl;					// This is the entire Taykor expansion of f (truncated at the 7th order)
	std::cout << f.derivative({4,3}) << std::endl;  // This is the value of the mixed derivative (d^7 / dx^4dy^3)
}