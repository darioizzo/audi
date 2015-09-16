#include <audi/audi.hpp>
#include <iostream>

using namespace audi;

int main() {
	// We want to compute the Taylor expansion of a function f (and thus all derivatives) at x=2, y=3
	// 1 - Define the generalized dual numbers (7 is the truncation order, i.e. the maximum order of derivation we will need)
	gdual x(2, "x", 7);
	gdual y(3, "y", 7);

	// 2 - Compute your function as usual
	gdual f = exp((x*x+cbrt(y)/log(x*y)));

	// 3 - Inspect the results
	std::cout << f << std::endl;						// This is the entire Taylor expansion of f (truncated at the 7th order)
	std::cout << f.get_derivative({4,3}) << std::endl;  // This is the value of the mixed derivative (d^7 / dx^4dy^3)
}