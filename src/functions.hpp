#ifndef AUDI_FUNCTIONS_HPP
#define AUDI_FUNCTIONS_HPP

#include "gdual.hpp"

namespace audi
{

	template <typename T>
	gdual pow(const gdual& d, const T& n)
	{
		gdual retval(d);
		for (auto i = 1u; i < n; ++i)
		{
			retval*=d;
		}
		return retval;
	}

} // end of namespace audi 
#endif