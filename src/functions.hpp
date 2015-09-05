#ifndef AUDI_FUNCTIONS_HPP
#define AUDI_FUNCTIONS_HPP

#include "gdual.hpp"
#include "piranha/binomial.hpp"

namespace audi
{

	gdual pow(const gdual& d, int n)
	{
		gdual retval(d);
		for (auto i = 1u; i < n; ++i)
		{
			retval*=d;
		}
		return retval;
	}

	gdual pow(const gdual& d, double alpha)
	{
		gdual retval(1, d.get_order());
		auto p0 = d.find_cf(std::vector<int>(d.get_n_variables(),0));
		double pow_p0 = std::pow(p0,alpha);
		if (!std::isfinite(p0) || p0==0) {
			throw std::domain_error("pow operation returns a non finite number");
		}
		auto phat = (d - p0);
		phat = phat/p0;
		gdual tmp(phat);

		retval+=alpha * phat;
		for (auto i = 2; i <= d.get_order(); ++i) {
			phat*=tmp;
			retval+=piranha::math::binomial(alpha,i) * phat;
		}
		retval*=pow_p0;
		return retval;
	}

} // end of namespace audi 
#endif