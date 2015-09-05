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
		if (!std::isfinite(pow_p0)) {
			throw std::domain_error("std::pow(double) returned a non finite number");
		}
		if (p0 == 0) {
			throw std::domain_error("exponentiation not defined for polynomial having zero constant term");
		}
		auto phat = d - p0;
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

	gdual sqrt(const gdual& d) 
	{
		return pow(d, 0.5);
	}

	gdual exp(const gdual& d)
	{
	    gdual retval(1, d.get_order());
	    double fact=1;
	    auto p0 = d.find_cf(std::vector<int>(d.get_n_variables(),0));
	    auto phat = d - p0;
		gdual tmp(phat);

		retval+=phat;
	    for (int i = 2; i <= d.get_order(); ++i) {
	    	phat*=tmp;
	        fact*=i;
	        retval+=phat / fact;
	    }
	    return retval * std::exp(p0);
	}

} // end of namespace audi 
#endif