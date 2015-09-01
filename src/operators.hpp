#ifndef AUDI_OPERATORS_HPP
#define AUDI_OPERATORS_HPP

namespace audi
{

gdual operator+(const gdual &d1, const gdual &d2)
{
	if (d1.get_order() != d2.get_order()) {
		throw std::invalid_argument("different truncation limit");
	}
	gdual retval(d1);
	retval.m_p += d2.m_p;
	return retval;
}

gdual operator-(const gdual &d1, const gdual &d2)
{
	if (d1.get_order() != d2.get_order()) {
		throw std::invalid_argument("different truncation limit");
	}
	gdual retval(d1);
	retval.m_p -= d2.m_p;
	return retval;
}

// This is for biscani boy
gdual operator*(const gdual &d1, const gdual &d2)
{
	if (d1.get_order() != d2.get_order()) {
		throw std::invalid_argument("different truncation limit");
	}
	gdual retval(d1);
	retval.m_p *= d2.m_p;
	return retval;
}

gdual operator/(const gdual& d1, const gdual& d2)
{
	if (d1.get_order() != d2.get_order()) {
		throw std::invalid_argument("different truncation limit");
	}
	gdual retval(1);

    double fatt = 1;
    auto m = d1.get_order();
    auto p0 = p2.find_cf(std::vector<int>(p2.get_symbol_set().size(),0));
    if (p0 == 0) {
        throw 20;
    }
    auto phat = (p2 - p0);
    phat = phat/p0;

    for (auto i = 1u; i <= m; ++i) {
        fatt*=-1;     
        retval+= fatt*math::pow(phat,i);
    }
    return (p1*retval)/p0;   
}

} // end of namespace audi 
#endif