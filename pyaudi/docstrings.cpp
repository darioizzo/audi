#include <string>

#include <pyaudi/docstrings.hpp>

namespace pyaudi
{

std::string gdual_docstring()
{
    return R"(The generalized dual number class

This class represents a generalized dual number, in a nutshell, a truncated multivariate Taylor polynomial.
A gdual is defined by its truncation order :math:`m` as well as on its expansion point :math:`\mathbf a`. All arithmetic
operators +,*,/,-,** are overloaded so that the Taylor expansion of arithmetic computations is obtained.

A basic example is:

gdual_double: :math:`m = 2`, :math:`\mathbf a = [1.2, -0.1]`, :math:`f = \frac{x1+x2}{x1-x2}`

>>> from pyaudi import gdual_double as gdual
>>> x1 = gdual(1.2, "x1", 2)
>>> x2 = gdual(-0.1, "x2", 2)
>>> f = (x1+x2) / (x1-x2)
>>> print(f) # doctest: +SKIP
1.42012*dx2+0.118343*dx1+1.0924*dx2**2+0.846154-1.00137*dx1*dx2-0.0910332*dx1**2

gdual_vdouble: :math:`m = 2`, :math:`\mathbf a = [[1.2, 1.1], [-0.1, -0.2]]`, :math:`f = \frac{x1+x2}{x1-x2}`

>>> from pyaudi import gdual_vdouble as gdual
>>> x1 = gdual([1.2, 1.1], "x1", 2)
>>> x2 = gdual([-0.1, -0.2], "x2", 2)
>>> f = (x1+x2) / (x1-x2)
>>> print(f) # doctest: +SKIP
[1.42012, 1.30178]*dx2+[0.118343, 0.236686]*dx1+[1.0924, 1.00137]*dx2**2+[0.846154, 0.692308]+[-1.00137, -0.819299]*dx1*dx2+[-0.0910332, -0.182066]*dx1**2


.. note::

   A gdual can operate on doubles (gdual_double) or on vectorized doubles (gdual_vdouble). The second case is
   immensely more efficient when applicable.

See also the docs of the C++ class :cpp:class:`gdual`.

)";
}

std::string gdual_degree_docstring()
{
    return R"(The degree of the generalized dual number

.. note::

   The degree and the order of a gdual can be different, the degree is always smaller or equal ot the order.

:return: The degree of the underlying Taylor polynomial

>>> from pyaudi import gdual_double as gdual
>>> x1 = gdual(1.2, "x1", 5)
>>> f = x1**2
>>> print(f.degree)
2

See also the docs of the C++ method of the gdual :cpp:func:`degree`.

)";
}


} // namespace
