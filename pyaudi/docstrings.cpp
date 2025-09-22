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

gdual_real128: :math:`m = 2`, :math:`\mathbf a = [1.2, -0.1]`, :math:`f = \frac{x1+x2}{x1-x2}`

>>> from pyaudi import gdual_real128 as gdual
>>> x1 = gdual("1.2", "x1", 2)
>>> x2 = gdual("-0.1", "x2", 2)
>>> f = (x1+x2) / (x1-x2)
>>> print(f) # doctest: +SKIP
1.42011834319526627218934911242603510e+00*dx2+1.18343195266272189349112426035503101e-01*dx1+1.09239872553482020937642239417387321e+00*dx2**2+8.46153846153846153846153846153845954e-01-1.00136549840691852526172052799271705e+00*dx1*dx2-9.10332271279016841147018661811561892e-02*dx1**2

.. note::

   A gdual can operate on doubles (gdual_double), on vectorized doubles (gdual_vdouble) or on quadruple precision doubles 
   (gdual_real128). The second case is immensely more efficient when applicable.

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

std::string taylor_model_docstring()
{
    return R"(
This class represents a Taylor model, containing a generalized dual number in combination with an
interval. This implementation is derived from Makino (1998) and is built on top of a gdual object
representing the Taylor polynomial (generalized dual number) part of the Taylor model.                             
                                                                                                                   }
A key concept here is Taylor's theorem (formulation used from Makino (1998) p.80) which allows for a
quantitative estimate of the error that is to be expected when approximating a function by its Taylor polynomial.
Furthermore it even offers a way to obtain bounds for the error in practice based on bounding the :math:`(n+1)` th 
derivative a method that has sometimes been employed in interval calculations.
                                                                                                                   
As a result, you get :math:`\forall \vec{x} \in [\vec{a}, \vec{b}]` , a given order  :math:`n` , and an expansion
point :math:`\vec{x_o}`:
                                                                                                                 
:math:`f(\vec{x}) \in P_{\alpha, f}(\vec{x} - \vec{x_0}) + I_{\alpha, f}`
                                                                                                                   
where f is the function you're creating a Taylor model for, P is the Taylor polynomial, and I is
the interval remainder.

A basic example would be:

>>> from pyaudi import gdual_double as gdual, taylor_model, int_d
>>> domain_size = 0.01
>>> exp_points = {"x": 1.1, "y": 1.2}
>>> dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
>>> rem = int_d(0.0, 0.0)
>>> tpol = gdual(exp_points["x"], "x", 6)
>>> x = taylor_model(tpol, rem, {"x": exp_points["x"]}, {"x": dom["x"]})
>>> y = taylor_model(tpol, rem, {"y": exp_points["y"]}, {"y": dom["y"]})
>>> test = 2 * x + 3 * y
>>> print(test)
  )";
}

} // namespace pyaudi
