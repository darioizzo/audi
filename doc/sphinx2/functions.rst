Mathematical functions
========================

Direct implementation
^^^^^^^^^^^^^^^^^^^^^

The following functions are implemented computing the truncated Taylor series of the argument:

 * exp
 * log
 * pow (4 overloads)

-----------------------------------------------------------------------

*#include <audi/functions.hpp>*

.. cpp:function:: template <typename T> inline T exp(const T &d)

   This templated function is enabled only if **T** is a :cpp:class:`gdual`.
   It implements the exponential performing the following computations in the :math:`\mathcal P_{n,m}` algebra:

   .. math::

      T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ... \right) 

   where :math:`T_f = f_0 + \hat f`.

   :param d: :cpp:class:`gdual` argument

   :return: the Taylor expansion of the exponential **d**

-----------------------------------------------------------------------

*#include <audi/functions.hpp>*

.. cpp:function:: template <typename T> inline T log(const T &d)

   This templated function is enabled only if **T** is a :cpp:class:`gdual`.
   It implements the logarithm performing the following computations in the :math:`\mathcal P_{n,m}` algebra:

   .. math::

      T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ... 

   where :math:`T_f = f_0 + \hat f`.

   :param d: :cpp:class:`gdual` argument

   :return: the Taylor expansion of the logarithm **d**

-----------------------------------------------------------------------

*#include <audi/functions.hpp>*

.. cpp:function:: template <typename T> inline T pow(double base, const T &d)

   This templated function is enabled only if **T** is a :cpp:class:`gdual`.
   It computes the exponentiation of a double to the power of a :cpp:class:`gdual`.
   If the exponent is a constant gdual, it calls the std::pow overload. Otherwise
   it converts :math:`a^{T_f}` to :math:`\exp(T_f \log(a))` and computes this
   last expression instead.

   :param base: base for the exponentiation
   :param d: :cpp:class:`gdual` exponent

   :return: the Taylor expansion *base* to the power of **d**

-----------------------------------------------------------------------

*#include <audi/functions.hpp>*

.. cpp:function:: template <typename T> inline T pow(const T &d, double alpha)

   This templated function is enabled only if **T** is a :cpp:class:`gdual`.
   It computes the exponentiation of a :cpp:class:`gdual` when the exponent is not an integer.

   .. math::
      
      T_{(f^\alpha)} = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\hat f / f_0\right)^k

   where,
   
   .. math:: 
   
      T_f = f_0 + \hat f

   :param d: a :cpp:class:`gdual` base for the exponentiation
   :param alpha: exponent

   :return: the Taylor expansion *base* to the power of **d**

-----------------------------------------------------------------------

Implementation from derivative
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions are implemented using the definition of their derivative:

 * erf: error function

-----------------------------------------------------------------------

*#include <audi/functions_from_d.hpp>*

.. cpp:function:: template <typename T> inline T erf(const T &d)

   This templated function is enabled only if **T** is a :cpp:class:`gdual`
   Essentially, it makes use of the definition:

   .. math::
   
      \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}}\exp(-x^2)
   
   where :math:`T_f = f_0 + \hat f`.
   
   :param d: :cpp:class:`gdual` argument
   
   :return: the Taylor expansion of the error function of **d**

