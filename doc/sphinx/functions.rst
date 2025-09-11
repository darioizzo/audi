Mathematical functions
========================

The following functions are implemented computing the truncated Taylor series of the argument:

 * :cpp:func:`exp`, :cpp:func:`log`, :cpp:func:`pow` (4 overloads), :cpp:func:`sqrt`, :cpp:func:`cbrt`
 * :cpp:func:`sin`, :cpp:func:`cos`, :cpp:func:`sin_and_cos`, :cpp:func:`tan`
 * :cpp:func:`atan`,  :cpp:func:`asin`, :cpp:func:`acos`
 * :cpp:func:`sinh`, :cpp:func:`cosh`, :cpp:func:`sinh_and_cosh`, :cpp:func:`tanh`
 * :cpp:func:`atanh`, :cpp:func:`asinh`, :cpp:func:`acosh`
 * :cpp:func:`abs`


The following functions are implemented using the definition of their derivative:

 * :cpp:func:`erf`: error function

-----------------------------------------------------------------------

*#include <audi/functions.hpp>*

.. doxygenfunction:: audi::exp(const gdual<T, M> &d)
.. doxygenfunction:: audi::log(const gdual<T, M> &d)
.. doxygenfunction:: audi::pow(U base, const gdual<T, M> &d)
.. doxygenfunction:: audi::binomial(T x, unsigned y)
.. doxygenfunction:: audi::pow(const gdual<T, M> &d, U alpha)
.. doxygenfunction:: audi::pow(const gdual<T, M> &d, int n)
.. doxygenfunction:: audi::pow(const gdual<T, M> &d1, const gdual<T, M> &d2)
.. doxygenfunction:: audi::sqrt(const gdual<T, M> &d)
.. doxygenfunction:: audi::cbrt(const gdual<T, M> &d)
.. doxygenfunction:: audi::sin(const gdual<T, M> &d)
.. doxygenfunction:: audi::cos(const gdual<T, M> &d)
.. doxygenfunction:: audi::sin_and_cos(const gdual<T, M> &d)
.. doxygenfunction:: audi::tan(const gdual<T, M> &d)
.. doxygenfunction:: audi::sinh(const gdual<T, M> &d)
.. doxygenfunction:: audi::cosh(const gdual<T, M> &d)
.. doxygenfunction:: audi::sinh_and_cosh(const gdual<T, M> &d)
.. doxygenfunction:: audi::tanh(const gdual<T, M> &d)
.. doxygenfunction:: audi::atanh(const gdual<T, M> &d)
.. doxygenfunction:: audi::atan(const gdual<T, M> &d)
.. doxygenfunction:: audi::asinh(const gdual<T, M> &d)
.. doxygenfunction:: audi::acosh(const gdual<T, M> &d)
.. doxygenfunction:: audi::asin(const gdual<T, M> &d)
.. doxygenfunction:: audi::acos(const gdual<T, M> &d)
.. doxygenfunction:: audi::abs(const gdual<T, M> &d)
.. doxygenfunction:: audi::abs(const gdual<vectorized<T>, M> &d)

.. .. cpp:function:: template <typename T> inline T exp(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It implements the exponential performing the following computations in the :math:`\mathcal P_{n,m}` algebra:
..
..    .. math::
..
..       T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ... \right) 
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the exponential **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T log(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It implements the logarithm performing the following computations in the :math:`\mathcal P_{n,m}` algebra:
..
..    .. math::
..
..       T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ... 
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the logarithm **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T pow(double base, const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It computes the exponentiation of a double to the power of a :cpp:class:`gdual`.
..    If the exponent is a constant gdual, it calls the std::pow overload. Otherwise
..    it converts :math:`a^{T_f}` to :math:`\exp(T_f \log(a))` and computes this
..    last expression instead.
..
..    :param base: base for the exponentiation
..    :param d: :cpp:class:`gdual` exponent
..
..    :return: the Taylor expansion *base* to the power of **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T pow(const T &d, double alpha)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It computes the exponentiation of a :cpp:class:`gdual` when the exponent is not an integer.
..
..    .. math::
..       
..       T_{(f^\alpha)} = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\hat f / f_0\right)^k
..
..    where,
..    
..    .. math:: 
..    
..       T_f = f_0 + \hat f
..
..    :param d: a :cpp:class:`gdual` base for the exponentiation
..    :param alpha: exponent
..
..    :return: the Taylor expansion **d** to the power of **alpha**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T pow(const T &d, int n)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It implements the integer exponentiation of a :cpp:class:`gdual`. Essentially,
..    it uses the :math:`\mathcal P_{n,m}` multiplication on **d** **n** times
..
..    :param d: a :cpp:class:`gdual` base for the exponentiation
..    :param n: integer exponent
..
..    :return: the Taylor expansion of **d** to the power of **n**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T (const T &d1, const T &d2)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    It implements the exponentiation of a :cpp:class:`gdual` when the exponent is also 
..    a :cpp:class:`gdual`. It computes the result as :math:`\exp(d2\log d1)`
..
..    :param d1: a :cpp:class:`gdual` base for the exponentiation
..    :param d2: a :cpp:class:`gdual` exponent
..
..    :return: the Taylor expansion of **d1** to the power of **d2**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T sqrt(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the square root of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{\sqrt{f}} = \sqrt{f_0} \sum_{k=0}^m {\frac 12 \choose k} \left(\hat f / f_0\right)^k
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the square root of **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T cbrt(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the cubic root of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{\sqrt[3]{f}} = \sqrt[3]{f_0} \sum_{k=0}^m {\frac 13 \choose k} \left(\hat f / f_0\right)^k
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the cubic root of **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T sin(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the sine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0 
..       \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the sine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T cos(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the cosine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\cos f)} = \cos f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0
..       \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the cosine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline std::array<T, 2> sin_and_cos(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the sine and cosine of a :cpp:class:`gdual`.
..    As most of the computations for the sine and cosine is the same, it is twice as fast
..    to get both sine and cosine at once rather than computing them in sequence.
..    Use this function when both sine and cosine are needed.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansions of sine and the cosine of **d** (first element, second element)
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T tan(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the tangent of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\tan f)} = \frac{\tan f_0 + \sum_{k=1}^{k \le 2k+1} B_{2k} \frac{(-4)^k(1-4^k)}{2k!}x^{2k - 1}}{1 - \tan f_0
..       \sum_{k=1}^{k \le 2k+1} \frac{B_{2k}(-4)^k(1-4^k)}{2k!}x^{2k - 1} }
..
..    where :math:`T_f = f_0 + \hat f` and :math:`B_{2k}\f` are the Bernoulli numbers.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the tangent of **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions.hpp>*
..
.. .. cpp:function:: template <typename T> inline T sinh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the hyperbolic sine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\sin f)} = \sinh f_0 \left(\sum_{i=0}^{2i\le m} \frac{\hat f^{2i}}{(2i)!}\right) + \cosh f_0
..       \left(\sum_{i=0}^{(2i+1)\le m} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the hyperbolic sine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T cosh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the hyperbolic cosine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\sin f)} = \cosh f_0 \left(\sum_{i=0}^{2i\le m} \frac{\hat f^{2i}}{(2i)!}\right) + \sinh f_0
..       \left(\sum_{i=0}^{(2i+1)\le m} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the hyperbolic cosine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline std::array<T, 2> sinh_and_cosh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the hyperbolic sine and cosine of a :cpp:class:`gdual`.
..    As most of the computations for the hyperbolic sine and cosine are the same, it is twice as fast
..    to get both the hyperbolic sine and cosine at once rather than computing them in sequence.
..    Use this function when both hyperbolic sine and cosine are needed.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansions of hyperbolic sine and the cosine of **d** (first element, second element)
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T tanh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the hyperbolic tangent of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\tan f)} = \frac{\tanh f_0 + \sum_{k=1}^{k \le 2k+1} B_{2k} \frac{4^k(4^k-1)}{2k!}x^{2k - 1}}{1 + \tanh f_0
..       \sum_{k=1}^{k \le 2k+1} \frac{B_{2k}4^k(4^k-1)}{2k!}x^{2k - 1} }
..
..    where :math:`T_f = f_0 + \hat f` and :math:`B_{2k}\f` are the Bernoulli numbers.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the hyperebolic tangent of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T atanh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse hyperbolic tangent of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{atanh} f)} =  \mbox{atanh} f_0 +\frac 12 \sum_{k=1}^m \left(\frac{1}{(1-f_0)^k} +
..       \frac{(-1)^{k+1}}{(1+f_0)^k}\right) \frac {\hat f^k}{k}
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse hyperebolic tangent of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T atan(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse tangent of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{atan} f)} =  \mbox{atan} f_0 + \sum_{k=1}^{2k-1\le m} \left(\frac{1 + \sum_{j=1}^{2j\le 2k-1} {2k-1 \choose
..       2j} f_0^{2j}(-1)^j}{(1+f_0^2)^{2k-1}}\right) \frac {\hat f^{2k-1}}{2k-1}(-1)^{k+1} + \\ + \sum_{k=1}^{2k\le m}
..       \left(\frac{\sum_{j=1}^{2j-1\le 2k} {2k \choose 2j-1} f_0^{2j-1}(-1)^{j+1}}{(1+f_0^2)^{2k}}\right) \frac {\hat
..       f^{2k}}{2k}(-1)^k 
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    This formula derives directly from the formula for :cpp:func:`atanh` noting that: :math:`\mbox{atan}(z) = i \mbox{atanh}(-iz)`
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse tangent of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T asinh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse hyperbolic sine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{asinh} f)} = T_{\left(\log\left(f + \sqrt{1 + f^2}\right)\right)}
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse hyperebolic sine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T acosh(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse hyperbolic cosine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{acosh} f)} = T_{\left(\log\left(f + \sqrt{f^2 - 1}\right)\right)}
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse hyperebolic cosine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T asin(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse sine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{asin} f)} = T_{\left(\mbox{atan} \left(f / \sqrt{1 - f^2}\right)\right)}
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse sine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T acos(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the inverse cosine of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{acos} f)} = T_{\left(\mbox{atan} \left(\sqrt{1 - f^2} / f\right)\right)}
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the inverse sine of **d**
..
.. -----------------------------------------------------------------------
..
.. .. cpp:function:: template <typename T> inline T abs(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`.
..    Implements the absolute value of a :cpp:class:`gdual`.
..    Essentially it performs the following computations in the :math:`\mathcal P_{n,m}`
..
..    .. math::
..
..       T_{(\mbox{abs} f)} = \left\{ \begin{array}{ll} T_f & f_0 \ge 0 \\ -T_f & f_0 < 0 \end{array} \right.
..
..    where :math:`T_f = f_0 + \hat f`.
..
..    .. note::
..    
..       If :math:`f_0` is zero, the right Taylor expansion will be returned rather than nans.
..
..    .. note::
..    
..       This operation is not availiable whtn **T** is std::complex.
..
..    :param d: :cpp:class:`gdual` argument
..
..    :return: the Taylor expansion of the absolute value of **d**
..
.. -----------------------------------------------------------------------
..
.. *#include <audi/functions_from_d.hpp>*
..
.. .. cpp:function:: template <typename T> inline T erf(const T &d)
..
..    This templated function is enabled only if **T** is a :cpp:class:`gdual`
..    Essentially, it makes use of the definition:
..
..    .. math::
..    
..       \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}}\exp(-x^2)
..    
..    where :math:`T_f = f_0 + \hat f`.
..    
..    :param d: :cpp:class:`gdual` argument
..    
..    :return: the Taylor expansion of the error function of **d**
..
