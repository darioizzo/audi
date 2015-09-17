.. contents::

The algebra of truncated polynomials 
====================================
Here we introduce, formally, a basic algebraic structure over the set of truncated polynomials and we show how such a structure allows to compute the partial derivatives of multivariate functions up to arbitrary order.

Formal definition 
-----------------
Consider the set :math:`\mathcal P_{n,m}` of all polynomials of order :math:`\le m` in \:math:`n` variables and having
coefficients in :math:`\mathbf K`. We indicate with the symbols :math:`T_f, T_g, T_h`, etc. the generic members of such a set. Such a set is an algebra over the field :math:`\mathbf K` if we introduce
the truncated multiplication as the standard polynomial multiplication truncated at order :math:`m`.
When needed, we will indicate such a multiplication with the symbol :math:`T_f \cdot T_g`.

This algebra is commonly referred to as the algebra of truncated polynomials. A first important
property of this algebra is that, under the multiplication, all polynomials having a zero constant
coefficient are nil-potent of order :math:`m+1`, as easily verified. We will indicate the generic
truncated polynomial :math:`\in \mathcal P_{n,m}` as :math:`T_f` and often we will consider its constant part
separated from the rest, writing :math:`T_f = f_0 + \hat f`.
It is worth noting at this point how such an algebra is unitary and associative.
The first property, in particular, deserves a few more words as it is a property that the
algebra of (non-truncated) polynomials does not possess. Formally
:math:`\forall T_f \in \mathcal P_{n,m}, \: \exists !\: T_g\in \mathcal P_{n,m}  \Rightarrow T_g\cdot T_f = 1`.
In practice:

.. math::
   T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
   :label: reciprocal

as its easily verified by performing the truncated multiplication :math:`T_g \cdot T_f` and accounting for the nilpotency of :math:`\hat f`. 

The link to Taylor expansions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We make use of the multi-index notation according to which

:math:`\alpha = (\alpha_1, ..,\alpha_n)` and :math:`\mathbf x = (x_1, .., x_n)`
are n-tuples and the Taylor expansion around the point :math:`\mathbf a` to order
:math:`m` of a multivariate function :math:`f` of :math:`\mathbf x` is written as:

.. math::
   T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)

where:

.. math::
   \partial^\alpha f = \frac{\partial^{|\alpha|} f}{\partial^{\alpha_1} x_1\partial^{\alpha_2} x_2\dots\partial^{\alpha_n} x_n}

.. math::
   \alpha ! = \prod_{i=j}^n \alpha_j

.. math::
   |\alpha| = \sum_{j=0}^n \alpha_j

The summation :math:`\sum_{|\alpha| = 0}^n` must then be taken over all possible
combinations of :math:`\alpha_j \in N` such that :math:`\sum_{j=1}^n \alpha_j = |\alpha|`. 
The expression above, i.e. the Taylor expansion truncated
at order :math:`m` of a generic function :math:`f`, is a polynomial 
:math:`\in \mathcal P_{n,m}` in the variables :math:`\mathbf{dx} = \mathbf x-\mathbf a`.
We now show that if :math:`T_f, T_g \in \mathcal P_{n,m}` are Taylor expansions
of two functions :math:`f, g` then the Taylor expansion of :math:`f\pm g, fg, f/g`
can be found operating on the algebra :math:`\mathcal P`, thus computing 
:math:`T_f\pm T_g, T_f\cdot T_g, T_f/T_g`. We may thus compute high order
derivatives of multivariate functions computing their Taylor expansions
and then extracting the desired coefficient.

Multiplication
^^^^^^^^^^^^^^

We here prove that the product of two truncated Taylor expansions is the
truncated Taylor expansion of the product. We perform the proof for
:math:`n=1` as the notation is there clearer. The multivariate case is
formally identical, requiring a more complex notation. The truncated
Taylor expansion of the product between two functions :math:`f` and :math:`g` is:

.. math::
   T_{(fg)} = \sum_{k=0}^m \frac{(x-a)^k}{k!}(fg)^{(k)}

where we indicate with the superscript :math:`(i)` the :math:`i`-th derivative with
respect to the independent variable.
We show how the same expression is derived by multiplying the Taylor
expansions of :math:`f` and :math:`g`:

.. math::
   T_f T_g = \sum_{k=0}^m \frac{(x-a)^k}{k!}(f)^{(k)}\sum_{k=0}^m \frac{(x-a)^k}{k!}(g)^{(k)} = \sum_{k=0}^m c_k (x-a)^k

The coefficients :math:`c_k` in the last power series are determined as the
Cauchy product of the two Taylor series (or discrete convolution) and are:

.. math::
   c_k = \sum_{n=0}^k \frac{f^{(n)}}{n!}  \frac{g^{(k-n)}}{(k-n)!}  = \frac{1}{k!}\sum_{n=0}^k {{k}\choose{n}}(f)^{(n)}(g)^{(k-n)}

applying now the general Leibniz rule to the last expression we get:

.. math::
   c_k =  \frac{1}{k!} (fg)^{(k)}

which allows us to conclude:

.. math::
   T_{(fg)} = T_f \cdot T_g.

Reciprocal
^^^^^^^^^^

We here prove that the reciprocal of a truncated Taylor expansion,
as defined in the algebra :math:`\mathcal P_{n,m}` is the Taylor expansion of
the reciprocal. Consider the generic
function :math:`f` and its truncated Taylor expansion :math:`T_f`.
We denote with :math:`T_{(1/f)}` the truncated Taylor expansion of the
reciprocal and apply the multiplication rule to derive that, necessarily,
:math:`T_f  T_{(1/f)} = 1`. We separate the constant part of :math:`T_f` from the
rest writing :math:`T_f = f_0 +\hat f` and we compute the product between :math:`T_f`
and the definition of reciprocal:

.. math::
   \left(f_0 + \hat f\right)\frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)= \frac 1f_0 \left(f_0 + \hat f\right)\left(1 - \frac{\hat f}{f_0} + \frac{\hat f^2}{f_0^2} - ... \right) = 1

which allows us to conclude:

.. math::
   T_{(1/f)} = \frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)

Elementary functions
====================

Consider the MacLaurin expansion of a generic function :math:`g(x) = \sum g_n x^n`. 
Consider now a multivariate function :math:`\hat f(\mathbf x) = \sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha` whose MacLaurin Taylor expansion does not have a constant term. The composition between these two functions will then be, trivially, 
:math:`(g \circ \hat f) (x) = \sum g_n (\sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha)^n`. 
If we now truncate such an expansion to order :math:`m`, we get 
:math:`T_{g\circ f}= \sum_{n=0}^m g_n (\sum_{|\alpha|=1}^m f_\alpha \mathbf x^\alpha)^n`, 
which can be written as:

.. math::
   T_{g\circ \hat f} = T_g\circ T_{\hat f}

The above equation is called the **composition rule** and is only valid for functions whose Taylor expansion 
does not have a constant term and, is thus nil-potent of order 
:math:`m+1` in :math:`\mathcal P_{n,m}`. In  general, we cannot compute the truncated 
Taylor expansion of a composition function directly composing the truncated 
Taylor expansions. For most elementary functions, though, we can consider 
:math:`T_f = f_0 + \hat f` and use some addition formula to be able to 
''extract`` :math:`\hat f` and thus exploit its nil-potency. The details 
on how this is done differ for each particular :math:`f` considered and are thus 
reported in the following subsections.

Exponential
-----------
Let us consider the case of the exponential:

.. math::
   g(x) = \exp(x) = \sum_{i=0} \frac{x^i}{i!} = 1 + x + \frac {x^2}{2} + ...

We want to compute the truncated Taylor expansion of \f$\exp(f(\mathbf x))\f$ 
starting from the truncated Taylor expansion \f$T_f = f_0 + \hat f\f$. 
We thus write:

.. math::
   (g \circ f) (\mathbf x) = \exp(f(\mathbf x)) =  \exp f_0 \exp (f(\mathbf x) - f_0)

note that, now, we can apply the **composition rule** to \f$\exp (f(\mathbf x) - f_0)\f$ 
since the MacLaurin Taylor expansion of \f$f(\mathbf x) - f_0\f$ does not have 
a constant term. Hence:

.. math::
   T_{g \circ f} = \exp f_0 T_g \circ T_{\hat f}
\f]
and, finally:

.. math::
   T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ... \right)
   :label: exp

Logarithm
---------
Let us consider the case of the natural logarithm:

.. math::
   g(x) = \log(x)

We want to compute the truncated Taylor expansion of 
\f$\log(f(\mathbf x))\f$ starting from the truncated Taylor expansion 
\f$T_f = f_0 + \hat f\f$. We thus write:

.. math::
   (g \circ f) (\mathbf x) = \log(f(\mathbf x)) =  \log (f_0 + (f(\mathbf x) - f_0)) = \log f_0 + \log(1 + \frac{f(\mathbf x) - f_0}{f_0})

We can now apply the **composition rule** to get:

.. math::
   T_{g \circ f} = \log f_0 + T_{\log(1+x)} \circ \frac{\hat f}{f_0}

and, using the known expression for MacLaurin expansion of \f$\log(1+x)\f$, we get:

.. math::
   T_{(\log f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} - \frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ...
   :label: log

Note that the above expression is only defined if \f$f_0 \ge 0\f$.

Sine and cosine
---------------
Let us consider the case of the sine and cosine functions:

.. math::
   g_1(x) = \sin(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i+1}}{(2i+1)!} = x - \frac{x^3}{3!} + \frac{x^5}{5!} - ... 

.. math::
   g_2(x) = \cos(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i}}{(2i)!} = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - ... 


We want to compute the truncated Taylor expansion of \f$\sin(f(\mathbf x))\f$, 
\f$\cos(f(\mathbf x))\f$ starting from the truncated Taylor expansion 
\f$T_f = f_0 + \hat f\f$. We thus write:

.. math::
   (g_1 \circ f) (\mathbf x) = \sin(f(\mathbf x)) =  \sin f_0 \cos(f(\mathbf x) - f_0) + \cos f_0 \sin(f(\mathbf x) - f_0) 

.. math::
   (g_2 \circ f) (\mathbf x) = \cos(f(\mathbf x)) =  \cos f_0 \cos(f(\mathbf x) - f_0) - \sin f_0 \sin(f(\mathbf x) - f_0) 

and, applying the **composition rule** to \f$\cos(f(\mathbf x) - f_0)\f$ and
\f$\sin(f(\mathbf x) - f_0)\f$, we get:

.. math::
   T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\
   T_{(\cos f)} = \cos f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right)
   :label: sinandcos


Exponentiation
--------------
Let us consider the case of the power function. 

.. math::
   g(x) = x^\alpha

We want to compute the truncated Taylor expansion of \f$f(\mathbf x)^\alpha\f$ 
assuming to have access to the truncated Taylor expansion of \f$f\f$, 
\f$T_f = f_0 + \hat f\f$. We thus write:

.. math::
   (g \circ f) (\mathbf x) = f(\mathbf x) ^ \alpha =  (f_0 + (f(\mathbf x) - f_0))^\alpha = f_0^\alpha \left( 1+ \frac{f(x) - f_0}{f_0}\right)^\alpha

We can now apply the **composition rule** to get:

.. math::
   T_{f(\mathbf x)^\alpha} = f_0^\alpha \left(T_{(1+x)^\alpha}\circ \frac{\hat f}{f_0}\right) = 

.. math::
   = f_0^\alpha \sum_{k=0}^m {\alpha \choose k} \left(\frac{\hat f}{f_0}\right)^k = f_0^\alpha\left(1 + \alpha \frac{\hat f}{f_0} + \frac{\alpha (\alpha - 1)}{2}\left(\frac{\hat f}{f_0}\right)^2 + ... \right)
   :label: pow
