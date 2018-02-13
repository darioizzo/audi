/*! \mainpage
\section z Listen (Introduction)
Audi (not the car, rather from latin: "listen!") is an open source, header only, C++ library that allows
for high order AUtomated DIfferentiation implementing the Taylor truncated polynomial
algebra (forward mode automated differentiation). It was created with the aim
to offer a solution to the user in need of a generic automated
differentiation system. While other automated differentiation codes may be more
efficient on some targeted application, AuDi was built to be fast and efficient
across all application ranges (low/high order of derivation, one variable,
many variables). Its use is straight forward and is illustarted
in the getting started example below.

\include getting_started.cpp

AuDi is internally powered by the truncated polynomial multiplication algorithm
of the open source project,  <a href="https://github.com/bluescarni/piranha">Piranha</a>
details of which are described in:

Biscani, Francesco. <a href="http://dl.acm.org/citation.cfm?id=2442845"> "Parallel sparse polynomial multiplication on
modern hardware architectures."</a> Proceedings of the 37th International Symposium on Symbolic and Algebraic
Computation. ACM, 2012.

Biscani, Francesco. <a href="http://arxiv.org/pdf/1004.4548v1.pdf">"Multiplication of sparse Laurent polynomials and
Poisson series on modern hardware architectures."</a> arXiv preprint arXiv:1004.4548 (2010).

<b>Piranha headers must be installed and accessible to the compiler for AuDi to function.</b>


In the following sections we present, briefly, the theory and some concrete examples
to help understanding the forward automated differentiation using
truncated Taylor polynomials implemented in AuDi.

\section a The algebra of truncated polynomials
Here we introduce formally a basic algebraic structure over the set of truncated polynomials and
we show how such a structure allows to compute the partial derivatives of multivariate functions
up to arbitrary order.

\subsection b Formal definition
Consider the set \f$ \mathcal P_{n,m}\f$ of all polynomials of order \f$\le m\f$ in \f$n\f$ variables and having
coefficients in \f$\mathbf K\f$. We indicate with the symbols \f$T_f, T_g, T_h\f$, etc. the generic members of such a
set. Such a set is an algebra over the field \f$\mathbf K\f$ if we introduce the truncated multiplication as the
standard polynomial multiplication truncated at order \f$m\f$. When needed, we will indicate such a multiplication with
the symbol \f$T_f \cdot T_g\f$.

This algebra is commonly referred to as the algebra of truncated polynomials. A first important
property of this algebra is that, under the multiplication, all polynomials having a zero constant
coefficient are nil-potent of order \f$m+1\f$, as easily verified. We will indicate the generic
truncated polynomial \f$\in \mathcal P_{n,m}\f$ as \f$T_f\f$ and often we will consider its constant part
separated from the rest, writing \f$T_f = f_0 + \hat f\f$.
It is worth noting at this point how such an algebra is unitary and associative.
The first property, in particular, deserves a few more words as it is a property that the
algebra of (non-truncated) polynomials does not possess. Formally
\f$\forall T_f \in \mathcal P_{n,m}, \: \exists !\: T_g\in \mathcal P_{n,m}  \Rightarrow T_g\cdot T_f = 1\f$.
In practice:
\f[
T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
\f]
as its easily verified by performing the truncated multiplication \f$T_g \cdot T_f\f$ and
accounting for the nilpotency of \f$\hat f\f$.

\subsection c The link to Taylor expansions
We make use of the multi-index notation according to which
\f$\alpha = (\alpha_1, ..,\alpha_n)\f$ and \f$\mathbf x = (x_1, .., x_n)\f$
are n-tuples and the Taylor expansion around the point \f$\mathbf a\f$ to order
\f$m\f$ of a multivariate function \f$f\f$ of \f$\mathbf x\f$ is written as:
\f[
T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)
\f]
where:
\f[
\partial^\alpha f = \frac{\partial^{|\alpha|} f}{\partial^{\alpha_1} x_1\partial^{\alpha_2} x_2\dots\partial^{\alpha_n}
x_n} \f] \f[ \alpha ! = \prod_{i=j}^n \alpha_j \f] and \f[
|\alpha| = \sum_{j=0}^n \alpha_j
\f]
The summation \f$ \sum_{|\alpha| = 0}^n\f$ must then be taken over all possible
combinations of \f$\alpha_j \in N\f$ such that \f$\sum_{j=1}^n \alpha_j = |\alpha|\f$.
The expression above, i.e. the Taylor expansion truncated
at order \f$m\f$ of a generic function \f$f\f$, is a polynomial
\f$\in \mathcal P_{n,m}\f$ in the variables \f$\mathbf{dx} = \mathbf x-\mathbf a\f$.
We now show that if \f$T_f, T_g \in \mathcal P_{n,m}\f$ are Taylor expansions
of two functions \f$f, g\f$ then the Taylor expansion of \f$f\pm g, fg, f/g\f$
can be found operating on the algebra \f$\mathcal P\f$, thus computing
\f$T_f\pm T_g, T_f\cdot T_g, T_f/T_g\f$. We may thus compute high order
derivatives of multivariate functions computing their Taylor expansions
and then extracting the desired coefficient.


\subsection d Multiplication
We here prove that the product of two truncated Taylor expansions is the
truncated Taylor expansion of the product. We perform the proof for
\f$n=1\f$ as the notation is there clearer. The multivariate case is
formally identical, requiring a more complex notation. The truncated
Taylor expansion of the product between two functions \f$f\f$ and \f$g\f$ is:
\f[
T_{(fg)} = \sum_{k=0}^m \frac{(x-a)^k}{k!}(fg)^{(k)}
\f]
where we indicate with the superscript \f$(i)\f$ the \f$i\f$-th derivative with
respect to the independent variable.
We show how the same expression is derived by multiplying the Taylor
 expansions of \f$f\f$ and \f$g\f$:
\f[
T_f T_g = \sum_{k=0}^m \frac{(x-a)^k}{k!}(f)^{(k)}\sum_{k=0}^m \frac{(x-a)^k}{k!}(g)^{(k)}
= \sum_{k=0}^m c_k (x-a)^k
\f]
The coefficients \f$c_k\f$ in the last power series are determined as the
Cauchy product of the two Taylor series (or discrete convolution) and are:
\f[
c_k = \sum_{n=0}^k \frac{f^{(n)}}{n!}  \frac{g^{(k-n)}}{(k-n)!}  = \frac{1}{k!}\sum_{n=0}^k
{{k}\choose{n}}(f)^{(n)}(g)^{(k-n)} \f] applying now the general Leibniz rule to the last expression we get: \f[ c_k =
\frac{1}{k!} (fg)^{(k)} \f] which allows us to conclude: \f[ T_{(fg)} = T_f \cdot T_g. \f]

\subsection e Reciprocal
We here prove that the reciprocal of a truncated Taylor expansion,
as defined in the algebra \f$\mathcal P_{n,m}\f$ is the Taylor expansion of
the reciprocal. Consider the generic
function \f$f\f$ and its truncated Taylor expansion \f$T_f\f$.
We denote with \f$T_{(1/f)}\f$ the truncated Taylor expansion of the
reciprocal and apply the multiplication rule to derive that, necessarily,
\f$T_f  T_{(1/f)} = 1\f$. We separate the constant part of \f$T_f\f$ from the
rest writing \f$T_f = f_0 +\hat f\f$ and we compute the product between \f$T_f\f$
and the definition of reciprocal:

\f[
\left(f_0 + \hat f\right)\frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)=
\f]

\f[
 = \frac 1f_0 \left(f_0 + \hat f\right)\left(1 - \frac{\hat f}{f_0} + \frac{\hat f^2}{f_0^2} - ... \right) = 1
\f]

which allows us to conclude:
\f[
T_{(1/f)} = \frac 1f_0 \left(1 +\sum_{j=1}^m (-1)^j (\hat f / f_0)^j\right)
\f]

\section f Elementary functions
Consider the MacLaurin expansion of a generic function \f$g(x) = \sum g_n x^n\f$.
Consider now a multivariate function \f$\hat f(\mathbf x) = \sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha\f$
whose MacLaurin Taylor expansion does not have a constant term. The
composition between these two functions will then be, trivially,
\f$(g \circ \hat f) (x) = \sum g_n (\sum_{|\alpha|=1} f_\alpha \mathbf x^\alpha)^n\f$.
If we now truncate such an expansion to order \f$m\f$, we get
\f$T_{g\circ f}= \sum_{n=0}^m g_n (\sum_{|\alpha|=1}^m f_\alpha \mathbf x^\alpha)^n \f$,
which can be written as:

\f[
T_{g\circ \hat f} = T_g\circ T_{\hat f}
\f]

The above equation is called the **composition rule** and is only valid for functions whose Taylor expansion
does not have a constant term and, is thus nil-potent of order
\f$m+1\f$ in \f$\mathcal P_{n,m}\f$. In  general, we cannot compute the truncated
Taylor expansion of a composition function directly composing the truncated
Taylor expansions. For most elementary functions, though, we can consider
\f$T_f = f_0 + \hat f\f$ and use some addition formula to be able to
''extract`` \f$\hat f\f$ and thus exploit its nil-potency. The details
on how this is done differ for each particular \f$g\f$ considered and are thus
reported in the following subsections.

\subsection g Exponential
Let us consider the case of the exponential:
\f[
g(x) = \exp(x) = \sum_{i=0} \frac{x^i}{i!} = 1 + x + \frac {x^2}{2} + ...
\f]
We want to compute the truncated Taylor expansion of \f$\exp(f(\mathbf x))\f$
starting from the truncated Taylor expansion \f$T_f = f_0 + \hat f\f$.
We thus write:
\f[
(g \circ f) (\mathbf x) = \exp(f(\mathbf x)) =  \exp f_0 \exp (f(\mathbf x) - f_0)
\f]
note that, now, we can apply the **composition rule** to \f$\exp (f(\mathbf x) - f_0)\f$
since the MacLaurin Taylor expansion of \f$f(\mathbf x) - f_0\f$ does not have
a constant term. Hence:
\f[
T_{g \circ f} = \exp f_0 T_g \circ T_{\hat f}
\f]
and, finally:
\f[
T_{(\exp f)} = \exp f_0 \sum_{i=0}^m \frac{\hat f^i}{i!} = \exp f_0 \left( 1 + \hat f + \frac {\hat f^2}{2!} + ...
\right) \f]

\subsection h Logarithm
Let us consider the case of the natural logarithm:
\f[
g(x) = \log(x)
\f]
We want to compute the truncated Taylor expansion of
\f$\log(f(\mathbf x))\f$ starting from the truncated Taylor expansion
\f$T_f = f_0 + \hat f\f$. We thus write:
\f[
(g \circ f) (\mathbf x) = \log(f(\mathbf x)) =  \log (f_0 + (f(\mathbf x) - f_0)) = \log f_0 + \log(1 + \frac{f(\mathbf
x) - f_0}{f_0}) \f] We can now apply the **composition rule** to get: \f[ T_{g \circ f} = \log f_0 + T_{\log(1+x)} \circ
\frac{\hat f}{f_0} \f] and, using the known expression for MacLaurin expansion of \f$\log(1+x)\f$, we get: \f[ T_{(\log
f)} = \log f_0 + \sum_{i=1}^m (-1)^{i+1} \frac 1i \left(\frac{\hat f}{f_0}\right)^i = \log f_0 + \frac{\hat f}{f_0} -
\frac 12 \left(\frac{\hat f}{f_0}\right)^2 + ... \f] Note that the above expression is only defined if \f$f_0 \ge 0\f$.

\subsection j Sine and cosine
Let us consider the case of the sine and cosine functions:
\f[
g_1(x) = \sin(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i+1}}{(2i+1)!} = x - \frac{x^3}{3!} + \frac{x^5}{5!} - ...
\f]
\f[
g_2(x) = \cos(x) = \sum_{i=0} (-1)^{i} \frac{x^{2i}}{(2i)!} = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - ...
\f]
We want to compute the truncated Taylor expansion of \f$\sin(f(\mathbf x))\f$,
\f$\cos(f(\mathbf x))\f$ starting from the truncated Taylor expansion
\f$T_f = f_0 + \hat f\f$. We thus write:
\f[
(g_1 \circ f) (\mathbf x) = \sin(f(\mathbf x)) =  \sin f_0 \cos(f(\mathbf x) - f_0) + \cos f_0 \sin(f(\mathbf x) - f_0)
\f]
\f[
(g_2 \circ f) (\mathbf x) = \cos(f(\mathbf x)) =  \cos f_0 \cos(f(\mathbf x) - f_0) - \sin f_0 \sin(f(\mathbf x) - f_0)
\f]
and, applying the **composition rule** to \f$\cos(f(\mathbf x) - f_0)\f$ and
\f$\sin(f(\mathbf x) - f_0)\f$, we get:
\f[
\begin{array}{l}
T_{(\sin f)} = \sin f_0 \left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) + \cos f_0
\left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i} \frac{\hat f^{2i+1}}{(2i+1)!}\right) \\ T_{(\cos f)} = \cos f_0
\left(\sum_{i=0}^{2i\le m} (-1)^{i} \frac{\hat f^{2i}}{(2i)!}\right) - \sin f_0 \left(\sum_{i=0}^{(2i+1)\le m} (-1)^{i}
\frac{\hat f^{2i+1}}{(2i+1)!}\right) \end{array} \f]

\subsection k Exponentiation
Let us consider the case of the power function.
\f[
g(x) = x^\alpha
\f]
We want to compute the truncated Taylor expansion of \f$f(\mathbf x)^\alpha\f$
assuming to have access to the truncated Taylor expansion of \f$f\f$,
\f$T_f = f_0 + \hat f\f$. We thus write:
\f[
(g \circ f) (\mathbf x) = f(\mathbf x) ^ \alpha =  (f_0 + (f(\mathbf x) - f_0))^\alpha = f_0^\alpha \left( 1+ \frac{f(x)
- f_0}{f_0}\right)^\alpha \f] We can now apply the **composition rule** to get: \f[ T_{f(\mathbf x)^\alpha} = f_0^\alpha
\left(T_{(1+x)^\alpha}\circ \frac{\hat f}{f_0}\right) = \f] \f[ = f_0^\alpha \sum_{k=0}^m {\alpha \choose k}
\left(\frac{\hat f}{f_0}\right)^k = f_0^\alpha\left(1 + \alpha \frac{\hat f}{f_0} + \frac{\alpha (\alpha -
1)}{2}\left(\frac{\hat f}{f_0}\right)^2 + ... \right) \f]

\section l Simple Examples (can be done by hand)
In the above sections we derived a number of results that allow operating
on simple Taylor expansions to compute Taylor expansions of increasingly
complex expressions. We summarize here those results (keep in mind that
\f$T_f = f_0 + \hat f\f$) :

\f[
\begin{array}{l}
T_{f\pm g} = T_f \pm T_g\\
T_{fg} = T_f \cdot T_g\\
T_{(1/f)} = \frac 1f_0 \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)\\
T_{(\exp f)} = \exp f_0 \sum_{k=0}^m \frac{\hat f^k}{k!} \\
T_{(\log f)} = \log f_0 - \sum_{k=1}^m \frac{(-1)^k}k \left(\hat f / f_0\right)^k \\
T_{(\sin f)} = \sin f_0 \left(\sum_{k=0}^{2k\le m} (-1)^{k} \frac{\hat f^{2k}}{(2k)!}\right) + \cos f_0
\left(\sum_{k=0}^{(2k+1)\le m} (-1)^k \frac{\hat f^{2k+1}}{(2k+1)!}\right) \\
T_{(\cos f)} = \cos f_0 \left(\sum_{k=0}^{2k\le m} (-1)^{k} \frac{\hat f^{2k}}{(2k)!}\right) - \sin f_0
\left(\sum_{k=0}^{(2k+1)\le m} (-1)^k \frac{\hat f^{2k+1}}{(2k+1)!}\right) \\ T_{(f^\alpha)} = f_0^\alpha \sum_{k=0}^m
{\alpha \choose k} \left(\hat f / f_0\right)^k \end{array} \f] It is worth mentioning here that other functions such as
the inverse functions, the hyperbolic functions etc. can also be treated in this way. The above equations can be used to
find Taylor expansions of increasingly complex functions by simply operating on the algebra \f$\mathcal P_{n,m}\f$. Once
a Taylor expansion is computed, its coefficients can be extracted to obtain the value of any desired derivative. We have
thus built an automated differentiation system. While the formalism presented can, at first, appear complex, the system
is rather simple as we hope will appear from the following examples.

\subsection cv Example 1 - A multiplication

Consider the simple function of two variables:
\f[
f(x,y) = x + 3xy + y^2
\f]
Its Taylor expansion \f$T_f \in \mathcal P_{2,2}\f$ can be computed as:
\f[
T_f = T_x + 3T_x \cdot T_y + T_y\cdot T_y
\f]
Let us explicitly compute such an expression at the point \f$x=3\f$, \f$y=7\f$. The exact sequence of computations to be
performed is: \f[ T_x = 3 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2 \f] \f[ T_y = 7 + 0 dx + 1 dy  + 0 dxdy + 0 dx^2 + 0
dy^2 \f] \f[ T_x \cdot T_y = 21 + 7 d x + 3 d y  + 1 dxdy + 0 dx^2 + 0 dy^2 \f] and \f[ T_y \cdot T_y = 49 + 0 dx + 14
dy  + 0 dxdy + 0 dx^2 + 1 dy^2 \f] We can then derive the final expression: \f[ T_f = 115 + 22 dx + 23 dy +3 dxdy + 0
dx^2 + 1 dy^2 \f] and we may easily extract the derivatives comparing this expression to the generic form of a Taylor
expansion: \f[ f = 115, \partial_x f = 22, \partial_y f = 23, \partial_{xy} f = 3, \partial_{xx} f = 0, \partial_{yy} f
= 2, \f]

\subsection x Example 2 - A division
Consider the simple function of two variables:
\f[
f = 1 / (x + 2xy + y^2) = 1 / p
\f]
Its Taylor expansion \f$T_f \in \mathcal P_{2,2}\f$ in (say) \f$x=0\f$, \f$y=1\f$
can be computed as follows:
\f[
T_x = 0 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2
\f]
\f[
T_y = 1 + 0 dx + 1 dy  + 0 dxdy + 0 dx^2 + 0 dy^2
\f]
\f[
T_p =  1 + 3 dx + 2 dy +2 dxdy + 0 dx^2 + 1 dy^2
\f]
and, applying the reciprocal rule, we conclude
\f[
T_f = ( 1 - \hat p + \hat p ^ 2 )
\f]
where \f$\hat p = 3 dx + 2 dy +2 dxdy + 0 dx^2 + 1 dy^2\f$, hence:
\f[
T_f = 1 -3 dx -2 dy + 10dxdy + 9dx^2 + 3dy^2
\f]
which allows, as in the previous example, to compute all derivatives up to order two:
\f[
f = 1,
\partial_x f = -3,
\partial_y f = -2,
\partial_{xy} f = 10,
\partial_{xx} f = 18,
\partial_{yy} f = 6,
\f]
\subsection cd Example 3 - An exponential

Consider the elementary function of two variables:
\f[
f = \exp(xy)
\f]
Its Taylor expansion \f$T_f \in \mathcal P_{2,2}\f$ in (say) \f$x=1\f$, \f$y=0\f$
can be computed as follows:
\f[
T_x = 1 + 1 dx + 0 dy  + 0 dxdy + 0 dx^2 + 0 dy^2
\f]
\f[
T_y = 0 + 0 d x + 1 dy  + 0 dxdy + 0 dx^2 + 0 dy^2
\f]
\f[
T_x \cdot T_y = 0 + 0 dx + 1 dy  + 1 dxdy + 0 dx^2 + 0dy^2
\f]
and, applying the rule for the exponential of Taylor series, we conclude:
\f[
T_f = 1 + dy  + dxdy + \frac 12 dy^2
\f]
and,
\f[
f = 1,
\partial_x f = 0,
\partial_y f = 1,
\partial_{xy} f = 1,
\partial_{xx} f = 0,
\partial_{yy} f = 1,
\f]
*/
