/*! \mainpage AuDi (Automated Differentiation)
\section a The algebra of truncated polynomials 
Here we introduce formally a basic algebraic structure over the set of truncated polynomials and
we show how such a structure allows to compute the partial derivatives of multivariate functions
up to arbitrary order.

\subsection b Formal definition 
Consider the set \f$ \mathcal P_{n,m}\f$ of all polynomials of order \f$\le m\f$ in \f$n\f$ variables and having
coefficients in \f$\mathbf K\f$. Such a set is an algebra over the field \f$\mathbf K\f$ if we introduce
the truncated multiplication as the standard polynomial multiplication truncated at order \f$m\f$.
When needed, we will indicate such a multiplication with the symbol \f$T_f \cdot T_g\f$.

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
\label{eq:div}
T_g = \frac 1{T_f} =  \frac 1{f_0} \left(1 +\sum_{k=1}^m (-1)^k (\hat f / f_0)^k\right)
\f]
as its easily verified by performing the truncated multiplication \f$T_g \cdot T_f\f$ and
accounting for the nilpotency of \f$\hat f\f$. 
*/
