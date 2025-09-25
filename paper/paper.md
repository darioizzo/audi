---
title: "pyaudi: A truncated Taylor polynomial algebra toolbox for differentiable intelligence, automatic differentiation, and verified integration applications."
tags:
  - Python
  - differentiable intelligence
  - automatic differentiation
  - verified integration
authors:
  - name: Dario Izzo
    affiliation: "1"
    orcid: 0000-0002-9846-8423
  - name: Francesco Biscani
    affiliation: "2"
  - name: Sean Cowan
    affiliation: "1"
    orcid: 0009-0000-9912-9354
affiliations:
  - index: 1
    name: Advanced Concepts Team, European Space Research and Technology Center (Noordwijk, NL)
  - index: 2
    name: Max Planck Institute for Astronomy (Heidelberg, DE)
date: 25 September 2025
bibliography: paper.bib
---

# Summary

<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

`pyaudi` is a Python toolbox developed at the [European Space Agency](https://www.esa.int) that implements the algebra of truncated Taylor
polynomials to achieve high-order order, forward mode, automatic differentiation in a multivarate setting. This form of forward mode automatic
differentiation is implemented via C++ class templates exposed to python using pybind11. This allows the generalized dual number type to
behave like a drop-in replacement for floats (or other scalar types), while operator overloading propagates derivatives automatically.

All standard mathematical functions are implemented exploiting the nil-potency property of exponentiation in the algebra of truncated Taylor polynomials.

On top of the algebra of truncated Taylor polynomials, `pyaudi` also offers an implementation of Taylor models [@makino1998rigorous], which combines truncated Taylor
polynomials with an interval bounding their truncation error as well as a number of miscellaneous algorithms useful for applications in differential intelligence,
high-order automatic differentiation, verified integration and more.

# Statement of need

<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. -->

`pyaudi` enables users to compute and manipulate order $n$ Taylor expansions of generic computational
graphs as well as bound precisely the truncation error introduced using
its corresponding Taylor model. The resulting representations
of program outputs can be used to perform fast Monte Carlo simulations, rigorous uncertainty analyses,
local inversions of output-input relations as well as high-order sensitivity
analysis. The package implements the approach to high-order automated differentiation perfectioned by Berz and Makino ([@berz2014introduction],  [@makino1998rigorous])
while introducing original implementation details aimed at increased efficiency in the
polynomial multiplication routines and bounding of Taylor models.

The C++/C library DACE [@massari2018differential] also offers an implementation of the differential algebra of truncated Taylor
polynomials, (same as `pyaudi` but with the addition of extra operators completing the algebra into a differential algebra). As opposed to `pyaudi`, DACE
makes use of a polynomial multiplication routine relying on static memory allocations for the storage of monomial coefficients. As outlined 
in the comparison reported , the difference makes DACE more efficient for single evaluations, an advantage that is lost when evaluating in batches.

The Julia library TaylorSeries.jl/TaylorModels.jl [@benet2019taylormodels] implements Taylor models to compute guaranteed bounds on
generic Taylor series, the approach used is very different from what implemented in `pyaudi` and a preliminary comparison, reported below, shows
how it can be greatly outpermed by `pyaudi`.



## Key aspects  

The main features of `pyaudi` are:  

- Efficient truncated polynomial arithmetic in arbitrary dimensions, built on top of [Obake](https://github.com/bluescarni/obake),
  a C++ library for symbolic manipulation of sparse multivariate polynomials, truncated power series, and Poisson series. 
  Unlike other packages, which often face severe memory bottlenecks as the polynomial order or number of variables grows,
  `pyaudi` avoids large static memory allocations and keeps computations memory-efficient, at the cost of some extra bookeeping.

- Vectorized generalized dual numbers, enabling simultaneous evaluation of identical computational graphs at
  multiple expansion points. This makes it possible to compute high-order tensors efficiently while amortizing
  the overhead of the extra bookeeping introduced by the use of [Obake](https://github.com/bluescarni/obake).

- Taylor models implemented using Bernstein polynomials for bounding the range of multivariate polynomials. A comparable open-source package, 
  called TaylorModels.jl, calculates bounds using Horner's scheme combined with interval arithmetic. A
  quick test in the next section shows the significant speedup introduced using `pyaudi` for even a relatively simple trivariate polynomial.

- Map inversion algorithm, implementing the algorithm described in [@berz2014introduction], thus allowing
  local inversion of input–output relations of generic computational graphs.

## Comparison with DACE

## Comparison with DACE

For the comparison with DACE, we take two cases: One with simple multiplication between Taylor polynomials. One with vectorization enabled (in audi since it does not exist in DACE). 

### Multiplication

We use two tenth-order polynomials of the form:

$$
\begin{array}{l}
\begin{aligned}
     p1 &= (1 + x1 + x2 + ... + xn)^{10} \\
     p2 &= (1 - x1 - x2 - ... - xn)^{10} \\
\end{aligned}
\\
\end{array}
$$

where $x1$ etc. are variables. These polynomials are then multiplied and timed. The speed up of pyaudi w.r.t. DACE is given below in seconds. It can be
seen that pyaudi is faster from nvars >= 10 and order >= 11. 

| nvars↓       Order→  |   6    |   7    |   8    |   9    |  10   |  11   |  12   |  13   |  14   |  15   |
|----------------------|--------|--------|--------|--------|-------|-------|-------|-------|-------|-------|
| 6                    | 0.0901 | 0.115  | 0.0978 | 0.125  | 0.166 | 0.264 | 0.369 | 0.426 | 0.966 | 1.04  |
| 8                    | 0.0791 | 0.226  | 0.166  | 0.475  | 0.635 | 0.977 | 1.42  | 1.72  | 1.98  | 2.14  |
| 10                   | 0.114  | 0.177  | 0.631  | 0.692  | 0.866 | 1.36  | 1.43  | 1.12  | 1.70  | 3.08  |
| 12                   | 0.138  | 0.422  | 0.544  | 0.723  | 0.914 | 1.32  | 1.73  | 2.18  | 2.84  | 4.08  |
                       
                       
### Vectorization-enabled multiplication

To showcase the vectorization, we multiply the following polynomial (where each variable has a
variable number of coefficients) five times with itself and time the operation:

$$
\begin{array}{l}
    p1 = \frac{c_v + x1 + x2 + ... + xn}{c_v - x1 - x2 - ... - xn}^5 \\
\end{array}
$$

where $c_v$ are coefficients. The results are displayed in three tables
per number of variables below. It can be seen that, from
64/256 points onwards, pyaudi is significantly faster than DACE. It should be noted that this test
is only indicative due to the limited study using an arbitrary choice of a polynomial both for the
single-thread multiplication as well as the vectorized version.

#### 2 variables
| points↓       Order →  |   1    |   2    |   3    |   4    |   5    |   6    |   7    |   8    |   9    |
|------------------------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| 16                     | 0.122  | 0.349  | 0.387  | 0.395  | 0.553  | 0.607  | 0.535  | 0.494  | 0.527  |
| 64                     | 1.54   | 1.42   | 1.17   | 1.10   | 1.52   | 1.80   | 1.50   | 1.42   | 0.637  |
| 256                    | 5.64   | 4.50   | 3.77   | 1.25   | 1.61   | 2.14   | 1.81   | 1.73   | 1.66   |
| 1024                   | 12.7   | 3.30   | 2.68   | 2.04   | 2.41   | 3.11   | 2.77   | 2.70   | 2.73   |
| 4096                   | 3.62   | 3.58   | 3.28   | 3.43   | 3.80   | 5.30   | 4.73   | 4.76   | 4.35   |
| 16384                  | 5.62   | 4.44   | 3.48   | 2.98   | 3.68   | 3.91   | 3.74   | 3.61   | 3.33   |


#### 5 variables
| points↓       Order →  |   1    |   2    |   3    |   4    |   5    |   6    |   7    |   8    |   9    |
|------------------------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| 16                     | 0.339  | 0.579  | 0.437  | 0.478  | 0.219  | 0.365  | 0.323  | 0.374  | 0.443  |
| 64                     | 1.42   | 1.87   | 0.591  | 0.470  | 0.507  | 0.799  | 0.703  | 0.921  | 1.11   |
| 256                    | 5.27   | 2.24   | 1.31   | 1.19   | 1.26   | 2.24   | 2.27   | 2.07   | 2.60   |
| 1024                   | 5.23   | 4.54   | 2.79   | 2.31   | 2.44   | 3.28   | 3.19   | 2.85   | 3.14   |
| 4096                   | 8.10   | 7.26   | 4.73   | 2.61   | 2.50   | 3.69   | 3.38   | 3.53   | 3.50   |
| 16384                  | 7.20   | 4.65   | 3.01   | 2.32   | 2.06   | 3.34   | 3.33   | 3.46   | 3.54   |


#### 10 variables
| points↓       Order →  |   1    |   2    |   3    |   4    |   5    |   6    |   7    |   8    |   9    |
|------------------------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| 16                     | 0.287  | 0.356  | 0.158  | 0.144  | 0.241  | 0.539  | 0.827  | 1.04   | 1.50   |
| 64                     | 1.19   | 0.461  | 0.394  | 0.387  | 0.594  | 1.20   | 1.55   | 1.97   | 2.89   |
| 256                    | 4.26   | 1.50   | 1.14   | 1.14   | 1.45   | 2.73   | 3.46   | 4.44   | 5.95   |
| 1024                   | 4.04   | 3.41   | 2.24   | 1.80   | 2.10   | 3.92   | 4.96   | 6.45   | 9.90   |
| 4096                   | 8.83   | 5.38   | 2.52   | 2.14   | 2.53   | 5.05   | 5.77   | 7.41   | 12.7   |
| 16384                  | 3.81   | 2.33   | 2.04   | 2.06   | 2.61   | 4.57   | 5.13   | 6.78   | 10.1   |


## Comparison with TaylorModels.jl

We test the performance of the implementation of Taylor models in `pyaudi` against the Julia package
TaylorModels.jl. To perform the comparison we use three functions $f,g,h$: one univariate, one bivariate and one trivariate defined
below. We then construct Taylor models of all the variables separately and time the evaluation of
the corresponding Taylor model. The comparison is made on a single CPU machine.

$$
\begin{array}{l}
\begin{aligned}
f(x, y, z) = {} & 
\frac{4 \tan(3y)}{3x + x \sqrt{\tfrac{6x}{-7(x-8)}}}
- 120 - 2x - 7z(1+2y) \\
& - \sinh\!\left(0.5 + \frac{6y}{8y+7}\right)
+ \frac{(3y+13)^2}{3z}
- 20z(2z-5) \\
& + \frac{5x \tanh(0.9z)}{\sqrt{5y}}
- 20y \sin(3z)
\end{aligned}
\\
g(x, y) = \sin(1.7x+0.5)(y+2)\sin(1.5y)
\\
h(x) = x(x-1.1)(x+2)(x+2.2)(x+2.5)(x+3)\sin(1.7x+0.5)
\end{array}
$$

| Dimension | Package         | Remainder Bound (Order 1) | Remainder Bound (Order 15) | Speed Comparison                  |
|-------------|-----------------|--------------------|--------------------|---------------------------------------|
| h(x)      | TaylorModels.jl | 1e+2                      | 1e-5                       | ~1–1.5× faster than pyaudi |
|           | pyaudi          | 1e+2                      | 1e-5                       | ~1–1.5× slower than TaylorModels.jl |
| g(x, y)   | TaylorModels.jl | 1e+1                      | 1e-6                       | Slower: pyaudi is 5× faster (order 3), 15× faster (order 15), 7800× faster (order 1, edge case) |
|           | pyaudi          | 1e+1                      | 1e-6                       | Faster (see above) |
| f(x, y, z)| TaylorModels.jl | 1e+0                      | 1e-11                      | Slower: pyaudi is 8× faster (order 3), 155× faster (order 15), 13000× faster (order 1, edge case) |
|           | pyaudi          | 1e-1                      | 1e-17                      | Faster (see above) |

In the table above, a clear trend can be seen both in terms of speed and accuracy. For univariate
Taylor models, TaylorModels.jl and `pyaudi` have similar performances. At two dimensions, while the
remainder bounds are comparable in size, `pyaudi` is significantly faster, with the speedup
increasing with the order of the polynomial. At three dimensions, `pyaudi` produces significantly
tighter bounds and is again significantly faster, with the speedup increasing with the order of
the polynomial.

# References

<!-- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline. -->

A number of references to relevant work and algorithms implemented in `pyaudi` are:

- [@obake2020biscani]
- [@makino1998rigorous]
- [@titi2019matrix]

Other software packages that do similar things are:

- JAX [@jax2018github]
- TensorFlow [@tensorflow2015-whitepaper]
- PyTorch [@paszke2019pytorch]
- COSY INFINITY [@makino2006cosy]
- DACE [@massari2018differential]
- TaylorSeries.jl/TaylorModels.jl [@benet2019taylormodels]
- CORA [@Althoff2015ARCH]

# Ongoing research

<!-- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. -->

- EclipseNET [@acciarini2024eclipsenets] [@acciarini2025eclipsenets]
- CR3BP stochastic continuation [@acciarini2024stochastic]
- Long-term propagation [@caleb2020can]
- Rapid nonlinear convex guidance [@burnett2025rapid]
- Differentiable genetic programming [@izzo2017differentiable]

# Acknowledgement of financial support

No financial support was provided for the development of this software.

# Bibliography
