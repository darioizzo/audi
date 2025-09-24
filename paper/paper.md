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
date: 22 September 2025
bibliography: paper.bib
---

# Summary

<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

`pyaudi` is a Python toolbox developed at the [European Space Agency](https://www.esa.int) that implements the algebra of Taylor
truncated polynomials to achieve high-order order, forward mode, automatic differentiation in a multivarate setting. The forward mode automatic
differentiation is implemented via C++ class templates exposed to python using pybind11. This allows the generalized dual number type to
behave like a drop-in replacement for floats (or other scalar types), while operator overloading propagates derivatives automatically.

All standard mathematical functions are all implemented exploiting the nihilpotency property of truncated polynomials lacking the constant term.

On top of the algebra of Taylor truncated polynomials, `pyaudi` also offers an implementation of Taylor models [@makino1998rigorous], which combine truncated Taylor
polynomials with an interval bounding its truncation error as well as a number of miscellaneous algorithms useful for applications in differential intelligence,
automatic differentiation, verified integration and more.



# Statement of need

<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. -->

`pyaudi` enables researchers to compute and manipulate order $n$ Taylor expansions of generic computational
trees as well as bound precisely the truncation error introduced using
its corresponding Taylor model. The resulting polynomial representations
of the program outputs can be used to perform fast Monte Carlo simulations, rigorous uncertainty analyses
local inversions of output-input relations as well as high-order sensitivity
analysis. The package implements the approach to high-order automated differentiation by [@berz2014introduction] and [@makino1998rigorous]
introducing original implementation details aimed at increased efficiency in the
polynomial multiplication routines and bounding of Taylor models.

## Key aspects  

The main features of `pyaudi` are:  

- Efficient truncated polynomial arithmetic in arbitrary dimensions, built on top of [Obake](https://github.com/bluescarni/obake),
  a C++ library for symbolic manipulation of sparse multivariate polynomials, truncated power series, and Poisson series. 
  Unlike other packages, which often face severe memory bottlenecks as the polynomial order or number of variables grows,
  `pyaudi` avoids large static memory allocations and keeps computations memory-efficient.  

- Vectorized generalized dual numbers, enabling simultaneous evaluation of identical computational graphs at
  multiple expansion points. This makes it possible to compute high-order tensors efficiently while amortizing
  the overhead of graph bookkeeping.  

- Taylor models implemented using Bernstein polynomials for bounding the range of multivariate polynomials. The only other comparable open-source package, 
  called TaylorModels.jl, calculates bounds using Horner's scheme combined with interval arithmetic. A
  quick test in the next section shows significant speedup for even a relatively simple trivariate polynomial.

- Map inversion algorithm, implementing the algorithm described in [@berz2014introduction], thus allowing
  local inversion of the input–output relation of generic computational graphs.

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
$$

| Dimension | Package         | Remainder Bound (Order 1) | Remainder Bound (Order 15) | Speed Comparison                  |
|-------------|-----------------|--------------------|--------------------|---------------------------------------|
| h(x)      | TaylorModels.jl | 1e-15                     | 1e-15                      | ~1.5–2× faster than pyaudi |
|           | pyaudi          | 1e+2                      | 1e-5                       | ~1.5–2× slower than TaylorModels.jl |
| g(x, y)   | TaylorModels.jl | 1e+1                      | 1e-6                       | Slower: pyaudi is 5× faster (order 3), 15× faster (order 15), 7800× faster (order 1, edge case) |
|           | pyaudi          | 1e+1                      | 1e-6                       | Faster (see above) |
| f(x, y, z)| TaylorModels.jl | 1e+0                      | 1e-11                      | Slower: pyaudi is 8× faster (order 3), 155× faster (order 15), 13000× faster (order 1, edge case) |
|           | pyaudi          | 1e-1                      | 1e-17                      | Faster (see above) |

In the table above, a clear trend can be seen both in terms of speed and accuracy. For univariate
Taylor models, TaylorModels.jl is both faster and produces tighter bounds. At two dimensions, the
remainder bounds are already of equal size, but `pyaudi` is significantly faster, with the speedup
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
