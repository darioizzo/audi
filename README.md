[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22009085.svg)](https://doi.org/10.5281/zenodo.22009085)
[![Build Status](https://travis-ci.org/darioizzo/audi.svg?branch=master)](https://travis-ci.org/darioizzo/audi)
[![PyPI](https://img.shields.io/pypi/v/pyaudi.svg)](https://pypi.python.org/pypi/pyaudi)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyaudi/badges/version.svg)](https://anaconda.org/conda-forge/pyaudi)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyaudi/badges/license.svg)](https://anaconda.org/conda-forge/pyaudi)

# Audi

Audi (not the car, rather from latin: “listen!”) is an open source, header only, C++ library (exposed to python in the pyaudi package) that implements the algebra of Taylor truncated polynomials, Taylor models and a few algorithms useful for its applications (differential intelligence, high-order automatic differentiation, verified integration, etc.)

The underlying truncated Taylor polynomial algebra is dealt with using [obake](https://github.com/bluescarni/obake) and can deal with high orders and many variables without eating up the whole system memory.

The polynomial multiplication algorithm used in obake (original with the software author [Francesco Biscani](https://github.com/bluescarni)) takes advantage of sparsity, multiple-threads and cache efficiency allowing a modest memory usage also at high orders.

AuDi computes with the full *truncated Taylor polynomial* (the jet) of a function around an expansion point, rather than propagating a single derivative or a directional derivative. All partial derivatives up to the truncation order and in an arbitrary number of variables are therefore available simultaneously as the coefficients of the polynomial. This makes AuDi a natural tool whenever the object of interest is not just a gradient but the local polynomial model of a function: high-order sensitivity analysis, Taylor-model based verified integration, perturbation theory, and differentiable genetic programming, among others.

In addition to the algebra of Taylor truncated polynomials, a novel Taylor model implementation is provided that is built on top of the Taylor truncated polynomials. This implementation exploits Bernstein polynomials for rapid multivariate polynomial bounding, which is crucial to the performance of Taylor model arithmetic.

AuDi also supports computations over complex numbers and a *vectorized* scalar type. The vectorized type evaluates the Taylor polynomial in many points at once (a form of SIMD/batched evaluation), which is useful in machine-learning style workloads where the same expansion has to be computed over a batch of inputs.

Documentation (preliminary) can be found [here](http://darioizzo.github.io/audi/)

## Relation to modern automatic differentiation libraries

When AuDi was first written (over a decade ago) general-purpose automatic differentiation was mostly the domain of specialized libraries and hand-written codes. The landscape has changed dramatically: modern machine-learning frameworks such as [JAX](https://github.com/jax-ml/jax), [PyTorch](https://pytorch.org/) and [TensorFlow](https://www.tensorflow.org/) now offer excellent, hardware-accelerated automatic differentiation, and JAX in particular can compose its transformations (`jacfwd`/`jacrev`, `jet`, `grad`) to obtain higher-order derivatives, internally relying on Faà di Bruno–style rules for the higher-order chain rule. For the vast majority of gradient-based optimization and deep-learning use cases these tools are the right choice: they are mature, GPU/TPU accelerated, and backed by large communities.

AuDi remains relevant and complementary in a more specific niche: the efficient manipulation of the *entire* truncated Taylor polynomial at high orders and in many variables. Rather than repeatedly applying a differentiation rule, AuDi treats the jet as a first-class algebraic object and multiplies polynomials using obake's sparse, cache-efficient, multi-threaded kernel. As a consequence:

- AuDi scales gracefully to high truncation orders and many variables (roughly `n >= 11`, `m >= 11`) where dense or repeated-differentiation approaches become memory- or time-prohibitive.
- It does not pre-allocate memory for dense polynomials, exploiting sparsity of the actual expansion.
- It provides Taylor-model arithmetic with verified (rigorous) enclosures, which general ML autodiff frameworks do not target.

Other projects in this more classical space are [libtaylor](https://code.google.com/p/libtaylor/), [DACE](https://github.com/dacelib/dace) and [COSY infinity](http://bt.pa.msu.edu/index_cosy.htm). Compared to libtaylor, AuDi can be used from a dynamic library and computes more efficiently at high orders. Compared to DACE, AuDi does not pre-allocate memory for dense polynomials and offers vectorization. Compared to COSY infinity, AuDi is entirely open source and offers vectorization. See the following paper for details on the speed-up w.r.t. COSY on a machine-learning style application:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable Genetic Programming." arXiv preprint arXiv:1611.04766 (2016).

# pyaudi

Pre-compiled pyaudi binaries are available both from the Pyhton Package Index (PyPi) and from conda-forge. Not all architectures are supported, namely only win64 (PyPi), linux 64 (PyPi and conda) and osx (only conda). The best is to try the following:

```
 conda config --add channels conda-forge
 conda install pyaudi
```

or

```
 pip install pyaudi --user
```

# Citing pyaudi

If you use pyaudi in your research, please cite the following paper:

Izzo, Dario, Francesco Biscani, and Sean Cowan. "pyaudi: A truncated Taylor polynomial algebra toolbox for differentiable intelligence, automatic differentiation, and verified integration applications." *Journal of Open Source Software* 11.124 (2026): 9905.

The corresponding BibTeX entry is:

```bibtex
@article{izzo2026pyaudi,
  title={pyaudi: A truncated Taylor polynomial algebra toolbox for differentiable intelligence, automatic differentiation, and verified integration applications.},
  author={Izzo, Dario and Biscani, Francesco and Cowan, Sean},
  journal={Journal of Open Source Software},
  volume={11},
  number={124},
  pages={9905},
  year={2026}
}
```
