[![DOI](https://zenodo.org/badge/41722955.svg)](https://zenodo.org/badge/latestdoi/41722955)
[![Build Status](https://travis-ci.org/darioizzo/audi.svg?branch=master)](https://travis-ci.org/darioizzo/audi)
[![PyPI](https://img.shields.io/pypi/v/pyaudi.svg)](https://pypi.python.org/pypi/pyaudi)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyaudi/badges/version.svg)](https://anaconda.org/conda-forge/pyaudi)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyaudi/badges/license.svg)](https://anaconda.org/conda-forge/pyaudi)

# Audi

Audi (not the car, rather from latin: “listen!”) is an open source, header only, C++ library (exposed to python in the pyaudi package) that implements the algebra of Taylor truncated polynomials, Taylor models and a few algorithms useful for its applications (differential intelligence, high-order automatic differentiation, verified integration, etc.)

The underlying truncated Taylor polynomial algebra (a differential algebra since integration and derivations are defined too) is dealt with using [obake](https://github.com/bluescarni/obake) and can deal with high orders and many variables without eating up the whole system memory.

The polynomial multiplication algorithm used in obake (original with the software author [Francesco Biscani](https://github.com/bluescarni)) takes advantage of sparsity, multiple-threads and cache efficiency allowing a modest memory usage also at high orders.

AuDi was developed with the goal to surpass the capabilities of [existing automated differentiation libraries](http://www.autodiff.org/?module=Tools) enabling high order differentiation in many variables.

While other automated differentiation codes may be more efficient on some targeted application requesting a specific order and sparsity, AuDi was built to be fast and efficient at very high orders and number of variables.

In addition to the algebra of Taylor truncated polynomials, a novel Taylor model implementation is provided that is built on top of the Taylor truncated polynomials. This implementation exploits Bernstein polynomials for rapid multivariate polynomial bounding, which is crucial to the performance of Taylor model arithmetic.

Documentation (preliminary) can be found [here](http://darioizzo.github.io/audi/)

## Comparison with existing code

Alternative projects that have similar capabilities to AuDi are [libtaylor](https://code.google.com/p/libtaylor/), [DACE](https://github.com/dacelib/dace) and [COSY infinity](http://bt.pa.msu.edu/index_cosy.htm). Unlike libtaylor AuDi can be used in a dynamic library and can compute at high orders with greater efficiency. Unlike DACE AuDi does not pre-allocate memory for dense polynomials and has vectorization capabilities. Unlike COSY infinity AuDi is entirely open source and has vectorization capabilities.

From the point of view of efficiency, the main difference of AuDi w.r.t. existing codes is in the polinomial multiplication algorithm. AuDi uses the third party [obake](https://github.com/bluescarni/obake) code and thus gets all the pros and cons of that particular algebraic manipulation system which is still actively developed and was born to deal with massively large polynomial manipulations typically encountered in celestial mechanics perturbation theory. To cut a long story short, AuDi will be "unbeatable" for high orders and many variables (n>=11, m>=11). Below this orders AuDi will still be incredibly memory efficient and fast when used in a machine where multiple threading capabilities are possible.

NEW: Audi also allows for making computations using complex numbers and a vectorized type. This last type allows to compute the derivatives (Taylor polynomial) in multiple points at once, making audi the fastest code of its kind for application such as Machine Learning where this has use. See the folloiwng paper for details on the speed-up w.r.t COSY:

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
