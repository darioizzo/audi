[![DOI](https://zenodo.org/badge/41722955.svg)](https://zenodo.org/badge/latestdoi/41722955)
[![Build status](https://ci.appveyor.com/api/projects/status/c8kjdcrxy52t8gy4?svg=true)](https://ci.appveyor.com/project/darioizzo/audi)
[![Build Status](https://travis-ci.org/darioizzo/audi.svg?branch=master)](https://travis-ci.org/darioizzo/audi)
[![PyPI](https://img.shields.io/pypi/v/pyaudi.svg)](https://pypi.python.org/pypi/pyaudi)

# AuDi (Automated Differentiation) 

Implementation of a high-order automated differentiation system using generalized dual numbers (i.e. truncated Taylor polynomials). The underlying truncated Taylor polynomial algebra (a differential algebra since integration and derivations are defined too) is dealt with using [Piranha](https://github.com/bluescarni/piranha) and can deal with high orders and many variables without eating up the whole system memory.

The polynomial multiplication algorithm used in piranha (original with the software author [Francesco Biscani](https://github.com/bluescarni)) takes advantage of sparsity, multiple-threads and cache efficiency allowing a modest memory usage also at high orders.

AuDi was developed with the goal to surpass the capabilities of [existing automated differentiation libraries](http://www.autodiff.org/?module=Tools) enabling high order differentiation in many variables.

While other automated differentiation codes may be more efficient on some targeted application requesting a specific order and sparsity, AuDi was built to be fast and efficient across all application ranges (low orders, high orders, one variable, many variables, sparse and dense). 

Documentation (preliminary) can be found [here](http://darioizzo.github.io/audi/)

## Comparison with existing code

Alternative projects that have similar capabilities to AuDi are [libtaylor](https://code.google.com/p/libtaylor/) and [COSY infinity](http://bt.pa.msu.edu/index_cosy.htm). Unlike libtaylor AuDi can be used in a dynamic library and can compute at high orders with greater efficiency. Unlike COSY infinity AuDi is entirely open source. 

From the point of view of efficiency, the main difference of AuDi w.r.t. existing codes is in the polinomial multiplication algorithm. AuDi uses the third party [Piranha](https://github.com/bluescarni/piranha) code and thus gets all the pros and cons of that particular algebraic manipulation system which is still actively developed and was born to deal with massively large polynomial manipulations typically encountered in celestial mechanics perturbation theory. To cut a long story short, AuDi will be "unbeatable" for high orders and many variables (n>=11, m>=11). Below this orders AuDi will still be incredibly memory efficient and fast when used in a machine where multiple threading capabilities are possible.

NEW: Audi also allows for making computations using complex numbers and a vectorized type. This last type allows to compute the derivatives (Taylor polynomial) in multiple points at once, making audi the fastest code of its kind for application such as Machine Learning where this has use. See the folloiwng paper for details on the speed-up w.r.t COSY:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable Genetic Programming." arXiv preprint arXiv:1611.04766 (2016).

# pyaudi
The python bindings are now available. You can build them activating the BUILD_PYAUDI option in CMake, or install the precompiled binaries (for the last release) typing:

 ```pip install pyaudi```

Only Windows 32 bits and OSX are not provided
