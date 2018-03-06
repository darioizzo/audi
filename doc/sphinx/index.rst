Audi and pyaudi
================================

Audi (not the car, rather from latin: "listen!") is an open source, header only, C++ library 
(exposed to python in the pyaudi package) that implements the algebra of Taylor truncated polynomials and
a few algorithms useful for its applications (Differential Intelligence, Taylor Models, etc.)

It was created with the aim to offer a solution to the user in need of a generic automated
differentiation system. While other automated differentiation codes may be more
efficient on some targeted application, Audi was built to be fast and efficient
across all application ranges (low/high order of derivation, one variable,
many variables, sparse/dense cases).

Audi is internally powered by the truncated polynomial multiplication algorithm
of the open source project `Piranha <https://github.com/bluescarni/piranha>`_ some details of which are described in:

Biscani, Francesco. `Parallel sparse polynomial multiplication on modern hardware architectures. <http://dl.acm.org/citation.cfm?id=2442845>`_  Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation. ACM, 2012.

Biscani, Francesco. `Multiplication of sparse Laurent polynomials and Poisson series on modern hardware architectures. <http://arxiv.org/pdf/1004.4548v1.pdf>`_ arXiv preprint arXiv:1004.4548 (2010).

.. note::

   Audi is thread-safe and, when possible, makes use of Piranha fine-grained parallelization of the truncated polynomial multiplication.
   The benefits of this fine grained parallelization are well visible for many variables and high differentiation orders.

Audi is open source (GPL3) and available in github <https://github.com/darioizzo/audi>`_


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   installation
   theory
   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
