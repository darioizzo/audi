.. AuDi documentation master file, created by
   sphinx-quickstart on Thu Sep 17 09:59:25 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/audi_logo.png

Not the car, rather from latin: "listen!") AuDI is an open source, header only, C++ library that allows
for AUtomated DIfferentiation implementing the Taylor truncated polynomial
algebra (forward mode automated differentiation). It was created with the aim
to offer a generic solution to the user in need of an automated
differentiation system. While other automated differentiation codes may be more
efficient on some targeted application requesting a specific order and sparsity, AuDi was built to be fast and efficient
across all application ranges (low orders, high orders, one variable, 
many variables, sparse and dense).

AuDi is internally powered by the truncated polynomial multiplication algorithm
of the open source project `Piranha <https://github.com/bluescarni/piranha>`_ details of which are described in:

Biscani, Francesco. `Parallel sparse polynomial multiplication on modern hardware architectures. <http://dl.acm.org/citation.cfm?id=2442845>`_  Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation. ACM, 2012.

Biscani, Francesco. `Multiplication of sparse Laurent polynomials and Poisson series on modern hardware architectures. <http://arxiv.org/pdf/1004.4548v1.pdf>`_ arXiv preprint arXiv:1004.4548 (2010).

AuDi is thread-safe and makes use of Piranha fine-grained parallelization of the truncated polynomial multiplication. A feature that starts to make a difference at high orders and many variables.