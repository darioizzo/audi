.. AuDi documentation master file, created by
   sphinx-quickstart on Thu Sep 17 09:59:25 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/audi_logo.png

Audi (not the car, rather from latin: "listen!") is an open source, header only, C++ library that allows
for high order AUtomated DIfferentiation implementing the Taylor truncated polynomial
algebra (forward mode automated differentiation). It was created with the aim
to offer a solution to the user in need of a generic automated
differentiation system. While other automated differentiation codes may be more
efficient on some targeted application, AuDi was built to be fast and efficient
across all application ranges (low/high order of derivation, one variable,
many variables, sparse/dense cases).

AuDi is internally powered by the truncated polynomial multiplication algorithm
of the open source project `Piranha <https://github.com/bluescarni/piranha>`_ some details of which are described in:

Biscani, Francesco. `Parallel sparse polynomial multiplication on modern hardware architectures. <http://dl.acm.org/citation.cfm?id=2442845>`_  Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation. ACM, 2012.

Biscani, Francesco. `Multiplication of sparse Laurent polynomials and Poisson series on modern hardware architectures. <http://arxiv.org/pdf/1004.4548v1.pdf>`_ arXiv preprint arXiv:1004.4548 (2010).

AuDi is thread-safe and, when possible, use of Piranha fine-grained parallelization of the truncated polynomial multiplication. The benefits of this fine grained parallelization are well visible for many variables and high differentiation orders.
