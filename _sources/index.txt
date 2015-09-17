.. AuDi documentation master file, created by
   sphinx-quickstart on Thu Sep 17 09:59:25 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome
================================

Audi (not the car, rather from latin: "listen!") is an open source, header only, C++ library that allows
for AUtomated DIfferentiation implementing the Taylor truncated polynomial
algebra (forward mode automated differentiation). It was created with the aim
to offer a generic solution to the user in need of an automated
differentiation system. While other automated differentiation codes may be more
efficient on some targeted application, AuDi was built to be fast and efficient
across all application ranges (low orders, high orders, one variable, 
many variables).

AuDi is internally powered by the truncated polynomial multiplication algorithm
of the open source project,  <a href="https://github.com/bluescarni/piranha">Piranha</a> details of which are described in:

Biscani, Francesco. <a href="http://dl.acm.org/citation.cfm?id=2442845"> "Parallel sparse polynomial multiplication on modern hardware architectures."</a> Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation. ACM, 2012.

Biscani, Francesco. <a href="http://arxiv.org/pdf/1004.4548v1.pdf">"Multiplication of sparse Laurent polynomials and Poisson series on modern hardware architectures."</a> arXiv preprint arXiv:1004.4548 (2010).

