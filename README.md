# AuDi (Automated Differentiation)

Implementation of an automated differentiation system using generalized dual numbers (i.e. truncated Taylor expansions). The underlying truncated polynomial algebra (differential algebra) is dealt with using [Piranha](https://github.com/bluescarni/piranha) and can deal with high orders and many variables. 

The polynomial multiplication algorithm used in piranha (original with the software author [Francesco Biscani](https://github.com/bluescarni)) takes advantage of sparsity, multiple-threads and cache efficiency. 

AuDi was developed with the goal to surpass the capabilities of [existing automated differentiation libraries](http://www.autodiff.org/?module=Tools) enabling high order differentiation in many variables.
