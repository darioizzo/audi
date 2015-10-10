# AuDi (Automated Differentiation)

Implementation of an automated differentiation system using generalized dual numbers (i.e. truncated Taylor expansions). The underlying truncated polynomial algebra (differential algebra) is dealt with using [Piranha](https://github.com/bluescarni/piranha) and can deal with high orders and many variables. 

The polynomial multiplication algorithm used in piranha (original with the software author [Francesco Biscani](https://github.com/bluescarni)) takes advantage of sparsity, multiple-threads and cache efficiency allowing a modest memory usage also at high orders.

AuDi was developed with the goal to surpass the capabilities of [existing automated differentiation libraries](http://www.autodiff.org/?module=Tools) enabling high order differentiation in many variables.

While other automated differentiation codes may be more efficient on some targeted application requesting a specific order and sparsity, AuDi was built to be fast and efficient across all application ranges (low orders, high orders, one variable, many variables, sparse and dense). 

Documentation (preliminary) can be found [here](http://darioizzo.github.io/audi/)
