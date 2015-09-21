.. cpp docs

C++ Documentation
=================

.. contents::

Classes
-------

gdual: Generalized dual number class.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: audi::gdual
   :project: AuDi
   :members:

Functions
---------

exp: Overload for the exponential.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::exp(const gdual&)
   :project: AuDi

----------------------------------------------------------

log: Overload for the logarithm.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::log(const gdual&)
   :project: AuDi

----------------------------------------------------------

pow: Overload for pow(const gdual&, double)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::pow(const gdual&, double)
   :project: AuDi

----------------------------------------------------------

pow: Overload for pow(const gdual&, int)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::pow(const gdual&, int)
   :project: AuDi

----------------------------------------------------------

pow: Overload for pow(double, const gdual&)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::pow(double, const gdual&)
   :project: AuDi

----------------------------------------------------------

pow: Overload for pow(const gdual&, const gdual&)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::pow(const gdual&, const gdual&)
   :project: AuDi

----------------------------------------------------------

sin: Overload for the sine.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::sin(const gdual&)
   :project: AuDi

----------------------------------------------------------

cos: Overload for the cosine.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::cos(const gdual&)
   :project: AuDi

----------------------------------------------------------

sin_and_cos: Computes both at once
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::sin_and_cos(const gdual&, gdual sine&, gdual& cosine)
   :project: AuDi

----------------------------------------------------------

tan: Overload for the tangent.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::tan(const gdual&)
   :project: AuDi

----------------------------------------------------------

sqrt: Overload for the square root.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::sqrt(const gdual&)
   :project: AuDi

----------------------------------------------------------

cbrt: Overload for the cubic root.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: audi::cbrt(const gdual&)
   :project: AuDi

