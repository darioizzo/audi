Generalized Dual Number
========================

*#include <audi/gdual.hpp>*

.. cpp:class:: template <typename Cf> gdual

This class represents an element of the algebra :math:`\mathcal P_{n,m}`, (see :ref:`formal_definition`) defined
over the field :math:`\mathbf K` (the field is represented by the template argument *Cf*). We call elements of this algebra
generalized dual numbers as, among other things, they generalize the dual numbers used for forward automatic differentiation.

Using the multi-index notation, a generalized dual number (over the field :math:`\mathbb R` represented by doubles) may be written as:

.. math::

   T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)

.. note::

   A generalized dual number is defined by its order :math:`m` as well as by its expansion point :math:`\mathbf a`. 
   All arithmetic operators +,*,/,- are overloaded so that computing with the type audi::gdual correspond to compute the
   Taylor expansion of the arithmetic computations.

A basic use case, where :math:`m = 2`, :math:`\mathbf a = [1.2, -0.1]` and :math:`f = \frac{x1+x2}{x1-x2}` is thus:

.. code::

   gdual<double> x1(1.2, "x1", 2);
   gdual<double> x2(-0.1, "x2", 2);
   std::cout << (x1+x2) / (x1-x2) << "\n";

resulting in the output:

.. code::

   0.118343*dx1+0.846154+1.42012*dx2+1.0924*dx2**2-0.0910332*dx1**2-1.00137*dx1*dx2

.. note::

  The class can be instantiated with any type that is suitable to be a coefficient in a piranha polynomial (
  piranha::is_cf<Cf>::value must be true). Classical examples would be double, float, std::complex<double>, and
  the audi::vectorized_double type. If piranha::is_differentiable<Cf>::value is also true then derivation
  and integration are also availiable and the resulting algebra is a differential algebra.


Constructors and destructors
------------------------------------------------------

.. cpp:function:: template <typename T> explicit gdual(const T &c)

   .. note::  

      This constructor is only enabled if *T* is not a gdual and if a piranha::polynomial can be constructed with coefficients in *T*.

   Constructs a gdual of order :math:`m = 0` with only a constant term :math:`c`. Essentially :math:`T_f = c` in :math:`\mathcal P_{0,0}`

   :param c: The constant term :math:`c`.

------------------------------------------------------

.. cpp:function:: template <typename T> explicit gdual(const T &c, const std::string &symbol, unsigned m)

   .. note::  

      This constructor is only enabled if *T* is not a gdual and if a piranha::polynomial can be constructed with coefficients in *T*.

   Constructs a gdual of order :math:`m` representing a variable. Essentially :math:`T_f = c + dx` in :math:`\mathcal P_{1,m}`

   :param c: The constant term :math:`c`.
   :param symbol: The symbol that will represent the variable (e.g. :math:`x`).
   :param m: The truncation order of the algebra.

   :exception: std::invalid_argument if the symbol names starts already with "d" a differential. This avoids symbols like ddx in the piranha::polynomial.

------------------------------------------------------

.. cpp:function:: gdual()

   Default constructor. Constructs a gdual of order :math:`m = 0` with only a zero constant term :math:`c`. Essentially :math:`T_f = 0` in :math:`\mathcal P_{0,0}`

------------------------------------------------------

.. cpp:function:: gdual(const gdual &) = default

   Defaulted copy constructor.

------------------------------------------------------

.. cpp:function:: gdual(gdual &&) = default

   Defaulted move constructor.

------------------------------------------------------

.. cpp:function:: ~gdual()

   Destructor. Performs a sanity check on the truncation order and degree of the gdual.

------------------------------------------------------


Methods
-------

Symbol set manipulation
^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: auto get_symbol_set_size() const

   Returns the size of the symbol set of the piranha::polynomial

------------------------------------------------------

.. cpp:function:: std::vector<std::string> get_symbol_set() const

   Returns the symbol set stripping the differentials.

------------------------------------------------------

.. cpp:function:: void extend_symbol_set(const std::vector<std::string> &sym_vars)

   Adds some symbolic variables to the current piranha::polynomial
   This is useful in situations where some differential :math:`dx` does not
   appear as its coefficient is zero but we still want to treat the gdual as a function of 
   :math:`x` too (for example when extracting the relative coefficient)

   :param sym_vars: list of symbolic names. It must contain all symbolic names of
     the current piranha::polynomial. It may contain more. All symbols must start with the letter "d".

   :exception: std::invalid_argument if any symbol in *sym_vars* does not start with the letter "d"
     or if *sym_vars* does not contain all current symbols too.

------------------------------------------------------

Differential algebra operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: template<> gdual integrate(const std::string &var_name)

   .. note::
      
      This template is only enabled if *Cf* satisfies piranha::is_differentiable, which is
      the case for float, double, std::complex and vectorized_double types. 

   Performs the integration of the gdual with respect to *var_name*. If the *var_name* differential is not in the symbol set
   of the piranha::polynomial it is added. 
   Note that Information may be lost as the truncation order is preserved.

   :param var_name: Symbol name (cannot start with "d").

   :exception: std::invalid_argument if *var_name* starts with the letter "d".

------------------------------------------------------

.. cpp:function:: template<> gdual partial(const std::string &var_name)

   .. note::
      
      This template is only enabled if *Cf* satisfies piranha::is_differentiable, which is
      the case for float, double, std::complex and vectorized_double types. 

   Performs the partial derivative of the gdual with respect to *var_name*. If the *var_name* differential is not in the symbol set
   of the piranha::polynomial it is added. 

   :param var_name: Symbol name (cannot start with "d").

   :exception: std::invalid_argument if *var_name* starts with the letter "d".

------------------------------------------------------

gdual manipulations
^^^^^^^^^^^^^^^^^^^

.. cpp:function:: template <typename T> gdual subs(const std::string &sym, const T &val) const

   Substitute the differential *sym* with *val*. The *Cf* type must be constructable from *val*.

   :param sym: The name of the differential to be substituted.
   :param val: The value to substitute *sym* with.
   :return: A new gdual with the substitution made.

------------------------------------------------------

.. cpp:function:: gdual subs(const std::string &sym, const gdual &val) const

   Substitute the differential *sym* with the gdual *val*

   :param sym: The name of the differential to be substituted.
   :param val: The value to substitute *sym* with.
   :return: A new gdual with the substitution made.

------------------------------------------------------

.. cpp:function:: gdual trim(double epsilon) const

   Sets to zero all coefficients of the gdual with absolute value smaller than *epsilon*.

   :param epsilon: Tolerance for the trim.
   :return: A new gdual without the trimmed coefficients.

------------------------------------------------------

.. cpp:function:: gdual extract_terms(unsigned degree) const

   Extracts all the monomials of a given *degree*.

   :param order: The monomials degree.
   :return:  A new gdual containing only the terms extracted, but preserving the truncation order of the original gdual.

   :exception: std::invalid_argument if the *degree* is higher than the gdual truncation order.