Generalized Dual Number
========================

*#include <audi/gdual.hpp>*

.. cpp:class:: template <typename Cf> gdual

This class represents an element of the algebra :math:`\mathcal P_{n,m}`, (see :ref:`formal_definition`) defined
over the field :math:`\mathbf K` (the field is represented by the template argument *Cf*). We call elements of this algebra
generalized dual numbers as, among other things and when :math:`\mathbf K` is :math:`\mathbb R`), they generalize the dual numbers used for forward automatic differentiation.

Using the multi-index notation, a generalized dual number (for example over the field :math:`\mathbb R` represented by doubles, i.e. :code:`gdual<double>`) 
may be written as:

.. math::

   T_f(\mathbf x) = \sum_{|\alpha| = 0}^m  \frac{(\mathbf x-\mathbf a)^\alpha}{\alpha!}(\partial^\alpha f)(\mathbf a)

.. note::

   A generalized dual number is defined by its truncation order :math:`m` as well as by its expansion point :math:`\mathbf a`. 
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

------------------------------------------------------

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

gdual manipulation
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

gdual inspection
^^^^^^^^^^^^^^^^^

.. cpp:function:: auto get_symbol_set_size() const

   Returns the size of the symbol set of the piranha::polynomial

------------------------------------------------------

.. cpp:function:: std::vector<std::string> get_symbol_set() const

   Returns the symbol set of the piranha::polynomial stripping the differentials (i.e. "dx" becomes "x")

------------------------------------------------------

.. cpp:function:: auto evaluate(const std::unordered_map<std::string, double> &dict) const

   Evaluates the Taylor polynomial using the values in *dict* for the differentials (variables variations)

   :param dict: a dictionary (unordered map) containing the values for the differentials.
   :return: the value of the Taylor polynomial.

   .. code::

      gdual<double> x1(1., "x1", 2);
      gdual<double> x2(1., "x2", 2);
      auto f = x1*x1 + x2;
      std::cout << f.evaluate({{"dx1", 1.}, {"dx2", 1.}}) << "\n";

------------------------------------------------------

.. cpp:function:: auto degree() const

   Returns the degree of the polynomial. This is necessarily smaller or equal the truncation order.

   :return: the polynomial degree.

------------------------------------------------------

.. cpp:function:: unsigned int get_order() const

   Returns the truncation order of the polynomial.

   :return: the polynomial truncation order.

------------------------------------------------------

.. cpp:function:: template <typename T> auto find_cf(const T &c) const

   Returns the coefficient of the monomial specified by the container *c*.
   The container contains the exponents of the requested monomial. In a three
   variable Taylor expansion with :math:`x, y, z` as symbols, [1, 3, 2] would denote
   the monomial :math:`dx dy^3 dz^2`.

   .. note::
   
      Alphabetical order is used to order the symbol set and thus specify
      the coefficients.

   .. warning::
     
     If the monomial requested is not found in the polynomial a zero coefficient will be returned.

   :return: the coefficient of the monomial, if present, zero otherwise.

   :exception: std::invalid_argument if the requested coefficient is beyond the truncation order
   
   .. code::

      gdual<double> x1(1.2, "x1", 2);
      gdual<double> x2(-0.2, "x2", 2);
      auto f = sin(x1*x1 + x2);
      std::cout << f.find_cf(std::vector<double>({1,1})) << "\n";

------------------------------------------------------

.. cpp:function:: template <typename T> auto find_cf(std::initializer_list<T> l) const

   Returns the coefficient of the monomial specified by the initializer list *l*.

   .. note::
   
      This method overloads the one above and is provided for convenience.

   :return: the coefficient of the monomial, if present, zero otherwise.

   :exception: std::invalid_argument if the requested coefficient is beyond the truncation order
   
   .. code::

      gdual<double> x1(1.2, "x1", 2);
      gdual<double> x2(-0.2, "x2", 2);
      auto f = sin(x1*x1 + x2);
      std::cout << f.find_cf({1,1}) << "\n";

------------------------------------------------------

.. cpp:function:: Cf constant_cf()

   Finds the constant coefficient of the Taylor polynomial. So that if :math:`T_{f} = f_0 + \hat f`, :math:`f_0` is returned

   :return: the constant coefficient.

------------------------------------------------------

.. cpp:function:: template <typename T> auto get_derivative(const T &c) const

   Returns the (mixed) derivative value of order specified by the container *c*

   .. note:: 
   
      The container describes the derivative requested. In a three
      variable polynomial with :math:`x, y, z` as symbols, [1, 3, 2] would denote
      the sixth order derivative :math:`\frac{d^6}{dxdy^3dz^2}`.

   .. note::

      No computations are made at this points as all derivatives are already
      contained in the Taylor expansion

   :return: the value of the derivative

   :exception: std::invalid_argument if the requested coefficient is beyond the truncation order

   .. code::

      gdual<double> x1(1.2, "x1", 2);
      gdual<double> x2(-0.2, "x2", 2);
      auto f = sin(x1*x1 + x2);
      // This streams the value of df^2/dx1/dx2 in x1=1.2, x2 = -0.2
      std::cout << f.get_derivative(std::vector<double>({1,1})) << "\n";

------------------------------------------------------

.. cpp:function:: template <typename T> auto get_derivative(std::initializer_list<T> l) const

   Returns the (mixed) derivative value of order specified by the initializer list *l*.

   .. note::
   
      This method overloads the one above and is provided for convenience.

   :return: the value of the derivative

   :exception: std::invalid_argument if the requested coefficient is beyond the truncation order

   .. code::

      gdual<double> x1(1.2, "x1", 2);
      gdual<double> x2(-0.2, "x2", 2);
      auto f = sin(x1*x1 + x2);
      // This streams the value of df^2/dx1/dx2 in x1=1.2, x2 = -0.2
      std::cout << f.get_derivative({1,1}) << "\n";

------------------------------------------------------

.. cpp:function:: template <typename T> auto get_derivative(const std::unordered_map<std::string, unsigned int> &dict) const

   Returns the (mixed) derivative value of the order specified in *dict*.

   :param dict: a dictionary (unordered map) containing the derivation order (assumes zero for symbols not present).

   :return: the value of the derivative

   :exception: std::invalid_argument if the requested derivative degree is beyond the truncation order

   .. code::

      gdual<double> x1(1.2, "x1", 2);
      gdual<double> x2(-0.2, "x2", 2);
      auto f = sin(x1*x1 + x2);
      // This streams the value of df^2/dx1/dx2 in x1=1.2, x2 = -0.2
      std::cout << f.get_derivative({{"dx1", 1u}, {"dx2", 1u}}) << "\n";

------------------------------------------------------

.. cpp:function:: bool is_zero(double tol) const

   Checks all coefficients of the gdual and returns true if all their absolute values are smaller
   than the defined tolerance *tol*.

   :return: whether the gdual is zero within a tolerance.

   .. code::

      gdual<double> x(0.1, "x", 6);
      auto f = 1 - sin(x)*sin(x) - cos(x)*cos(x);
      if (f.is_zero(1e-13)) {
        std::cout << "The trigonomoetric identity holds!!" << std::endl;
      }

------------------------------------------------------

Operators
---------

The following operators are implemented: 

  * <<, streaming 
  * ==, equal to 
  * =,  assignement
  * !=, not equal to 
  * +=, addition assignment
  * -=, subtraction assignment
  * \*=, multiplication assignment
  * /=, division assignment
  * -, unary minus
  * +, unary plus
  * +, addition
  * -, subtraction
  * \*, multiplication
  * /, division

They allow to compute with the type gdual as you would operate with a basic type.

.. note:: 

   When relevant, the operators implement order promotion so that, for example, if a gdual of order 2 is added to a
   gdual of order 3 the resulting gdual will have order three.

We specify the documentation of a few operators with non trivial meaning.

------------------------------------------------------

.. cpp:function:: friend std::ostream &operator<<(std::ostream &os, const gdual &d)

   Will direct to stream a human-readable representation of the generalized dual number.

   .. note::
      
      The print order of the terms will be undefined. At most piranha::settings::get_max_term_output() terms
      are printed, and terms in excess are represented with ellipsis "..."
      at the end of the output; if piranha::settings::get_max_term_output()
      is zero, all the terms will be printed. piranha::settings::set_max_term_output()
      is used to set this parameter.

   :param os: target stream.
   :param d: gdual argument.

   :return: reference to *os*

------------------------------------------------------

.. cpp:function:: friend bool operator==(const gdual &d1, const gdual &d2)

   Equality operator. Two gduals are considered equal if all their coefficients are equal.

   .. note:: 
   
      The truncation order of *d1* and *d2* may be different

   :param d1: first argument.
   :param d2: second argument.

   :return: The result of the comparison.

------------------------------------------------------

.. cpp:function:: friend bool operator!=(const gdual &d1, const gdual &d2)

   Non equality operator.

   .. note:: 
   
      The truncation order of *d1* and *d2* may be different

   :param d1: first argument.
   :param d2: second argument.

   :return: The result of the comparison.

