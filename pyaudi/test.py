from __future__ import absolute_import as _ai

import unittest as _ut


def some_complex_rational_f(x, y, z):
    return (x * x * x + x * y * z + z * x * y) * (x * x * x + x * y * z + z * x * y) / (x * x * x + x * y * z + z * x * y) * (x * x * x + x * y * z + z * x * y)


def some_complex_irrational_f(x, y, z):
    from pyaudi import exp, log, cos, sin, tan, sqrt, cbrt, cos, sin, tan, acos, asin, atan, cosh, sinh, tanh, acosh, asinh, atanh
    from pyaudi import abs as gd_abs
    from pyaudi import sin_and_cos, sinh_and_cosh
    f = (x + y + z) / 10.
    retval = exp(f) + log(f) + f**2 + sqrt(f) + cbrt(f) + cos(f) + sin(f)
    retval += tan(f) + acos(f) + asin(f) + atan(f) + cosh(f) + sinh(f)
    retval += tanh(f) + acosh(f) + asinh(f) + atanh(f)
    a = sin_and_cos(f)
    b = sinh_and_cosh(f)
    retval += a[0] + a[1] + b[0] + b[1]
    return retval


class test_gdual_double(_ut.TestCase):

    def test_construction(self):
        from pyaudi import gdual_double as gdual
        x = gdual(2., "x", 3)
        c = gdual(-2.)
        with self.assertRaises(ValueError):
            x = gdual(2., "dx", 3)
        self.assertEqual(x.constant_cf, 2)
        self.assertEqual(x.find_cf([0]), 2)
        self.assertEqual(c.constant_cf, -2)
        with self.assertRaises(ValueError):
            x.find_cf([1, 3])
            c.find_cf([1, 3])
            c.find_cf([0])  # no symbols in constants
        with self.assertRaises(ValueError):
            x.find_cf([5])
            c.find_cf([1])
        self.assertEqual(x.degree, 1)
        self.assertEqual(x.order, 3)
        self.assertEqual(x.symbol_set, ["x"])
        self.assertEqual(x.symbol_set_size, 1)
        self.assertEqual(c.degree, 0)
        self.assertEqual(c.order, 0)
        self.assertEqual(c.symbol_set, [])
        self.assertEqual(c.symbol_set_size, 0)
        with self.assertRaises(TypeError):
            gdual([1.], "x", 2)
            gdual([1.])

    def test_order_promotion(self):
        from pyaudi import gdual_double as gdual
        c = gdual(1)
        x = gdual(2, "x", 4)
        y = gdual(-2, "y", 2)
        self.assertEqual(c.order, 0)
        self.assertEqual((x + y).order, 4)
        self.assertEqual((x - y).order, 4)
        self.assertEqual((x * y).order, 4)
        self.assertEqual((x / y).order, 4)
        self.assertEqual((x + c).order, 4)
        self.assertEqual((x - c).order, 4)
        self.assertEqual((x * c).order, 4)
        self.assertEqual((x / c).order, 4)
        self.assertEqual((y + x).order, 4)
        self.assertEqual((y - x).order, 4)
        self.assertEqual((y * x).order, 4)
        self.assertEqual((y / x).order, 4)
        self.assertEqual((c + x).order, 4)
        self.assertEqual((c - x).order, 4)
        self.assertEqual((c * x).order, 4)
        self.assertEqual((c / x).order, 4)

    def test_addition(self):
        from pyaudi import gdual_double as gdual
        p1 = gdual(1)
        p2 = gdual(0, "x", 4)
        self.assertEqual(p1 + p2, p2 + p1)
        self.assertEqual(1 + p2, p2 + 1)
        self.assertEqual(1. + p2, p2 + 1.)

    def test_subtraction(self):
        from pyaudi import gdual_double as gdual
        p1 = gdual(1)
        p2 = gdual(0, "x", 4)
        self.assertEqual(p1 - p1, gdual(0))
        self.assertEqual(p2 - p2, gdual(0))
        self.assertEqual(p1 - p2, - (p2 - p1))
        self.assertEqual(1 + p1 - p1, gdual(1))
        self.assertEqual(1 + p2 - p2, gdual(1))
        self.assertEqual(1. + p1 - p1, gdual(1))
        self.assertEqual(1. + p2 - p2, gdual(1))
        self.assertEqual((1 - p1) + (p1 + p2), 1 + p2)

    def test_multiplication(self):
        from pyaudi import gdual_double as gdual
        x = gdual(3, "x", 2)
        y = gdual(7, "y", 2)
        p1 = x + 3 * x * y + y * y
        self.assertEqual(p1.find_cf([0, 0]), 115)
        self.assertEqual(p1.find_cf([1, 0]), 22)
        self.assertEqual(p1.find_cf([0, 1]), 23)
        self.assertEqual(p1.find_cf([1, 1]), 3)
        self.assertEqual(p1.find_cf([2, 0]), 0)
        self.assertEqual(p1.find_cf([0, 2]), 1)

    def test_division(self):
        from pyaudi import gdual_double as gdual
        x = gdual(0, "x", 2)
        y = gdual(1, "y", 2)
        p1 = 1 / (x + 2 * x * y + y * y)
        self.assertEqual(p1.find_cf([0, 0]), 1)
        self.assertEqual(p1.find_cf([1, 0]), -3)
        self.assertEqual(p1.find_cf([0, 1]), -2)
        self.assertEqual(p1.find_cf([1, 1]), 10)
        self.assertEqual(p1.find_cf([2, 0]), 9)
        self.assertEqual(p1.find_cf([0, 2]), 3)

    def test_identities(self):
        from pyaudi import gdual_double as gdual
        x = gdual(2, "x", 3)
        y = gdual(3, "y", 3)

        p1 = x * x + y - x * x * x * x * y - y * x * x
        p2 = y * y - x + y * y * y * y * x - 2 * x
        self.assertEqual((x + y) * (x + y), x * x + y * y + 2 * x * y)
        self.assertEqual((p1 + p2) * (p1 + p2), p1 *
                         p1 + p2 * p2 + 2 * p1 * p2)
        self.assertEqual(x * x * x * x - y * y * y * y,
                         (x - y) * (x + y) * (x * x + y * y))
        self.assertEqual(p1 * p1 * p1 * p1 - p2 * p2 * p2 * p2,
                         (p1 - p2) * (p1 + p2) * (p1 * p1 + p2 * p2))
        self.assertTrue(((p1 / p2) * (p2 / p1) - 1).is_zero(1e-12))
        self.assertTrue(((p1 / p2) * p2 - p1).is_zero(1e-12))

    def test_derivatives(self):
        from pyaudi import gdual_double as gdual

        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = x * x * x * x * x + x * y * z * x * x + z * x * y * y * y

        self.assertEqual(f.get_derivative([1, 1, 1]), 6.)
        self.assertEqual(f.get_derivative([2, 1, 1]), 6.)
        self.assertEqual(f.get_derivative([1, 2, 1]), 6.)
        self.assertEqual(f.get_derivative([1, 1, 2]), 0.)
        self.assertEqual(f.get_derivative([4, 0, 0]), 120.)
        with self.assertRaises(ValueError):
            f.get_derivative([4, 1, 1])

        x = gdual(1, "x", 8)
        y = gdual(1, "y", 8)
        f = (x * y + 2 * x * x * y) / (1 + x + y)
        self.assertEqual(f.get_derivative([1, 1]), 1.)
        self.assertAlmostEqual(f.get_derivative([2, 2]), 0, delta=1e-12)
        self.assertAlmostEqual(f.get_derivative([3, 3]), 0, delta=1e-12)
        self.assertAlmostEqual(f.get_derivative([4, 4]), 0, delta=1e-12)

        # we test the dictionary interface
        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = x * x * x * x * x + x * y * z * x * x + z * x * y * y * y
        dictionary = {"dx": 1, "dy": 1, "dz": 1}
        self.assertEqual(f.get_derivative(dictionary), 6.)
        with self.assertRaises(ValueError):
            f.get_derivative([4, 1, 1])
            dictionary = {"dx": 4, "dy": 1, "dz": 1}
            f.get_derivative(dictionary)
        dictionary = {"dx": 1, "dr": 1, "dz": 1}
        self.assertEqual(f.get_derivative(dictionary), 0.)

    def test_integrate_partial(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = x * x * x + x * y * z + z * x * y
        fx = 3 * x * x + y * z + z * y
        fy = x * z + z * x
        fz = y * x + x * y

        self.assertEqual(f.partial("x"), fx)
        self.assertEqual(f.partial("y"), fy)
        self.assertEqual(f.partial("z"), fz)

    def test_subs(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1, "x", 1)
        y = gdual(-1, "y", 1)
        f = x * y * x / (x - y)  # -0.75*dx-0.5+0.25*dy
        res2 = f.subs("dx", 1)
        res3 = f.subs("dy", 1)
        self.assertEqual(res2.constant_cf, -1.25)
        self.assertEqual(res2.constant_cf, -1.25)
        self.assertEqual(res2.get_derivative({"dy": 1}), 0.25)
        self.assertEqual(res2.get_derivative({"dx": 1}), 0.)
        self.assertEqual(res3.constant_cf, -0.25)
        self.assertEqual(res3.constant_cf, -0.25)
        self.assertEqual(res3.get_derivative({"dx": 1}), -0.75)
        self.assertEqual(res3.get_derivative({"dy": 1}), 0.)

        x = gdual(1, "x", 4)
        x = x**2
        y = gdual(0., "y", 1)
        newx = x.subs("dx", y)
        self.assertEqual(newx, gdual(1, "y", 4)*gdual(1, "y", 4))
        self.assertEqual(newx.order, 4)


    def test_trim(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1e-4, "x", 1)
        self.assertEqual(x, x.trim(1e-5))
        self.assertEqual(x.trim(1e-3), gdual(0,"x",1))

    def test_extract_terms(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sin
        x = gdual(1e-4, "x", 3)
        sinx = sin(x)
        self.assertEqual(sinx, sinx.extract_terms(0) + sinx.extract_terms(1)+ sinx.extract_terms(2)+ sinx.extract_terms(3))
        self.assertEqual(x.trim(1e-3), gdual(0,"x",1))

    def test_serialization(self):
        from pyaudi import gdual_double as gdual
        import pickle as pk
        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = (x * x * x + x * y * z + z * x * y) * (x * x * x + x * y * z + z * x * y) * \
            (x * x * x + x * y * z + z * x * y) * \
            (x * x * x + x * y * z + z * x * y)
        pk.dump(f, open("tmp.pk", "wb"))
        new_f = pk.load(open("tmp.pk", "rb"))
        self.assertEqual(f, new_f)
        self.assertEqual(f.order, new_f.order)


class test_function_calls(_ut.TestCase):

    def runTest(self):
            self.test_exponentiation()
            self.test_sqrt()
            self.test_cbrt()
            self.test_exp_log()
            self.test_sin_and_cos()
            self.test_tan()
            self.test_sinh_and_cosh()
            self.test_tanh()
            self.test_inverse_functions()
            self.test_abs()

    def test_exponentiation(self):
        from pyaudi import gdual_double as gdual

        x = gdual(0., "x", 3)
        y = gdual(0., "y", 3)

        p1 = x * x * y + x * y * x * x * x - 3 * y * y * y * y * x * y * x + 3.2
        self.assertEqual(p1**3, p1 * p1 * p1)
        self.assertEqual(p1**3., p1 * p1 * p1)
        self.assertTrue((p1**gdual(3) - p1 * p1 * p1).is_zero(1e-12))
        self.assertTrue((p1**gdual(3.1) - p1**3.1).is_zero(1e-12))

        x = gdual(0., "x", 3)
        y = gdual(0., "y", 3)
        p1 = x + y - 3 * x * y + y * y
        p2 = p1 - 3.5
        self.assertTrue((3**gdual(3.2) - 3**3.2).is_zero(1e-12))
        self.assertTrue((p2**3 - p2**3.).is_zero(1e-12))
        self.assertTrue((p2**-1 - 1 / p2).is_zero(1e-12))
        self.assertTrue((p2**-1. - 1 / p2).is_zero(1e-12))
        self.assertTrue(((p1 + 3.5)**gdual(-1.1) -
                         (p1 + 3.5)**-1.1).is_zero(1e-12))

    def test_sqrt(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sqrt
        x = gdual(2.3, "x", 3)
        y = gdual(1.5, "y", 3)

        p1 = x * x * y - x * y * x * x * x + 3 * \
            y * y * y * y * x * y * x  # positive p0
        p2 = x * x * y - x * y * x * x * x - 3 * y * \
            y * y * y * x * y * x  # negative coefficient
        self.assertTrue((sqrt(p1) * sqrt(p1) - p1).is_zero(1e-12))
        self.assertTrue((sqrt(p2) * sqrt(p2) - p2).is_zero(1e-12))

    def test_cbrt(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import cbrt
        x = gdual(2.3, "x", 3)
        y = gdual(1.5, "y", 3)

        p1 = x * x * y - x * y * x * x * x + 3 * \
            y * y * y * y * x * y * x  # positive p0
        p2 = x * x * y - x * y * x * x * x - 3 * y * \
            y * y * y * x * y * x  # negative coefficient
        self.assertTrue((cbrt(p1) * cbrt(p1) * cbrt(p1) - p1).is_zero(1e-12))
        self.assertTrue((cbrt(p2) * cbrt(p2) * cbrt(p2) - p2).is_zero(1e-12))

    def test_exp_log(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import exp, log
        x = gdual(2.3, "x", 5)
        y = gdual(1.5, "y", 5)

        p1 = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x
        self.assertTrue((exp(log(p1)) - p1).is_zero(1e-10))

    def test_sin_and_cos(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sin, cos, sin_and_cos
        x = gdual(2.3, "x", 8)
        y = gdual(1.5, "y", 8)

        p1 = x + y

        self.assertTrue((sin(2 * p1) - 2 * sin(p1) * cos(p1)).is_zero(1e-12))
        self.assertTrue((cos(2 * p1) - 1 + 2 * sin(p1)
                         * sin(p1)).is_zero(1e-12))
        self.assertTrue((cos(2 * p1) + 1 - 2 * cos(p1)
                         * cos(p1)).is_zero(1e-12))
        self.assertTrue((cos(2 * p1) - cos(p1) * cos(p1) +
                         sin(p1) * sin(p1)).is_zero(1e-12))
        self.assertTrue((sin(p1) * sin(p1) + cos(p1)
                         * cos(p1) - 1).is_zero(1e-12))

        res = sin_and_cos(p1)
        self.assertTrue((res[0] - sin(p1)).is_zero(1e-12))
        self.assertTrue((res[1] - cos(p1)).is_zero(1e-12))

    def test_tan(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sin, cos, tan

        x = gdual(2.3, "x", 10)
        y = gdual(1.5, "y", 10)

        p1 = x + y
        self.assertTrue((tan(p1) - sin(p1) / cos(p1)).is_zero(1e-10))

    def test_sinh_and_cosh(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sinh, cosh, sinh_and_cosh
        x = gdual(2.3, "x", 8)
        y = gdual(1.5, "y", 8)

        p1 = x + y

        self.assertTrue((sinh(2 * p1) - 2 * sinh(p1)
                         * cosh(p1)).is_zero(1e-12))
        self.assertTrue((cosh(2 * p1) - 1 - 2 * sinh(p1)
                         * sinh(p1)).is_zero(1e-12))
        self.assertTrue((cosh(2 * p1) + 1 - 2 * cosh(p1)
                         * cosh(p1)).is_zero(1e-12))
        self.assertTrue((cosh(2 * p1) - cosh(p1) * cosh(p1) -
                         sinh(p1) * sinh(p1)).is_zero(1e-12))
        self.assertTrue((- sinh(p1) * sinh(p1) + cosh(p1)
                         * cosh(p1) - 1).is_zero(1e-12))

        res = sinh_and_cosh(p1)
        self.assertTrue((res[0] - sinh(p1)).is_zero(1e-12))
        self.assertTrue((res[1] - cosh(p1)).is_zero(1e-12))

    def test_tanh(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sinh, cosh, tanh

        x = gdual(2.3, "x", 10)
        y = gdual(1.5, "y", 10)

        p1 = x + y
        self.assertTrue((tanh(p1) - sinh(p1) / cosh(p1)).is_zero(1e-12))

    def test_inverse_functions(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sinh, cosh, tanh
        from pyaudi import asinh, acosh, atanh
        from pyaudi import sin, cos, tan
        from pyaudi import asin, acos, atan

        x = gdual(1.1, "x", 6)
        y = gdual(1.2, "y", 6)
        p1 = 1. / (x + y)

        self.assertTrue((cos(acos(p1)) - p1).is_zero(1e-12))
        self.assertTrue((acos(cos(p1)) - p1).is_zero(1e-12))

        self.assertTrue((sin(asin(p1)) - p1).is_zero(1e-12))
        self.assertTrue((asin(sin(p1)) - p1).is_zero(1e-12))

        self.assertTrue((tan(atan(p1)) - p1).is_zero(1e-12))
        self.assertTrue((atan(tan(p1)) - p1).is_zero(1e-12))

        self.assertTrue((cosh(acosh(p1)) - p1).is_zero(1e-12))
        self.assertTrue((acosh(cosh(p1)) - p1).is_zero(1e-12))

        self.assertTrue((sinh(asinh(p1)) - p1).is_zero(1e-12))
        self.assertTrue((asinh(sinh(p1)) - p1).is_zero(1e-12))

        self.assertTrue((tanh(atanh(p1)) - p1).is_zero(1e-12))
        self.assertTrue((atanh(tanh(p1)) - p1).is_zero(1e-12))

    def test_abs(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import abs as gd_abs

        x = gdual(1.1, "x", 6)
        y = gdual(1.2, "y", 6)

        self.assertEqual(x + y, gd_abs(x + y))
        self.assertEqual(-(x + y), - gd_abs(x + y))


class test_gdual_vdouble(_ut.TestCase):

    def runTest(self):
            self.test_construction()
            self.test_consistency()
            self.test_consistency_functions()
            self.test_trim()
            self.test_extract_terms()
            self.test_subs()

    def test_construction(self):
        from pyaudi import gdual_vdouble as gdual
        x = gdual([2., 3.], "x", 3)
        c = gdual([-2., -3])
        with self.assertRaises(ValueError):
            x = gdual([2., 3.], "dx", 3)
        self.assertEqual(x.constant_cf, [2, 3])
        self.assertEqual(x.find_cf([0]), [2, 3])
        self.assertEqual(c.constant_cf, [-2, -3])
        with self.assertRaises(ValueError):
            x.find_cf([1, 3])
            c.find_cf([1, 3])
            c.find_cf([0])  # no symbols in constants
        with self.assertRaises(TypeError):
            c.find_cf(["x"])
        with self.assertRaises(TypeError):
            x.find_cf(5)
        self.assertEqual(x.degree, 1)
        self.assertEqual(x.order, 3)
        self.assertEqual(x.symbol_set, ["x"])
        self.assertEqual(x.symbol_set_size, 1)
        self.assertEqual(c.degree, 0)
        self.assertEqual(c.order, 0)
        self.assertEqual(c.symbol_set, [])
        self.assertEqual(c.symbol_set_size, 0)
        with self.assertRaises(TypeError):
            gdual(1, "x", 2)
            gdual(1)

    def test_consistency(self):
        from pyaudi import gdual_vdouble as gdual_v
        from pyaudi import gdual_double as gdual_d
        x1 = gdual_d(1, "x", 4)
        y1 = gdual_d(0.1, "y", 4)
        z1 = gdual_d(-0.3, "z", 4)
        f1 = some_complex_rational_f(x1, y1, z1)
        x2 = gdual_d(-0.5, "x", 4)
        y2 = gdual_d(0.03, "y", 4)
        z2 = gdual_d(0.23, "z", 4)
        f2 = some_complex_rational_f(x2, y2, z2)
        x3 = gdual_d(0.441413, "x", 4)
        y3 = gdual_d(-0.2341243241, "y", 4)
        z3 = gdual_d(0.2421413, "z", 4)
        f3 = some_complex_rational_f(x3, y3, z3)
        xv = gdual_v([1, -0.5, 0.441413], "x", 4)
        yv = gdual_v([0.1, 0.03, -0.2341243241], "y", 4)
        zv = gdual_v([-0.3, 0.23, 0.2421413], "z", 4)
        fv = some_complex_rational_f(xv, yv, zv)
        for dx in range(0, 4):
            for dy in range(0, 4):
                for dz in range(0, 4):
                    if dx + dy + dz <= 4:
                        self.assertAlmostEqual(f1.get_derivative(
                            [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[0], delta=1e-12)
                        self.assertAlmostEqual(f2.get_derivative(
                            [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[1], delta=1e-12)
                        self.assertAlmostEqual(f3.get_derivative(
                            [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[2], delta=1e-12)

    def test_consistency_functions(self):
        from pyaudi import gdual_vdouble as gdual_v
        from pyaudi import gdual_double as gdual_d
        from numpy import isnan, nan
        x1 = gdual_d(0.23, "x", 4)
        y1 = gdual_d(0.1, "y", 4)
        z1 = gdual_d(-0.3, "z", 4)
        f1 = some_complex_irrational_f(x1, y1, z1)
        x2 = gdual_d(-0.5, "x", 4)
        y2 = gdual_d(0.03, "y", 4)
        z2 = gdual_d(0.23, "z", 4)
        f2 = some_complex_irrational_f(x2, y2, z2)
        x3 = gdual_d(0.441413, "x", 4)
        y3 = gdual_d(-0.2341243241, "y", 4)
        z3 = gdual_d(0.2421413, "z", 4)
        f3 = some_complex_irrational_f(x3, y3, z3)
        xv = gdual_v([0.23, -0.5, 0.441413], "x", 4)
        yv = gdual_v([0.1, 0.03, -0.2341243241], "y", 4)
        zv = gdual_v([-0.3, 0.23, 0.2421413], "z", 4)
        fv = some_complex_irrational_f(xv, yv, zv)
        for dx in range(0, 4):
            for dy in range(0, 4):
                for dz in range(0, 4):
                    if dx + dy + dz <= 4:
                        if isnan(f1.get_derivative([dx, dy, dz])):
                            self.assertTrue(
                                isnan(fv.get_derivative([dx, dy, dz])[0]))
                        else:
                            self.assertAlmostEqual(f1.get_derivative(
                                [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[0], delta=1e-12)
                        if isnan(f2.get_derivative([dx, dy, dz])):
                            self.assertTrue(
                                isnan(fv.get_derivative([dx, dy, dz])[1]))
                        else:
                            self.assertAlmostEqual(f2.get_derivative(
                                [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[1], delta=1e-12)
                        if isnan(f3.get_derivative([dx, dy, dz])):
                            self.assertTrue(
                                isnan(fv.get_derivative([dx, dy, dz])[2]))
                        else:
                            self.assertAlmostEqual(f3.get_derivative(
                                [dx, dy, dz]), fv.get_derivative([dx, dy, dz])[2], delta=1e-12)

    def test_trim(self):
        from pyaudi import gdual_vdouble as gdual
        x = gdual([1e-4, -1e-4], "x", 1)
        self.assertTrue(x == x.trim(1e-5))
        self.assertTrue(x.trim(1e-3) == gdual([0, 0],"x",1))

    def test_extract_terms(self):
        from pyaudi import gdual_vdouble as gdual
        from pyaudi import sin
        x = gdual([1e-4, -1e-4], "x", 3)
        sinx = sin(x)
        self.assertTrue(sinx == sinx.extract_terms(0) + sinx.extract_terms(1)+ sinx.extract_terms(2)+ sinx.extract_terms(3))
        self.assertTrue(x.trim(1e-3) == gdual([0, 0],"x",1))

    def test_subs(self):
        from pyaudi import gdual_vdouble as gdual
        x = gdual([1], "x", 1)
        y = gdual([-1], "y", 1)
        f = x * y * x / (x - y)  # -0.75*dx-0.5+0.25*dy
        res2 = f.subs("dx", [1])
        res3 = f.subs("dy", [1])
        self.assertEqual(res2.constant_cf, [-1.25])
        self.assertEqual(res2.constant_cf, [-1.25])
        self.assertEqual(res2.get_derivative({"dy": 1}), [0.25])
        self.assertEqual(res2.get_derivative({"dx": 1}), [0.])
        self.assertEqual(res3.constant_cf, [-0.25])
        self.assertEqual(res3.constant_cf, [-0.25])
        self.assertEqual(res3.get_derivative({"dx": 1}), [-0.75])
        self.assertEqual(res3.get_derivative({"dy": 1}), [0.])

        x = gdual([1], "x", 4)
        x = x**2
        y = gdual([0.], "y", 1)
        newx = x.subs("dx", y)
        self.assertEqual(newx, gdual([1], "y", 4)*gdual([1], "y", 4))
        self.assertEqual(newx.order, 4)


class test_utilities(_ut.TestCase):

    def runTest(self):
            self.test_map_inversion()

    def test_map_inversion(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import sin, exp, invert_map, cos

        x = gdual(0., "x", 11)
        f = 1./(1+exp(sin(x)+1./(x+1))) + x
        g = invert_map([f], False)[0]
        dx = g.evaluate({"dp0": 0.01})
        newf = f.evaluate({"dx": dx})
        self.assertAlmostEqual(newf, f.constant_cf + 0.01, delta=1e-10)

        x = gdual(0, "x", 4)
        y = gdual(0, "y", 4)
        f0 = 1./(1+exp(sin(x*y)+1./(x+1))) + x - y
        f1 = 1./(1+exp(cos(x*y)+1./(y+1))) + x + y
        g0, g1 = invert_map([f0, f1], False)
        dx = g0.evaluate({"dp0": 0.01, "dp1":-0.02})
        dy = g1.evaluate({"dp0": 0.01, "dp1":-0.02})
        newf0 = f0.evaluate({"dx": dx, "dy": dy}) 
        newf1 = f1.evaluate({"dx": dx, "dy": dy})
        self.assertAlmostEqual(newf0, f0.constant_cf + 0.01, delta=1e-6)
        self.assertAlmostEqual(newf1, f1.constant_cf - 0.02, delta=1e-6)

        # We test the API
        g0, g1 = invert_map(map = [f0, f1], verbose = False)
        g0, g1 = invert_map(map = [f0, f1])
        g0, g1 = invert_map([f0, f1])

class test_numpy_integration(_ut.TestCase):

    def runTest(self):
            self.test_function_methods()
            self.test_numpy_calls()

    def test_function_methods(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import exp, log, sin ,cos
        x = gdual(0.5, "x", 11)
        self.assertEqual(exp(x), x.exp())
        self.assertEqual(log(x), x.log())
        self.assertEqual(sin(x), x.sin())
        self.assertEqual(cos(x), x.cos())

    def test_numpy_calls(self):
        from pyaudi import gdual_double as gdual
        from pyaudi import exp, log, sin ,cos
        import numpy as np
        x = gdual(0.5, "x", 11)
        self.assertEqual(np.exp(x), x.exp())
        self.assertEqual(np.log(x), x.log())
        self.assertEqual(np.sin(x), x.sin())
        self.assertEqual(np.cos(x), x.cos())


class test_taylor_model(_ut.TestCase):

    def runTest(self):
        self.test_int_d()
        self.test_comparison_methods()
        self.test_validity_check()
        self.test_construction_and_getters_univariate()
        self.test_construction_and_getters_multivariate()
        self.test_construction_and_getters_identity()
        self.test_comparison_of_construction_order()
        self.test_get_bounds()
        self.test_addition()
        self.test_subtraction()
        self.test_multiplication()
        self.test_division()
        self.test_division_2()
        self.test_power()
        self.test_makino1998_simpleexample()

    def test_int_d(self):
        from pyaudi import int_d

        # construct interval
        interval = int_d(1.0, 3.0)
        self.assertEqual(interval.lower, 1.0)
        self.assertEqual(interval.upper, 3.0)

    def test_comparison_methods(self):
        from pyaudi import taylor_model, int_d

        exp1 = {"x": 0.0, "y": 1.0}
        exp2 = {"x": 0.0}
        self.assertFalse(taylor_model.map_equal(exp1, exp2))

        exp_2_1 = {"x": 0.0, "y": 1.0}
        exp_2_2 = {"x": 0.0, "y": 1.0}
        self.assertTrue(taylor_model.map_equal(exp_2_1, exp_2_2))

        dom1 = {"x": int_d(2.0, 3.0), "y": int_d(2.0, 3.0)}
        dom2 = {"x": int_d(2.0, 3.0)}
        self.assertFalse(taylor_model.map_interval_equal(dom1, dom2))

        dom_2_1 = {"x": int_d(2.0, 3.0), "y": int_d(2.0, 3.0)}
        dom_2_2 = {"x": int_d(2.0, 3.0), "y": int_d(2.0, 3.0)}
        self.assertTrue(taylor_model.map_interval_equal(dom_2_1, dom_2_2))

    def test_validity_check(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        x = gdual(0.0, "x", order)
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0}
        dom = {"y": int_d(0.0, 1.0)}
        with self.assertRaises(ValueError):  # maps to std::invalid_argument
            taylor_model(x, rem, exp, dom)

    def test_construction_and_getters_univariate(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        x = gdual(exp["x"], "x", order)
        tm_x = taylor_model(x, rem, exp, dom)

        self.assertEqual(tm_x.tpol, x)
        self.assertEqual(tm_x.ndim, 1)
        self.assertEqual(tm_x.order, order)
        self.assertTrue(taylor_model.interval_equal(tm_x.rem_bound, rem))

        self.assertEqual(len(tm_x.exp_point), len(exp))
        self.assertTrue(taylor_model.map_equal(tm_x.exp_point, exp))
        self.assertTrue(taylor_model.map_interval_equal(tm_x.domain, dom))

    def test_construction_and_getters_multivariate(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0, "y": 1.0}
        dom = {"x": int_d(0.0, 1.0), "y": int_d(1.0, 2.0)}
        x = gdual(exp["x"], "x", order)
        y = gdual(exp["y"], "y", order)
        f_xy = x * x + y + 20
        tm_xy = taylor_model(f_xy, rem, exp, dom)

        self.assertEqual(tm_xy.tpol, f_xy)
        self.assertEqual(tm_xy.ndim, 2)
        self.assertEqual(tm_xy.order, order)
        self.assertTrue(taylor_model.interval_equal(tm_xy.rem_bound, rem))

        self.assertEqual(len(tm_xy.exp_point), len(exp))
        self.assertTrue(taylor_model.map_equal(tm_xy.exp_point, exp))
        self.assertTrue(taylor_model.map_interval_equal(tm_xy.domain, dom))

    def test_construction_and_getters_identity(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        tm_id = taylor_model.identity()
        self.assertEqual(tm_id.tpol, gdual(1.0))
        self.assertEqual(tm_id.ndim, 0)
        self.assertEqual(tm_id.order, 0)
        self.assertTrue(taylor_model.interval_equal(tm_id.rem_bound, int_d(0.0, 0.0)))

        order = 5
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0, "y": 1.0}
        dom = {"x": int_d(0.0, 1.0), "y": int_d(1.0, 2.0)}
        x = gdual(exp["x"], "x", order)
        y = gdual(exp["y"], "y", order)
        f_xy = x * x + y + 20
        tm_xy = taylor_model(f_xy, rem, exp, dom)

        tm_id_2 = taylor_model.identity(rem, exp, dom)
        self.assertEqual(tm_id_2.tpol, gdual(1.0))
        self.assertEqual(tm_id_2.ndim, 0)
        self.assertEqual(tm_id_2.order, 0)
        self.assertTrue(taylor_model.interval_equal(tm_id_2.rem_bound, rem))
        self.assertTrue(taylor_model.map_equal(tm_id_2.exp_point, exp))
        self.assertTrue(taylor_model.map_interval_equal(tm_id_2.domain, dom))

    def test_comparison_of_construction_order(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        x = gdual(exp["x"], "x", order)
        tm_x = taylor_model(x, rem, exp, dom)

        gx = 2
        prod_ans = tm_x + gx
        exp_x = x + 2
        exp_ans = taylor_model(exp_x, rem, exp, dom)

        self.assertEqual(prod_ans.tpol, exp_ans.tpol)
        self.assertTrue(taylor_model.map_equal(prod_ans.exp_point, exp_ans.exp_point))
        self.assertTrue(taylor_model.map_interval_equal(prod_ans.domain, exp_ans.domain))
        self.assertEqual(prod_ans, exp_ans)

    def test_get_bounds(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0, "y": 1.0}
        dom = {"x": int_d(0.0, 1.0), "y": int_d(1.0, 2.0)}
        x = gdual(exp["x"], "x", order)
        y = gdual(exp["y"], "y", order)

        f_xy = x * x + y + 20
        tm_xy = taylor_model(f_xy, rem, exp, dom)

        self.assertTrue(taylor_model.interval_equal(tm_xy.get_bounds(), int_d(21.0, 23.0)))

    def test_addition(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        x = gdual(0.0, "x", order)
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        tm_fx = taylor_model(x, rem, exp, dom)

        gx = 2
        prod_ans = tm_fx + gx
        exp_x = x + 2
        exp_ans = taylor_model(exp_x, rem, exp, dom)

        self.assertEqual(prod_ans, exp_ans)

    def test_subtraction(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        x = gdual(0.0, "x", order)
        rem = int_d(0.0, 2.0)
        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        tm_fx = taylor_model(x, rem, exp, dom)

        gx = -5
        prod_ans = tm_fx + gx   # equivalent to tm_fx - 5

        exp_x = x - 5
        exp_ans = taylor_model(exp_x, rem, exp, dom)

        self.assertEqual(prod_ans, exp_ans)
        
    def test_multiplication(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}

        x = gdual(exp["x"], "x", 5)
        f_x = x

        tpol_range = int_d(0.0, 2.0)
        tm_fx = taylor_model(f_x, tpol_range, exp, dom)

        # multiplication by int
        g_x = 2
        prod_ans = tm_fx * g_x
        exp_ans = taylor_model(f_x * 2.0, int_d(0.0, 4.0), exp, dom)
        self.assertEqual(prod_ans, exp_ans)

        # right-hand side multiplication
        prod_ans = g_x * tm_fx
        self.assertEqual(prod_ans, exp_ans)

        # multiplication by double
        h_x = 2.0
        prod_ans = tm_fx * h_x
        exp_ans = taylor_model(f_x * 2.0, int_d(0.0, 4.0), exp, dom)
        self.assertEqual(prod_ans, exp_ans)

        # multiplication by interval
        i_x = int_d(0.0, 2.0)
        prod_ans = tm_fx * i_x
        exp_ans = taylor_model(f_x, int_d(0.0, 8.0), exp, dom)  # [0,8]
        self.assertEqual(prod_ans, exp_ans)

        # multiplication by another TaylorModel
        prod_ans = tm_fx * tm_fx
        exp_ans = taylor_model(f_x * f_x, int_d(0.0, 8.0), exp, dom)
        self.assertEqual(prod_ans, exp_ans)

    def test_division(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        exp = {"x": 1.0}
        dom = {"x": int_d(0.5, 1.5)}
        x = gdual(exp["x"], "x", 5)

        rem = int_d(0.0, 0.0)
        tm_fx = taylor_model(x, rem, exp, dom)

        prod_ans = 1 / tm_fx
        exp_ans = taylor_model(1 / x, int_d(0.0, 2.0), {"x": 1.0}, {"x": int_d(0.5, 1.5)})
        self.assertEqual(prod_ans, exp_ans)

    def test_division_2(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        x = gdual(exp["x"], "x", 5)

        rem = int_d(0.0, 2.0)
        tm_fx = taylor_model(x, rem, exp, dom)

        prod_ans = tm_fx / 2
        exp_ans = taylor_model(x / 2, int_d(0.0, 1.0), {"x": 0.0}, {"x": int_d(0.0, 1.0)})
        self.assertEqual(prod_ans, exp_ans)

        prod_ans_2 = tm_fx / 2.0
        exp_ans_2 = taylor_model(x / 2.0, int_d(0.0, 1.0), {"x": 0.0}, {"x": int_d(0.0, 1.0)})
        self.assertEqual(prod_ans_2, exp_ans_2)

    def test_power(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        order = 5
        rem = int_d(0.0, 0.0)
        exp = {"x": 0.0}
        dom = {"x": int_d(0.0, 1.0)}
        x = gdual(exp["x"], "x", order)

        fx = x**3
        tm_fx_exp = taylor_model(fx, rem, exp, dom)
        tm_x = taylor_model(x, rem, exp, dom)

        tm_fx = tm_x**3
        self.assertEqual(tm_fx, tm_fx_exp)

    def test_makino1998_simpleexample(self):
        from pyaudi import taylor_model, int_d
        from pyaudi import gdual_double as gdual

        rem_bound = int_d(0.0, 0.0)
        exp_points = {"x": 2.0}
        domain = {"x": int_d(1.9, 2.1)}

        # Expected answers (Table 4.2 p.95)
        exp_ans = [
            (0.0, 1.4579384e-3), (-7.6733603e-5, 7.6733603e-5), (0.0, 4.0386107e-6),
            (-2.1255845e-7, 2.1255845e-7), (0.0, 1.1187287e-8), (-5.8880459e-10, 5.8880459e-10),
            (0.0, 3.0989715e-11), (-1.6310376e-12, 1.6310376e-12), (0.0, 8.5844087e-14),
            (-4.5181098e-15, 4.5181098e-15), (0.0, 2.3779525e-16), (-1.2515539e-17, 1.2515539e-17),
            (0.0, 6.5871262e-19), (-3.4669086e-20, 3.4669085e-20), (0.0, 1.8246887e-21)
        ]

        # Loop over orders
        for order in range(1, 16):
            x = gdual(exp_points["x"], "x", order)
            const_fx = taylor_model(x, rem_bound, exp_points, domain)
            const_fx_2 = 1 / const_fx
            const_fx_3 = const_fx_2 + const_fx

            self.assertAlmostEqual(const_fx_3.rem_bound.lower, exp_ans[order - 1][0], delta=1e-7)
            self.assertAlmostEqual(const_fx_3.rem_bound.upper, exp_ans[order - 1][1], delta=1e-7)

        # Construction test without build-up
        x = gdual(exp_points["x"], "x", 3)
        tm_fx = taylor_model(x, rem_bound, exp_points, domain)
        tm_fx_2 = 1 / tm_fx + tm_fx

        exp_points_2 = {"y": 2.0}
        domain_2 = {"y": int_d(1.9, 2.1)}
        y = gdual(exp_points_2["y"], "y", 3)
        const_gy = taylor_model(y, rem_bound, exp_points_2, domain_2)
        const_gy_2 = 1 / const_gy
        const_gy_3 = const_gy_2 + const_gy

        self.assertEqual(float(tm_fx_2.rem_bound.lower), float(const_gy_3.rem_bound.lower))
        self.assertEqual(float(tm_fx_2.rem_bound.upper), float(const_gy_3.rem_bound.upper))

        # Verify bound values
        exp_bounds = int_d(2.42631, 2.57618)
        self.assertAlmostEqual(tm_fx_2.get_bounds().lower, exp_bounds.lower, delta=1e-5)
        self.assertAlmostEqual(tm_fx_2.get_bounds().upper, exp_bounds.upper, delta=1e-5)

        # Verify bound interval (Taylor model order)
        exp_taymodorder_ans = [0.15145793, 0.15015346, 0.14987903, 0.14987542, 0.14987469, 0.14987468]

        for order in range(5, 7):
            it = order - 1
            x = gdual(exp_points["x"], "x", order)
            T_x = taylor_model(x, rem_bound, exp_points, domain)
            T_fx = 1 / T_x + T_x

            interval_size = float(T_fx.get_bounds().upper) - float(T_fx.get_bounds().lower)
            self.assertAlmostEqual(interval_size, exp_taymodorder_ans[it], delta=1e-7)

class test_taylor_model_function_calls(_ut.TestCase):

    def runTest(self):
            self.test_exponentiation()
            self.test_sqrt()
            self.test_exp_log()
            self.test_sin_and_cos()
            self.test_tan()
            self.test_sinh_and_cosh()
            self.test_tanh()
            self.test_inverse_functions()
            self.test_abs()

    def test_exponentiation(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d

        exp_points = {"x": 2.0, "y": 2.0}
        dom = {"x": int_d(1.0, 3.0), "y": int_d(1.0, 3.0)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 3), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 3), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        test1 = y * x * x * x
        temp = x * test1
        p1 = x * x * y + x * y * x * x * x + 3 * y * y * y * y * x * y * x + 3.2
        self.assertEqual(p1**3, p1 * p1 * p1)
        self.assertEqual(p1**3., p1 * p1 * p1)
        self.assertTrue((p1**3 - p1 * p1 * p1).tpol == gdual(0.0, "irrelevant", 0))
        self.assertTrue((p1**3. - p1 * p1 * p1).tpol == gdual(0.0, "irrelevant", 0))

        p2 = x + y + 3 * x * y + y * y + 3.5
        self.assertTrue((p2**3 - p2**3.).tpol == gdual(0.0, "irrelevant", 0))
        self.assertTrue((p2**-1 - 1 / p2).tpol == gdual(0.0, "irrelevant", 0))

    def test_sqrt(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, sqrt

        domain_size = 0.2
        exp_points = {"x": 0.0, "y": 0.5}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 3), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 3), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x * x * y + x * y * x * x * x + 3 * y * y * y * y * x * y * x + 1
        p2 = x * x * y + x * y * x * x * x - 3 * y * y * y * y * x * y * x + 1

        test1 = sqrt(p1) * sqrt(p1)
        test2 = sqrt(p2) * sqrt(p2)
        self.assertTrue((test1.tpol - p1.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0) and taylor_model.map_equal(test1.exp_point, p1.exp_point) and taylor_model.map_interval_equal(test1.domain, p1.domain))
        self.assertTrue((test1.tpol - p1.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0) and taylor_model.map_equal(test2.exp_point, p2.exp_point) and taylor_model.map_interval_equal(test2.domain, p2.domain))

    def test_exp_log(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, exp, log

        domain_size = 0.2
        exp_points = {"x": 0.3, "y": 0.2}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 3), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 3), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x * x * y - x * y * x * x * x + 3 * y * y * y * y * x * y * x + 1

        test1 = exp(log(p1))
        self.assertTrue((test1.tpol - p1.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0) and taylor_model.map_equal(test1.exp_point, p1.exp_point) and taylor_model.map_interval_equal(test1.domain, p1.domain))

    def test_sin_and_cos(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, sin, cos, sin_and_cos

        exp_points = {"x": 2.3, "y": 1.5}
        dom = {"x": int_d(2.0, 3.0), "y": int_d(1.0, 2.0)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 8), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 8), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x + y

        test_pairs = [
            (sin(2 * p1), 2 * sin(p1) * cos(p1)),
            (cos(2 * p1), 1 - 2 * sin(p1) * sin(p1)),
            (cos(2 * p1), 2 * cos(p1) * cos(p1) - 1),
            (cos(2 * p1), cos(p1) * cos(p1) - sin(p1) * sin(p1)),
            (sin(p1) * sin(p1) + cos(p1) * cos(p1), taylor_model.identity(rem, exp_points, dom)),
        ]

        for test1, test2 in test_pairs:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-15) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

        res = sin_and_cos(p1)
        test_pairs_res = [(res[0], sin(p1)), (res[1], cos(p1))]
        for test1, test2 in test_pairs_res:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

    def test_tan(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, sin, cos, tan

        domain_size = 0.2
        exp_points = {"x": 2.3, "y": 1.5}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 10), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 10), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x + y
        test1, test2 = tan(p1), sin(p1) / cos(p1)
        self.assertTrue(
            (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
            and taylor_model.map_equal(test1.exp_point, test2.exp_point)
            and taylor_model.map_interval_equal(test1.domain, test2.domain)
        )

    def test_sinh_and_cosh(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, sinh, cosh, sinh_and_cosh

        domain_size = 0.2
        exp_points = {"x": 2.3, "y": 1.5}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 8), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 8), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x + y

        test_pairs = [
            (sinh(2 * p1), 2 * sinh(p1) * cosh(p1)),
            (cosh(2 * p1), 1 + 2 * sinh(p1) * sinh(p1)),
            (cosh(2 * p1), 2 * cosh(p1) * cosh(p1) - 1),
            (cosh(2 * p1), cosh(p1) * cosh(p1) + sinh(p1) * sinh(p1)),
            (cosh(p1) * cosh(p1) - sinh(p1) * sinh(p1), taylor_model.identity(rem, exp_points, dom)),
        ]

        # Lower precicision than expected: 1e-12 instead of 1e-16
        for test1, test2 in test_pairs:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-12) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

        res = sinh_and_cosh(p1)
        test_pairs_res = [(res[0], sinh(p1)), (res[1], cosh(p1))]
        for test1, test2 in test_pairs_res:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

    def test_tanh(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d, sinh, cosh, tanh

        domain_size = 0.2
        exp_points = {"x": 2.3, "y": 1.5}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 10), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 10), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = x + y
        test1, test2 = tanh(p1), sinh(p1) / cosh(p1)
        self.assertTrue(
            (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
            and taylor_model.map_equal(test1.exp_point, test2.exp_point)
            and taylor_model.map_interval_equal(test1.domain, test2.domain)
        )
    def test_inverse_functions(self):
        from pyaudi import (
            taylor_model,
            gdual_double as gdual,
            int_d,
            sinh, cosh, tanh,
            asinh, acosh, atanh,
            sin, cos, tan,
            asin, acos, atan,
        )

        domain_size = 0.1
        exp_points = {"x": 1.2, "y": 1.3}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 6), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 6), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = 1.0 / (x + y)

        test_pairs = [
            (cos(acos(p1)), p1),
            (acos(cos(p1)), p1),
            (sin(asin(p1)), p1),
            (asin(sin(p1)), p1),
            (tan(atan(p1)), p1),
            (atan(tan(p1)), p1),
        ]
        for test1, test2 in test_pairs:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-15) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

        exp_points = {"x": 1.1, "y": 1.2}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 6), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 6), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = 1.0 / (x + y) + 10.0
        test_pairs_2 = [
            (cosh(acosh(p1)), p1),
            (acosh(cosh(p1)), p1),
            (sinh(asinh(p1)), p1),
            (asinh(sinh(p1)), p1),
        ]
        for test1, test2 in test_pairs_2:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-14) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

        exp_points = {"x": 2.0, "y": 2.0}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 6), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 6), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p1 = 1.0 / (x + y)
        test_pairs_3 = [
            (tanh(atanh(p1)), p1),
            (atanh(tanh(p1)), p1),
        ]

        for test1, test2 in test_pairs_3:
            self.assertTrue(
                (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
                and taylor_model.map_equal(test1.exp_point, test2.exp_point)
                and taylor_model.map_interval_equal(test1.domain, test2.domain)
            )

    def test_abs(self):
        from pyaudi import taylor_model, gdual_double as gdual, int_d
        from pyaudi import abs as tm_abs

        domain_size = 0.2
        exp_points = {"x": 1.1, "y": 1.2}
        dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
        rem = int_d(0.0, 0.0)

        x = taylor_model(gdual(exp_points["x"], "x", 6), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
        y = taylor_model(gdual(exp_points["y"], "y", 6), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

        p = x + y

        # Check abs(p) ≈ p
        test1, test2 = tm_abs(p), p
        self.assertTrue(
            (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
            and taylor_model.map_equal(test1.exp_point, test2.exp_point)
            and taylor_model.map_interval_equal(test1.domain, test2.domain)
        )

        # Check -abs(p) ≈ -p
        test1, test2 = -tm_abs(p), -p
        self.assertTrue(
            (test1.tpol - test2.tpol).trim(1e-16) == gdual(0.0, "irrelevant", 0)
            and taylor_model.map_equal(test1.exp_point, test2.exp_point)
            and taylor_model.map_interval_equal(test1.domain, test2.domain)
        )

def run_test_suite():
    """Run the full test suite.
    This function must raise an exception if at least one test fails.
    """
    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(test_gdual_double)
    suite.addTest(test_function_calls())
    suite.addTest(test_gdual_vdouble())
    suite.addTest(test_utilities())
    suite.addTest(test_numpy_integration())
    suite.addTest(test_taylor_model())
    suite.addTest(test_taylor_model_function_calls())


    test_result = _ut.TextTestRunner(verbosity=2).run(suite)
    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')

if __name__ == "__main__":
    from pyaudi import (
        taylor_model,
        gdual_double as gdual,
        int_d,
        sinh, cosh, tanh,
        asinh, acosh, atanh,
        sin, cos, tan,
        asin, acos, atan,
    )


    domain_size = 0.01
    exp_points = {"x": 1.1, "y": 1.2}
    dom = {"x": int_d(exp_points["x"] - domain_size, exp_points["x"] + domain_size), "y": int_d(exp_points["y"] - domain_size, exp_points["y"] + domain_size)}
    rem = int_d(0.0, 0.0)

    x = taylor_model(gdual(exp_points["x"], "x", 6), rem, {"x": exp_points["x"]}, {"x": dom["x"]})
    y = taylor_model(gdual(exp_points["y"], "y", 6), rem, {"y": exp_points["y"]}, {"y": dom["y"]})

    p1 = 1.0 / (x + y) + 10.0
    test_pairs_2 = [
        (cosh(acosh(p1)), p1),
        (acosh(cosh(p1)), p1),
        (sinh(asinh(p1)), p1),
        (asinh(sinh(p1)), p1),
    ]
    for test1, test2 in test_pairs_2:
        x = (
            (test1.tpol - test2.tpol).trim(1e-14) == gdual(0.0, "irrelevant", 0)
            and taylor_model.map_equal(test1.exp_point, test2.exp_point)
            and taylor_model.map_interval_equal(test1.domain, test2.domain)
        )
