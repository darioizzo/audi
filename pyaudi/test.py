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
        self.assertTrue((tan(p1) - sin(p1) / cos(p1)).is_zero(1e-12))

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
        with self.assertRaises(AttributeError):
            x.find_cf(5)
        self.assertEqual(x.degree, 1)
        self.assertEqual(x.order, 3)
        self.assertEqual(x.symbol_set, ["x"])
        self.assertEqual(x.symbol_set_size, 1)
        self.assertEqual(c.degree, 0)
        self.assertEqual(c.order, 0)
        self.assertEqual(c.symbol_set, [])
        self.assertEqual(c.symbol_set_size, 0)
        with self.assertRaises(AttributeError):
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

class test_utilities(_ut.TestCase):

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




def run_test_suite():
    """Run the full test suite.
    This function will raise an exception if at least one test fails.
    """
    retval = 0
    suite_gdouble = _ut.TestLoader().loadTestsFromTestCase(test_gdual_double)
    suite_function_calls = _ut.TestLoader().loadTestsFromTestCase(test_function_calls)
    suite_gvdouble = _ut.TestLoader().loadTestsFromTestCase(test_gdual_vdouble)
    suite_utilities = _ut.TestLoader().loadTestsFromTestCase(test_utilities)

    print("\nRunning gdual_double tests")
    test_result = _ut.TextTestRunner(verbosity=2).run(suite_gdouble)
    print("\nRunning function calls tests")
    test_result = _ut.TextTestRunner(verbosity=2).run(suite_function_calls)
    print("\nRunning gdual_vdouble tests")
    test_result = _ut.TextTestRunner(verbosity=2).run(suite_gvdouble)
    print("\nRunning utilities tests")
    test_result = _ut.TextTestRunner(verbosity=2).run(suite_utilities)
