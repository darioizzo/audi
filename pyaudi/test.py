import unittest as _ut

def gdual_norm(p):
    order = p.order
    res = [p.find_cf([x,y]) for x in range(0,order) for y in range(0,order) if (x+y) <= order]
    return sum(res) / len(res)


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
            x.find_cf([1,3])
            c.find_cf([1,3])
            c.find_cf([0]) # no symbols in constants
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
        p1 = gdual(1);
        p2 = gdual(0, "x", 4);
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
        x = gdual(3, "x",2)
        y = gdual(7, "y",2)
        p1 = x + 3*x*y + y*y
        self.assertEqual(p1.find_cf([0,0]), 115)
        self.assertEqual(p1.find_cf([1,0]), 22)
        self.assertEqual(p1.find_cf([0,1]), 23)
        self.assertEqual(p1.find_cf([1,1]), 3)
        self.assertEqual(p1.find_cf([2,0]), 0)
        self.assertEqual(p1.find_cf([0,2]), 1)

    def test_division(self):
        from pyaudi import gdual_double as gdual
        x = gdual(0, "x",2)
        y = gdual(1, "y",2)
        p1 = 1 / (x + 2*x*y + y*y)
        self.assertEqual(p1.find_cf([0,0]), 1)
        self.assertEqual(p1.find_cf([1,0]), -3)
        self.assertEqual(p1.find_cf([0,1]), -2)
        self.assertEqual(p1.find_cf([1,1]), 10)
        self.assertEqual(p1.find_cf([2,0]), 9)
        self.assertEqual(p1.find_cf([0,2]), 3)

    def test_identities(self):
        from pyaudi import gdual_double as gdual
        x = gdual(2, "x",3);
        y = gdual(3, "y",3);

        p1 = x*x+y-x*x*x*x*y-y*x*x
        p2 = y*y-x+y*y*y*y*x-2*x
        self.assertEqual((x + y)*(x + y), x*x + y*y + 2*x*y)
        self.assertEqual((p1 + p2)*(p1 + p2), p1*p1 + p2*p2 + 2*p1*p2)
        self.assertEqual(x*x*x*x-y*y*y*y, (x-y)*(x+y)*(x*x+y*y))
        self.assertEqual(p1*p1*p1*p1-p2*p2*p2*p2, (p1-p2)*(p1+p2)*(p1*p1+p2*p2))
        self.assertAlmostEqual(gdual_norm((p1/p2) * (p2/p1) - 1), 0, 10)
        self.assertAlmostEqual(gdual_norm((p1/p2) * p2 - p1), 0, 10)

    def test_derivatives(self):
        from pyaudi import gdual_double as gdual

        x = gdual(1, "x",4)
        y = gdual(1, "y",4)
        z = gdual(1, "z",4)
        f = x*x*x*x*x + x*y*z*x*x + z*x*y*y*y

        self.assertEqual(f.get_derivative([1,1,1]), 6.)
        self.assertEqual(f.get_derivative([2,1,1]), 6.)
        self.assertEqual(f.get_derivative([1,2,1]), 6.)
        self.assertEqual(f.get_derivative([1,1,2]), 0.)
        self.assertEqual(f.get_derivative([4,0,0]),120.)
        with self.assertRaises(ValueError):
            f.get_derivative([4,1,1])

        x = gdual(1, "x",8)
        y = gdual(1, "y",8)
        f = (x * y + 2 * x * x * y) / (1 + x + y)
        self.assertEqual(f.get_derivative([1,1]), 1.)
        self.assertAlmostEqual(f.get_derivative([2,2]), 0, 10)
        self.assertAlmostEqual(f.get_derivative([3,3]), 0, 10)
        self.assertAlmostEqual(f.get_derivative([4,4]), 0, 10)

        # we test the dictionary interface
        x = gdual(1, "x",4)
        y = gdual(1, "y",4)
        z = gdual(1, "z",4)
        f = x*x*x*x*x + x*y*z*x*x + z*x*y*y*y
        dictionary =  {"dx": 1, "dy": 1, "dz": 1}
        self.assertEqual(f.get_derivative(dictionary), 6.)
        with self.assertRaises(ValueError):
            f.get_derivative([4,1,1])
            dictionary = {"dx": 4, "dy": 1, "dz": 1}
            f.get_derivative(dictionary)
        dictionary = {"dx": 1, "dr": 1, "dz": 1}
        self.assertEqual(f.get_derivative(dictionary), 0.)

    def test_integrate_partial(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = x*x*x + x*y*z + z*x*y
        fx = 3*x*x + y*z + z*y
        fy = x*z + z*x
        fz = y*x + x*y

        self.assertEqual(f.partial("x"), fx)
        self.assertEqual(f.partial("y"), fy)
        self.assertEqual(f.partial("z"), fz)

    def test_subs(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1, "x", 1)
        y = gdual(-1, "y", 1)
        f = x*y*x / (x-y) # -0.75*dx-0.5+0.25*dy
        res2 = f.subs("dx", 1)
        res3 = f.subs("dy", 1)
        self.assertEqual(res2.constant_cf, -1.25)
        self.assertEqual(res2.constant_cf, -1.25)
        self.assertEqual(res2.get_derivative({"dy":1}), 0.25)
        self.assertEqual(res2.get_derivative({"dx":1}), 0.)
        self.assertEqual(res3.constant_cf, -0.25)
        self.assertEqual(res3.constant_cf, -0.25)
        self.assertEqual(res3.get_derivative({"dx":1}), -0.75)
        self.assertEqual(res3.get_derivative({"dy":1}), 0.)

    def test_serialization(self):
        from pyaudi import gdual_double as gdual
        x = gdual(1, "x", 4)
        y = gdual(1, "y", 4)
        z = gdual(1, "z", 4)
        f = (x*x*x + x*y*z + z*x*y)*(x*x*x + x*y*z + z*x*y)*(x*x*x + x*y*z + z*x*y)*(x*x*x + x*y*z + z*x*y)
        new_f = gdual(0)
        new_f.__setstate__(f.__getstate__())
        self.assertEqual(f, new_f)
        self.assertEqual(f.order, new_f.order)

def run_test_suite():
    """Run the full test suite.
    This function will raise an exception if at least one test fails.
    """
    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(test_gdual_double)
    test_result = _ut.TextTestRunner(verbosity=2).run(suite)
