import unittest
import numpy as np
from repro_numerator_discrepancy import setup, interval_diagnostic


class OriginalReproducerRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, cls.r, _, cls.g, cls.F, _ = setup()

    def test_original_number_and_both_endpoint_repairs(self):
        for lo, hi in ((40.,50.), (50.,60.), (60.,70.), (70.,80.)):
            with self.subTest(interval=(lo,hi)):
                row = interval_diagnostic(self.r,self.g,self.F,lo,hi)
                self.assertLess(abs(abs(row['original'])-0.999835), 4e-6)
                self.assertGreater(abs(row['original']-1), 1e-3)
                self.assertLess(abs(row['include_endpoint']-1), 5e-7)
                self.assertLess(abs(row['move_primitive_endpoint']-1), 5e-7)

    def test_endpoint_accounting_on_exact_polynomial(self):
        r = np.linspace(0,2,21)
        row = interval_diagnostic(r,r**2,r**3/3,0.5,1.5)
        self.assertAlmostEqual(row['include_endpoint'],1.,places=13)
        self.assertAlmostEqual(row['move_primitive_endpoint'],1.,places=13)
        self.assertGreater(abs(row['original']-1),0.1)


if __name__ == '__main__':
    unittest.main()
