import unittest
import numpy as np
import sympy as sp
from kerr_madelung_analytic import geometric_term


class GeometricTerm(unittest.TestCase):
    def test_schwarzschild_limit(self):
        r=np.array([2.5,3.,5.,20.])
        np.testing.assert_allclose(geometric_term(r,0),(1-2/r)*2/r**3,rtol=1e-14)

    def test_symbolic_tortoise_derivative(self):
        r,a=sp.symbols('r a',positive=True)
        h=sp.sqrt(r*r+a*a)
        speed=(r*r-2*r+a*a)/(r*r+a*a)
        expression=speed*sp.diff(speed*sp.diff(h,r),r)/h
        ref=sp.lambdify((r,a),expression,'numpy')
        for spin in (.3,.6,.9):
            x=np.array([2.,3.,5.,10.])
            np.testing.assert_allclose(geometric_term(x,spin),ref(x,spin),rtol=1e-13)


if __name__=='__main__':
    unittest.main()
