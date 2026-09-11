import unittest
from vaidya_transport_denominator import compute,outgoing_Domega
from outgoing_asymptotics import log_derivative_omega


class TransportDenominator(unittest.TestCase):
    def test_norm_identity_and_cutoff(self):
        rows=compute(uppers=(40.,50.,60.))
        base=complex(*rows[0]['P'])
        for row in rows:
            self.assertLess(row['identity_error'],1e-7)
            self.assertLess(abs(complex(*row['P'])/base-1),1e-7)

    def test_horizon_width(self):
        a=compute(width=.2,uppers=(40.,))[0]
        b=compute(width=1.,uppers=(40.,))[0]
        self.assertLess(abs(complex(*a['P'])/complex(*b['P'])-1),1e-7)

    def test_analytic_outgoing_frequency_derivative(self):
        w=.483643872210713-.09675877597828786j
        a=outgoing_Domega(2,0,w,60.,20)
        b=log_derivative_omega(2,0,w,60.,step=1e-5)
        self.assertLess(abs(a-b),1e-8)


if __name__=='__main__':
    unittest.main()
