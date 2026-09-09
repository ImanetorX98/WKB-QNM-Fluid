"""Separate endpoint, quadrature, profile and primitive errors on identical data.

This reconstructs possible mechanisms, not Claude's unavailable original driver.
"""
import json
import numpy as np
from scipy.integrate import simpson, solve_ivp
from vaidya_solvability import radial_mode, tortoise, leaver_qnm, source_integrand
from asymptotic_series import integrand_series, antiderivative_series
from outgoing_asymptotics import series_coefficients


def frobenius_mode(om, offset, points, order=2, rtol=2e-11):
    """Independent integration with local ingoing h=1+h1*d+h2*d^2.

    R=exp(-i*om*r*)h, ell=2, spin=0. Diagnostic only.
    """
    start = 2+offset
    t = np.linspace(float(tortoise(start)), float(tortoise(70.)), points)
    u0, u1 = 7/4, -15/8
    h1 = u0/(0.5-2j*om) if order >= 1 else 0j
    h2 = ((u0+0.5)*h1+u1)/(2-4j*om) if order >= 2 else 0j
    h, hp = 1+h1*offset+h2*offset**2, h1+2*h2*offset
    phase = np.exp(-1j*om*t[0])
    initial = [phase*h, phase*(-1j*om*h+(1-2/start)*hp), start+0j]
    def rhs(_, y):
        R, dR, rc = y
        r = rc.real
        f = 1-2/r
        return [dR, (f*(6/r**2+2/r**3)-om**2)*R, f]
    sol = solve_ivp(rhs, (t[0],t[-1]), initial, t_eval=t,
                    method='DOP853',rtol=rtol,atol=rtol/100)
    if not sol.success:
        raise RuntimeError(sol.message)
    return dict(x=sol.y[2].real, xstar=t, R=sol.y[0], dR=sol.y[1])


def budget(points=100000, lo=40., hi=60., offset=1e-4, horizon_order=None, rtol=2e-11):
    om = leaver_qnm(2, 0, 0)
    data = (radial_mode(2, 0, om, 2+offset, 70., points) if horizon_order is None
            else frobenius_mode(om, offset, points, horizon_order, rtol))
    mask = (data['x'] >= lo) & (data['x'] <= hi)
    x, t, R = (data[k][mask] for k in ('x', 'xstar', 'R'))
    w = integrand_series(2, 0, om, 24)
    v = antiderivative_series(w, om, 28, 24)
    def primitive(r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        return v.evaluate(r)*np.exp(2j*om*tortoise(r))
    predicted = w.evaluate(x)*np.exp(2j*om*tortoise(x))
    delta = np.diff(primitive(x[[0, -1]]))[0]
    nominal = np.diff(primitive([lo, hi]))[0]
    q = simpson(predicted, x=x)
    u = np.polynomial.polynomial.polyval(1/x, series_coefficients(2,0,om,24))
    c = (R*np.exp(-1j*om*t)/u)[-1]
    numeric = source_integrand(data, om, 2)[mask]/c**2
    q_num = simpson(numeric, x=x)
    # Discrete weighted L1 bound, valid even for signed Simpson weights.
    # For this grid the weights are positive; use trapezoid for an independent
    # rigorous discrete inequality that does not require forming Simpson weights.
    trap = np.trapezoid(predicted, x)
    trap_num = np.trapezoid(numeric, x)
    eps = float(np.max(abs(numeric/predicted-1)))
    kappa = float(np.trapezoid(abs(predicted),x)/abs(trap))
    residual = v.derivative().evaluate(x)+(2j*om/(1-2/x))*v.evaluate(x)-w.evaluate(x)
    endpoint_linear = (predicted[-1]*(x[-1]-hi)-predicted[0]*(x[0]-lo))/nominal
    wrong = q/nominal-1
    return {
        'points':points, 'window':[lo,hi], 'horizon_offset':offset,
        'horizon_order':horizon_order, 'rtol':rtol,
        'endpoint_offsets':[float(x[0]-lo),float(x[-1]-hi)],
        'trap_error':float(abs(trap/delta-1)),
        'simpson_error':float(abs(q/delta-1)),
        'wrong_endpoints_error':float(abs(wrong)),
        'wrong_endpoints_complex':[float(wrong.real),float(wrong.imag)],
        'endpoint_linear_prediction':float(abs(endpoint_linear)),
        'endpoint_linear_remainder':float(abs(wrong-endpoint_linear)),
        'max_pointwise_profile_error':eps,
        'quadrature_condition':kappa,
        'discrete_profile_error_bound':eps*kappa,
        'discrete_profile_error':float(abs((trap_num-trap)/trap)),
        'corrected_total_error':float(abs(q_num/delta-1)),
        'tortoise_consistency':float(np.max(abs(t-tortoise(x)))),
        'primitive_differential_residual':float(np.max(abs(residual/w.evaluate(x)))),
    }


if __name__ == '__main__':
    for n in (24000, 100000, 200000):
        for lo,hi in ((40.,50.), (40.,60.), (50.,70.)):
            print(json.dumps(budget(n,lo,hi)))
    for order in (0,1,2):
        for offset in (1e-3,1e-4,1e-5,1e-6):
            print(json.dumps(budget(24000, offset=offset, horizon_order=order)))
    for tol in (2e-9,2e-11,2e-13):
        print(json.dumps(budget(100000, horizon_order=2, rtol=tol)))
