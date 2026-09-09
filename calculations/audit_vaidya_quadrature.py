"""Reproducible exterior-interval audit; NOT a regularized QNM overlap."""
import numpy as np
from scipy.integrate import simpson
from vaidya_solvability import radial_mode, leaver_qnm, source_integrand
from asymptotic_series import integrand_series, antiderivative_series
from outgoing_asymptotics import series_coefficients


def audit(points, lo=40., hi=60.):
    omega = leaver_qnm(2, 0, 0)
    data = radial_mode(2, 0, omega, 2.0001, 70., points)
    r, t, R = (data[k] for k in ('x', 'xstar', 'R'))
    mask = (r >= lo) & (r <= hi)
    x = r[mask]
    w = integrand_series(2, 0, omega, 24)
    v = antiderivative_series(w, omega, 28, 24)
    pred = w.evaluate(x)*np.exp(2j*omega*t[mask])
    boundary = np.diff(v.evaluate(x[[0,-1]])*np.exp(2j*omega*t[mask][[0,-1]]))[0]
    u = np.polynomial.polynomial.polyval(1/x, series_coefficients(2,0,omega,24))
    cs = R[mask]*np.exp(-1j*omega*t[mask])/u
    c = cs[-1]
    numeric = source_integrand(data, omega, 2)[mask]/c**2
    G = np.exp(1j*omega*t)*R
    dG = np.exp(1j*omega*t)*(1j*omega*R+data['dR'])/(1-2/r)
    gradient = (2*np.exp(-2j*omega*t)*G*np.gradient(r*dG,r))[mask]/c**2
    return dict(points=points, dr=float(np.max(np.diff(x))),
        trap_error=abs(np.trapezoid(pred,x)/boundary-1),
        simpson_error=abs(simpson(pred,x=x)/boundary-1),
        ode_error=abs(simpson(numeric,x=x)/boundary-1),
        gradient_error=abs(simpson(gradient,x=x)/boundary-1),
        normalization_spread=float(np.max(abs(cs/c-1))))


if __name__ == '__main__':
    for points in (1000, 2000, 4000, 8000, 24000):
        print(audit(points))
