import sys; sys.path.insert(0,'calculations')
import sympy as sp, cmath
from claude_kerr_derivative_validation import spheroidal

spin, ell, m_int = sp.Rational(3,5), 28, 19
omega = sp.Rational(3,1) - sp.Rational(1,4)*sp.I
eig,_,_ = spheroidal(ell, m_int, complex(spin*omega))
lam = sp.nsimplify(complex(eig), rational=False) + spin**2*omega**2 - 2*spin*m_int*omega
rp = 1+sp.sqrt(1-spin**2); rm = 1-sp.sqrt(1-spin**2); b = rm-rp
r,z = sp.symbols('r z')
Delta=(r-rp)*(r-rm); K=(r**2+spin**2)*omega-spin*m_int
Pz = sp.simplify(((2*r-rm-rp)/Delta).subs(r,rp+b*z)*b)
Qz = sp.simplify(((K**2-lam*Delta)/Delta**2).subs(r,rp+b*z)*b**2)

kappa=(float(rp)**2+float(spin)**2)*complex(omega)-float(spin)*m_int
p  = -1j*kappa/float(rp-rm)            # ramo ENTRANTE all'orizzonte
qe = cmath.sqrt(-complex(sp.limit(Qz*(z-1)**2, z, 1)))
s  = -1j*complex(omega)*float(b)       # uscente all'infinito
for name,val in (('p',p),('q_e',qe),('s',s)): print(f'  {name} = {val:+.10f}')

P_,Q_ = sp.nsimplify(Pz,rational=False), sp.nsimplify(Qz,rational=False)
w = z**sp.nsimplify(p,rational=False)*(z-1)**sp.nsimplify(qe,rational=False)*sp.exp(sp.nsimplify(s,rational=False)*z)
Gam = sp.simplify(2*sp.diff(w,z)/w + P_)
Lam = sp.simplify(sp.diff(w,z,2)/w + P_*sp.diff(w,z)/w + Q_)
print('\n  Gamma(z) deve essere  gamma/z + delta/(z-1) + eps :')
print('    gamma =', complex(sp.limit(Gam*z,z,0)), '  atteso 2p+1 =', 2*p+1)
print('    delta =', complex(sp.limit(Gam*(z-1),z,1)), '  atteso 2q_e+1 =', 2*qe+1)
print('    eps   =', complex(sp.limit(Gam,z,sp.oo)), '  atteso 2s =', 2*s)
prod = sp.simplify(sp.expand(Lam*z*(z-1)))
alpha = complex(sp.limit(prod/z, z, sp.oo)); qH = -complex(prod.subs(z,0))
print('\n  Lambda(z) z(z-1) deve essere  alpha z - q :')
print(f'    alpha = {alpha:+.10f}')
print(f'    q     = {qH:+.10f}')
resid = sp.simplify(prod - (sp.nsimplify(alpha,rational=False)*z - sp.nsimplify(qH,rational=False)))
print(f'    residuo del matching (deve essere 0): {complex(sp.simplify(resid.subs(z,sp.Rational(1,3)))):.3e}')
print(f'\nHeunC[q, alpha, gamma, delta, epsilon, z] con')
print(f'  q={qH!r}\n  alpha={alpha!r}\n  gamma={2*p+1!r}\n  delta={2*qe+1!r}\n  epsilon={2*s!r}')
