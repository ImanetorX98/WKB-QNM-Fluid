"""Bounded §9.6 audit. Prints Markdown; no claim of exact QNM frequencies."""
import numpy as np
from kerr_madelung_analytic import analytic_profile
from kerr_wkb3_selfconsistent import selfconsistent_frequency


def run():
    print('| a | ell,m | frequency | BC order | offset | median Q_M | r_end |')
    print('|---|---|---|---|---|---|---|', flush=True)
    for spin, ell in ((0.,28),(.6,28),(.6,61)):
        mu=2/3
        m=(2*ell+1)//3
        assert abs(m/(ell+.5)-mu)<1e-15
        result=selfconsistent_frequency(ell,spin,mu)
        if not result['converged'] or result['residual']>1e-8:
            raise RuntimeError(result)
        for name,omega in (('eikonal',result['omega_eikonal']),('WKB3',result['omega'])):
            for order,offset,points in ((0,1e-4,12001),(0,1e-6,12001),
                                        (1,1e-4,12001),(1,1e-6,12001),
                                        (1,1e-6,24001)):
                data=analytic_profile(ell,spin,mu,omega=omega,
                    horizon_offset=offset,points=points,boundary_order=order)
                mask=(data['r']>20)&(data['r']<50)
                if not mask.any():
                    raise RuntimeError('Empty measurement window')
                print(f'| {spin} | {ell},{m} | {name} ({points}) | {order} | {offset:g} '
                      f'| {np.median(data["q_madelung"][mask]):.12g} | {data["r"][-1]:.9f} |',flush=True)
        print(f'\nMode {ell},{m}: omega_WKB3={result["omega"]}, '
              f'residual={result["residual"]:.3g}.\n',flush=True)


if __name__=='__main__':
    run()
