# Rigenerazione del riferimento Kerr — 12 settembre 2026

Comando: `python3.13 calculations/kerr_corrected_profile_audit.py`.
30 integrazioni, completate con exit 0; M=1, n=0, mu=2/3 esatto.
La mediana è su 20<r<50; Q_M=-epsilon^2 |psi|''/|psi|.
Derivate ricavate dall'ODE, geometria corretta, span tortoise analitico.

## Sintesi WKB3

| a | ell,m | mediana Q_M (BC1, offset 1e-6, 24001) | variazione relativa offset BC1 | variazione griglia |
|---|---|---|---|---|
| 0 | 28,19 | -1.15374320942e-5 | 5.27e-7 | 8.91e-7 |
| 0.6 | 28,19 | -1.02841097333e-5 | 5.19e-7 | sotto cifre stampate |
| 0.6 | 61,41 | -2.20839990718e-6 | 2.59e-6 | sotto cifre stampate |

Questi numeri sono stabilità della mediana, non errore certificato del profilo
né della frequenza. La griglia è di campionamento t_eval: DOP853 sceglie
indipendentemente i passi interni. Raddoppiarla non verifica la tolleranza ODE.
L'offset varia fra 1e-4 e 1e-6, le griglie fra 12001 e 24001.
A ell61 il bordo di ordine zero a offset1e-6 differisce da BC1 di circa0.214%:
la sensibilità al bordo è ancora materialmente superiore a molte cifre delle
vecchie tabelle. La frequenza eikonale presenta maggiore sensibilità del
campionamento e un warning del vecchio solutore (riportato sotto).

Non è stata riprodotta l'intera tabella d'errore §9.6: manca il confronto
con una previsione corretta e frequenze QNM indipendenti. Non interpretare
il residuo WKB3 come accuratezza spettrale. Nessuna rivendicazione di novità
fisica deriva da questo sweep.

## Output completo (incluso warning)

```text
| a | ell,m | frequency | BC order | offset | median Q_M | r_end |
|---|---|---|---|---|---|---|
| 0.0 | 28,19 | eikonal (12001) | 0 | 0.0001 | -1.14939808586e-05 | 60.000000000 |
| 0.0 | 28,19 | eikonal (12001) | 0 | 1e-06 | -1.15047535863e-05 | 60.000000001 |
| 0.0 | 28,19 | eikonal (12001) | 1 | 0.0001 | -1.15070943049e-05 | 60.000000000 |
| 0.0 | 28,19 | eikonal (12001) | 1 | 1e-06 | -1.15070810162e-05 | 59.999999999 |
| 0.0 | 28,19 | eikonal (24001) | 1 | 1e-06 | -1.15071094469e-05 | 59.999999999 |
| 0.0 | 28,19 | WKB3 (12001) | 0 | 0.0001 | -1.14958843377e-05 | 60.000000000 |
| 0.0 | 28,19 | WKB3 (12001) | 0 | 1e-06 | -1.15349208089e-05 | 60.000000001 |
| 0.0 | 28,19 | WKB3 (12001) | 1 | 0.0001 | -1.15374484531e-05 | 60.000000000 |
| 0.0 | 28,19 | WKB3 (12001) | 1 | 1e-06 | -1.15374423752e-05 | 60.000000001 |
| 0.0 | 28,19 | WKB3 (24001) | 1 | 1e-06 | -1.15374320942e-05 | 60.000000001 |

Mode 28,19: omega_WKB3=(5.485046110897005-0.09622922179660573j), residual=2.91e-15.

/Users/iman.rosignoli/Desktop/ISTRUZIONE/PROGRAMMAZIONE/WKBQNMMadelungFluid/calculations/kerr_radial_order_profile.py:91: RuntimeWarning: The iteration is not making good progress, as measured by the
 improvement from the last ten iterations.
  current = fsolve(residuals, current, args=(current_spin,), xtol=1.0e-12)
| 0.6 | 28,19 | eikonal (12001) | 0 | 0.0001 | -1.17369435839e-05 | 60.000000000 |
| 0.6 | 28,19 | eikonal (12001) | 0 | 1e-06 | -1.17362015002e-05 | 59.999999998 |
| 0.6 | 28,19 | eikonal (12001) | 1 | 0.0001 | -1.17366003229e-05 | 60.000000000 |
| 0.6 | 28,19 | eikonal (12001) | 1 | 1e-06 | -1.17359903853e-05 | 60.000000003 |
| 0.6 | 28,19 | eikonal (24001) | 1 | 1e-06 | -1.17358320116e-05 | 60.000000003 |
| 0.6 | 28,19 | WKB3 (12001) | 0 | 0.0001 | -1.0262340527e-05 | 60.000000000 |
| 0.6 | 28,19 | WKB3 (12001) | 0 | 1e-06 | -1.02704520915e-05 | 59.999999998 |
| 0.6 | 28,19 | WKB3 (12001) | 1 | 0.0001 | -1.02841150675e-05 | 60.000000000 |
| 0.6 | 28,19 | WKB3 (12001) | 1 | 1e-06 | -1.02841097333e-05 | 59.999999999 |
| 0.6 | 28,19 | WKB3 (24001) | 1 | 1e-06 | -1.02841097333e-05 | 59.999999999 |

Mode 28,19: omega_WKB3=(6.726476235883466-0.09105585810999496j), residual=9.63e-14.

| 0.6 | 61,41 | eikonal (12001) | 0 | 0.0001 | -2.52004038681e-06 | 60.000000000 |
| 0.6 | 61,41 | eikonal (12001) | 0 | 1e-06 | -2.51891442801e-06 | 59.999999999 |
| 0.6 | 61,41 | eikonal (12001) | 1 | 0.0001 | -2.51960740828e-06 | 60.000000000 |
| 0.6 | 61,41 | eikonal (12001) | 1 | 1e-06 | -2.51872088225e-06 | 60.000000002 |
| 0.6 | 61,41 | eikonal (24001) | 1 | 1e-06 | -2.51909223788e-06 | 60.000000002 |
| 0.6 | 61,41 | WKB3 (12001) | 0 | 0.0001 | -2.2032584435e-06 | 60.000000000 |
| 0.6 | 61,41 | WKB3 (12001) | 0 | 1e-06 | -2.20366727394e-06 | 60.000000002 |
| 0.6 | 61,41 | WKB3 (12001) | 1 | 0.0001 | -2.20839418193e-06 | 60.000000000 |
| 0.6 | 61,41 | WKB3 (12001) | 1 | 1e-06 | -2.20839990718e-06 | 59.999999999 |
| 0.6 | 61,41 | WKB3 (24001) | 1 | 1e-06 | -2.20839990718e-06 | 59.999999999 |

Mode 61,41: omega_WKB3=(14.514470147120745-0.09105324749728791j), residual=1.11e-14.


```

## Controlli aggiuntivi

`verify_diagnostic_scope.wl`, eseguito con WolframKernel con licenza:
3/3 identità True (riscalatura ampiezza, linearizzazione Riccati, coefficiente
al bordo). I due test geometrici precedenti sono stati rieseguiti: 2/2 OK.
