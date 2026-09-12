# Controllo G — secondo ordine e zeri mobili

**Esito: chiuso, su tutti e sei i punti.** Il termine degli zeri mobili non solo
è necessario: in questo benchmark **domina**, e ometterlo cambia il segno di
\(E_2\).

Script: `calculations/claude_compact_second_order_response.py` — 9.4 s.
Nessun commit oltre i file miei.

---

## 1. Il controesempio esatto, che non coinvolge QNM

\(Q=x-\eta\), \(w=1\), \(J=[-1,1]\), dunque \(U(\eta)=1+\eta^2\) e \(U_2=2\).
Qui \(Q_2=0\) identicamente: **tutto** il risultato viene dallo zero mobile.

| | |
|---|---|
| termine integrale \(\int w\,\mathrm{sgn}(Q)\,Q_2\) | 0 |
| termine degli zeri \(2\,w(x_j)Q_1^2/\lvert Q_x\rvert\) | **2.0** |
| somma | 2.0 |
| differenze centrali seconde | 1.9999999994 |

Omettendo il termine si otterrebbe \(0\) contro un valore vero di \(2\).

---

## 2. \(d=\omega'\) ed \(e=\omega''\)

Sono **derivate**, non coefficienti: \(\omega=\omega_0+\eta d+\eta^2 e/2\).
Benchmark esterno, centro 2, larghezza 0.4.

| passo \(\eta\) | \(d\) | \(e\) |
|---|---|---|
| 1e−2 | \(+0.0103755185-0.1760166061\,i\) | \(+0.11077958-0.30985192\,i\) |
| 3e−3 | \(+0.0103630505-0.1760007385\,i\) | \(+0.11072617-0.30981624\,i\) |
| 1e−3 | \(+0.0103619525-0.1759993440\,i\) | \(+0.11072272-0.30981268\,i\) |

\(d\) converge al valore analitico \(B/N=0.010361818212-0.175999169467\,i\) del
gate A, con lo scarto \(O(\eta^2)\) atteso.

---

## 3. \(j\) dalle due vie

ODE \(j'+2zj=-2-2k^2\) con \(j(a)=0\), contro la forma integrale
\(j=\psi^{-2}\big[-\int_a^x(2+2k^2)\psi^2\big]\):

| \(x\) | \(j\) | scarto |
|---|---|---|
| −0.8 | \(-0.000618410-0.000330570\,i\) | 6.0e−16 |
| −0.4 | \(-0.025532787-0.293806032\,i\) | 3.2e−15 |
| 0.0 | \(-1.264479740-2.501678400\,i\) | 5.5e−14 |
| +0.4 | \(-17.161830547+7.704515409\,i\) | 1.9e−12 |
| +0.8 | \(+39.169341827-19.406350685\,i\) | 2.7e−13 |

---

## 4. Gli zeri e il loro moto

Due zeri semplici in \(J\), simmetrici. \(Q_x\) da derivate controllate —
\(p_x=-\operatorname{Im}(z^2)-\operatorname{Im}(\omega^2)\) dalla Riccati,
\(Q_x=-V_0'-2p\,p_x\) — non da una differenza grossolana.

| \(x_j\) | \(Q_x\) | \(Q_1(x_j)\) | \(-Q_1/Q_x\) | \(x_j'\) misurata | scarto |
|---|---|---|---|---|---|
| \(-0.3586793241\) | \(+9.52235125\) | \(+0.25287581\) | \(-0.02655603\) | \(-0.02655601\) | **1.3e−8** |
| \(+0.3586793241\) | \(-9.52235125\) | \(-3.31012729\) | \(-0.34761659\) | \(-0.34761624\) | **3.5e−7** |

Gli zeri sono stati **seguiti** a \(\pm\eta\), non ricercati da capo: restano
semplici e lontani dagli estremi, quindi le ipotesi del §4 della nota reggono.

Nota: i due zeri sono simmetrici in posizione e hanno \(Q_x\) opposto, ma
\(Q_1\) **non** è simmetrico — 0.253 contro −3.310. La condizione uscente rompe
la simmetria del profilo anche dove il potenziale la conserva.

---

## 5. \(U_2\) ed \(E_2\)

| | |
|---|---|
| \(U_1\) | \(+0.1279372872\) |
| \(D_1\) | \(+1.1246457120\) |
| \(E_1\) | \(+0.0018943203\) |
| \(U_2\), termine integrale | \(+0.1496918781\) |
| \(U_2\), **contributo degli zeri** | \(+2.0353076450\) |
| \(U_2\) totale | \(+2.1849995231\) |
| \(D_2\) | \(+3.7658660751\) |
| **\(E_2\) completa** | \(+0.0621331857\) |
| \(E_2\) senza il termine degli zeri | \(-0.0026657930\) |

\(E_1=0.0018943203\) coincide con la derivata prima del gate B
(\(0.001894279539\)) a \(2\times10^{-5}\), che è il troncamento della
differenza finita usata qui per \(d\).

### Contro le differenze centrali seconde

\(E\) calcolata a quadratura **spezzata sugli zeri di ciascun \(\eta\)**:

| passo \(\eta\) | \(E_2\) numerica | scarto con il termine degli zeri | scarto senza |
|---|---|---|---|
| 1e−2 | \(+0.0621511766\) | 2.90e−4 | 24.31 |
| 3e−3 | \(+0.0621342466\) | 1.71e−5 | 24.31 |
| 1e−3 | \(+0.0621331535\) | **5.19e−7** | 24.31 |

Target del brief \(10^{-3}\): superato di **tre ordini e mezzo**. Nessun plateau
osservato nell'intervallo esplorato — la discesa è \(O(\eta^2)\) pulita.

---

## 6. Controllo negativo

Il brief chiedeva di non presupporre che il termine degli zeri domini. **In
questo benchmark domina**, e nettamente:

* vale \(2.035\) su un \(U_2\) totale di \(2.185\), cioè il **93%**;
* ometterlo sposta \(E_2\) da \(+0.062133\) a \(-0.002666\), un cambio del
  **104%** che ne inverte il segno;
* lo scarto con le differenze finite passa da \(5\times10^{-7}\) a **24.3**,
  cioè un fattore \(4.7\times10^7\) — ben oltre il fattore 10 richiesto.

Nessuno zero è diventato multiplo né ha attraversato un estremo nell'intervallo
di \(\eta\) esplorato, quindi non si è dovuto interrompere.

---

## Limiti dichiarati

* Vale per **questo** benchmark compatto con \(q=\omega^2-V_0\), dove
  \(q_{\omega\omega}=2\) e \(j(a)=0\). **Per Kerr non si può porre
  \(q_{\omega\omega}=2\)**: servirebbe anche la seconda derivata dell'autovalore
  sferoidale, che non è stata calcolata.
* Le formule valgono per il **rapporto integrale**, non per la mediana pesata,
  che non è un integrale di \(\lvert Q\rvert\) e non è differenziabile in
  presenza di plateaux.
* Il rango al più cinque al secondo ordine non è stato testato: il brief non lo
  chiedeva in G, e un rango superiore a due a \(\eta\) finito **non** confuta il
  teorema lineare.

## Comando

```sh
python3.13 calculations/claude_compact_second_order_response.py
```
