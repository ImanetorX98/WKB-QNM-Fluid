# Barriera esatta, residui per ordine e zeri del fluido di Madelung

**Nota di calcolo — 7 settembre 2026.** Questa nota registra sia un controllo
positivo fuori dalla famiglia di Schwarzschild sia un risultato negativo che
obbliga a modificare il programma di ricerca.

## 1. Perché Pöschl–Teller

La barriera

\[
V(y)=L^2\operatorname{sech}^2y
\tag{1}
\]

ha frequenze QNM esatte

\[
\omega_n=\sqrt{L^2-\frac14}-i\left(n+\frac12\right).
\tag{2}
\]

Posto \(z=(1+\tanh y)/2\) e
\(\lambda=\sqrt{L^2-1/4}\), una soluzione con le condizioni QNM è

\[
\psi_n=[z(1-z)]^{-i\omega_n/2}
{}_2F_1\!\left(-n,-n-2i\lambda;
\frac12-n-i\lambda;z\right).
\tag{3}
\]

L'ipergeometrico termina dopo \(n\) termini. Il benchmark non usa quindi né
Leaver né un'integrazione approssimata per costruire frequenza e autofunzione.
Le frequenze e le torri QNM di questa barriera sono risultati standard; si
vedano Cardona e Molina (2017) e Beyer (1999). L'uso della curvatura di
Madelung come diagnostico dell'errore WKB è la parte sottoposta a test, non lo
spettro (2).

## 2. Controllo del diagnostico completo

Si applica lo stesso funzionale del benchmark Schwarzschild,

\[
\mathcal E_M=
\frac{\int w|Q_M|dy}
{\int w\left(|\omega/L|^2+\operatorname{sech}^2y+|P|^2\right)dy},
\quad
Q_M=-L^{-2}\frac{A''}{A},
\quad
P=L^{-1}\operatorname{Im}\frac{\psi'}\psi.
\tag{4}
\]

Per \(L=1.5,2,3,4,6,8,12,16\):

| overtone | dominio | Pearson con errore WKB1 | Pearson con errore WKB3 | Spearman |
|---:|---|---:|---:|---:|
| 0 | regolare | 1.0000 | 1.0000 | 1.0000 |
| 1 | nodale | non definito | non definito | non definito |
| 2 | regolare | 0.9994 | 0.9989 | 1.0000 |

Il controllo a differenze finite di \(-L^{-2}A''/A\) ha residuo relativo
massimo inferiore a \(9\times10^{-9}\) sui modi regolari. Questo conferma che
la correlazione di Schwarzschild non dipende dal metodo di Leaver o dalla
forma specifica del potenziale di Regge–Wheeler. Essendo Pöschl–Teller una
famiglia a un solo parametro, il test non separa però il contenuto di
\(\mathcal E_M\) dallo scaling in \(1/L\): è un controllo esterno di
consistenza, non ancora una dimostrazione di potere predittivo aggiuntivo.

## 3. Risultato nodale

Per \(n=1\), il polinomio della (3) è

\[
{}_2F_1(-1,-1-2i\lambda;-1/2-i\lambda;z)=1-2z=-\tanh y.
\]

Quindi \(\psi_1(0)=0\). L'ampiezza reale \(A=|\psi|\) ha uno zero e
\(Q_M=-A''/A\) non è definito in quel punto. Il programma ora marca questi
casi come *nodali* invece di produrre un numero finito spurio.

Questo non invalida la formulazione ampiezza–fase locale: ne rende esplicito il
dominio. Un eventuale diagnostico globale dovrà essere costruito per carte
prive di zeri, deflazionare i fattori nodali oppure abbandonare la densità reale
positiva in favore di una fase complessa. Le tre scelte non sono equivalenti e
vanno confrontate, non nascoste dentro una regolarizzazione numerica.

## 4. Fallimento del residuo locale ingenuo

Per la chiusura

\[
q=u^2-\varepsilon^2
\frac{(u^{-1/2})''}{u^{-1/2}},
\]

sono stati costruiti

\[
u^{(0)}=u_0,\qquad
u^{(1)}=u_0+\varepsilon^2u_2,\qquad
u^{(2)}=u_0+\varepsilon^2u_2+\varepsilon^4u_4
\]

e i residui di chiusura \(\mathfrak R_N\). Mathematica verifica che essi
iniziano formalmente a \(\varepsilon^2,\varepsilon^4,\varepsilon^6\).
Nonostante ciò, integrandoli sulla barriera dei fondamentali Schwarzschild si
osserva tipicamente

\[
\mathfrak R_0\sim0.06{-}0.17,\qquad
\mathfrak R_1\sim0.18{-}0.50,\qquad
\mathfrak R_2\sim0.80{-}0.98.
\]

La gerarchia locale peggiora invece di convergere. Il risultato è coerente con
la sua singolarità quando i turning point coalescono: l'ordine formale è valido
a \(q\) fissato e non uniformemente nella regione dove \(q\) diventa piccolo.

Di conseguenza il candidato
\(Q_M-Q_M^{(N)}\) costruito con la serie locale **non** è il diagnostico giusto
per l'errore della WKB di barriera. La prossima costruzione deve partire dalla
forma normale parabolico-cilindrica usata nel matching uniforme. Un candidato
più fedele è la differenza fra la curvatura di ampiezza esatta e quella della
soluzione uniforme,

\[
\Delta Q_M^{\rm unif}
=-\varepsilon^2\left(
\frac{A''}{A}-\frac{A_{\rm PC}''}{A_{\rm PC}}
\right),
\tag{5}
\]

con una prescrizione separata per gli zeri.

## 5. Riproduzione e verifica indipendente

```text
python3.13 calculations/poschl_teller_madelung_benchmark.py --n 0 --n 1 --n 2
python3.13 calculations/order_resolved_madelung_benchmark.py --n 0
python3.13 -m unittest calculations/test_poschl_teller_madelung.py
python3.13 -m unittest calculations/test_order_resolved_madelung.py
```

`calculations/verify_formalism.wl` controlla con Mathematica sia gli ordini
formali dei tre residui sia le autofunzioni esatte \(n=0,1\) e lo zero nodale.

## Riferimenti essenziali

- A. F. Cardona e C. Molina, *Quasinormal modes of generalized
  Pöschl–Teller potentials*, Class. Quantum Grav. 34 (2017) 245002,
  arXiv:1711.00479.
- H. R. Beyer, *On the Completeness of the Quasinormal Modes of the
  Poeschl–Teller Potential*, Commun. Math. Phys. 204 (1999) 397–423,
  arXiv:gr-qc/9803034.
