# Audit dell'archivio PDF e controlli di priorità

**Data:** 10 settembre 2026. **Bibliografia:** 35 voci. **PDF in `papers/`:** 30.

## 1. Stato dell'archivio

| stato | voci |
|---|---|
| **presenti** (24) | 7, 8, 9, 11, 12, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35 |
| **mancanti, scaricabili** (0) | — |
| **mancanti, paywall** (8) | 1 Madelung, 2 Bohm, 3 Regge–Wheeler, 4 Zerilli, 6 Schutz–Will, 10 Konoplya JPS, 13 Ferrari–Mashhoon, 21 Unruh |
| **monografie** (2) | 5 Chandrasekhar, 20 Bender–Orszag |

Scaricati oggi da arXiv, con identificativo verificato sulla pagina `abs`:
`1007.5097` (Dolan), `1207.4253` (Yang *et al.*), `1711.00321` (KMM PNAS),
`1904.10333` (Konoplya–Zhidenko–Zinhailo), `2605.28887` (Meza-Domínguez–Matos).

Le otto voci mancanti sono **classici citati per attribuzione**: nessuna
affermazione del manoscritto dipende dal loro contenuto oltre a ciò che è
universalmente noto. Restano nel manifest per Codex.

## 2. Controllo di priorità sulla §9 (Kerr)

**Esito: la §9.5 del manoscritto è corretta così com'è, e la mia affermazione in
chat del 9 settembre era sovrastimata.** Avevo detto che la §9 «non risulta
preceduta» avendo controllato solo Seidel–Iyer. Letti ora Dolan e Yang:

- **Dolan 2010** sviluppa esplicitamente un metodo di espansione in $1/l$ e
  ottiene «a spin-dependent correction to the frequency at order $1/l$ for
  equatorial modes». Per l'autovalore angolare però **cita** l'espansione in
  $a\omega$ di altri, non la ricostruisce.
- **Yang *et al.* 2012**, Appendice A, espande $A_{\ell m}$ in serie di
  $a\omega/l$, ottiene $A^0_{\ell m}=(l+\tfrac12)^2$ — **il valore di Langer** —
  dalla condizione di Bohr–Sommerfeld, e la correzione
  $A_{\ell m}=l(l+1)-\tfrac{a^2\omega^2}{2}\big(1-\tfrac{m^2}{l(l+1)}\big)$.

Yang dichiara però che «the error in this approximation scales as $1/l$»: il
coefficiente di ordine $1/L$ è **esattamente ciò che il loro schema scarta**. E
il loro parametro $a\omega/l$ non è il nostro $\varepsilon=1/L$ — a $\omega\sim L$
esso è $O(1)$, numericamente piccolo ma non infinitesimo nell'eikonale.

Il manoscritto dice già la cosa giusta: «che $A_1\neq0$ è **noto**, contenuto in
[23, 24]; non è questa la rivendicazione». La rivendicazione è la collocazione
di quel termine nella stessa casella della connessione di spin di Dirac
attraverso una misura operativa comune. **Non va rafforzata.**

## 3. Tre lacune bibliografiche, in ordine di gravità

### (a) Letteratura sullo pseudospettro — la più seria

`2004.06434` (Jaramillo, Panosso Macedo, Al Sheikh) e `2107.12865`
(Gasperín–Jaramillo) sono in `papers/` da prima, **non citati**.

Il risultato di testa del manoscritto è un enunciato di **condizionamento**, e lo
pseudospettro dei QNM *è* la letteratura consolidata sul condizionamento di
questo problema. Sono oggetti diversi — noi trattiamo la sensibilità di $Q_M$ a
$\omega$, loro la sensibilità dello spettro a perturbazioni del potenziale — ma
un referee di CQG si aspetta che la distinzione sia dichiarata, non ignorata.

### (b) Vaidya: un lavoro del 2026 sullo stesso sistema

`2510.25062` — Yoo, Kimura, Ishibashi, Ohashi, *Ringdown in Vaidya spacetimes:
time-dependent frequencies, Penrose limit and time-domain analyses*. In
`papers/`, non citato, e la nostra §11 è nuova.

Metodo distinto: limite di Penrose attorno alla fotosfera dinamica, confronto
con simulazioni nel dominio del tempo. Nessuna occorrenza di «solvability»,
«Madelung», «quantum potential». Ma osservano che «the adiabatic limit of the
PL geometry reflects more details of the dynamical geometry compared to the
naive adiabatic limit of the original geometry» — affermazione **adiacente** alla
nostra (11.2). Va citato e distinto.

Anche `2104.06631` (Lin, Sun, Zhang, *Quasinormal Modes for Dynamical Black
Holes*) è in archivio e non citato.

### (c) Scelta del prodotto scalare, §11.2

`2107.12865` sostiene un prodotto scalare **di energia** per i QNM. Noi usiamo
il prodotto **bilineare** senza coniugazione di Leung *et al.*, appropriato ai
problemi di risonanza. È una scelta, e come tale va motivata contro
l'alternativa esistente invece che presentata come l'unica.

## 4. Da fare

1. Citare (a), (b), (c) con le distinzioni dichiarate sopra.
2. `1512.04611` (Fusca, *The Madelung Transform as a Momentum Map*) è in
   archivio e discusso in sessione: verificare se sia citato.
3. DOI di [10] Konoplya, J. Phys. Stud. 8, 93: da verificare.
