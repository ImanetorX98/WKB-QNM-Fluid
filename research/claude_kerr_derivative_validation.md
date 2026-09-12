# Estensione E e ricorrenze ai bordi — validazione

**Esito: E1–E4 chiusi. Le ricorrenze dei bordi funzionano, ma quella
all'orizzonte ha un raggio di convergenza di \(1.7\times10^{-3}\)**, che era la
riserva della nota sugli «zeri del profilo», qui quantificata.

Script: `claude_kerr_derivative_validation.py` (6.1 s),
`claude_kerr_derivative_boundary.py` (≈40 s).

---

## E1–E2 — \(A_c\) bilineare

Nella base armonica sferica, ortonormale in \(\int\!\cdot\,\sin\theta\,d\theta\),
con \(S=\sum_j v_j Y_j\):

\[
\int S^2\sin\theta\,d\theta=v^{T}v,\qquad
\int S^2\cos^2\theta\,\sin\theta\,d\theta=v^{T}M_2 v,
\]

**trasposto, non coniugato**, quindi \(A_c=-2c\,(v^{T}M_2v)/(v^{T}v)\).

Contro differenze centrali, per \((\ell,m)=(28,19)\) e \((61,41)\) — entrambi con
\(m/L=2/3\) **esatto** — e \(c/L=0.2,\;0.4,\;0.3+0.1i\):

| \((\ell,m)\) | \(c/L\) | passo ottimo | dir. reale | dir. immaginaria | \(v^{T}v\) |
|---|---|---|---|---|---|
| (28,19) | 0.2 | 1e−4 | **4.2e−11** | 5.9e−11 | 1.0000 |
| (28,19) | 0.4 | 1e−4 | 2.7e−10 | 1.9e−10 | 1.0000 |
| (28,19) | 0.3+0.1i | 1e−4 | 3.1e−10 | 1.6e−10 | **0.9852** |
| (61,41) | 0.2 | 1e−4 | **4.5e−11** | 6.5e−11 | 1.0000 |
| (61,41) | 0.4 | 1e−4 | 1.9e−10 | 2.0e−10 | 1.0000 |
| (61,41) | 0.3+0.1i | 1e−4 | 3.2e−10 | 1.9e−10 | **0.9332** |

Target del brief \(10^{-6}\): superato di **quattro ordini**.

Il passo \(10^{-3}\) è dominato dal troncamento (\(\sim2\times10^{-8}\)),
\(10^{-5}\) dall'arrotondamento: l'ottimo è \(10^{-4}\), come atteso per una
differenza centrale.

**Il denominatore bilineare va registrato, e qui si vede perché**: per \(c\)
reale \(v^{T}v=1\) esatto, ma per \(c\) complessa vale 0.9852 e 0.9332. È la
misura di quanto \(v^{T}v\) differisca da \(v^{\dagger}v\), cioè esattamente
l'errore che si commetterebbe sostituendo \(S^2\) con \(\lvert S\rvert^2\).
Entrambi restano lontani da zero: nessuna auto-ortogonalità in vista.

**Dimensione della base**: `pad` 40 contro 80 dà \(2.4\times10^{-15}\) e
\(1.6\times10^{-15}\). Convergente.

---

## E3 — \(q_\omega\) a tre raggi

Formula della nota contro differenziazione dell'**intero** potenziale, con \(A\)
ricalcolata a ogni \(\omega\). \(\omega=3-0.25i\), \((\ell,m)=(28,19)\).

| \(a\) | \(r\) | scarto formula/FD | **\(A\) congelata** (controllo negativo) |
|---|---|---|---|
| 0.0 | 3.20 | 7.3e−12 | 7.3e−12 |
| 0.0 | 3.00 | 7.3e−12 | 7.3e−12 |
| 0.0 | 6.00 | 7.3e−12 | 7.3e−12 |
| 0.6 | 2.88 | **3.7e−10** | **1.95e−01** |
| 0.6 | 3.00 | 3.5e−10 | 1.85e−01 |
| 0.6 | 6.00 | 1.3e−10 | 6.76e−02 |

Il controllo negativo è netto: congelare \(A\) sbaglia dal **7% al 19.5%** a
\(a=0.6\). A \(a=0\) non cambia nulla, ed è corretto che sia così — il
coefficiente di \(A_c\) è proporzionale ad \(a\).

## E4 — limite \(a=0\)

\(q_\omega=2\omega\) con scarto **0.00e+00**. Va però detto che a \(a=0\) il
termine in \(A_c\) ha prefattore nullo: il limite è esatto **per costruzione**, e
non è un controllo forte. Il controllo che pesa è E3 a \(a=0.6\).

---

## Ricorrenze ai bordi

Implementate le ricorrenze §2 e §3 della nota. I coefficienti di \(v\) e \(q\)
non sono trascritti a mano: \(q(r)\) è costruita in forma simbolica dalle
definizioni e sviluppata per **divisione polinomiale esatta**. (`sympy.series`
fallisce con `NotImplementedError` su queste razionali a coefficienti complessi;
la divisione numeratore/denominatore è esatta e non ha quel problema.)

Controllo di consistenza: la serie di \(q\) riproduce \(q\) valutata
direttamente a \(3\times10^{-15}\), e \(q_0=-h_0^2\) come richiede l'equazione
indiciale.

### All'infinito: si comporta come deve

| ordine \ \(r\) | 40 | 80 | 160 |
|---|---|---|---|
| 4 | 4.7e−4 | 1.8e−5 | 6.1e−7 |
| 6 | 2.1e−5 | 1.9e−7 | 1.6e−9 |
| 8 | 1.1e−6 | 2.4e−9 | **5.0e−12** |

Il residuo cala sia con l'ordine sia con il raggio, come per una serie
asintotica ben posta.

### All'orizzonte: raggio \(1.7\times10^{-3}\)

Il primo tentativo a \(\rho=0.02\) dava residui fino a \(3\times10^{16}\), **in
crescita con l'ordine**. Non è un errore della ricorrenza: è divergenza.

\[
\left|\frac{h_n}{h_{n-1}}\right|\longrightarrow 591.0
\quad\text{(stabile a quattro cifre da }n=6\text{)}
\;\Longrightarrow\;
R\simeq 1.69\times10^{-3}.
\]

| ordine \ \(\rho\) | 2e−3 (fuori) | 5e−4 | 2e−4 |
|---|---|---|---|
| 4 | 2.95 | 2.0e−3 | 1.9e−5 |
| 6 | 3.10 | 2.6e−4 | 4.0e−7 |
| 8 | **7.51** | 3.0e−5 | **7.4e−9** |

Dentro il raggio il residuo cala con l'ordine; fuori **cresce**. Il denominatore
di Frobenius non c'entra: \(\min\lvert n v_1+2h_0\rvert=0.338\), lontano da zero.

> La causa è un **polo di \(D\)**, cioè uno zero di \(\psi\), a distanza
> \(1.7\times10^{-3}\) dall'orizzonte. È la riserva della nota — «lontano da
> risonanze della ricorrenza e da **zeri del profilo**» — qui quantificata per
> \(a=0.6\), \((\ell,m)=(28,19)\), \(\omega=3-0.25i\).

**Conseguenza operativa.** La serie all'orizzonte va inizializzata a
\(\rho\lesssim10^{-4}\) e propagata con l'ODE, esattamente la «procedura
operativa» della §4 della nota. Non è utilizzabile a un raggio di matching
ordinario.

### \(K=D_\omega\)

| bordo | posizione | \(K\) | \(\lvert\)reale − immag.\(\rvert\) | \(\lvert K-\text{limite}\rvert\) |
|---|---|---|---|---|
| orizzonte | \(r_+ +2\!\cdot\!10^{-4}\) | \(-0.031438-1.367861\,i\) | 2.4e−8 | 3.7e−1 |
| orizzonte | \(r_+ +5\!\cdot\!10^{-4}\) | \(+0.302697-2.037990\,i\) | 7.2e−8 | 1.1e+0 |
| infinito | 60 | \(-0.001935+1.012070\,i\) | 2.6e−11 | 1.2e−2 |
| infinito | 120 | \(-0.000495+1.003026\,i\) | 8.7e−12 | 3.1e−3 |

Le due direzioni di differenziazione concordano: \(K\) è olomorfa, e \(A\) è
stata ricalcolata a ogni \(\omega\) — nessun congelamento.

All'infinito \(K\to i\) con scarto che cala come \(1/r^2\) (1.2e−2 → 3.1e−3 al
raddoppiare di \(r\)), coerente con \(D_{+,\omega}=i+t_{2,\omega}/r^2+\dots\).

All'orizzonte lo scarto da \(-i\) **non** è piccolo, ed è corretto che non lo
sia: con \(h_1\approx-39-273i\) il termine \(h_{1,\omega}\rho\) domina già a
\(\rho=2\times10^{-4}\). Il limite \(-i\) si raggiunge solo molto più vicino.

---

## Che cosa resta non validato

Il punto 4 del brief resta **aperto**, e non per mancanza di formule:

* si è verificata la **consistenza** delle serie, non la **selezione del ramo**;
* per \(\operatorname{Im}\omega<0\) il ramo uscente cresce lungo l'asse reale e
  una contaminazione entrante è esponenzialmente subdominante: la serie in
  \(1/r\) da sola non lo seleziona, e la stabilità in \(R\) e nell'ordine non lo
  dimostra;
* serve la continuazione dal semipiano \(\operatorname{Im}\omega>0\), oppure
  Leaver/Jost, oppure un contorno giustificato.

**La norma globale Kerr resta quindi non validata**, come la nota prescrive di
dichiarare.

## Comandi

```sh
python3.13 calculations/claude_kerr_derivative_validation.py
python3.13 calculations/claude_kerr_derivative_boundary.py
```
