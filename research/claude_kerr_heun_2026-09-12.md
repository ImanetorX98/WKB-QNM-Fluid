# N1 — `HeunC` per i bordi di Kerr: **il bordo interno è risolto**

Esito in una riga: **sì per l'orizzonte, no per l'infinito.** Il compito analitico
A1 di Codex resta necessario, ma solo per il bordo esterno e solo se il bersaglio
di accuratezza lo richiede.

## 1. La riduzione a Heun confluente

La radiale di Teukolsky con \(s=0\) è confluente Heun. Con
\(z=(r-r_+)/(r_--r_+)\), che manda \(r_+\mapsto0\) e \(r_-\mapsto1\), la forma
normale ha

$$P_z=\frac1z+\frac1{z-1},$$

e posto \(R=z^{p}(z-1)^{q_e}e^{sz}H(z)\) la funzione \(H\) soddisfae l'equazione
di Heun confluente nella convenzione di *Mathematica*,
\(H''+(\gamma/z+\delta/(z-1)+\epsilon)H'+\frac{\alpha z-q}{z(z-1)}H=0\), con

$$\gamma=2p+1,\qquad \delta=2q_e+1,\qquad \epsilon=2s .$$

Gli esponenti non sono trascritti: escono dai poli doppi di \(Q_z\) e
**coincidono con le formule attese**,

| | calcolato | atteso |
|---|---|---|
| \(p\) | \(0.5625-0.3750\,i\) | \(iK(r_+)/(r_+-r_-)\) ✓ |
| \(s\) | \(0.4000+4.8000\,i\) | \(-i\omega(r_--r_+)\) ✓ |

Il matching dei restanti parametri ha **residuo esattamente nullo**:
`calculations/claude_kerr_heun_parameters.py`. Per \(a=3/5\),
\((\ell,m)=(28,19)\), \(\omega=3-\tfrac{i}{4}\):

\[
q=746.823871198005+16.110026570174\,i,\qquad \alpha=0.8+9.6\,i .
\]

## 2. \(D_-\) senza serie asintotiche

Il prefattore multivalore **non serve**: la derivata logaritmica ne conserva solo
\(p/z+q_e/(z-1)+s\), che non ha tagli. Quindi

$$D_-(r)=v(r)\left[\frac{r}{\Sigma}+\frac1b\left(\frac pz+\frac{q_e}{z-1}+s
+\frac{H'}{H}\right)\right],\qquad \Sigma=r^2+a^2,$$

con \(H=\texttt{HeunC}\), \(H'=\texttt{HeunCPrime}\).

### Confronto dove entrambe valgono

A \(\rho=10^{-4}\), dentro il raggio della serie di Frobenius:

| | \(D_-\) |
|---|---|
| Frobenius (Python, ordine 8) | \(-0.252977188245+0.138107963531\,i\) |
| `HeunC` (Mathematica) | \(-0.252977188242+0.138107963532\,i\) |

**Undici cifre.** Le due costruzioni non condividono nulla: convenzioni ricavate
separatamente, un potenziale costruito per divisione razionale e uno simbolico,
linguaggi diversi.

### E dove la serie muore

Residuo della Riccati \(vD'+D^2+q\):

| \(\rho\) | \(10^{-4}\) | \(10^{-3}\) | \(2\times10^{-2}\) | \(0.1\) | \(0.5\) | \(1.0\) |
|---|---|---|---|---|---|---|
| residuo | 4.5e−11 | 6.7e−10 | **8.9e−13** | **4.8e−14** | **1.1e−13** | **1.4e−13** |

La serie di Frobenius aveva raggio \(1.7\times10^{-3}\) e a \(\rho=2\times10^{-2}\)
divergeva con residuo \(7.5\). `HeunC` dà \(9\times10^{-13}\): **il raggio
utilizzabile passa da \(1.7\times10^{-3}\) a \(\ge1\)**, cioè un fattore
\(\sim600\), e in linea di principio fino a \(|b|=1.6\).

*Nota sui due valori peggiori.* A \(\rho=10^{-4}\) e \(10^{-3}\) il residuo è
1e−11–1e−10 invece di 1e−13: è l'arrotondamento della differenza centrale con
passo \(10^{-7}\) là dove \(D\) varia rapidamente, non un difetto della
rappresentazione.

## 3. Un difetto mio, trovato dal residuo

La prima esecuzione dava residui che **crescevano con \(\rho\)** — fino a
\(3\times10^{-2}\). Non era Heun: nel trascrivere il potenziale avevo scritto
`D[D0/(x^2+a^2), x]` con `D0` già valutato in \(r\), quindi la derivata toccava
solo il denominatore. Derivando \(v\) **simbolicamente** il residuo crolla a
\(10^{-13}\).

È stato il confronto con la serie di Frobenius a dire che Heun era giusto e il
potenziale no: senza un secondo calcolo indipendente avrei accusato Heun.

## 4. Che cosa **non** risolve

`HeunC` è la soluzione locale al punto singolare **regolare**. Per il ramo
uscente all'infinito, che è un punto singolare **irregolare**, servirebbero i
coefficienti di connessione, e *Mathematica* non li espone: i built-in
disponibili — `HeunB`, `HeunC`, `HeunD`, `HeunG`, `HeunT` e le rispettive
primate — sono tutte soluzioni a \(z=0\). Nessun pacchetto Teukolsky installato.

Quindi:

* **bordo interno**: risolto, convergente, nessuna ambiguità di ramo — e non ce
  n'era mai stata, perché il punto regolare seleziona l'esponente entrante; il
  problema era solo il raggio;
* **bordo esterno**: resta la serie asintotica, con il budget
  \(|2\omega|e^{2\operatorname{Im}(\omega)r}\) già quantificato.

## 5. Conseguenza per la sincronizzazione

Il compito **A1 di Codex non è annullato, ma è dimezzato**: serve solo la
rappresentazione convergente al bordo **esterno**, e solo se il bersaglio di
accuratezza supera quanto il budget concede. A \(r=120\) l'ambiguità vale
\(1.4\times10^{-9}\) per il fondamentale.

## Comandi

```sh
python3.13 calculations/claude_kerr_heun_parameters.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt \
  -script calculations/claude_kerr_heun_boundary.wl
```
