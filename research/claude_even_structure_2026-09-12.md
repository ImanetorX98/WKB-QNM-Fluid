# Controlli numerici sulla quantizzazione angolare di Codex

Verifica di `research/bohr_sommerfeld_action_2026-09-12.md`. **Tutto confermato,
compreso il difetto che Codex ha trovato nel mio numero.**

## 1. Il mio \(A_0\) era estrapolato male, e la deriva era quello

Codex calcola \(A_0\) **dall'azione**, senza solutore di autovalori:

\[
A_0 = 0.897478121218106117251403583,
\]

stabile a precisione 35 e 55, residuo dell'azione \(4.5\times10^{-36}\).

Io lo avevo estrapolato da due \(L\) finiti ottenendo \(0.8974781250072938\):
**più grande di \(3.79\times10^{-9}\)**.

E quella differenza spiega la deriva che avevo attribuito a \(A_3/L^3\). Con
\(\delta=3.79\times10^{-9}\), il termine spurio \(-L^2\delta\) vale
\(+6.2\times10^{-6}\) a \(L=40.5\) e \(+1.17\times10^{-5}\) a \(L=55.5\): una
deriva prevista di \(+5.5\times10^{-6}\) contro gli \(+8.6\times10^{-6}\)
osservati. Stesso ordine, stesso segno.

> Avevo scritto «costante a cinque cifre, con deriva compatibile con
> \(A_3/L^3\)». La deriva era **mia**. Codex ha ragione a chiedere di non
> estrapolare \(A_0\) assumendo già assente \(A_1\).

## 2. Il test vincolante a \(\hat c=0\): esatto

A \(\hat c=0\) vale \(A=\ell(\ell+1)=L^2-\tfrac14\), quindi
\(\hat A=1-\tfrac{1}{4L^2}\) e \(A_2=-\tfrac14\) **esattamente**.

| \(L\) | 40.5 | 55.5 | 70.5 |
|---|---|---|---|
| \(L^2(\hat A-1)\) | −0.25 | −0.25 | −0.25 |

Sedici cifre. La catena — convenzione \(\lambda_{\rm MMA}-c^2\), \(m\) intero
esatto, estrazione — è sana.

## 3. Solo potenze pari: confermato, e in modo discriminante

Codex deriva che la quantizzazione formale è
\(I(\hat A)+\varepsilon^2I_2(\hat A)+\dots=\pi(1-|\mu|)\), **senza potenze
dispari**, da cui \(A_1I_A=0\) e quindi \(A_1=0\), e altrettanto \(A_3=0\).

Adattando \(L^2(\hat A-A_0)\) su \(\ell=40,43,\dots,100\) con \(A_0\) esatto:

| modello | parametri | residuo massimo |
|---|---|---|
| \(A_2+b/L+c/L^2\) (\(b\) libero) | 3 | 4.64e−9 |
| \(A_2+b/L+c/L^2+d/L^3\) | 4 | 1.89e−10 |
| \(A_2+c/L^2\) (solo pari) | 2 | 2.39e−8 |
| **\(A_2+c/L^2+d/L^4\)** (solo pari) | **3** | **1.18e−12** |

Il modello **a sole potenze pari con tre parametri batte di quasi tre ordini
quello con l'esponente dispari libero e lo stesso numero di parametri**. I
coefficienti dispari escono a \(10^{-5}\)–\(10^{-6}\), cioè al livello
dell'errore di modello: stanno adattando il residuo, non un segnale.

\[
A_2(\mu=2/3,\;\hat c=3/5) = -0.2310651878 .
\]

Il mio valore precedente, \(-0.2311\), portava dentro l'errore su \(A_0\).

## 4. \(\partial_{\hat c}A_0\): due strade, stesso limite

Formula bilineare \(A_c=-2c\,(v^TM_2v)/(v^Tv)\) divisa per \(L\), contro il
valore dall'azione \(-0.34922245866756854380202\):

| \(\ell\) | 40 | 61 | 82 | 121 | 160 | 241 |
|---|---|---|---|---|---|---|
| scarto rel. | 1.54e−5 | 6.81e−6 | 3.81e−6 | 1.76e−6 | 1.01e−6 | **4.48e−7** |

Converge come \(1/L^2\) (rapporto 34 su un rapporto di \(L\) di 6, contro 36
atteso): ancora la struttura a sole potenze pari, vista da un'altra direzione.

## 5. Stato

* \(A_1=0\) ha ora **una dimostrazione** (Codex, dalla struttura pari della
  quantizzazione) e **due conferme numeriche indipendenti** — la nostra matrice
  armonica e `SpheroidalEigenvalue` di Mathematica.
* Il compito (b) del `CODEX_HANDOFF_KERR_A1.md` è quindi chiuso per \(A_1\).
* Resta aperto \(I_2\), cioè \(A_2\) in forma chiusa: richiede il periodo su
  ciclo complesso, e Codex avverte giustamente di non darlo a un integratore
  reale con estremi tagliati a mano. Il vincolo \(A_2(\hat c=0)=-\tfrac14\) è
  verificato sopra ed è il test di segni e normalizzazione per quando ci si
  arriverà.

## Comando

```sh
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt \
  -script calculations/claude_verify_even_structure.wl
```
