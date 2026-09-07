# Benchmark statico del diagnostico di Madelung

**Nota di calcolo — 7 settembre 2026.** Questo è un risultato pilota e non
una rivendicazione di universalità. La specifica del diagnostico è stata
fissata prima della scansione numerica.

## 1. Domanda falsificabile

Si vuole stabilire se la curvatura dell'ampiezza nella regione di barriera
contenga più informazione sull'errore WKB del solo parametro eikonale
\(\varepsilon=(\ell+1/2)^{-1}\). Per una soluzione QNM esatta nel senso delle
frazioni continue di Leaver definiamo

\[
P=\varepsilon\,\operatorname{Im}\frac{\psi'}\psi,
\qquad
Q_M=-\varepsilon^2\frac{A''}{A},
\qquad A=|\psi|,
\]

e il numero adimensionale

\[
\mathcal E_M=
\frac{\int w\,|Q_M|\,dx_*}
{\int w\left(|\varepsilon\omega|^2+
\varepsilon^2|V|+|P|^2\right)dx_*}.
\tag{1}
\]

La finestra è gaussiana,
\(w=\exp[-(x_*-x_{*0})^2/(2\sigma^2)]\), centrata al massimo di \(V\), con
\(\sigma=\sqrt{V_0/(-V_0'')}\). L'integrale è troncato a
\(|x_*-x_{*0}|\leq2.25\sigma\). Questa prescrizione evita di adattare la
finestra dopo aver visto le correlazioni.

## 2. Indipendenza del benchmark

Il profilo impiegato nella (1) non viene generato dalla WKB di barriera.

1. La frequenza complessa è trovata con la frazione continua di Leaver.
2. La serie di Frobenius di Leaver fissa il rapporto \(\psi'/\psi\) sul bordo
   vicino all'orizzonte.
3. L'ODE di Regge–Wheeler viene integrata attraverso la barriera.
4. Soltanto alla fine si confrontano WKB1 e WKB3 con la frequenza di Leaver.

Il valore analitico di \(A''/A\), ottenuto dalla Riccati esatta, è confrontato
con differenze finite su \(|\psi|\). Il massimo residuo relativo è
\(4.99\times10^{-6}\) per \(n=0\), \(3.13\times10^{-4}\) per \(n=1\) e
\(4.74\times10^{-5}\) per \(n=2\). Dimezzare la profondità della frazione
continua cambia le frequenze di meno di \(2.2\times10^{-16}\) in relativo.

## 3. Risultati

La scansione comprende i settori scalare, elettromagnetico e gravitazionale
assiale, con \(1\leq\ell\leq8\) per \(s=0,1\) e
\(2\leq\ell\leq8\) per \(s=2\): 23 modi per ogni overtone.

| \(n\) | errore | Pearson \(r\) | Spearman \(\rho\) | \(r\) parziale dato \(\varepsilon^2\) | RMSE LOOCV con \(\mathcal E_M\) [dex] | RMSE LOOCV con \(\varepsilon^2\) [dex] |
|---:|---|---:|---:|---:|---:|---:|
| 0 | WKB1 | 1.0000 | 1.0000 | 0.9997 | 0.0014 | 0.0387 |
| 0 | WKB3 | 0.9998 | 1.0000 | 0.9855 | 0.0234 | 0.0608 |
| 1 | WKB1 | 0.9978 | 0.9872 | 0.9641 | 0.0305 | 0.0517 |
| 1 | WKB3 | 0.9987 | 0.9872 | 0.9214 | 0.0527 | 0.0953 |
| 2 | WKB1 | 0.9999 | 1.0000 | 0.9959 | 0.0059 | 0.0670 |
| 2 | WKB3 | 0.9989 | 1.0000 | 0.9925 | 0.0421 | 0.1129 |

Qui il coefficiente parziale è calcolato sui residui dopo regressione lineare
in \(\log\varepsilon^2\). La LOOCV rifà la regressione lasciando fuori un modo
alla volta. Per il fondamentale, variare la larghezza gaussiana fra
\(0.75\sigma\) e \(1.25\sigma\) lascia Spearman esattamente uguale a 1 e porta
la RMSE WKB1 soltanto da 0.0013 a 0.0015 dex. Il segnale non è quindi prodotto
da una scelta molto fine della finestra.

Se i tre overtone vengono mescolati senza condizionare su \(n\), la correlazione
globale scende a \(r=0.921\) per WKB1 e \(r=0.673\) per WKB3. Questo è un limite
utile: \(\mathcal E_M\) completo misura soprattutto quanto contenuto viene
omesso dall'eikonale. Dopo che WKB3 ha già cancellato parte di quel contenuto,
il suo errore residuo dovrebbe essere confrontato con un **resto di Madelung
dipendente dall'ordine**, non con \(Q_M\) completo.

## 4. Interpretazione prudente

Il test sostiene una versione precisa e limitata dell'ipotesi:

> A overtone fissato, il funzionale integrato di Madelung ordina gli errori WKB
> in Schwarzschild meglio del solo scaling \(1/L^2\), includendo la dipendenza
> dallo spin del potenziale.

Non è ancora un predittore economico: per calcolare \(\mathcal E_M\) in questo
benchmark si è usato il modo di Leaver. Il risultato dimostra che il profilo
contiene l'informazione cercata, non che la si possa estrarre senza risolvere il
problema esatto. I prossimi tentativi che possono falsificare l'idea sono:

1. sostituire il profilo esatto con una forma uniforme parabolico-cilindrica;
2. costruire il resto \(Q_M-Q_M^{(N)}\) appropriato all'ordine WKB;
3. ripetere il protocollo, senza riadattarlo, su Reissner–Nordström o
   Schwarzschild–de Sitter;
4. verificare la dipendenza dalla variabile master e dalla finestra.

## 5. Riproduzione

```text
python3.13 calculations/static_madelung_benchmark.py
python3.13 calculations/static_madelung_benchmark.py --n 1
python3.13 calculations/static_madelung_benchmark.py --n 2
python3.13 calculations/static_madelung_benchmark.py --n 0 --n 1 --n 2
python3.13 -m unittest calculations/test_static_madelung_benchmark.py
```

Le identità simboliche indipendenti sono controllate dal kernel Mathematica
13.3.1 con `calculations/verify_formalism.wl`. Il controllo comprende anche la
riduzione radiale scalare esatta di Kerr, il limite Schwarzschild e il termine
lineare di frame dragging slow Kerr.
