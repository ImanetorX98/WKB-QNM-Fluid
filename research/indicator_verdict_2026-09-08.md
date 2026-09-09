# Verdetto sull'indicatore E_M: robustificazione, confronti e ridondanza

**Data:** 8 settembre 2026.
**Natura:** risultati **misurati**, con derivazioni analitiche indicate come tali.
**Rapporto con l'audit:** esegue i passi richiesti dal protocollo
[novelty_publication_assessment](novelty_publication_assessment_2026-09-08.md) §5,
punti 4 e 5. Non contraddice l'audit: ne conferma il verdetto e lo rende definitivo.
**Riproduzione:** `python3.13 calculations/robust_indicator_test.py`.

## 0. Sintesi

L'indicatore E_M non va rivendicato. Non perché non sia validato, ma perché è
**dimostrato ridondante**: a overtone fissato è una riscalatura costante di
\(|\Lambda_3|\), che la formula WKB3 calcola già gratuitamente dal getto del
potenziale al massimo. Il criterio d'arresto fissato dall'audit è scattato.

## 1. Verifiche indipendenti dell'audit precedente

Tutte le formule chiuse dell'audit sono state ricontrollate e sono corrette.

| Enunciato | Metodo di verifica | Esito |
|---|---|---|
| \(\omega_{\rm ex}=\sqrt{L^2-1/4}-iN\) | condizione \(C_{\rm in}=0\iff a=-n\) | esatto, scarto 0.0 |
| \(Q_M=-\tfrac{\varepsilon^2}{4}(1+\operatorname{sech}^2y)\) | simbolico da \(A_0\propto\sqrt{\cosh y}\) | differenza identicamente 0 |
| \(\mathcal E_M=\varepsilon^2(2-c_w)/(8-\varepsilon^2c_w)\) | contro il codice, \(L=2\ldots64\) | \(4\times10^{-17}\) |
| \(\nu-n=N(\sqrt{1-\varepsilon^2/4}-1)-iB\varepsilon/2\) | derivazione indipendente | coincide |
| \(E_1\sim B/(2L^2)\), \(E_3\sim N/(128L^5)\) | numerico \(L=20\ldots320\) | rapporti \(\to1.0000\) |
| 32 test (12 calculations + 20 core) | rilanciati | superati |

## 2. Perché la degenerazione Pöschl–Teller è strutturale (derivato)

L'audit constata che \(\mathcal E_M\) non aggiunge informazione indipendente da
\(\varepsilon\) nella famiglia fondamentale. Il motivo è più forte di così. Per
Pöschl–Teller con \(n=0\) valgono **due identità esatte**, qui derivate e
verificate numericamente (scarti \(0.0\) e \(2\times10^{-16}\)):

\[
\left|\frac{\omega}{L}\right|^2=1,
\qquad
P^2=\left(1-\frac{\varepsilon^2}{4}\right)\tanh^2y .
\]

Da queste il denominatore dell'indicatore vale esattamente
\(2-\tfrac{\varepsilon^2}{4}\tanh^2y\), e la formula chiusa dell'audit segue in
due righe. Non è che \(\mathcal E_M\) *risulti* funzione del solo \(\varepsilon\):
**ogni ingrediente dell'indicatore è una forma fissa moltiplicata per una
funzione di \(\varepsilon\)**. Non esiste margine geometrico perché trasporti
informazione in quella famiglia.

## 3. Il segnale apparente a n=1 è il nodo (misurato)

Su Schwarzschild lo spread di \(\mathcal E_M/\varepsilon^2\) su \(\ell=3\ldots8\)
sembra distinguere gli overtone:

| | n=0 | n=1 |
|---|---|---|
| s=0 | 1.6% | 85.3% |
| s=2 | 13.0% | 83.9% |

Ma l'integrale a \(n=1\) è dominato dal picco di \(|Q_M|\) al nodo, e la
dominanza **cresce con \(\ell\)** (79% a \(\ell=4\), 89% a \(\ell=8\)): è questo
a produrre la crescita apparente. Sostituendo l'integrale con statistiche robuste
al polo, la dipendenza da \(\ell\) collassa sul livello di \(n=0\):

| variante | s=0, n=1 | s=2, n=1 |
|---|---|---|
| integrale | 85.3% | 83.9% |
| **mediana pesata** | **1.9%** | **14.6%** |
| troncato al 90° percentile | 16.8% | 33.0% |
| escisso \(\pm0.25\sigma\) dal nodo | 17.9% | 35.6% |

Il difetto del nodo segnalato dall'audit al §5.3 e la correlazione riportata al
§5.1 non sono due problemi distinti: **sono lo stesso problema.**

## 4. Test a ε fissato con forma variata (misurato)

Il protocollo richiede una validazione su un parametro di forma indipendente da
\(\varepsilon\). A \(\ell\) fissato lo spin \(s=0,1,2\) cambia la barriera
tramite \(2(1-s^2)/x^3\) lasciando \(\varepsilon=1/(\ell+\tfrac12)\) identico, e
le frequenze esatte restano disponibili da Leaver. Su \(\ell=3\ldots8\), 18
coppie:

| predittore | concordanza | errore medio sul rapporto di errore |
|---|---|---|
| \(\mathcal E_M\) (mediana) | 18/18 | 0.0179 |
| \(\Delta_2=|\omega_3-\omega_1|/2\) | 18/18 | 0.0185 |
| \(|\Lambda_3|\) | 18/18 | 0.0187 |

A \(\varepsilon\) fissato un predittore basato su \(\varepsilon^2\) non discrimina
per costruzione, quindi il 18/18 di \(\mathcal E_M\) è un risultato reale — ma
**non è distinguibile** da quello dei due controlli, che sono gratuiti.

## 5. Il risultato decisivo: E_M è Λ₃ riscalato (misurato)

| overtone | \(\mathcal E_M/|\Lambda_3|\), media su \(s=0,2\) e \(\ell=3\ldots8\) | dispersione relativa |
|---|---|---|
| n=0 | 1.0426 | **0.15%** |
| n=1 | 3.95 | 1.6% |

A overtone fissato l'indicatore è proporzionale a \(|\Lambda_3|\) con costante
indipendente da \(\ell\) e dallo spin. La costante cambia con \(n\), ma \(n\) è
un input noto e non costituisce informazione.

Conseguenza sui costi: \(\Lambda_3\) è un sottoprodotto algebrico della formula
WKB3, calcolato dal getto di \(V\) nel massimo. \(\mathcal E_M\) richiede la
frequenza esatta di Leaver, i dati iniziali di Frobenius e l'integrazione
dell'ODE attraverso la barriera. Stessa informazione, costo incomparabile.

## 6. Nota su Δ_k

Data la definizione citata dall'audit, \(\Delta_k=|\omega_{k+1}-\omega_{k-1}|/2\),
la quantità \(|\omega_3-\omega_1|/2\) **è** un \(\Delta_k\) legittimo con \(k=2\),
non un surrogato. Il repository può quindi già eseguire il confronto con
\(\Delta_2\) senza implementare nulla; solo \(\Delta_3\) richiede WKB2 e WKB4.
L'audit è corretto nel vietare l'etichetta \(\Delta_3\), ma più pessimista del
necessario su ciò che è già disponibile.

## 7. Conseguenze

1. **Non rivendicare \(\mathcal E_M\)** in nessuna forma, nemmeno robustificata.
   La robustificazione elimina l'artefatto ma con esso l'apparente segnale.
2. **Il risultato negativo è però pubblicabile come tale**, ed è più utile della
   rivendicazione che sostituisce: mostra che la riscrittura di Madelung
   dell'ampiezza WKB riproduce la correzione \(\Lambda_3\) invece di estenderla.
   È un enunciato preciso sul contenuto informativo del formalismo.
3. Kerr e Vaidya non vanno affrontati per salvare l'indicatore: non c'è
   indicatore da salvare.
4. Resta in piedi il criterio di ordinamento in \(\varepsilon\) (settore bosonico
   pendenza 2, fermionico 1), che è indipendente da \(\mathcal E_M\) e non è
   toccato da questo verdetto.
