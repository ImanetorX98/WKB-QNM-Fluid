# Previsione analitica di Q_M, e il suo condizionamento

**Data:** 9 settembre 2026.
**Natura:** teoria **derivata**, validazione e sensibilità **misurate**.
**Origine:** richiesta di approfondire l'analisi teorica invece di limitarsi a
misurare pendenze, per capire quanto potere predittivo abbia il formalismo.
**Riproduzione:** `python3.13 calculations/madelung_wkb_prediction.py`.

---

## 1. Il Teorema 2 non copre Kerr

Il Teorema 2 del manoscritto — la serie WKB non contiene ordini dispari — è
dimostrato assumendo \(q\) indipendente da \(\varepsilon\). In Kerr **il
potenziale riscalato ha una sua espansione**,
\(\hat q=q_0+\varepsilon q_1+\varepsilon^2q_2\) con \(q_1=-\Delta A_1/H^2\).
Rifacendo l'espansione del momento \(u=u_0+\varepsilon u_1+\varepsilon^2u_2\):

\[
u_0=\sqrt{q_0},\qquad
\boxed{u_1=\frac{q_1}{2\sqrt{q_0}}\neq0},\qquad
u_2=\frac{q_2-u_1^2+Q[u_0]}{2u_0}.
\]

**La parità si rompe.** Questo va dichiarato nel manoscritto: l'enunciato del
Teorema 2 vale per potenziali \(\varepsilon\)-indipendenti, e Kerr non lo è.

Conviene allora **non separare** \(\hat q\) in ordini, e applicare la serie WKB
standard al \(\hat q\) complesso completo.

---

## 2. La previsione, e due errori che costa evitare

\[
u=\sqrt{q}+\varepsilon^2\,\frac{5q'^2-4qq''}{32\,q^{5/2}}+O(\varepsilon^4),
\qquad {}'=\frac{d}{dx_*} .
\]

**Errore 1 — il fattore esponenziale.** Con \(\omega\) complessa l'ampiezza reale
della decomposizione di Madelung **non** è \(|u|^{-1/2}\). Da
\(\psi=u^{-1/2}e^{i\int u/\varepsilon}\) segue

\[
\ln A=-\tfrac12\ln|u|-\frac1\varepsilon\int\operatorname{Im}u,
\qquad
Q_M=-\varepsilon^2\left[(\ln A)''+\big((\ln A)'\big)^2\right].
\]

Omettere il secondo termine sbaglia \(Q_M\) di un fattore ~50. Il termine
\(-\operatorname{Im}(u)/\varepsilon\) è \(O(1)\) perché per un QNM
\(\operatorname{Im}\hat\Omega=O(\varepsilon)\); il suo quadrato contribuisce a
\(Q_M\) allo stesso ordine del funzionale di Bohm.

**Errore 2 — il ramo.** \(u\) e \(-u\) sono entrambe soluzioni. Usare lo stesso
ramo sui due lati della barriera produce un errore del 10% **che non svanisce
al crescere di \(L\)** — quindi indistinguibile da un difetto di teoria se non
si guarda la convergenza.

Entrambi gli errori sono stati commessi e corretti in questa analisi.

---

## 3. Validazione su un caso esattamente risolubile

Pöschl–Teller \(n=0\), dove \(Q_M=-\tfrac{\varepsilon^2}{4}(1+\operatorname{sech}^2y)\)
è esatto. Errore relativo mediano su \(1<y<2.5\), ramo uscente:

| \(L\) | \(u=\sqrt q\) | \(u=\sqrt q+\varepsilon^2u_2\) |
|---|---|---|
| 8 | \(3.9\times10^{-3}\) | \(6.9\times10^{-4}\) |
| 20 | \(6.3\times10^{-4}\) | \(4.0\times10^{-5}\) |
| 50 | \(1.0\times10^{-4}\) | \(5.1\times10^{-6}\) |
| 120 | \(1.8\times10^{-5}\) | \(8.9\times10^{-7}\) |
| 300 | \(2.8\times10^{-6}\) | \(1.7\times10^{-7}\) |

L'errore di testa va come \(\varepsilon^2\); il termine \(u_2\) guadagna un
ordine. **La previsione non richiede di integrare l'ODE**: usa solo il
potenziale complesso e la frequenza.

Controllo analitico indipendente: per PT si verifica a mano che
\(A_{\rm WKB}=|u|^{-1/2}\exp(-\varepsilon^{-1}\!\int\operatorname{Im}u)=\sqrt{\cosh y}\),
che è l'ampiezza esatta.

---

## 4. Il risultato che conta: Q_M è mal condizionato nella frequenza

Perturbando \(\omega\) di una frazione relativa e misurando lo scarto fra
previsione e numerica (Schwarzschild, \(\ell=70\), finestra \(20<r<50\)):

| perturbazione relativa su \(\omega\) | errore su \(Q_M\) | amplificazione |
|---|---|---|
| 0 | \(5.5\times10^{-4}\) | — |
| \(10^{-5}\) | \(1.0\times10^{-3}\) | ~100 |
| \(10^{-4}\) | \(1.2\times10^{-2}\) | ~120 |
| \(10^{-3}\) | \(1.4\times10^{-1}\) | ~140 |
| \(3\times10^{-3}\) | \(4.8\times10^{-1}\) | ~160 |

**Un errore relativo di \(10^{-3}\) sulla frequenza produce un errore del 14% su
\(Q_M\).** L'amplificazione è di due ordini di grandezza ed è stabile.

### 4.1 Conseguenza sul potere predittivo

Questo chiude, dal lato teorico, la questione aperta dal verdetto su
\(\mathcal E_M\). Un funzionale che amplifica di \(\sim10^2\) l'errore sulla
frequenza **non può essere un predittore della frequenza**: per calcolarlo con
tre cifre servono cinque cifre di ciò che si vorrebbe prevedere. Non è una
difficoltà numerica da superare con più risoluzione, è il condizionamento del
problema.

Il potere predittivo del formalismo di Madelung sta quindi **nell'ordinamento**,
che è robusto — le pendenze in \(\varepsilon\) sono sopravvissute a errori su
\(\omega\) di ordine \(10^{-3}\) — e **non nei valori**, che non lo sono.

---

## 5. Diagnosi del comportamento erratico su Kerr

Il \(Q_M\) integrato su Kerr risultava irregolare e non monotono in \(\ell\).
Tre ipotesi sono state testate e **due scartate**:

| ipotesi | esito |
|---|---|
| risoluzione insufficiente della griglia | **scartata**: converge al 2% raddoppiando i punti |
| quasi-zeri di \(\psi\) nella finestra | **scartata**: zero minimi di \(|A|\) in \(20<r<50\) |
| autovalore sferoidale troncato alla parte reale | **scartata**: correggendolo il risultato non cambia |
| **errore sulla frequenza eikonale** | **confermata** |

Lo scarto osservato su Kerr (\(a=0.9\), \(\ell=70\): \(2.3\times10^{-1}\))
corrisponde, sulla tabella del §4, a un errore relativo su \(\omega\) di
\(\approx1.5\times10^{-3}\). È l'ordine atteso: `eikonal_solution` stima
\(A_0\) a \(\ell_{\rm rif}=400\), e per \(a\neq0\) l'errore residuo è
\(A_1/400\approx5\times10^{-4}\), contro \(1.6\times10^{-6}\) per \(a=0\).
La non monotonia in \(\ell\) si spiega con la fase dell'ammistione di ramo, che
dipende da \(\ell\).

Il bordo all'orizzonte era invece un problema reale e già corretto: la pendenza
di \(Q_M\) passa da 0.569 a **2.026** spostando l'offset da \(10^{-3}\) a
\(10^{-6}\); sotto \(10^{-7}\) l'integratore cede.

---

## 6. Stato del passo B e cosa serve per chiuderlo

**Stabilito.** Il termine di rotazione ha pendenza 0.89–0.95 con rotazione e
2.04 senza, coerente con il §9 del manoscritto. Su Schwarzschild
(\(a=0\)) previsione e numerica concordano a \(1.9\times10^{-4}\) a
\(\ell=100\), migliorando con \(\ell\).

**Non stabilito.** Il confronto su Kerr resta limitato dall'errore sulla
frequenza. Per chiuderlo serve una frequenza **autoconsistente al potenziale
effettivamente usato a quel \(\ell\)** — non la frequenza eikonale. Basta un
WKB3 di barriera costruito sul potenziale esatto, con iterazione sull'autovalore
sferoidale: molto meno di un Leaver per Kerr.

Finché non c'è, i numeri di \(Q_M\) su Kerr **non vanno riportati come
risultato**. Quelli sul potenziale (§9 del manoscritto) sono invece indipendenti
da questo problema, perché non richiedono di integrare l'ODE.

---

## ADDENDUM del 10 settembre — la tabella del §4 non è riproducibile

Tentando di disegnarne la figura è emerso che **il codice che genera la tabella
di amplificazione qui sopra non è nel repository.** `madelung_wkb_prediction.py`
contiene solo la validazione su Pöschl–Teller e rimanda a questo file. La
tabella «Schwarzschild, ℓ=70, finestra 20<r<50» non ha uno script.

Due ricostruzioni indipendenti, entrambe negative.

**Su Pöschl–Teller**, dove $Q_M$ è noto in forma chiusa, l'amplificazione vale
**1.41 e non dipende da $L$** (misurata a $L=8,20,50,120,300,700$). Il caso è
però speciale: $Q_M=-\tfrac{\varepsilon^2}{4}(1+\operatorname{sech}^2y)$ **non
dipende dalla frequenza**, quindi non può esibire il fenomeno. Test inconcludente
per costruzione.

**Su Schwarzschild** (`calculations/madelung_conditioning_schwarzschild.py`, nuovo)
l'amplificazione risulta $\lesssim2$ a ogni $\ell$ provato — ma la ricostruzione
è essa stessa difettosa: lo scarto a frequenza esatta **cresce come $L^2$**
(3.1e-3, 1.2e-2, 3.7e-2, 7.8e-2 a $\ell=20,40,70,100$) invece di calare. Il
pavimento è indipendente dalla risoluzione a tre valori di griglia, quindi è uno
scarto di modello, non rumore. La sostituzione di Langer, primo sospetto, non lo
cambia di una cifra.

### Stato

L'enunciato del §5 del manoscritto — l'ostruzione di condizionamento, che due
giorni fa ho proposto come **risultato di testa** — non è al momento verificabile:
il codice originale è assente e la ricostruzione non converge, quindi non può né
confermarlo né smentirlo.

Va risolto prima di qualunque invio. Le opzioni sono ricostruire la misura
correttamente, oppure ridurre l'enunciato a ciò che è dimostrabile.

**La figura 2 non va prodotta finché il numero non è verificato.**
