# Validazione dei gate A–D del brief del 12 settembre

**Esito: A, B, C e D chiusi.** Nessun gate fallito, quindi nessun arresto.
Dati completi in `claude_compact_data.json`. Nessun commit eseguito.

Convenzioni comuni: \(V_0=10\,b(x)\) con \(b\) il bump liscio a supporto
\([-1,1]\); perturbazione \(\eta\,b((x-c)/w)\); \(e^{-i\omega t}\); bordi uscenti
esatti \(\psi'(a)=-i\omega\psi(a)\), \(\psi'(c)=+i\omega\psi(c)\), leciti perché
\(V=0\) fuori dal supporto; \(\varepsilon=1\); \(J=[-0.8,0.8]\);
\(w=e^{-x^2}\); \(\psi(a)=1\).

---

## Gate B — la quadratura (eseguito per primo, come prescritto)

`calculations/claude_compact_quadrature.py` — **4.5 s**

### Diagnosi

Il problema non era la precisione ma la **regolarità**. \(U=\int_J w|Q_M|\,dx\)
ha integrando \(C^0\): \(|Q_M|\) ha angoli dove \(Q_M\) cambia segno, e
\(\dot U=\int_J w\,\mathrm{sgn}(Q_M)\,\dot Q_M\,dx\) ha **salti** negli stessi
punti. Una quadratura globale liscia su un integrando \(C^0\) converge a una
potenza bassa del passo: nessun raffinamento uniforme la salva. La cura è
spezzare sugli zeri e usare quadratura adattiva in ciascun tratto, dove
l'integrando è analitico.

### Risultati

**Zeri di \(Q_M\) in \(J\): esattamente 2**, in \(\pm0.3586793241\), identici su
griglie di ricerca da 1001 a 8001 punti. Il conteggio è stabile.

**\(\psi\) non ha zeri**: \(\min|\psi|=8.41\times10^{-1}\) su \(J\), a
\(x=+0.61\). Gli zeri trovati sono di \(Q_M\), non di \(\psi\): la distinzione
che il brief chiede di registrare.

**\(k\) con due metodi**, ODE \(k'+2zk=-2\omega\) contro formula integrale:
accordo fra \(5\times10^{-13}\) e \(7\times10^{-12}\) su cinque punti di prova.

**Derivata analitica con quadratura spezzata: 0.001894279539**, stabile alla
dodicesima cifra fra tolleranze ODE \(10^{-12}\) e \(10^{-14}\).

### Confronto con le differenze centrali

A tolleranza \(10^{-14}\), confrontando con \(E\) calcolata anch'essa a
quadratura spezzata:

| passo \(\eta\) | FD centrale | scarto relativo |
|---|---|---|
| 1e−3 | 0.001894335099 | 2.93e−5 |
| 3e−4 | 0.001894287496 | 4.20e−6 |
| **1e−4** | 0.001894279553 | **7.43e−9** |
| 3e−5 | 0.001894278775 | 4.04e−7 |

L'andamento è quello atteso: troncamento \(O(\eta^2)\) in discesa, poi
arrotondamento in risalita a \(3\times10^{-5}\). Target del brief \(10^{-5}\):
**superato di tre ordini**.

### Il plateau, diagnosticato

A tolleranza ODE \(10^{-12}\) l'accordo si fermava a \(1.2\times10^{-5}\).
Stringendo a \(10^{-14}\) scende a \(7.4\times10^{-9}\): **era la tolleranza
della risonanza, non la quadratura**. Il brief chiedeva di non attribuire il
plateau a una sola fonte a priori; la fonte è identificata variandola.

### Sul valore discreto

Il continuo dà \(0.0018942795\), il funzionale discreto di Codex
\(0.00187362\): differenza 1.1%. Sono livelli di approssimazione diversi e non
sono stati forzati a coincidere, come il brief prescrive.

---

## Gate A — secondo e terzo metodo spettrale

`calculations/claude_compact_spectrum.py` — **11.8 s**

Il brief esclude che «cambiare linguaggio/CAS allo stesso shooting» conti come
indipendenza. Sono stati implementati due metodi con discretizzazione diversa:

* **matching dai due lati** con Wronskiano normalizzato in un punto interno
  \(x_m=0.7\): stessa ODE, ma residuo e propagazione degli errori diversi;
* **collocazione di Chebyshev**, problema agli autovalori **quadratico** in
  \(\omega\) (i bordi dipendono da \(\omega\)), linearizzato in forma compagna e
  risolto in blocco. Nessuna iterazione, nessun valore iniziale.

### Classificazione dello spettro

Il terzo metodo restituisce lo spettro senza inseguire una radice. Nella finestra
\(\operatorname{Re}\omega\in(0.5,12)\), \(\operatorname{Im}\omega>-3\):

| | |
|---|---|
| \(3.4232955-0.6073876\,i\) | la meno smorzata |
| \(4.5194795-1.6474855\,i\) | |
| \(6.0830050-2.3869388\,i\) | |
| \(7.7700169-2.9310218\,i\) | |

Quattro risonanze, conteggio stabile a \(N=140,180,220\). La prima è quella del
brief. **«Fondamentale» resta una classificazione relativa alla finestra**, non
un'etichetta intrinseca.

### Accordo fra metodi

| \(\eta\) | shooting − matching | shooting − collocazione (\(N=220\)) |
|---|---|---|
| 0 | 1.28e−14 | 1.22e−08 |
| +0.001 | 8.89e−12 | 1.23e−08 |
| −0.001 | 7.45e−12 | 1.26e−08 |
| +0.01 | 6.64e−13 | 1.19e−08 |
| −0.01 | 2.94e−13 | 2.81e−08 |

La collocazione era al limite del gate a \(N=220\). È **troncamento di
Chebyshev**, non un disaccordo: il bump è \(C^\infty\) ma **non analitico**
(singolarità essenziale a \(y=\pm1\)), quindi la convergenza spettrale è
sub-geometrica. Raffinando:

| \(N\) | 140 | 180 | 220 | 260 | 300 | 340 |
|---|---|---|---|---|---|---|
| scarto | 2.4e−7 | 5.5e−8 | 1.2e−8 | 3.1e−9 | 1.1e−9 | **2.1e−10** |

**Gate superato**: \(2.1\times10^{-10}\) contro la soglia \(10^{-8}\).

### Variazioni risolte

| \(\eta\) | \(|\Delta\omega|\) | errore fra metodi | rapporto |
|---|---|---|---|
| 1e−3 | 1.765e−4 | 8.9e−12 | **2.0e+7** |
| 1e−2 | 1.779e−3 | 6.6e−13 | **2.7e+9** |

Il brief chiedeva almeno \(10\times\).

### Norma generalizzata

\(d=B/N\) con termini di superficie, traslando i bordi in regioni dove \(V=0\):

| \([a,c]\) | \(d\) | scarto dal riferimento |
|---|---|---|
| \([-1.0,\,2.4]\) | \(0.010361818212-0.175999169467\,i\) | — |
| \([-1.5,\,2.4]\) | idem | 5.1e−15 |
| \([-1.0,\,3.0]\) | idem | 2.4e−15 |
| \([-2.0,\,3.5]\) | idem | 9.9e−13 |
| \([-3.0,\,4.0]\) | idem | 3.0e−14 |

**Invariante a \(10^{-13}\) su cinque domini.** È la verifica che i termini di
superficie sono quelli giusti: senza, spostare il bordo cambierebbe \(d\).

Contro differenze centrali di \(\omega\):

| \(\eta\) | 1e−2 | 3e−3 | 1e−3 | 3e−4 | 1e−4 |
|---|---|---|---|---|---|
| scarto rel. | 1.26e−4 | 1.13e−5 | 1.25e−6 | 1.23e−7 | **7.79e−9** |

Discesa \(O(\eta^2)\) pulita. Target \(10^{-5}\): superato di tre ordini.

---

## Gate C — fattorizzazione e rango

`calculations/claude_compact_rank.py` — **15.9 s**

Dominio esteso a \([-1,3.2]\) per contenere il supporto fino a 2.8, con il bordo
libero fuori da ogni bump. Massimo dominante \(=10.0000000000\), invariato.

### \(h=d\,k\) sulle sei perturbazioni esterne

| centro | largh. | \(d\) | L2 rel. | err. ass. max |
|---|---|---|---|---|
| 1.6 | 0.20 | \(-0.022327+0.084303\,i\) | 1.20e−9 | 9.9e−10 |
| 1.6 | 0.40 | \(-0.048318+0.097093\,i\) | 1.77e−9 | 1.8e−9 |
| 2.0 | 0.20 | \(-0.020353-0.140303\,i\) | 3.64e−9 | 4.9e−9 |
| 2.0 | 0.40 | \(+0.010362-0.175999\,i\) | 6.10e−9 | 1.0e−8 |
| 2.4 | 0.20 | \(+0.119877+0.196842\,i\) | 1.08e−8 | 2.4e−8 |
| 2.4 | 0.40 | \(+0.096701+0.269803\,i\) | 1.68e−8 | 4.6e−8 |

### Rango

Valori singolari della matrice \(6\) colonne \(\times\) \(801\) punti:

\[
1.529\times10^{2},\;\; 5.354\times10^{1},\;\; \mathbf{6.52\times10^{-7}},\;\;
2.04\times10^{-7},\;\; 1.53\times10^{-8},\;\; 1.41\times10^{-9}
\]

\(\sigma_3/\sigma_1=4.26\times10^{-9}\).

**Rumore stimato**, non soglia arbitraria: ricalcolando la stessa colonna con
passo \(\eta\) triplicato si ottiene uno scarto di \(9.68\times10^{-6}\). Quindi
\(\sigma_3\) sta **quindici volte sotto** il rumore di differenziazione. Il brief
chiedeva esattamente questo confronto, e non chiedeva che \(\sigma_2\) fosse non
nullo — non lo è: \(\sigma_2/\sigma_1=0.35\), le due direzioni sono entrambe
attive.

### Predizione con \(\{f_1,f_2\}\), senza rifittare \(d\)

\(f_1=2[\operatorname{Re}\omega-p\operatorname{Im}k]\),
\(f_2=-2[\operatorname{Im}\omega+p\operatorname{Re}k]\), calcolate **una volta**
sul problema imperturbato. L2 relativo: 2.24e−9, 3.12e−9, 4.91e−9, 8.97e−9,
1.88e−8, 2.45e−8.

### Controllo negativo

Perturbazione **interna** (centro 0.3, larghezza 0.1), che cambia l'ODE dentro
\(J\): \(h\) contro \(d\,k\) dà L2 relativo **8.28e−1**, e la colonna contro la
base **1.24**. L'ipotesi fallisce dove deve fallire — sette ordini di grandezza
di separazione dai casi esterni.

---

## Gate D — utilità predittiva

`calculations/claude_compact_predictive.py` — **18.1 s**, 30 righe

Il test è netto per costruzione: le bump sono nulle su \([-1,1]\), quindi il
getto centrale è invariato **a tutti gli ordini** e

\[
\Lambda_2=0.059292706128,\quad \Lambda_3=0.001562500000,\quad
\omega_{\rm WKB3}=3.257504253653-0.486141916079\,i
\]

sono **costanti** su tutte le 30 righe. L'errore della WKB3 varia solo perché
varia il bersaglio.

### Definizioni, come richiesto

* \(E_{\rm int}\): rapporto degli integrali pesati, quadratura spezzata sugli
  zeri (gate B);
* \(E_{\rm med}\): rapporto delle **mediane pesate**, dove \(\mathrm{med}_w(g)\)
  è il valore di \(g\) nel punto in cui la somma cumulata dei pesi, ordinata per
  \(g\) crescente, raggiunge metà del peso totale. Quantile **discreto**, nessuna
  interpolazione, griglia uniforme di 4001 punti.
* La mediana è usata **solo come valore, mai derivata**: la formula del §4 non
  vale per essa.

### Esito

| | correlazione di Pearson con l'errore WKB3 | escursione relativa |
|---|---|---|
| \(E_{\rm int}\) | **−0.2151** | 6.45e−2 |
| \(E_{\rm med}\) | **−0.4247** | 7.36e−2 |
| errore WKB3 | — | 5.37e−2 |

Le escursioni sono **comparabili** — il diagnostico si muove quanto l'errore —
ma la correlazione è debole **e di segno sbagliato**.

> Il diagnostico **risponde** alle perturbazioni esterne, con la struttura di
> rango due del gate C, ma **non segue** l'errore della WKB3. È una diagnosi a
> posteriori della soluzione, non una previsione.

Nessun fit è stato eseguito, quindi la questione training/test non si pone:
si confrontano andamenti su un campione unico e dichiarato.

---

## Costi, fallimenti, comandi

| script | tempo | esito |
|---|---|---|
| `claude_compact_quadrature.py` | 4.5 s | gate B |
| `claude_compact_spectrum.py` | 11.8 s | gate A |
| `claude_compact_rank.py` | 15.9 s | gate C |
| `claude_compact_predictive.py` | 18.1 s | gate D |

**Fallimenti da registrare.** Due, entrambi miei e corretti:
`eikonal_frequency` chiamata con la firma sbagliata, e `ndarray.ptp` rimosso in
NumPy 2. Nessun fallimento fisico o di convergenza.

```sh
python3.13 calculations/claude_compact_quadrature.py
python3.13 calculations/claude_compact_spectrum.py
python3.13 calculations/claude_compact_rank.py
python3.13 calculations/claude_compact_predictive.py
```

`WolframKernel` è presente e con licenza; non è servito, perché il plateau si è
risolto identificandone la causa invece di chiedere più cifre.

## Limiti dichiarati

* Tutto vale per **questa** barriera compatta, non per Schwarzschild o Kerr.
* Il rango due è verificato su sei perturbazioni esterne di una sola famiglia;
  è coerente con la proposizione, non una sua dimostrazione.
* Non è stata toccata la questione dell'originalità, né il manoscritto, né gli
  handoff precedenti.
* L'estensione E sulla derivata spettrale di Kerr **non è stata avviata**: i
  gate A–B che la abilitano sono ora chiusi.
