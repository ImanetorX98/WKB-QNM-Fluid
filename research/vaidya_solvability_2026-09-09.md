# Vaidya: la condizione di solvibilità, e dove si ferma

**Data:** 9 settembre 2026.
**Natura:** **misurato**, con esito **parziale** e ostruzione identificata.
**Riproduzione:** `python3.13 calculations/vaidya_solvability.py`.

---

## 1. La domanda

A ordine \(\dot M\) la correzione \(Z_1\) risolve \(L_M Z_1=-2\partial_r\partial_M Z\).
Ma alla frequenza QNM l'operatore congelato è **singolare**: \(Z\) sta nel suo
nucleo. La forzatura è risonante, e la componente della sorgente parallela al
modo non produce \(Z_1\) ma uno spostamento di frequenza. Solo la componente
ortogonale genera una correzione vera — ed è lì che potrebbe entrare una
dipendenza dalla storia, perché invertire l'operatore fuori dal nucleo richiede
una funzione di Green, quindi un integrale su \(v'<v\).

Se la proiezione fosse **esattamente nulla**, non ci sarebbe \(Z_1\) e non ci
sarebbe memoria a nessun ordine. È l'esito banale da escludere per primo.

## 2. Impostazione

Da \(2\psi_{vr}+\partial_r(f\psi_r)-U_\ell\psi=0\) con
\(\psi=Z(r;M)e^{-i\int\omega\,dv}\), l'ordine \(\dot M^0\) dà

\[
L_M Z\equiv fZ''+(f'-2i\omega)Z'-U_\ell Z=0,\qquad {}'=\partial_r .
\]

Il fattore integrante che porta \(L_M\) in forma di Sturm–Liouville soddisfa
\(\mu'/\mu=-2i\omega/f\), cioè \(\mu=e^{-2i\omega r_*}\). Con quel peso \(L_M\) è
autoaggiunto rispetto al prodotto **bilineare senza coniugazione**
\(\langle u,v\rangle=\int\mu\,uv\,dr\), che è quello corretto per problemi di
risonanza non autoaggiunti.

**Semplificazione di scala.** Schwarzschild ha un solo parametro, quindi
\(\partial_M\) è una derivata di scala: con \(x=r/M\), \(\Omega=M\omega\),
\(Z(r;M)=G(x)\) e \(G=e^{i\Omega x_*}R\),

\[
\partial_M Z=-\frac{x}{M}G',\qquad
\partial_r\partial_M Z=-\frac{1}{M^2}\left(xG'\right)' .
\]

La proiezione da valutare è dunque, a meno di potenze di \(M\),

\[
I=\int e^{-2i\Omega x_*}\,G\,(xG')'\,dx .
\]

## 3. Risultato

**(a) La cancellazione esatta è esclusa.** Il rapporto fra il modulo
dell'integrale e l'integrale del modulo dell'integrando è di ordine \(10^{-1}\),
non \(10^{-3}\) o meno:

| \(\ell\) | \(M\omega\) | \(|I|/\int|{\rm integrando}|\) |
|---|---|---|
| 2 | \(0.4836-0.0968i\) | 0.244 |
| 3 | \(0.6754-0.0965i\) | 0.176 |
| 4 | \(0.8674-0.0964i\) | 0.138 |
| 6 | \(1.2519-0.0963i\) | 0.096 |

Non c'è ortogonalità accidentale: **la sorgente ha sovrapposizione genuina con
il modo**. Ne segue che \(Z_1\) esiste ed è non nulla, e che la via alla memoria
non è chiusa dall'esito banale.

**(b) Ma l'integrale non converge.** Allargando la finestra, \(|I|\) diverge come
atteso per un QNM, che cresce a entrambi i bordi:

| finestra (\(\ell=2\)) | \(4<x<10\) | \(<20\) | \(<30\) | \(<40\) |
|---|---|---|---|---|
| \(|I|\) | \(1.3\times10^{2}\) | \(2.8\times10^{3}\) | \(3.4\times10^{4}\) | \(3.5\times10^{5}\) |

## 4. La regolarizzazione di Leung, e dove si ferma davvero

**Aggiornamento del 9 settembre, dopo aver letto Leung *et al.* 1998**
(J. Phys. A **31**, 3271, doi:10.1088/0305-4470/31/14/013 — nota: gli autori sono
**cinque**, Leung, Liu, Suen, Tam, Young, non tre come citato finora).

Lo strumento c'è ed è la loro eq. (2.16): la **norma generalizzata**

\[
\langle\phi_0|\phi_0\rangle=\int_{L_-}^{L_+}\phi_0^2\,dx
+\frac{1}{2\omega_0}\Big[D'_{+0}\phi_0^2(L_+)-D'_{-0}\phi_0^2(L_-)\Big],
\]

con \(D_{\pm0}(\omega)\) la derivata logaritmica della soluzione asintotica ai
punti di raccordo. La proprietà decisiva è dichiarata da loro: numeratore e
denominatore della (2.14) devono essere **separatamente indipendenti da
\(L_\pm\)**. È esattamente ciò che serve al nostro integrale divergente.

**Test 1 — la forma semplice non basta.** La (2.19), con \(D_{\pm0}=\pm i\omega\),
vale per potenziali **senza code**, e Leung lo dice esplicitamente. Applicandola
a Schwarzschild, che ha code di potenza, i termini di bordo riducono la
divergenza di circa tre ordini ma **non la eliminano**:

| \(L_+\) (in \(x\)) | integrale nudo | con bordi (2.19) |
|---|---|---|
| 20 | \(1.3\times10^{2}\) | \(7.4\times10^{0}\) |
| 60 | \(4.9\times10^{5}\) | \(1.7\times10^{3}\) |
| 100 | \(1.4\times10^{9}\) | \(1.7\times10^{6}\) |

Cresce ancora. Serve la (2.16) con \(D'\) vero.

**Test 2 — e \(D'\) non si ottiene per integrazione all'indietro.** Calcolando la
derivata logaritmica della soluzione uscente integrando da \(r\) grande verso
l'interno si ottiene \(D_+\approx-i\omega\) invece di \(+i\omega\): segno
opposto. Il motivo è strutturale: con \(\operatorname{Im}\omega<0\) la soluzione
uscente \(e^{i\omega x_*}\) **cresce** con \(x_*\), quindi integrando verso
l'interno è il ramo **recessivo**, e la contaminazione dal ramo entrante lo
sovrasta. Il numero misura il ramo sbagliato.

La via corretta è la **serie asintotica** della soluzione uscente a grande
\(r\), valutata in \(L_+\) — che è ciò che Leung intende con "assume that the
asymptotic regions have been solved".

**Test 3 — implementata, e funziona.** Con \(R=e^{i\omega r_*}u(r)\) l'equazione
radiale diventa

\[
\left(1-\frac2r\right)u''+\left(\frac{2}{r^2}+2i\omega\right)u'
-\left[\frac{\ell(\ell+1)}{r^2}+\frac{2(1-s^2)}{r^3}\right]u=0,
\]

e ponendo \(u=\sum_k a_k r^{-k}\), \(a_0=1\):

\[
a_k=\frac{[k(k-1)-\ell(\ell+1)]\,a_{k-1}-[2k(k-2)+2(1-s^2)]\,a_{k-2}}{2i\omega k},
\]

troncata al termine minimo, essendo asintotica. Il residuo dell'equazione scende
a \(10^{-12}\) con dodici termini, e **\(D_+\) ha ora il segno corretto**:
\(+0.0968+0.4832i\) contro \(i\omega=+0.0968+0.4836i\) a \(r=120\), con
\(D'_+\to i\) come deve. Per il lato orizzonte si usa la serie di Frobenius già
implementata in `static_madelung_benchmark.frobenius_log_derivative`.

**Test 4 — la norma regolarizzata è indipendente da \(L_+\).** Con i \(D'\) veri:

| \(L_+\) | integrale nudo | norma regolarizzata |
|---|---|---|
| 40 | \(8.5\times10^{3}\) | **10.21** |
| 70 | \(3.6\times10^{6}\) | **10.19** |

Accordo allo **0.2%** fra due punti di raccordo che differiscono di tre ordini di
grandezza nell'integrale nudo. La regolarizzazione di Leung è validata nella
nostra implementazione.

Il valore a \(L_+=70\) converge solo stringendo la tolleranza dell'ODE
(11.10 → 10.36 → 10.19 per rtol \(10^{-9},10^{-11},10^{-13}\)): **il limite è la
doppia precisione**, non la formula. A \(L_+=100\) la cancellazione richiede nove
cifre e non converge. Per \(L_+\) grandi servirebbe `mpmath`.

## 5. Cosa manca ora, con precisione

### 5.1 Tentativo sul numeratore, e perché è stato abbandonato

L'asintotica della sorgente è stata ricavata: con \(Z=e^{2i\omega r_*}u\) nella
regione esterna (il fattore è doppio perché in coordinate entranti un'onda
uscente ha fase raddoppiata),

\[
S=-2\partial_r\partial_M Z\;\sim\;2(2i\omega)^2\,r\,e^{2i\omega r_*},
\]

cioè l'integrando del numeratore ha **lo stesso tasso esponenziale della norma
più un fattore \(r\)**. Si è quindi tentato di costruire l'antiderivata
asintotica per parti ripetute e sottrarla, che è equivalente al termine di
superficie ma più diretto.

**Non converge, e il modo in cui fallisce è una trappola.** Al crescere dei
termini nell'antiderivata:

| termini | \(L_+=30\) | \(L_+=45\) | scarto relativo |
|---|---|---|---|
| 2 | \(8.9\times10^{2}\) | \(3.4\times10^{3}\) | 4.4 |
| 3 | \(4.0\times10^{4}\) | \(3.6\times10^{4}\) | 0.10 |
| 4 | \(1.1\times10^{5}\) | \(1.1\times10^{5}\) | 0.036 |
| 5 | \(3.0\times10^{9}\) | \(3.0\times10^{9}\) | **\(1.1\times10^{-5}\)** |
| 8 | \(2.1\times10^{15}\) | \(2.1\times10^{15}\) | \(2.1\times10^{-3}\) |

Lo **scarto** migliora fino a \(10^{-5}\), ma il **valore** diverge di undici
ordini di grandezza. L'indipendenza da \(L_+\) è spuria: due numeri enormi che
coincidono perché dominati dallo stesso termine divergente. Guardando la sola
colonna dello scarto si concluderebbe "converge" e si riporterebbe
\(3\times10^{9}\).

**Regola che ne segue, da applicare sempre:** l'indipendenza dal parametro di
regolarizzazione **non è sufficiente**. Il valore deve essere stabile
separatamente. Sono due controlli, non uno.

La causa tecnica è probabilmente la differenziazione **numerica** ripetuta
(`np.gradient`) nella ricorsione dell'antiderivata: ogni applicazione amplifica
il rumore, e dopo cinque iterazioni domina. La via corretta è differenziare
**analiticamente** la serie \(\sum a_k r^{-k}\), non numericamente il profilo.

### 5.2 Seconda iterazione: algebra di serie esatta

Il difetto diagnosticato al §5.1 — differenziazione numerica ripetuta — è stato
corretto scrivendo `calculations/asymptotic_series.py`: rappresentazione di ogni
funzione come serie troncata in \(1/r\), con derivate **esatte** sui
coefficienti.

**Cosa funziona ora, verificato:**

| controllo | esito |
|---|---|
| algebra di base: \(f\cdot(1/f)=1\), derivate | esatte |
| ricorsione \(V\leftarrow(f/2i\omega)(W-V')\) soddisfa \(V'+(2i\omega/f)V=W\) | residuo \(1.8\times10^{-15}\) |
| **valore stabile** al variare delle iterazioni | identico da 6 iterazioni in poi |
| serie dell'integrando contro numerica, rapporto **complesso** | \(0.9999998\) |
| costante di normalizzazione \(C\) (l'integrando è quadratico in \(Z\)) | convergente a 12 cifre |
| antiderivata contro integrale della serie, griglia fine | \(1.000000000\) |

L'asintotica della sorgente, ora esplicita: con \(Z=Ce^{2i\omega r_*}u\),

\[
P=r\left(\frac{2i\omega u}{f}+u'\right),\qquad
W=2u\left[\frac{2i\omega}{f}P+P'\right],
\]

e l'integrando del numeratore è \(C^2W(r)e^{2i\omega r_*}\), con \(W\) di
grado principale \(r^{+1}\).

**Cosa resta irrisolto.** Con tutti gli ingredienti accurati a \(10^{-7}\) o
meglio, l'integrale definito calcolato numericamente e la differenza
dell'antiderivata differiscono di \(1.65\times10^{-4}\), **in modo costante
sull'intervallo**. Un errore *relativo* costante su un integrale che cresce
esponenzialmente produce un residuo che cresce esponenzialmente: è questo a
impedire l'indipendenza da \(L_+\).

I pezzi tornano singolarmente e non tornano insieme. Il sospetto residuo è la
quadratura sulla griglia dell'ODE, ma la stima dell'errore di trapezio a quella
risoluzione è \(10^{-7}\), non \(10^{-4}\). **Senza un'ipotesi precisa non si
prosegue**: è il modo in cui, in questa sessione, sono stati prodotti sei numeri
plausibili e sbagliati.

### 5.3 Stato

Il **denominatore** della (2.14) è fatto e validato. Manca il **numeratore**:
la (2.15) richiede i termini \(\Delta_{\pm}\), che dipendono dall'asintotica
della **nostra sorgente specifica** \(-2\partial_r\partial_M Z\) nelle due
regioni. Le serie per costruirli ci sono ora entrambe — uscente e di Frobenius —
quindi è un passo delimitato, non una direzione da esplorare.

Senza quella regolarizzazione:

- si sa che \(Z_1\neq0\) e che il canale non adiabatico è aperto;
- **non** si sa come la sorgente si ripartisca fra spostamento di frequenza
  (componente parallela) e correzione genuina (ortogonale);
- e quindi **non** si può ancora dire se \(\mathcal M_{vr}\) contenga memoria
  vera, cioè dipendenza da \(M(v')\) non catturata dagli istantanei.

## 6. Osservazione da non sovrainterpretare

Il rapporto decresce con \(\ell\) (0.244, 0.176, 0.138, 0.096) in modo
compatibile con \(\sim1/\ell\). Sarebbe coerente con un'adiabaticità migliore nel
limite eikonale, ma è un rapporto **dipendente dalla finestra** su un integrale
divergente: non ci costruirei sopra nulla finché non c'è la norma regolarizzata.

## 7. Stato del programma Vaidya

| domanda | esito |
|---|---|
| \(Q_\times\) ridondante a ordine \(\dot M\)? | **sì**, residuo \(=2\partial_r\partial_M Z\) |
| \(\mathcal M_{vr}\) di quale ordine? | \(O(\dot M^2)\), per costruzione e confermato |
| la sorgente è ortogonale al modo? | **no**, sovrapposizione genuina |
| quanto vale la proiezione? | **aperto**: serve la norma radiativa |
| c'è memoria? | **aperto**, ma con percorso identificato |

Per il manoscritto Vaidya resta il terzo caso della tesi «riproduce, non
estende», in forma analitica. Il termine di memoria è lavoro dichiarato, ora con
una definizione precisa, un percorso e un'ostruzione nominata.
