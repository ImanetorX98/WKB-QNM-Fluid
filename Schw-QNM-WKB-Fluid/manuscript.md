# Il potenziale di Madelung come chiusura della gerarchia WKB per le perturbazioni di Schwarzschild

**Bozza di manoscritto — destinazione: rivista di fisica matematica**

---

## Abstract

Per l'equazione master di tipo Schrödinger che governa le perturbazioni di un
buco nero di Schwarzschild mostriamo che l'intera serie WKB, riorganizzata
attorno al momento di fase pari, si chiude su un unico funzionale: il potenziale
quantistico di Bohm–Madelung valutato sull'ampiezza WKB completa. La chiusura

$$q = u^2 - \varepsilon^2\,\frac{(u^{-1/2})''}{u^{-1/2}}$$

è esatta, equivalente all'equazione differenziale di partenza, e non
un'approssimazione. Ne deriviamo tre conseguenze strutturali. Primo: la serie
non contiene ordini dispari, e il funzionale di Madelung compare per la prima
volta al **terzo termine della serie di fase** — dopo l'eikonale e dopo il
trasporto — generando da solo tutti gli ordini successivi per valutazione
autoconsistente. Secondo: nella riscalatura eikonale del potenziale di
Regge–Wheeler lo spin del campo **non entra nel termine di testa**; tutta la sua
traccia vive allo stesso ordine $\varepsilon^2$ del termine di Madelung, con cui
diventa quindi indistinguibile per scaling e separabile solo per forma
funzionale. Il caso gravitazionale assiale è quello in cui i due competono più
da vicino. Terzo: la gerarchia si rompe nella regione di coalescenza dei turning
point — precisamente dove la WKB di barriera costruisce le correzioni di
Iyer–Will — e questo esclude una corrispondenza canonica fra i coefficienti
locali di Madelung e quelle correzioni.

Documentiamo inoltre due limiti negativi, entrambi misurati e non solo
argomentati: per una frequenza quasi-normale complessa nessuna decomposizione a
densità reale positiva ammette una legge di continuità conservata — il flusso
obbedisce a $(\rho v)' = -\Gamma\rho$ con $\Gamma=\operatorname{Im}\omega^2\neq0$
— e la struttura non è universale, poiché per un campo di Dirac massless compare
un termine di connessione di spin all'ordine $\varepsilon^1$, che precede il
funzionale di Madelung e non ne fa parte.

Il risultato è un'identità esatta con un dominio di validità netto, non una
proposta di nuova fisica; non formuliamo previsioni osservative.

---

## 1. Introduzione

Il metodo WKB applicato al massimo della barriera è, dagli anni Ottanta, lo
strumento semianalitico standard per i modi quasi-normali (QNM) dei buchi neri.
Schutz e Will [6] adattarono la formula di trasmissione di barriera al problema
di scattering di Regge–Wheeler; Iyer e Will [7] costruirono la sistematica di
ordine superiore mediante matching uniforme, applicata a Schwarzschild nel
lavoro successivo [8]; Konoplya [9] la estese al sesto ordine — nella sua
formulazione, "al sesto ordine **oltre l'approssimazione eikonale**" —, e
Matyjasek e Opala [11] al tredicesimo con risommazione di Padé. La struttura di
quei risultati è però quella di una condizione di quantizzazione: correzioni
scalari $\Lambda_j$ costruite dal getto del potenziale nel punto di massimo.

In parallelo, la trasformazione di Madelung riscrive un'equazione di tipo
Schrödinger come sistema idrodinamico, in cui l'unico termine estraneo alla
meccanica dei fluidi classica è il potenziale quantistico di Bohm. È naturale
chiedersi se la gerarchia WKB dei QNM ammetta una lettura in quei termini: se
cioè le correzioni di ordine crescente siano, in qualche senso preciso, una
famiglia di potenziali di Madelung.

Questo lavoro risponde separando la domanda in tre affermazioni di forza
decrescente, e stabilendo quale sia vera.

- **H1, chiusura formale.** La serie WKB può essere espressa attraverso un
  funzionale di tipo Madelung e una ricorsione compatta. → **Vera**, e nel
  §3 la dimostriamo come identità esatta valida lontano da zeri e turning point.
- **H2, potenziali per ordine.** A ogni correzione WKB corrisponde un potenziale
  locale autonomo e canonico, in corrispondenza con i $\Lambda_j$. → **Falsa in
  generale**: §5 e §7.
- **H3, fluido fisico.** La gerarchia descrive correzioni dispersive di un fluido
  reale effettivo. → **Vincolata al punto di non essere sostenibile** nella forma
  ingenua: §8.

La distinzione serve a impedire che un'identità algebrica utile venga promossa a
dichiarazione ontologica.

### 1.1 Contributo

1. Dimostrazione della chiusura esatta e della sua equivalenza con l'ODE
   (Teorema 1), con verifica simbolica indipendente.
2. Dimostrazione dell'assenza di ordini dispari e della ricorsione generata da un
   solo funzionale (Teorema 2, Corollario 3), con $u_2$ e $u_4$ espliciti.
3. Individuazione precisa dell'ordine WKB di prima comparsa del termine di
   Madelung (Proposizione 4).
4. Osservazione, per Schwarzschild, che lo spin del campo è un effetto di ordine
   $\varepsilon^2$ e quindi degenere in scaling con il termine di Madelung
   (Proposizione 5), con tabella numerica per $s=0,1,2$.
5. Misura quantitativa della rottura della gerarchia alla coalescenza dei turning
   point (§7).
6. Due risultati negativi netti: fluido aperto (§8) e non universalità rispetto
   allo spin del campo, via Dirac (§9).

### 1.2 Stato della letteratura

Una ricerca bibliografica mirata non ha identificato un precedente diretto che
reinterpreti la gerarchia WKB dei QNM gravitazionali come gerarchia di potenziali
di Madelung. Esistono tre filoni adiacenti ma distinti: (i) WKB di ordine elevato
e Padé per QNM; (ii) equazione di Hamilton–Jacobi quantistica ed exact WKB, dove
il momento quantistico e i suoi periodi sono oggetti centrali; (iii) buchi neri
acustici in condensati di Bose–Einstein, dove il termine di pressione quantistica
è fisicamente microscopico e non un artificio di riscrittura. L'assenza
riscontrata è un esito di ricerca bibliografica, non una dimostrazione di
inesistenza, e va letta come tale.

---

## 2. Impostazione, convenzioni e i tre conteggi di "ordine"

Adottiamo $G=c=1$, $x=r/M$, coordinata tortoise adimensionale
$x_*=r_*/M$ con

$$x_* = x + 2\log\!\left(\frac{x}{2}-1\right),\qquad \frac{d}{dx_*}=f\frac{d}{dx},\qquad f=1-\frac{2}{x},$$

e dipendenza temporale $e^{-i\omega t}$, per cui un modo smorzato ha
$\operatorname{Im}\omega<0$. Le frequenze sono adimensionali, $\Omega=M\omega$.

L'equazione master è

$$\frac{d^2\psi}{dx_*^2}+\left[\Omega^2-V_s(x)\right]\psi=0,\qquad
V_s(x)=f\left[\frac{\ell(\ell+1)}{x^2}+\frac{2(1-s^2)}{x^3}\right],$$

con $s=0$ scalare, $s=1$ elettromagnetico, $s=2$ gravitazionale assiale
(Regge–Wheeler). In forma canonica scriviamo

$$\varepsilon^2\psi''+q(x_*)\psi=0,\qquad q=\Omega^2-V,$$

dove $\varepsilon$ è un parametro formale di conteggio e l'apice indica $d/dx_*$.

### 2.1 Avvertenza sui conteggi

Nella letteratura convivono almeno tre nozioni di "ordine", e la loro confusione
è la sorgente principale di affermazioni scorrette in questo ambito. Le fissiamo
qui una volta per tutte.

| Conteggio | Significato | Indicizzazione |
|---|---|---|
| **(O1) Ordine della serie di fase** | Numero del termine in $S=S_0+S_1+S_2+\dots$ | primo = eikonale $S_0$; secondo = trasporto $S_1$; terzo = $S_2$ |
| **(O2) Potenza di $\varepsilon$** | Grado nel parametro formale | $u=u_0+\varepsilon^2u_2+\varepsilon^4u_4+\dots$ |
| **(O3) Ordine WKB di barriera** | Troncamento di Iyer–Will nella condizione di quantizzazione, contato *oltre l'eikonale* | "WKB3" = $\Lambda_2,\Lambda_3$ inclusi |

Le tre non sono in corrispondenza banale. In particolare, l'affermazione centrale
di questo lavoro — *il potenziale di Madelung si riconosce dal terzo ordine in
poi* — è formulata in **(O1)**: il funzionale compare per la prima volta nel
terzo termine della serie di fase, che in **(O2)** è l'ordine $\varepsilon^2$.
In **(O3)** non ha invece alcun analogo diretto, ed è esattamente questo il
contenuto negativo della Proposizione 6.

---

## 3. La chiusura esatta

Sia $q$ regolare e non nulla su un dominio semplicemente connesso $D$, con una
scelta di ramo di $\sqrt q$ e senza zeri di $\psi$.

> **Teorema 1 (chiusura di Madelung).** Sia $u$ una funzione non nulla su $D$ e
> si ponga
> $$\psi = u^{-1/2}\exp\left(\frac{i}{\varepsilon}\int^{x} u\right).$$
> Allora $\psi$ risolve $\varepsilon^2\psi''+q\psi=0$ **se e solo se**
> $$q = u^2 - \varepsilon^2\,\frac{(u^{-1/2})''}{u^{-1/2}}. \tag{3.1}$$

*Dimostrazione.* Dalla definizione,
$\psi'/\psi = \tfrac{i}{\varepsilon}u - \tfrac{u'}{2u}$, e quindi

$$\frac{\psi''}{\psi} = \left(\frac{\psi'}{\psi}\right)' + \left(\frac{\psi'}{\psi}\right)^2
= -\frac{u^2}{\varepsilon^2} - \frac{u''}{2u} + \frac{3u'^2}{4u^2},$$

dove i due termini in $\tfrac{i}{\varepsilon}u'$ si cancellano identicamente.
Posto $A=u^{-1/2}$ si ha $A''/A = \tfrac{3u'^2}{4u^2}-\tfrac{u''}{2u}$, da cui
$\varepsilon^2\psi''/\psi = -u^2+\varepsilon^2 A''/A$, e la richiesta
$\varepsilon^2\psi''/\psi+q=0$ è la (3.1). $\square$

Il secondo termine della (3.1) è, a meno del fattore $-\tfrac12$ della
normalizzazione di Schrödinger, esattamente il potenziale quantistico di
Bohm–Madelung

$$Q_{\rm M} = -\frac{\varepsilon^2}{2}\frac{A''}{A}$$

valutato non sull'ampiezza di ordine zero, ma sull'**ampiezza WKB completa**
$A=u^{-1/2}$. Questa è la differenza sostanziale rispetto all'uso consueto del
potenziale di Bohm come termine correttivo: qui è il funzionale che chiude la
serie, non un contributo aggiunto.

La (3.1) è un'identità, non un'approssimazione: non abbiamo troncato nulla. Il
contenuto approssimativo entra solo quando si risolve la (3.1) per serie.

---

## 4. Assenza di ordini dispari e ricorsione

> **Teorema 2 (parità della serie).** Cercando una soluzione della (3.1) nella
> forma $u=\sqrt q+\varepsilon u_1+O(\varepsilon^2)$, il coefficiente di
> $\varepsilon^1$ impone $u_1\equiv0$. La serie è quindi in potenze di
> $\varepsilon^2$:
> $$u=u_0+\varepsilon^2u_2+\varepsilon^4u_4+\dots,\qquad u_0=\sqrt q .$$

Il contenuto dell'ordine dispari non è perduto: è stato assorbito nel prefattore
$u^{-1/2}$, cioè nell'equazione di trasporto dell'ampiezza. Il momento $u$ è per
costruzione la parte pari della derivata logaritmica.

Risolvendo la (3.1) ordine per ordine si ottiene:

$$u_2=\frac{5q'^2-4qq''}{32\,q^{5/2}}, \tag{4.1}$$

$$u_4=\frac{64q^3q^{(4)}-448q^2q'q'''-304q^2q''^2+1768q\,q'^2q''-1105\,q'^4}{2048\,q^{11/2}} . \tag{4.2}$$

La (4.1) coincide con la classica correzione WKB di secondo ordine
$5q'^2/(32q^{5/2})-q''/(8q^{3/2})$, il che fornisce un controllo esterno della
normalizzazione adottata.

> **Corollario 3 (un solo funzionale).** I coefficienti $u_{2j}$ per $j\geq2$ non
> richiedono ingredienti nuovi: si ottengono reinserendo il troncamento
> $u_0+\varepsilon^2u_2+\dots+\varepsilon^{2j-2}u_{2j-2}$ nello stesso funzionale
> $A''/A$ della (3.1) e raccogliendo. In particolare
> $$u_4=\frac{1}{2u_0}\left[\left.\frac{(u^{-1/2})''}{u^{-1/2}}\right|_{u_0+\varepsilon^2u_2}\right]_{\varepsilon^2}-\frac{u_2^2}{2u_0}.$$

Questa è la formulazione precisa — e delimitata — dell'intuizione di partenza.
Non c'è una *famiglia* di potenziali di Madelung: c'è **un** potenziale di
Madelung, la cui valutazione autoconsistente genera l'intera gerarchia. La
differenza non è terminologica. Una famiglia di potenziali indipendenti
suggerirebbe una successione di gradi di libertà fisici; un solo funzionale
iterato non suggerisce nulla del genere.

> **Proposizione 4 (ordine di prima comparsa).** Nel conteggio (O1) il funzionale
> di Madelung compare per la prima volta al **terzo termine** della serie di
> fase: $S_0$ è eikonale e non lo contiene, $S_1$ è il trasporto ed è assorbito
> nell'ampiezza, e la prima occorrenza è in $S_2$, cioè in $u_2$, di ordine
> $\varepsilon^2$ nel conteggio (O2). Da lì in poi il funzionale controlla ogni
> ordine successivo, per il Corollario 3.

Le verifiche simboliche di tutti gli enunciati di questa sezione sono in
`verification/madelung_recursion.py` (cinque controlli indipendenti, tutti
superati).

---

## 5. Schwarzschild: lo spin è un effetto del terzo ordine

Passiamo al caso concreto. Introduciamo il parametro eikonale di Langer

$$L=\ell+\tfrac12,\qquad \varepsilon=\frac1L,\qquad \Omega = L\,\hat\Omega,$$

e usiamo $\ell(\ell+1)=L^2-\tfrac14$. Il potenziale di Regge–Wheeler si riscrive
allora esattamente come

$$\boxed{\;\frac{V_s}{L^2}=h^2+\varepsilon^2 v_2\;},\qquad
h^2=\frac{f}{x^2},\qquad
v_2=f\left[\frac{2(1-s^2)}{x^3}-\frac{1}{4x^2}\right]. \tag{5.1}$$

L'indipendenza dallo spin del termine di testa non è un fatto nuovo, ed è
importante attribuirlo correttamente. Cardoso, Miranda, Berti, Witek e Zanchin
hanno mostrato che nel limite eikonale il potenziale di una vasta classe di
perturbazioni massless è "universale" e il suo estremo coincide con la geodetica
nulla circolare, da cui la corrispondenza
$\omega_{\rm QNM}=\Omega_c\ell-i(n+\tfrac12)\lvert\lambda\rvert$ fra frequenza
e esponente di Lyapunov dell'orbita instabile. Nella nostra notazione il loro
$Q_0\simeq\omega^2-f\ell^2/r^2$ è precisamente $L^2h^2$ a meno della
sostituzione di Langer. Ciò che aggiungiamo è quanto segue.

> **Proposizione 5 (degenerazione di ordine spin–Madelung).** La (5.1) è una
> riscrittura **esatta**, non asintotica, del potenziale di Regge–Wheeler: con
> la sostituzione di Langer $L=\ell+\tfrac12$ la serie in $\varepsilon$ termina
> al secondo ordine e non contiene termini di ordine $\varepsilon^1$. Tutta la
> dipendenza da $s$ è quindi confinata in $v_2$, cioè all'ordine
> $\varepsilon^2$ — lo stesso ordine a cui compare il termine di Madelung. Per
> scaling in $\varepsilon$ i due contributi sono **degeneri**: si distinguono
> solo per forma funzionale, non per ordine.

Il contenuto nuovo non è dunque l'universalità del termine di testa, che è
nota, ma il fatto che il primo termine **non** universale cada esattamente
all'ordine del funzionale di Madelung, rendendo i due indistinguibili con il
solo criterio di scaling. È questa coincidenza di ordini a rendere il caso
bosonico strutturalmente diverso da quello fermionico del §9.

Il caso elettromagnetico è degenere in modo ancora più stretto: per $s=1$ si ha
$2(1-s^2)=0$ e $v_2=-f/(4x^2)$ si riduce al **puro termine di Langer**, senza
alcun contributo di curvatura. Il caso gravitazionale assiale ha invece il
coefficiente di modulo massimo, $2(1-s^2)=-6$, ed è quindi quello in cui il
contributo di spin compete più da vicino con quello di Madelung. In questo senso
preciso il caso $s=2$ è il più interessante dei tre.

### 5.1 Chiusura numerica e misura degli ordini

Integrando il modo con condizione entrante all'orizzonte e decomponendo
$\psi=Ae^{iS/\varepsilon}$, $P=S'$, la chiusura assume la forma

$$\operatorname{Re}\hat\Omega^2=P^2+h^2+\varepsilon^2v_2+Q_{\rm M},\qquad
Q_{\rm M}=-\varepsilon^2\frac{A''}{A}. \tag{5.2}$$

La (5.2) è un'identità algebrica: è soddisfatta a precisione macchina per
costruzione e **non costituisce una verifica**. Il controllo indipendente è
ricostruire $A''/A$ per differenze finite dall'ampiezza integrata, il che
richiede una griglia equispaziata in $x_*$; adottiamo perciò $x_*$ come variabile
di integrazione, trasportando $x$ tramite $dx/dx_*=f$. Il residuo così misurato
resta fra $10^{-8}$ e $10^{-6}$ su tutto l'intervallo di $\ell$ esaminato.

Misura delle pendenze in $\varepsilon$, con $n=0$, finestra $20<x<50$ scelta
lontano dai turning point, $\ell\in\{2,4,8,16,32\}$:

| $s$ | $v_2$ | pendenza di $\varepsilon^2v_2$ | pendenza di $Q_{\rm M}$ | $\varepsilon^2v_2/Q_{\rm M}$ a $\ell=2$ |
|---|---|---|---|---|
| 0 | $f[+2/x^3-1/(4x^2)]$ | 2.0000 | 1.907 | 0.035 |
| 1 | $f[\;\;0\;\;-1/(4x^2)]$ | 2.0000 | 1.894 | 0.061 |
| 2 | $f[-6/x^3-1/(4x^2)]$ | 2.0000 | 1.854 | **0.150** |

Entrambe le pendenze sono 2, come richiesto dalla Proposizione 5, e il rapporto
fra i due contributi cresce monotonamente con lo spin, culminando nel caso
gravitazionale. Il difetto residuo della pendenza di $Q_{\rm M}$ rispetto al
valore esatto 2 non è fisico ma numerico: la condizione entrante
$\psi\sim e^{-i\Omega x_*}$ è esatta solo per $V\to0$, e a raggio iniziale finito
il ramo riflesso contamina $Q_{\rm M}$ con un termine che non scala in
$\varepsilon$. Spostando il bordo verso l'orizzonte la pendenza sale
monotonamente verso 2:

| $x_{\min}$ | 2.02 | 2.001 | 2.0001 | 2.00005 |
|---|---|---|---|---|
| pendenza di $Q_{\rm M}$ | 1.557 | 1.749 | 1.839 | 1.865 |

Riproducibile con `verification/scalar_eikonal_scaling.py`.

---

## 6. Assenza di corrispondenza con le correzioni di Iyer–Will

> **Proposizione 6.** Non esiste corrispondenza canonica
> $Q_{2j}(x)\leftrightarrow\Lambda_j$.

Le ragioni sono tre, di natura diversa.

**(i) Tipo di oggetto.** I $\Lambda_j$ sono scalari: numeri costruiti dal getto
del potenziale nel singolo punto di massimo e dall'indice di overtone
$\alpha=n+\tfrac12$. I $Q_{2j}$ sono funzioni della coordinata. Nessuna mappa
naturale porta le une negli altri senza una scelta arbitraria di funzionale di
valutazione.

**(ii) Non unicità.** I singoli $Q_{2j}$ dipendono dalla coordinata scelta, dalla
variabile master, dalla normalizzazione e dalla convenzione con cui si introduce
$\varepsilon$. L'exact WKB moderno insiste su questo punto: l'introduzione del
parametro formale non è unica, e qualunque lettura "fisica" dei coefficienti deve
dichiarare per intero quelle scelte. Il nostro §5 ne è un esempio: la scelta di
Langer $L=\ell+\tfrac12$ produce la (5.1) senza termine $\varepsilon^1$, mentre
la scelta $L=\ell$ ne produrrebbe uno spurio.

**(iii) Rottura nel punto sbagliato.** Le correzioni $\Lambda_j$ sono costruite
proprio nella regione di coalescenza dei turning point, dove — come mostriamo
quantitativamente nel §7 — la serie locale in potenze inverse di $q$ è singolare
e va sostituita dalla forma normale parabolico-cilindrica. La gerarchia di
Madelung e la gerarchia di Iyer–Will vivono in regioni di validità disgiunte.

---

## 7. Dove la gerarchia si rompe, misurato

Il punto (iii) è verificabile. Al massimo della barriera un QNM ha
$q\to0$, i due turning point coalescono e $P\to0$: la (3.1), che divide per $u$,
degenera. Misuriamo l'effetto confrontando le pendenze in $\varepsilon$ del
termine di Madelung in due finestre, una lontana dai turning point e una centrata
sul massimo. Per il campo di Dirac del §9, dove il confronto è più netto per la
presenza di due termini di ordine diverso:

| Regione | $\min\lvert P\rvert$ | pendenza del termine $\varepsilon^1$ | pendenza di $Q_{\rm M}$ |
|---|---|---|---|
| lontano dai turning point, $20<x<50$ | $\approx0.187$ | **1.0000** | 1.84 |
| al massimo di barriera | $2\times10^{-5}$ – $2\times10^{-4}$ | 1.14 | **1.40** |

Nella prima riga la gerarchia è quella prevista, con il termine analitico di
ordine $\varepsilon$ riprodotto a precisione macchina. Nella seconda, dove
$\lvert P\rvert$ crolla di tre-quattro ordini di grandezza, entrambe le pendenze si
degradano e $Q_{\rm M}$ cessa di seguire $\varepsilon^2$. La rottura non è un
artefatto numerico: è il comportamento atteso di una serie asintotica valutata
nel proprio punto singolare, e la sua posizione coincide con il punto di
costruzione dei $\Lambda_j$.

---

## 8. Il fluido di un QNM è aperto

Il passaggio dalla chiusura formale a un'interpretazione idrodinamica richiede
una decomposizione a densità reale e positiva. Per una frequenza complessa questo
è ostacolato in modo strutturale.

Poniamo $\psi=Ae^{iS}$ con $A,S$ reali e $\Omega^2=E+i\Gamma$. Separando parte
reale e immaginaria dell'equazione master si ottiene il sistema

$$A''-A\,S'^2+(E-V)A=0, \tag{8.1}$$
$$2A'S'+AS''+\Gamma A=0 \iff (\rho v)'=-\Gamma\rho,\qquad \rho=A^2,\;v=S'. \tag{8.2}$$

La (8.1) è l'equazione di Hamilton–Jacobi con il termine di Madelung; la (8.2)
è una legge di continuità **con sorgente**. Per un modo smorzato,
$\operatorname{Re}\Omega>0$ e $\operatorname{Im}\Omega<0$, si ha
$\Gamma=2\operatorname{Re}\Omega\operatorname{Im}\Omega<0$ e la sorgente
$-\Gamma\rho$ è ovunque positiva.

Verifica numerica ($s=2$, $n=0$), con residuo relativo alla scala
$\lvert\Gamma\rvert\rho$:

| $\ell$ | $\Gamma$ | residuo della legge con sorgente | violazione della legge conservata |
|---|---|---|---|
| 2 | $-0.0666$ | $8.5\times10^{-7}$ | 1.000 |
| 4 | $-0.1524$ | $9.4\times10^{-7}$ | 1.000 |
| 8 | $-0.3073$ | $1.2\times10^{-6}$ | 1.000 |
| 16 | $-0.6073$ | $2.2\times10^{-6}$ | 1.000 |

La legge con sorgente è soddisfatta a precisione numerica; la conservazione del
flusso è violata di ordine uno rispetto alla scala naturale del problema, in modo
uniforme in $\ell$. Restano due sole opzioni, entrambe con un costo. Mantenere
densità e velocità reali richiede di accettare un fluido aperto, con creazione
distribuita di materia. Mantenere la forma WKB esatta della (3.1) preserva la
chiusura ma rende $u$ e $A$ complessi, e con essi si perde il significato
ordinario di densità e velocità. Non è disponibile una terza via che conservi
entrambe le proprietà. Riproducibile con `verification/open_continuity.py`.

---

## 9. Non universalità: il contrappunto fermionico

La struttura dei §§3–5 è specifica del settore bosonico. Per un campo di Dirac
massless su Schwarzschild il sistema radiale è esattamente un sistema di Dirac
unidimensionale con massa spaziale $W=\kappa\sqrt f/x$, e le componenti scalari
obbediscono a equazioni di Schrödinger con potenziali partner di Darboux
$V_\pm=W^2\pm W'$. Inserendo il parametro semiclassico **prima** di disaccoppiare
il sistema — con $K=\lvert\kappa\rvert$, $\varepsilon=1/K$, $\omega=K\hat\Omega$,
$h=\sqrt f/x$, $\tau=\sigma\operatorname{sgn}\kappa$ — si ottiene

$$\frac{V_\tau}{K^2}=h^2+\tau\varepsilon h', \tag{9.1}$$

da confrontare direttamente con la (5.1). Il termine di testa è **lo stesso**
$h^2$ del caso bosonico, come deve essere: è il potenziale delle geodetiche nulle
e non conosce lo spin del campo. Ma la prima correzione è ora di ordine
$\varepsilon^1$, non $\varepsilon^2$. La chiusura diventa

$$\operatorname{Re}\hat\Omega^2=P^2+h^2+\underbrace{\tau\varepsilon h'}_{O(\varepsilon)}+\underbrace{Q_{\rm M}}_{O(\varepsilon^2)} . \tag{9.2}$$

Il termine di ordine $\varepsilon$ è la connessione di spin, coincide con metà
della separazione fra i partner, $V_+-V_-=2Kh'$, e **precede** il funzionale di
Madelung senza farne parte. La misura riportata nel §7 conferma la pendenza 1 a
precisione macchina.

La conclusione è netta: per Dirac la famiglia efficace esiste localmente, ma non
è una successione di soli potenziali di Madelung. La chiusura del Teorema 1
sopravvive — è un enunciato su equazioni scalari di secondo ordine, e le
componenti disaccoppiate lo sono — ma la sua interpretazione come *unico*
contenuto subprincipale della gerarchia è falsa non appena il campo ha spin
semintero. Ne segue un criterio di riconoscimento operativo, che è forse il
risultato più utile del lavoro:

> **Criterio.** Si misuri la pendenza in $\varepsilon$ del primo termine
> subprincipale della chiusura. Pendenza 2: il regime è puramente di Madelung, e
> spin e curvatura entrano degeneri con esso (settore bosonico, §5). Pendenza 1:
> esiste un contributo geometrico che precede il Madelung e va trattato a parte
> (settore fermionico, §9).

---

## 10. Discussione

Il verdetto complessivo è che l'ipotesi iniziale è **confermata in senso
formale-locale e complessificato**, e **non confermata in senso forte, canonico e
fisico**. Il candidato più robusto per una formulazione invariante non è una
collezione di potenziali locali, bensì il momento quantistico visto come 1-forma,
con i suoi periodi e i simboli di Voros — oggetti che sopravvivono alla
continuazione analitica e alle linee di Stokes, cosa che i singoli $Q_{2j}$ non
fanno.

**Sulla verificabilità sperimentale.** Non formuliamo previsioni osservative, e
riteniamo che non ve ne siano di accessibili. La riscrittura di Madelung è una
trasformazione esatta di variabili: non modifica lo spettro dei QNM, e dunque non
modifica nessuna forma d'onda. Ogni contenuto è di struttura matematica, non di
fenomenologia. L'unico contesto in cui un termine di questa forma è fisicamente
microscopico — e quindi in linea di principio misurabile — è quello dei buchi neri
acustici in condensati di Bose–Einstein, dove la pressione quantistica è un
termine reale della dinamica del condensato e non un artificio di riscrittura.
Anche lì, però, l'analogia riguarda il termine, non la gerarchia: la
corrispondenza fra i due contesti non è stabilita da questo lavoro e non ne
rivendichiamo alcuna.

**Direzioni aperte.** (a) Sostituire la serie locale con la forma normale
parabolico-cilindrica nella regione di coalescenza, e verificare se la chiusura
di Madelung ammetta lì un analogo uniforme. (b) Estendere il §9 a Kerr, dove la
ricorsione proiettiva del sistema di Chandrasekhar è già stata verificata
simbolicamente fino a $O(\varepsilon^2)$, e stabilire se la connessione di spin
mantenga l'ordine $\varepsilon^1$ in presenza di rotazione. (c) Chiarire se il
caso $s=1$, in cui $v_2$ si riduce al puro termine di Langer, ammetta una
caratterizzazione invariante — è l'unico dei tre spin bosonici in cui il termine
subprincipale non contiene curvatura.

---

## Appendice A. Riproducibilità

Interprete: `python3.13` (numpy, scipy, sympy, matplotlib). Tutti gli script sono
deterministici e girano in meno di un minuto ciascuno.

| Enunciato | Script | Esito |
|---|---|---|
| Teorema 1, chiusura esatta | `verification/madelung_recursion.py` (a) | verificato |
| Teorema 2, assenza ordini dispari | `verification/madelung_recursion.py` (b) | verificato |
| Eq. (4.1)–(4.2), ricorsione | `verification/madelung_recursion.py` (c,d) | verificato |
| Corollario 3, singolo funzionale | `verification/madelung_recursion.py` (e) | verificato |
| Proposizione 5 e tabella §5.1 | `verification/scalar_eikonal_scaling.py` | pendenze 2.0000 / ≈1.87 |
| Legge con sorgente, §8 | `verification/open_continuity.py` | residuo $\sim10^{-6}$ |
| Eq. (9.1)–(9.2), pendenze §7 | `../core/dirac_madelung_profile.py --scaling` | pendenza spin 1.0000 |
| QNM di riferimento | `../core/test_wkb.py`, `../core/test_dirac_madelung.py` | 14 test superati |

I QNM di riferimento sono confrontati con Iyer–Will ($s=2$, $\ell=2$, $n=0$:
$0.3732-0.0892i$) e con Cho ($\kappa=1,2$).

## Appendice B. Nota sulla scelta di Langer

L'uso di $L=\ell+\tfrac12$ anziché $L=\ell$ nella (5.1) non è cosmetico. Con
$L=\ell$ si otterrebbe $\ell(\ell+1)=L^2+L$, che genera un termine di ordine
$\varepsilon^1$ nel potenziale riscalato. Quel termine è però un artefatto della
parametrizzazione, non una connessione di spin: scompare con la sostituzione di
Langer, mentre il termine $\tau\varepsilon h'$ del caso di Dirac non scompare per
nessuna riparametrizzazione di $\kappa$, essendo fissato dalla separazione dei
partner di Darboux. La distinzione fra i due è essenziale perché il criterio del
§9 sia ben posto, e va dichiarata esplicitamente in ogni applicazione.

---

## Riferimenti

DOI risolti via Crossref e verificati per autori, volume e pagine; identificativi
arXiv via Semantic Scholar. Versione BibTeX in `references.bib`; gli otto PDF ad
accesso libero sono in `papers/`, i restanti sono elencati in
`codex-download-manifest.tsv`.

1. E. Madelung, *Quantentheorie in hydrodynamischer Form*, Z. Phys. **40**, 322–326 (1927). doi:10.1007/BF01400372
2. D. Bohm, *A suggested interpretation of the quantum theory in terms of "hidden" variables. I*, Phys. Rev. **85**, 166–179 (1952). doi:10.1103/PhysRev.85.166 — Parte II: doi:10.1103/PhysRev.85.180
3. T. Regge, J. A. Wheeler, *Stability of a Schwarzschild singularity*, Phys. Rev. **108**, 1063–1069 (1957). doi:10.1103/PhysRev.108.1063
4. F. J. Zerilli, *Effective potential for even-parity Regge–Wheeler gravitational perturbation equations*, Phys. Rev. Lett. **24**, 737–738 (1970). doi:10.1103/PhysRevLett.24.737
5. S. Chandrasekhar, *The Mathematical Theory of Black Holes*, Oxford University Press (1983). ISBN 978-0198512912
6. B. F. Schutz, C. M. Will, *Black hole normal modes: a semianalytic approach*, Astrophys. J. Lett. **291**, L33 (1985). doi:10.1086/184453
7. S. Iyer, C. M. Will, *Black-hole normal modes: a WKB approach. I. Foundations and application of a higher-order WKB analysis of potential-barrier scattering*, Phys. Rev. D **35**, 3621–3631 (1987). doi:10.1103/PhysRevD.35.3621
8. S. Iyer, *Black-hole normal modes: a WKB approach. II. Schwarzschild black holes*, Phys. Rev. D **35**, 3632–3636 (1987). doi:10.1103/PhysRevD.35.3632
9. R. A. Konoplya, *Quasinormal behavior of the $D$-dimensional Schwarzschild black hole and the higher order WKB approach*, Phys. Rev. D **68**, 024018 (2003). doi:10.1103/PhysRevD.68.024018, arXiv:gr-qc/0303052
10. R. A. Konoplya, *Quasinormal modes of the Schwarzschild black hole and higher order WKB approach*, J. Phys. Stud. **8**, 93–100 (2004). doi:10.30970/jps.08.93
11. J. Matyjasek, M. Opala, *Quasinormal modes of black holes: the improved semianalytic approach*, Phys. Rev. D **96**, 024011 (2017). doi:10.1103/PhysRevD.96.024011, arXiv:1704.00361
12. H. T. Cho, *Dirac quasinormal modes in Schwarzschild black hole spacetimes*, Phys. Rev. D **68**, 024003 (2003). doi:10.1103/PhysRevD.68.024003, arXiv:gr-qc/0303078
13. V. Ferrari, B. Mashhoon, *New approach to the quasinormal modes of a black hole*, Phys. Rev. D **30**, 295–304 (1984). doi:10.1103/PhysRevD.30.295
14. V. Cardoso, A. S. Miranda, E. Berti, H. Witek, V. T. Zanchin, *Geodesic stability, Lyapunov exponents, and quasinormal modes*, Phys. Rev. D **79**, 064016 (2009). doi:10.1103/PhysRevD.79.064016, arXiv:0812.1806
15. E. Berti, V. Cardoso, A. O. Starinets, *Quasinormal modes of black holes and black branes*, Class. Quantum Grav. **26**, 163001 (2009). doi:10.1088/0264-9381/26/16/163001, arXiv:0905.2975
16. R. A. Konoplya, A. Zhidenko, *Quasinormal modes of black holes: from astrophysics to string theory*, Rev. Mod. Phys. **83**, 793–836 (2011). doi:10.1103/RevModPhys.83.793, arXiv:1102.4014
17. M. V. Berry, K. E. Mount, *Semiclassical approximations in wave mechanics*, Rep. Prog. Phys. **35**, 315–397 (1972). doi:10.1088/0034-4885/35/1/306
18. A. Voros, *The return of the quartic oscillator. The complex WKB method*, Ann. Inst. H. Poincaré A **39**(3), 211–338 (1983). NUMDAM: AIHPA_1983__39_3_211_0 (nessun DOI assegnato)
19. E. Delabaere, H. Dillinger, F. Pham, *Exact semiclassical expansions for one-dimensional quantum oscillators*, J. Math. Phys. **38**, 6126–6184 (1997). doi:10.1063/1.532206
20. C. M. Bender, S. A. Orszag, *Advanced Mathematical Methods for Scientists and Engineers I*, Springer (1999). doi:10.1007/978-1-4757-3069-2
21. W. G. Unruh, *Experimental black-hole evaporation?*, Phys. Rev. Lett. **46**, 1351–1353 (1981). doi:10.1103/PhysRevLett.46.1351
22. C. Barceló, S. Liberati, M. Visser, *Analogue gravity*, Living Rev. Relativity **29** (2026). doi:10.1007/s41114-026-00064-9, arXiv:gr-qc/0505065 — sostituisce l'edizione **14** (2011), doi:10.12942/lrr-2011-3
