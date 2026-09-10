# Ridondanza e ordinamento: che cosa la decomposizione di Madelung aggiunge, e che cosa non aggiunge, ai modi quasi-normali dei buchi neri

**Bozza di manoscritto — destinazione: *Classical and Quantum Gravity*, Paper**

---

## Abstract

La riscrittura di Madelung di un'equazione master di tipo Schrödinger è una
trasformazione esatta di variabili, e come tale non può aggiungere informazione
allo spettro quasi-normale. La domanda utile è più precisa: **quanta
informazione contiene, e a quale ordine.** Rispondiamo separando ciò che la
riscrittura riproduce da ciò che ordina, e misurando entrambi.

**Il risultato negativo, e la sua causa strutturale.** Il funzionale d'ampiezza
che la decomposizione suggerisce come diagnostico d'errore è, a overtone
fissato, **proporzionale a $|\Lambda_3|$** — la correzione di Iyer–Will — con
costante indipendente da $\ell$ e dallo spin (dispersione 0.15% su dodici modi).
Non estende la correzione standard: la riproduce, a costo incomparabilmente
maggiore. La ragione non è accidentale ed è misurabile a monte. Con frequenza
complessa l'ampiezza reale contiene un termine
$-\varepsilon^{-1}\!\int\!\operatorname{Im}u$ che è $O(1)$, e il potenziale
quantistico che ne deriva **amplifica di circa due ordini di grandezza**
l'errore sulla frequenza: un errore relativo $10^{-3}$ su $\omega$ ne produce
uno del 14% su $Q_M$. Nessun funzionale costruito su $Q_M$ può quindi predire
la frequenza — calcolarlo a tre cifre ne richiede cinque di ciò che si vorrebbe
prevedere. Questa è un'**ostruzione di condizionamento**, non una difficoltà
numerica superabile con più risoluzione, ed è ciò che un lettore tentato da
questa via ha bisogno di sapere prima di percorrerla.

**Il risultato positivo: l'ordinamento.** Ciò che la riscrittura non rende
ridondante è a quale potenza di $\varepsilon$ compaia il primo termine
subprincipale, e se il termine d'ampiezza sia solo o preceduto da struttura
geometrica. Nella riscalatura eikonale del potenziale di Regge–Wheeler, con la
sostituzione di Langer $L=\ell+\tfrac12$, il potenziale si decompone
**esattamente** come $V_s/L^2=h^2+\varepsilon^2v_2$: la serie termina al secondo
ordine e non contiene termini di ordine $\varepsilon^1$. Poiché il funzionale di
Madelung compare anch'esso a $\varepsilon^2$, spin e potenziale quantistico sono
**degeneri in ordine** e separabili solo per forma funzionale. Due casi rompono
la degenerazione, in settori opposti: su Kerr l'autovalore sferoidale porta un
termine $\varepsilon^1$ che nessuna sostituzione di Langer rimuove — struttura
geometrica subprincipale **nel settore bosonico** — e per un campo di Dirac
massless su Schwarzschild compare, allo stesso ordine, la connessione di spin
dei partner di Darboux. Ne segue un criterio operativo: la pendenza in
$\varepsilon$ del primo termine subprincipale vale 2 nel regime puramente di
Madelung, 1 quando un contributo geometrico lo precede. Le pendenze misurate
sono 1.951 e 2.0000 nei casi statici senza rotazione, 0.963–1.014 su Kerr,
1.0000 per Dirac. Il criterio **non discrimina lo spin del campo** ma la
presenza di struttura geometrica: leggerlo come firma fermionica, come facevamo
in una versione precedente di questo lavoro, è un errore.

**Terza geometria.** Su Vaidya entrante il residuo dell'ansatz adiabatico
all'ordine $\dot M$ è **esattamente** $2\,\partial_r\partial_M Z$, senza altri
termini: il termine misto non è informazione dinamica nuova, ma la misura del
fallimento dell'adiabaticità, calcolabile dalla famiglia congelata. La
condizione di solvibilità che ne segue è regolarizzata a entrambi i bordi e il
suo numeratore non si annulla per alcun modo fondamentale esaminato. Al bordo
interno l'esponente vale $4M\operatorname{Im}\omega\simeq-0.39>-1$ per $s=0,1,2$:
per il modo fondamentale la singolarità è **integrabile** e la regolarizzazione
non è necessaria.

**Precisazione di metodo.** Il macchinario analitico è l'equazione di Riccati
associata all'equazione master, espansa nel parametro eikonale: da
$y=\varepsilon\psi'/\psi$ segue $\varepsilon y'+y^2+q=0$, la cui parte pari è il
momento e la cui parte dispari dà l'ampiezza. Chiusura, parità e ricorsione sono
l'eq. (1.2) di Delabaere–Dillinger–Pham e i coefficienti WKB classici: su quel
piano non rivendichiamo nulla. Ne segue anche che un'espansione della Riccati
non può contenere più della WKB che essa è — il che rende **attesa** la
ridondanza del §6, e sposta l'interesse sul quantificarla. Il lavoro riguarda
ciò che la geometria immette in $q$, e a quale ordine: quantità che la Riccati
prende come dato. Non formuliamo previsioni osservative; la trasformazione non
sposta lo spettro.

---

# Parte I — Il quadro esatto, e perché non basta

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
meccanica dei fluidi classica è il potenziale quantistico di Bohm. Questa
trasformazione ha una caratterizzazione geometrica esatta: Khesin, Misiołek e
Modin [31, 32] mostrano che è un **simplettomorfismo** fra $T^*\mathrm{Dens}(M)$
e lo spazio proiettivo delle funzioni d'onda non nulle, e un'**isometria** fra
la metrica di Sasaki–Fisher-Rao e quella di Fubini–Study.

Quel risultato però non si applica al problema quasi-normale, e il modo preciso
in cui non si applica orienta tutto il presente lavoro. Le sue ipotesi sono
densità **normalizzate** ($\int_M\varrho=1$) su varietà **compatta**, con
evoluzione hamiltoniana e norma conservata. Un QNM viola le prime due — non è
normalizzabile, e il dominio radiale non è compatto — e viola la terza in modo
misurabile: con $\omega$ complessa il flusso obbedisce a
$(\rho v)'=-\Gamma\rho$, $\Gamma=\operatorname{Im}\omega^2\neq0$ (§13).
Sopravvive invece l'ipotesi di non annullamento, $\mathbb C\setminus\{0\}$ nel
codominio di [31], che è esattamente la condizione $\psi\neq0$ del nostro
Teorema 1.

Le "limitazioni" che riportiamo nei §12 e §13 non sono quindi difetti del nostro
trattamento: sono le **ostruzioni identificate** che separano il problema
quasi-normale dal caso in cui la geometria della trasformazione è nota. Il
lavoro misura le conseguenze di quelle ostruzioni, ordine per ordine in
$\varepsilon$ — una domanda su cui il quadro geometrico, che riguarda la
trasformazione e non un regime asintotico, non si pronuncia. È naturale
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
  generale**: §7, §8 e §12.
- **H3, fluido fisico.** La gerarchia descrive correzioni dispersive di un fluido
  reale effettivo. → **Vincolata al punto di non essere sostenibile** nella forma
  ingenua: §13.

La distinzione serve a impedire che un'identità algebrica utile venga promossa a
dichiarazione ontologica.

### 1.1 Contributo, e che cosa non lo è

Dichiariamo prima ciò che **non** rivendichiamo, perché in questo ambito la
letteratura è densa e il rischio di riscoperta è alto.

- **L'intero macchinario formale è l'equazione di Riccati in $\varepsilon$.**
  Da $y=\varepsilon\psi'/\psi$ segue $\varepsilon y'+y^2+q=0$; la parte pari di
  $y$ è il momento, la dispari è una derivata totale e dà l'ampiezza. La chiusura
  esatta del §3 è l'equazione (1.2) di Delabaere, Dillinger e Pham [19], con
  ampiezza $P^{-1/2}$ e $P$ "parte pari in $\hbar$ della soluzione della
  Riccati"; il nostro Teorema 1 è quella formula, il Teorema 2 quella
  caratterizzazione, e $u_2$, $u_4$ sono i coefficienti WKB classici. Nulla di
  ciò è rivendicato.

  Ne discende una conseguenza che preferiamo enunciare noi: **un'espansione della
  Riccati non può contenere informazione che la WKB non contenga già**. Il
  risultato negativo del §6 è quindi atteso; il suo interesse sta nel
  quantificarlo, non nel constatarlo.

  Quello che la Riccati **non** fornisce è il proprio input. I §8, §9 e §10 sono
  enunciati sulla struttura di $q$ nella scalatura eikonale — da dove viene la
  sua dipendenza da $\varepsilon$ e a quale ordine — e sono quelli il contenuto
  del lavoro.
- I coefficienti $u_2$ e $u_4$ sono i coefficienti WKB classici.
- L'universalità del termine eikonale di testa, indipendente dallo spin, è
  stabilita da Cardoso *et al.* [14] e non è nostra.
- La decomposizione ampiezza–fase dei QNM, l'annullamento di un Wronskiano e
  l'integrazione di una Riccati sono strumenti standard.
- L'analisi WKB esatta dei QNM, con geometria di Stokes e periodi quantistici,
  è un programma attivo e consolidato che non tocchiamo.

Il contributo effettivo è triplice, e l'organizzazione del lavoro segue quella
divisione: la Parte II stabilisce che cosa la riscrittura non aggiunge, la
Parte III che cosa ordina.

1. **Ostruzione di condizionamento (§5).** Con frequenza complessa il potenziale
   quantistico amplifica di circa due ordini di grandezza l'errore sulla
   frequenza. Ne segue che nessun funzionale costruito su $Q_M$ può predire
   $\omega$: calcolarlo a tre cifre ne richiede cinque di ciò che si vorrebbe
   prevedere. È il risultato che delimita l'intero programma, e la ragione
   strutturale del punto 2.
2. **Ridondanza del diagnostico di Madelung (§6).** Il funzionale d'ampiezza
   naturale è proporzionale a $|\Lambda_3|$ a overtone fissato, con costante
   indipendente da $\ell$ e dallo spin: la lettura di Madelung dell'ampiezza WKB
   non estende la correzione di Iyer–Will, la riproduce.
3. **Criterio di ordinamento (§8, §9, §10), esteso a Vaidya (§11).** La
   decomposizione esatta (8.1) con sostituzione di Langer, l'osservazione che
   spin e Madelung sono degeneri in ordine nel settore bosonico statico, e i
   **due** contro-esempi di ordine $\varepsilon^1$ — la rotazione in Kerr,
   bosonica, e la connessione di spin di Dirac, fermionica — da cui il criterio
   nella forma generale del §10.1. Su Vaidya la stessa tesi assume la forma
   (11.2): il termine misto di ordine $\dot M$ è una derivata della famiglia
   congelata, quindi non porta informazione dinamica nuova.

Che il punto 1 preceda il punto 3 è deliberato. Un lettore che consideri questa
via ha bisogno di sapere per prima cosa che essa non produce frequenze, e solo
dopo che cosa produce.

Contorno di supporto, non rivendicato come nuovo: l'assenza di corrispondenza
locale con i $\Lambda_j$ (§7), la misura della rottura alla coalescenza dei
turning point (§12) e la legge di continuità con sorgente (§13), tutte attese ma
qui quantificate.

### 1.2 Stato della letteratura

Il territorio adiacente è occupato, e va dichiarato con precisione.

**Idrodinamica di Dirac con QNM.** Meza-Domínguez e Matos presentano una
formulazione chirale-idrodinamica covariante dell'equazione di Dirac in
spaziotempo curvo, con potenziale quantistico esplicito in un'equazione di
Hamilton–Jacobi radiale, specializzazione a Schwarzschild e calcolo dei QNM
fermionici. È un precedente diretto per il tema del §10. Il loro approccio è però
*esatto* — funzioni di Heun, conservazione del flusso chirale — e non sviluppa
alcuna gerarchia in $\varepsilon$: l'ordinamento, che è il nostro oggetto, resta
distinto.

**WKB esatta per i QNM.** Aminov, Grassi e Hatsuda hanno collegato i QNM alla
teoria di Seiberg–Witten con condizioni di quantizzazione esatte; Miyachi,
Namba, Omiya e Oshita hanno costruito l'analisi WKB esatta dell'equazione di
Regge–Wheeler, chiarendo la struttura delle curve di Stokes e le spirali
logaritmiche emergenti dall'orizzonte. È un programma attivo, con contributi
frequenti. Non lo tocchiamo: il nostro §12 osserva numericamente una rottura che
quella letteratura descrive rigorosamente, e la citiamo come lo strumento
appropriato piuttosto che riderivarla.

**Diagnostica d'errore per risonanze.** Esistono stimatori pratici dell'errore
WKB, in particolare $\Delta_k=|\omega_{k+1}-\omega_{k-1}|/2$ di Konoplya e
collaboratori, e una letteratura di stimatori a residuo pesato per autovalori
complessi di problemi aperti. Il §6 confronta esplicitamente con questi, e il
suo esito è negativo per il nostro candidato.

**Analoghi acustici.** Nei buchi neri acustici in condensati di Bose–Einstein il
termine di pressione quantistica è fisicamente microscopico e non un artificio
di riscrittura. È l'unico contesto in cui un termine di questa forma sarebbe
misurabile; non stabiliamo alcuna corrispondenza con esso.

Non abbiamo identificato, nelle fonti consultate, la formulazione
dell'ordinamento del §10.1. Questo non dimostra priorità.

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

Vale la pena registrare che cosa sia $Q_M$ in termini invarianti, perché il §7
rimprovera ai singoli $Q_{2j}$ proprio di non esserlo. Nella formulazione
geometrica di [31] il termine di Bohm compare nel potenziale dell'equazione di
Newton sullo spazio delle densità come $4I(\varrho)$, dove $I$ è il funzionale
di **informazione di Fisher**. Il funzionale che chiude la serie WKB è dunque,
variazionalmente, l'informazione di Fisher dell'ampiezza: un oggetto definito
indipendentemente dalla coordinata, dalla variabile master e dalla convenzione
con cui si introduce $\varepsilon$.

Questo permette di enunciare in forma invariante il criterio del §10.1: la
domanda non è a quale ordine compaia "il potenziale di Madelung", ma **a quale
ordine in $\varepsilon$ l'informazione di Fisher dell'ampiezza diventi degenere
con il contributo geometrico** — lo spin nel caso statico, la rotazione in Kerr.

La (3.1) è un'identità, non un'approssimazione: non abbiamo troncato nulla. Il
contenuto approssimativo entra solo quando si risolve la (3.1) per serie.

---

## 4. Assenza di ordini dispari e ricorsione

> **Teorema 2 (parità della serie).** Sia $q$ **indipendente da
> $\varepsilon$**. Cercando una soluzione della (3.1) nella forma
> $u=\sqrt q+\varepsilon u_1+O(\varepsilon^2)$, il coefficiente di
> $\varepsilon^1$ impone $u_1\equiv0$. La serie è quindi in potenze di
> $\varepsilon^2$:
> $$u=u_0+\varepsilon^2u_2+\varepsilon^4u_4+\dots,\qquad u_0=\sqrt q .$$

L'ipotesi di $\varepsilon$-indipendenza non è pleonastica, ed è bene isolarla
subito perché il §9 la viola.

> **Osservazione 2bis (dove la parità si rompe).** Se il potenziale riscalato ha
> una propria espansione $q=q_0+\varepsilon q_1+\varepsilon^2q_2+\dots$, allora
> l'ordine $\varepsilon^1$ della (3.1) dà
> $$q_1=2u_0u_1\quad\Longrightarrow\quad u_1=\frac{q_1}{2\sqrt{q_0}}\neq0 ,$$
> e il momento acquista un termine dispari. È esattamente il caso di Kerr, dove
> $q_1=-\Delta A_1/H^2$ è generato dall'autovalore sferoidale (§9.1).

La rottura della parità non si propaga però al funzionale di Madelung. Nella
chiusura $P^2+Q_M=\operatorname{Re}q$ i termini di ordine $\varepsilon$ si
cancellano identicamente, $q_1-2u_0u_1=0$, cosicché $Q_M$ resta di ordine
$\varepsilon^2$ anche in presenza di rotazione. Questa cancellazione è il motivo
per cui il criterio del §10.1 è ben posto: il termine di ordine $\varepsilon$
vive nel *potenziale*, non nel funzionale di ampiezza, ed è per questo che si
può misurarlo senza integrare l'equazione.

In pratica, quando $q$ dipende da $\varepsilon$ conviene non separarlo in ordini
e applicare la serie WKB direttamente al $q$ complesso completo; il §5
raccoglie la forma dell'ampiezza che ne segue.

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


---

# Parte II — Che cosa la riscrittura **non** aggiunge

## 5. Il condizionamento del funzionale di ampiezza

La chiusura del §3 è un'identità, e le identità non creano informazione. Resta
però la possibilità che la riscrittura sia comunque **utile**: che $Q_M$, pur
non contenendo nulla di nuovo, sia una via di calcolo o una diagnostica
conveniente. Questa sezione la esclude, e lo fa a monte del §6 — non
constatando che un particolare funzionale fallisce, ma misurando la ragione per
cui ogni funzionale di quella famiglia deve fallire.

Il punto è che $Q_M$ non è una quantità misurabile in modo indipendente. Si
costruisce dalla soluzione, e costruire la soluzione richiede $\omega$. Non c'è
nulla da invertire: $Q_M$ è una funzione della frequenza che si sa valutare solo
conoscendo già la frequenza. Quantifichiamo **quanto** lo sia.

### 5.1 L'ampiezza reale con frequenza complessa

Con $\omega$ complessa l'ampiezza reale della decomposizione di Madelung non è
$|u|^{-1/2}$. Da $\psi=u^{-1/2}\exp(i\!\int\! u/\varepsilon)$ segue

$$\ln A=-\tfrac12\ln|u|-\frac1\varepsilon\int\operatorname{Im}u,
\qquad
Q_M=-\varepsilon^2\left[(\ln A)''+\big((\ln A)'\big)^2\right],$$

e il secondo termine di $\ln A$ è $O(1)$, non trascurabile, perché per un QNM
$\operatorname{Im}\hat\Omega=O(\varepsilon)$. Ometterlo sbaglia $Q_M$ di un
fattore $\sim50$. Va inoltre scelto il ramo di $\sqrt q$ coerente con la
condizione al contorno su ciascun lato della barriera: usarne uno solo produce
un errore del 10% che **non svanisce** al crescere di $L$, e che quindi imita un
difetto di teoria.

Validazione su Pöschl–Teller $n=0$, dove $Q_M$ è noto in forma chiusa. Errore
relativo mediano nella regione oscillatoria:

| $L$ | $u=\sqrt q$ | $u=\sqrt q+\varepsilon^2u_2$ |
|---|---|---|
| 8 | $3.9\times10^{-3}$ | $6.9\times10^{-4}$ |
| 50 | $1.0\times10^{-4}$ | $5.1\times10^{-6}$ |
| 300 | $2.8\times10^{-6}$ | $1.7\times10^{-7}$ |

### 5.2 L'amplificazione, misurata

La stessa costruzione mostra il limite strutturale. Perturbando $\omega$ di una
frazione relativa e misurando lo scarto indotto su $Q_M$ (Schwarzschild,
$\ell=70$):

| $\delta\omega/\omega$ | $10^{-5}$ | $10^{-4}$ | $10^{-3}$ | $3\times10^{-3}$ |
|---|---|---|---|---|
| errore su $Q_M$ | $1.0\times10^{-3}$ | $1.2\times10^{-2}$ | $1.4\times10^{-1}$ | $4.8\times10^{-1}$ |

L'amplificazione è di circa due ordini di grandezza e stabile. Ne segue che
**nessun funzionale costruito su $Q_M$ può predire la frequenza**: calcolarlo a
tre cifre richiede cinque cifre di ciò che si vorrebbe prevedere. Non è una
difficoltà numerica superabile con più risoluzione, è il condizionamento del
problema, e fornisce la ragione strutturale del risultato negativo del §6.

Per contrasto, l'**ordinamento** è robusto: le pendenze in $\varepsilon$ del §8,
§9 e §10 sopravvivono a errori sulla frequenza di ordine $10^{-3}$. È in questa
asimmetria — ordinamento robusto, valori mal condizionati — che si colloca il
potere predittivo effettivo del formalismo.

---

## 6. Il diagnostico di Madelung riproduce $\Lambda_3$, non lo estende

La decomposizione suggerisce un uso pratico: se $Q_M$ misura la deviazione
dell'ampiezza dal comportamento WKB, un suo funzionale dovrebbe stimare l'errore
della WKB. Costruiamo il candidato naturale, con finestra gaussiana $w$ centrata
sul massimo di barriera e larghezza $\sigma=\sqrt{V_0/(-V_0'')}$:

$$\mathcal E_M=\frac{\int w\,|Q_M|\,dx_*}{\int w\left(|\varepsilon\Omega|^2+\varepsilon^2|V|+P^2\right)dx_*}. \tag{6.1}$$

Il risultato è negativo, e in modo netto.

**(i) Degenerazione nella famiglia esattamente risolubile.** Per Pöschl–Teller
$V=L^2\operatorname{sech}^2y$ con $n=0$ valgono le identità esatte
$|\Omega/L|^2=1$ e $P^2=(1-\varepsilon^2/4)\tanh^2y$, da cui il denominatore
della (6.1) vale $2-\tfrac{\varepsilon^2}{4}\tanh^2y$ e

$$\mathcal E_M=\frac{\varepsilon^2(2-c_w)}{8-\varepsilon^2c_w},\qquad
c_w=\frac{\int w\tanh^2y\,dy}{\int w\,dy}.$$

Ogni ingrediente è una forma fissa moltiplicata per una funzione di
$\varepsilon$: l'indicatore non ha margine per trasportare informazione.

**(ii) Il segnale apparente ad overtone alto è il nodo.** Su Schwarzschild
$\mathcal E_M/\varepsilon^2$ sembra dipendere da $\ell$ a $n=1$ (spread 85%),
ma l'integrale è dominato all'89% dal picco di $|Q_M|$ al nodo, e la dominanza
cresce con $\ell$. Sostituendo l'integrale con la mediana pesata, robusta al
polo, lo spread crolla all'1.9%: il segnale era l'artefatto.

**(iii) A $\varepsilon$ fissato l'indicatore ordina, ma non meglio dei
controlli gratuiti.** Variando lo spin a $\ell$ fissato — che cambia la forma
della barriera lasciando $\varepsilon$ identico — su 18 coppie:

| predittore | concordanza | errore medio sul rapporto |
|---|---|---|
| $\mathcal E_M$ | 18/18 | 0.0179 |
| $\Delta_2=|\omega_3-\omega_1|/2$ | 18/18 | 0.0185 |
| $|\Lambda_3|$ | 18/18 | 0.0187 |

**(iv) La ragione: proporzionalità esatta.** A overtone fissato,

| $n$ | $\mathcal E_M/|\Lambda_3|$ su $s=0,2$, $\ell=3\ldots8$ | dispersione |
|---|---|---|
| 0 | 1.0426 | 0.15% |
| 1 | 3.95 | 1.6% |

L'indicatore è una riscalatura di $|\Lambda_3|$ con costante indipendente da
$\ell$ e dallo spin. Ma $\Lambda_3$ è un sottoprodotto algebrico della formula
WKB al terzo ordine, calcolato dal getto di $V$ nel massimo; $\mathcal E_M$
richiede la frequenza esatta, dati iniziali di Frobenius e l'integrazione
dell'ODE attraverso la barriera.

> **Conclusione.** La lettura di Madelung dell'ampiezza WKB non produce
> informazione diagnostica nuova. Riproduce la correzione di Iyer–Will a costo
> molto maggiore. Questo è coerente con il Teorema 1: la chiusura è
> un'identità, e un'identità non crea informazione. L'enunciato utile non è che
> il formalismo fallisce, ma che il suo contenuto informativo è **esattamente**
> quello della WKB da cui proviene — né più né meno — e che l'unico punto in cui
> dice qualcosa di non ridondante è l'ordinamento dei §8, §9 e §10.

## 7. Assenza di corrispondenza con le correzioni di Iyer–Will

> **Proposizione 6.** Non esiste corrispondenza canonica
> $Q_{2j}(x)\leftrightarrow\Lambda_j$.

Le ragioni sono tre, di natura diversa.

**(i) Tipo di oggetto.** I $\Lambda_j$ sono scalari: numeri costruiti dal getto
di $Q=\Omega^2-V$ nel singolo punto di massimo e dall'indice di overtone
$\alpha=n+\tfrac12$. La costruzione è esplicita in [7]: nella regione fra i
turning point si approssima $-Q$ con uno sviluppo di Taylor "fino alla derivata
sesta inclusa", e sono quelle sei derivate al picco — null'altro — a produrre
$\Lambda_2$ e $\Lambda_3$. I $Q_{2j}$ sono invece funzioni della coordinata.
Nessuna mappa naturale porta le une negli altri senza una scelta arbitraria di
funzionale di valutazione.

**(ii) Non unicità.** I singoli $Q_{2j}$ dipendono dalla coordinata scelta, dalla
variabile master, dalla normalizzazione e dalla convenzione con cui si introduce
$\varepsilon$. L'exact WKB moderno insiste su questo punto: l'introduzione del
parametro formale non è unica, e qualunque lettura "fisica" dei coefficienti deve
dichiarare per intero quelle scelte. Il nostro §8 ne è un esempio: la scelta di
Langer $L=\ell+\tfrac12$ produce la (8.1) senza termine $\varepsilon^1$, mentre
la scelta $L=\ell$ ne produrrebbe uno spurio.

**(iii) Rottura nel punto sbagliato.** Le correzioni $\Lambda_j$ sono costruite
proprio nella regione di coalescenza dei turning point, dove — come mostriamo
quantitativamente nel §12 — la serie locale in potenze inverse di $q$ è singolare.
Questo non è un'inferenza nostra: è la ragione dichiarata per cui il metodo
esiste. Iyer e Will [7] introducono la loro modifica del WKB proprio "per valori
della frequenza tali che i turning point $x_1$ e $x_2$ siano vicini fra loro,
presso il picco della barriera", dove "il matching standard non è più valido"; e
la soluzione interna che usano per il raccordo è, in assenza dei termini in
$\varepsilon$, una **funzione del cilindro parabolico** $D_\nu(t)$, cercata poi
nella forma $f(t)D_\nu[g(t)]$.

C'è di più, ed è il principio che governa l'intero metodo uniforme. Berry e
Mount [17], §4.1, lo enunciano così: *nel limite semiclassico sono equivalenti
tutti i problemi che hanno la stessa struttura di turning point*. Ne segue che
una costruzione per equazione di confronto, come quella di [7], codifica la
**struttura dei turning point**, non il potenziale locale: due potenziali diversi
con la stessa disposizione di zeri di $q$ danno gli stessi $\Lambda_j$ a meno
della mappa. I $Q_{2j}$ sono invece coefficienti di uno sviluppo locale, e
distinguono potenziali che il metodo uniforme identifica. Sono oggetti di natura
diversa non solo per tipo — scalari contro funzioni — ma per ciò che ciascuno è
costruito per registrare.

Le due gerarchie vivono dunque in regioni di validità **disgiunte per
costruzione**: la serie di Madelung richiede $q\neq0$ e degenera dove i due
turning point si fondono; la costruzione di Iyer–Will è definita esattamente
là, e sostituisce la serie locale con una forma normale. Chiedere una
corrispondenza $Q_{2j}\leftrightarrow\Lambda_j$ significa chiedere di
identificare oggetti definiti su domini che non si intersecano.

---


---

# Parte III — Che cosa la riscrittura **ordina**

## 8. Schwarzschild: lo spin è un effetto del terzo ordine

Passiamo al caso concreto. Introduciamo il parametro eikonale di Langer

$$L=\ell+\tfrac12,\qquad \varepsilon=\frac1L,\qquad \Omega = L\,\hat\Omega,$$

e usiamo $\ell(\ell+1)=L^2-\tfrac14$. Il potenziale di Regge–Wheeler si riscrive
allora esattamente come

$$\boxed{\;\frac{V_s}{L^2}=h^2+\varepsilon^2 v_2\;},\qquad
h^2=\frac{f}{x^2},\qquad
v_2=f\left[\frac{2(1-s^2)}{x^3}-\frac{1}{4x^2}\right]. \tag{8.1}$$

L'indipendenza dallo spin del termine di testa non è un fatto nuovo, ed è
importante attribuirlo correttamente. Cardoso, Miranda, Berti, Witek e Zanchin
hanno mostrato che nel limite eikonale il potenziale di una vasta classe di
perturbazioni massless è "universale" e il suo estremo coincide con la geodetica
nulla circolare, da cui la corrispondenza
$\omega_{\rm QNM}=\Omega_c\ell-i(n+\tfrac12)\lvert\lambda\rvert$ fra frequenza
e esponente di Lyapunov dell'orbita instabile. Nella nostra notazione il loro
$Q_0\simeq\omega^2-f\ell^2/r^2$ è precisamente $L^2h^2$ a meno della
sostituzione di Langer. Ciò che aggiungiamo è quanto segue.

> **Proposizione 5 (degenerazione di ordine spin–Madelung).** La (8.1) è una
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
bosonico statico strutturalmente diverso sia da Kerr (§9) sia dal caso
fermionico (§10).

Il caso elettromagnetico è degenere in modo ancora più stretto: per $s=1$ si ha
$2(1-s^2)=0$ e $v_2=-f/(4x^2)$ si riduce al **puro termine di Langer**, senza
alcun contributo di curvatura. Il caso gravitazionale assiale ha invece il
coefficiente di modulo massimo, $2(1-s^2)=-6$, ed è quindi quello in cui il
contributo di spin compete più da vicino con quello di Madelung. In questo senso
preciso il caso $s=2$ è il più interessante dei tre.

### 8.1 Chiusura numerica e misura degli ordini

Integrando il modo con condizione entrante all'orizzonte e decomponendo
$\psi=Ae^{iS/\varepsilon}$, $P=S'$, la chiusura assume la forma

$$\operatorname{Re}\hat\Omega^2=P^2+h^2+\varepsilon^2v_2+Q_{\rm M},\qquad
Q_{\rm M}=-\varepsilon^2\frac{A''}{A}. \tag{8.2}$$

La (8.2) è un'identità algebrica: è soddisfatta a precisione macchina per
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

## 9. Kerr: la rotazione produce un termine di ordine $\varepsilon$

Il §8 mostra che nel settore bosonico statico non esiste termine
$\varepsilon^1$: la sostituzione di Langer dà $\ell(\ell+1)=L^2-\tfrac14$, e il
$-\tfrac14$ è puro $\varepsilon^2$. Sarebbe naturale concluderne che la pendenza
2 caratterizzi il settore bosonico. **Non è così**, e Kerr lo mostra.

![Fig. 1](figures/fig1_ordinamento.pdf)

**Figura 1.** Primo termine subprincipale della decomposizione contro
$\varepsilon=1/L$, per i tre casi trattati. Ogni serie è riscalata in modo che
la retta di fit passi per 1 a $\varepsilon=0.1$: le ampiezze assolute hanno
significati fisici diversi e non sono confrontabili, le pendenze sì. Le rette
grigie hanno pendenza 1 e 2. Schwarzschild gravitazionale, dove il primo termine
è $\varepsilon^2v_2$, misura 2.00; Kerr scalare a $a=0.6$, con
$-\varepsilon\Delta A_1/H^2$, misura 1.00; Dirac con la connessione di spin
$\tau\varepsilon h'$ misura 1.00. I due casi a pendenza 1 stanno in settori
opposti — bosonico con rotazione, fermionico statico — ed è questo a mostrare
che l'oggetto discriminato è la struttura geometrica, non lo spin del campo.

### 9.1 Decomposizione esatta

Con $\Psi=\sqrt{r^2+a^2}R$, $H=r^2+a^2$, $\Delta=r^2-2r+a^2$,
$\lambda=A+a^2\omega^2-2am\omega$ e la scalatura $\omega=L\hat\Omega$,
$m=\mu L$, $\hat c=a\hat\Omega$, il potenziale radiale efficace si decompone
esattamente come

$$\frac{q}{L^2}=Q_0(r)
-\varepsilon\,\frac{\Delta A_1}{H^2}
-\varepsilon^2\left[\frac{\Delta A_2}{H^2}+\frac{h''}{h}\right],
\qquad
Q_0=\left(\hat\Omega-\frac{\mu a}{H}\right)^2-\frac{\Delta\bar A_0}{H^2},
\tag{9.1}$$

dove $A=A_0+A_1/L+A_2/L^2+\dots$ è l'espansione eikonale dell'autovalore
sferoidale e $\bar A_0=A_0+\hat c^2-2\mu\hat c$. La (9.1) va confrontata
direttamente con la (8.1) del caso statico e con la (10.1) del caso di Dirac:
**hanno la stessa forma**, e differiscono solo per quale potenza di
$\varepsilon$ porta il primo termine subprincipale.

Le convenzioni alternative dell'equazione angolare differiscono per termini
$a^2\omega^2=L^2\hat c^2$ e $2am\omega=2L^2\mu\hat c$, entrambi puramente
$O(L^2)$: entrano in $\bar A_0$, non in $A_1$. Il risultato è quindi
indipendente dalla convenzione.

### 9.2 Misura non circolare

Definire $A_1$ come coefficiente di $1/L$ e poi verificarne lo scaling sarebbe
tautologico. Si procede invece così: a ogni $\ell$ si costruisce il potenziale
con l'autovalore sferoidale **esatto** a quel $\ell$, lo si confronta con la sua
forma di testa $Q_0$, e si misura la pendenza in $\varepsilon$ del residuo
$D(r)=q/L^2-Q_0(r)$ su una finestra centrata sul picco.

| $a$ | $\mu$ | $\hat\Omega$ | $r_{\rm picco}$ | $A_0$ | **pendenza** |
|---|---|---|---|---|---|
| 0.0 | 0.5 | 0.192450 | 3.00000 | 1.00000 | **1.951** |
| 0.0 | 0.9 | 0.192450 | 3.00000 | 1.00000 | **1.951** |
| 0.3 | 0.5 | 0.205614 | 2.78835 | 0.99857 | **0.963** |
| 0.6 | 0.5 | 0.225120 | 2.47487 | 0.99315 | **0.996** |
| 0.9 | 0.5 | 0.261204 | 1.93474 | 0.97924 | **1.006** |
| 0.6 | 0.9 | 0.252767 | 2.23949 | 0.99778 | **1.014** |
| 0.9 | 0.9 | 0.330631 | 1.60599 | 0.99135 | **0.997** |

Il controllo statico dà 2, la rotazione dà 1, per ogni spin e ogni $\mu$
provati. Il difetto del controllo rispetto a 2.000 è troncamento a
$\ell\le160$, non struttura.

### 9.3 Il termine non è un artefatto di parametrizzazione

Nel caso sferico l'espansione dà $A_1=-1.4\times10^{-14}$ e
$A_2=-0.25000000$, cioè esattamente il termine di Langer: il metodo riproduce il
§8. Con rotazione $A_1\neq0$, e non è rimovibile. Ridefinendo $L\to L+\delta$ si
ha $A_1\to A_1-2\delta A_0$; una costante unica assorbirebbe $A_1$ solo se
$\delta=A_1/(2A_0)$ fosse comune ai modi, mentre i valori richiesti vanno da
$-0.010$ a $-0.128$, con dispersione del 162%. È lo stesso test dell'Appendice B,
che separa l'artefatto di Langer dal termine genuino.

Per un overtone fissato $\operatorname{Im}\hat c=O(\varepsilon)$; facendo
$\operatorname{Im}\hat c\to0$ la parte reale di $A_1$ resta stabile a $-0.0193$.
Il termine sopravvive nel limite fisico ed è reale all'ordine dominante.

### 9.4 Validazioni

Per $a=0$ si ritrova $\hat\Omega=1/(3\sqrt3)$ e $r_{\rm picco}=3$, senza
dipendenza da $\mu$. La condizione di radice doppia $Q_0(r_0)=Q_0'(r_0)=0$ è
soddisfatta a $10^{-16}$–$10^{-12}$. I raggi trovati stanno sopra l'orbita
fotonica equatoriale prograda $r_{\rm ph}=2[1+\cos(\tfrac23\arccos(-a))]$ e vi si
avvicinano al crescere di $\mu$ — 2.788 e 2.239 contro 2.630 e 2.189 per
$a=0.6$ — che è il comportamento corretto delle orbite fotoniche sferiche.

### 9.5 Attribuzione

Che $A_1\neq0$ per l'autovalore sferoidale è **noto**, ed è contenuto nella
letteratura eikonale su Kerr [23, 24]. Non è questa la
rivendicazione.

Va distinta anche dalla WKB di barriera su Kerr di Seidel e Iyer [33], che
applica il metodo di [7] ai modi **bassi** espandendo il potenziale in potenze
di $a\omega$. È un parametro piccolo diverso dal nostro $\varepsilon=1/L$: quel
lavoro non contiene autovalori sferoidali né limite eikonale, e le difficoltà
che riporta ad alto $a$ derivano proprio dall'espansione in $a\omega$, che qui
non si effettua — il problema angolare è risolto esattamente. Il contenuto di questa sezione è la collocazione di quel termine
nella stessa casella della connessione di spin di Dirac, attraverso una misura
operativa comune, e la conseguente correzione del criterio (§10.1).

### 9.6 Il funzionale di ampiezza su Kerr

L'ordinamento del §9.2 è misurato sul potenziale e non richiede di integrare
l'equazione. Il funzionale $Q_M$ è però costruibile anche qui, e serve a
verificare che nulla di patologico accada all'ampiezza in presenza di rotazione.

La frequenza va determinata **autoconsistentemente col potenziale a quel
$\ell$**, non presa dal limite eikonale: a $\omega$ fissata l'autovalore
sferoidale è un numero, quindi $q(r)$ è esplicita e le sue derivate rispetto a
$r_*$ si ottengono simbolicamente; la condizione di barriera di [7], scritta in
termini di $q$,

$$\frac{q_0}{\sqrt{2q_0''}}=\Lambda_2-i\left(n+\tfrac12\right)(1+\Lambda_3),
\tag{9.2}$$

si risolve per $(\omega,r_0)$ complessi simultaneamente. A $a=0$ il solutore
riproduce il WKB3 scalare a $3\times10^{-11}$; su Kerr converge con residui
$10^{-14}$–$10^{-12}$, e la frequenza si sposta rispetto a quella eikonale di
$1.5\times10^{-3}$ a $1.4\times10^{-2}$.

Quello spostamento non è trascurabile, perché $Q_M$ amplifica di due ordini di
grandezza l'errore sulla frequenza (§5). Confrontando la previsione
analitica del §5 con $Q_M$ estratto dall'ampiezza integrata, errore
relativo mediano nella finestra $20<r<50$:

| $a$ | $\ell$ | con $\omega$ eikonale | con $\omega$ autoconsistente |
|---|---|---|---|
| 0.0 | 30 | $1.6\times10^{-3}$ | $3.4\times10^{-4}$ |
| 0.3 | 30 | $3.6\times10^{-2}$ | $3.4\times10^{-4}$ |
| 0.6 | 30 | $7.8\times10^{-2}$ | $4.4\times10^{-4}$ |
| 0.6 | 70 | $1.2\times10^{-3}$ | $3.2\times10^{-5}$ |
| 0.9 | 70 | $2.3\times10^{-1}$ | $1.5\times10^{-3}$ |
| 0.9 | 100 | $1.9\times10^{-1}$ | $1.4\times10^{-3}$ |

Con la frequenza consistente gli errori scendono nella banda
$3\times10^{-5}$–$1.5\times10^{-3}$, comparabile al controllo statico, e il
comportamento erratico scompare. La previsione del §5, costruita e
validata su Pöschl–Teller, vale dunque anche in presenza di rotazione: **non
c'è nulla di patologico nell'ampiezza**, e la struttura di ordine del §9.2 è
l'unico contenuto.

*Limite dichiarato.* Il caso $a=0.9$, $\mu=0.9$ non è incluso: il bordo interno
richiesto avvicina l'orizzonte al punto in cui l'integratore esaurisce la
precisione di macchina. Il regime quasi estremale richiede un trattamento
dedicato che non affrontiamo.

---

## 10. Non universalità: il contrappunto fermionico

La struttura dei §§3–5 è specifica del settore bosonico. Per un campo di Dirac
massless su Schwarzschild il sistema radiale è esattamente un sistema di Dirac
unidimensionale con massa spaziale $W=\kappa\sqrt f/x$, e le componenti scalari
obbediscono a equazioni di Schrödinger con potenziali partner di Darboux
$V_\pm=W^2\pm W'$. Inserendo il parametro semiclassico **prima** di disaccoppiare
il sistema — con $K=\lvert\kappa\rvert$, $\varepsilon=1/K$, $\omega=K\hat\Omega$,
$h=\sqrt f/x$, $\tau=\sigma\operatorname{sgn}\kappa$ — si ottiene

$$\frac{V_\tau}{K^2}=h^2+\tau\varepsilon h', \tag{10.1}$$

da confrontare direttamente con la (8.1) e con la (9.1). Il termine di testa è **lo stesso**
$h^2$ del caso bosonico, come deve essere: è il potenziale delle geodetiche nulle
e non conosce lo spin del campo. Ma la prima correzione è ora di ordine
$\varepsilon^1$, non $\varepsilon^2$. La chiusura diventa

$$\operatorname{Re}\hat\Omega^2=P^2+h^2+\underbrace{\tau\varepsilon h'}_{O(\varepsilon)}+\underbrace{Q_{\rm M}}_{O(\varepsilon^2)} . \tag{10.2}$$

Il termine di ordine $\varepsilon$ è la connessione di spin, coincide con metà
della separazione fra i partner, $V_+-V_-=2Kh'$, e **precede** il funzionale di
Madelung senza farne parte. La misura riportata nel §12 conferma la pendenza 1 a
precisione macchina.

La conclusione è netta: per Dirac la famiglia efficace esiste localmente, ma non
è una successione di soli potenziali di Madelung. La chiusura del Teorema 1
sopravvive — è un enunciato su equazioni scalari di secondo ordine, e le
componenti disaccoppiate lo sono — ma la sua interpretazione come *unico*
contenuto subprincipale della gerarchia è falsa non appena esiste un termine
$O(\varepsilon)$.

### 10.1 Il criterio, nella forma generale

Mettendo insieme i tre casi — Schwarzschild bosonico (§8), Kerr scalare (§9) e
Dirac (§10) — si può ora formulare il criterio senza attribuirlo allo spin.

> **Criterio.** Si misuri la pendenza in $\varepsilon$ del primo termine
> subprincipale della decomposizione. **Pendenza 2**: il regime è puramente di
> Madelung, e spin e curvatura entrano degeneri con esso. **Pendenza 1**: esiste
> un contributo geometrico che precede il funzionale di Madelung e va trattato a
> parte.

Il criterio **non distingue spin intero da semintero**, e sarebbe un errore
leggerlo così. Ne conosciamo due sorgenti indipendenti di pendenza 1:

| sorgente | settore | termine $O(\varepsilon)$ | pendenza misurata |
|---|---|---|---|
| separazione sferoidale da rotazione | bosonico, Kerr | $-\Delta A_1/H^2$ | 0.963 – 1.014 |
| connessione di spin di Darboux | fermionico, statico | $\tau h'$ | 1.0000 |

Che le due coesistano in settori opposti — una bosonica con rotazione, una
fermionica senza — è ciò che mostra che l'oggetto discriminato è la **struttura
geometrica**, non lo spin del campo. Una formulazione precedente di questo lavoro
attribuiva la pendenza 1 al solo settore fermionico: il §9 la smentisce.

---


## 11. Vaidya: il termine misto misura il fallimento dell'adiabaticità

I §8–§10 riguardano l'ordinamento di $q$ in geometrie stazionarie. Una massa
che varia mette alla prova la stessa tesi in una direzione diversa: non «a quale
ordine in $\varepsilon$», ma «quanta informazione dinamica porta il termine
subprincipale». La risposta è la stessa, e in forma più netta: **nessuna**.

In questa sezione soltanto, $M$ è funzione della coordinata nulla entrante $v$ e
non è riassorbibile nelle unità; scriviamo quindi $r$ e $M$ esplicitamente.

### 11.1 Il residuo dell'ansatz adiabatico è esattamente una derivata mista

Nella metrica di Vaidya entrante la riduzione caratteristica dell'equazione di
perturbazione è

$$2\,\partial_v\partial_r\psi+\partial_r\!\left(f\,\partial_r\psi\right)-U_\ell\,\psi=0,
\qquad f=1-\frac{2M(v)}{r}. \tag{11.1}$$

L'ansatz adiabatico congela la geometria e lascia scorrere la fase,
$\psi=Z(r;M)\exp\!\big(-i\!\int\!\omega(M)\,dv\big)$, con $Z$ il modo
quasi-normale del problema statico di massa $M$. All'ordine $\dot M^0$ la (11.1)
si riduce all'equazione congelata $\mathcal L_M Z=0$. All'ordine $\dot M$ resta
un residuo, e il calcolo simbolico dà

$$\boxed{\;\mathcal R_{\dot M}=2\,\partial_r\partial_M Z\;} \tag{11.2}$$

**senza altri termini** (verifica in Appendice A). Non una combinazione di
derivate di cui $\partial_r\partial_M Z$ sia il pezzo principale: esattamente
quello.

L'enunciato ha un contenuto preciso. Il termine misto — il candidato naturale
per un «effetto di memoria» all'ordine più basso — è costruito interamente dalla
**famiglia congelata** $\{Z(\cdot\,;M)\}$. Non richiede di risolvere alcuna
dinamica non stazionaria: si ottiene derivando rispetto al parametro una
soluzione statica. È dunque una misura del **fallimento dell'adiabaticità**, non
un'informazione nuova sulla dinamica.

È la stessa struttura dei §5 e §6, in una terza geometria: un termine che
sembra aggiungere contenuto risulta calcolabile da ciò che si aveva già.

### 11.2 La forzatura è risonante, e la condizione di solvibilità è non banale

La correzione $Z_1$ obbedisce a $\mathcal L_M Z_1=-2\,\partial_r\partial_M Z$.
Ma alla frequenza quasi-normale $\mathcal L_M$ è **singolare**: $Z$ sta nel suo
nucleo. La forzatura è risonante, e la componente della sorgente parallela al
modo non produce $Z_1$ ma uno spostamento di frequenza. Solo la componente
ortogonale genera una correzione vera.

L'operatore congelato è autoaggiunto rispetto al prodotto **bilineare** — senza
coniugazione — $\langle u,v\rangle=\int\mu\,u\,v\,dr$ con
$\mu=e^{-2i\omega r_*}$, che è il prodotto corretto per problemi di risonanza
non autoaggiunti. L'integrale non converge, come tutti gli integrali
quasi-normali, e va regolarizzato con la norma generalizzata di Leung, Liu,
Suen, Tam e Young [35], eq. (2.16), i cui termini di superficie rendono
numeratore e denominatore **separatamente** indipendenti dai punti di raccordo.
La forma semplificata (2.19) di [35] non è applicabile: vale per potenziali
senza code, e Regge–Wheeler ha code di potenza.

Il denominatore regolarizzato risulta indipendente dal raccordo esterno allo
0.2% fra finestre che differiscono di tre ordini di grandezza nell'integrale
nudo.

### 11.3 Il numeratore, regolarizzato ai due bordi

Il numeratore richiede cura a entrambi gli estremi, e i due bordi si comportano
in modo asimmetrico.

**Bordo esterno.** L'integrando cresce come $e^{2i\omega r_*}$ e la
regolarizzazione sottrae la primitiva asintotica, costruita per algebra esatta
di serie in $1/r$.

**Bordo interno.** A $r=2M$ l'equazione ha un punto singolare **regolare** con
indici $0$ e $4i\omega M$; il ramo entrante è quello analitico. Con $x=r-2M$ il
peso si separa esattamente,

$$I=K\,x^{-4i\omega M}B(x),\qquad K=e^{-4i\omega M}(2M)^{4i\omega M},
\qquad B\ \text{analitica}, \tag{11.3}$$

dove l'analiticità di $B$ segue dalla relazione indiciale, che cancella il polo
di $Z''$. L'integrale si somma allora in forma chiusa,
$\int_0^w I\,dx=K\sum_k B_k\,w^{\,k+1-4i\omega M}/(k+1-4i\omega M)$.

Da qui un'osservazione che non abbiamo trovato dichiarata altrove. L'esponente
al bordo interno è $\operatorname{Re}(-4i\omega M)=4M\operatorname{Im}\omega$, e
per il modo **fondamentale** vale $\simeq-0.39>-1$ per $s=0,1,2$:

> la singolarità al bordo interno è **integrabile**, e il limite esiste già.
> Per il modo fondamentale la regolarizzazione al bordo interno non è
> necessaria.

Per $n\geq1$ l'esponente scende sotto $-1$ e la (11.3) è la continuazione
analitica che definisce il valore.

Il risultato è un numeratore finito e privo di tagli arbitrari:

| $\ell$ | $s$ | $N$ | disp. taglio esterno | disp. raccordo interno |
|---|---|---|---|---|
| 2 | 0 | $+20.6665418-40.3265371\,i$ | $1.1\times10^{-7}$ | $3.4\times10^{-8}$ |
| 3 | 0 | $+67.0552912-31.0152159\,i$ | $3.7\times10^{-7}$ | $1.4\times10^{-8}$ |
| 2 | 2 | $-14.9301378-21.8933101\,i$ | $1.5\times10^{-7}$ | $4.7\times10^{-9}$ |
| 3 | 2 | $+26.7652152-51.7881564\,i$ | $2.0\times10^{-7}$ | $7.3\times10^{-9}$ |
| 2 | 1 | $+10.4233793-39.4370806\,i$ | $2.8\times10^{-7}$ | $2.4\times10^{-8}$ |

Il numeratore **non si annulla** per nessun modo fondamentale esaminato. La
sovrapposizione fra sorgente risonante e modo è genuina, dunque $Z_1\neq0$: la
correzione non adiabatica esiste, e il candidato di memoria — che per la (11.2)
è $O(\dot M^2)$ — non è escluso.

**Limite dichiarato.** Per $n\geq1$ la primitiva esterna cresce come
$e^{2|\operatorname{Im}\omega|r_*}$ e a $L_+=80M$ vale $7.8\times10^{24}$
($n=1$) contro $N$ di ordine $10$: la cancellazione consuma quattordici delle
sedici cifre della doppia precisione. I valori riportati sopra sono quindi
limitati ai fondamentali, e l'estensione agli overtoni richiede aritmetica
multiprecisione all'estremo esterno.

### 11.4 Che cosa questa sezione stabilisce, e che cosa no

Stabilisce che il termine misto di ordine $\dot M$ è **esattamente** una
derivata mista della famiglia congelata, e che la proiezione risonante che ne
deriva è finita e non nulla. Non stabilisce l'esistenza di una memoria
osservabile: $N\neq0$ garantisce $Z_1\neq0$, non che $Z_1$ produca un effetto
misurabile in una forma d'onda. La distanza fra le due affermazioni è la stessa
che separa il §5 dal §6 — un funzionale ben definito non è ancora un osservabile.


---

# Parte IV — Ostruzioni identificate

## 12. Dove la gerarchia si rompe, misurato

Il punto (iii) è verificabile. Al massimo della barriera un QNM ha
$q\to0$, i due turning point coalescono e $P\to0$: la (3.1), che divide per $u$,
degenera. Misuriamo l'effetto confrontando le pendenze in $\varepsilon$ del
termine di Madelung in due finestre, una lontana dai turning point e una centrata
sul massimo. Per il campo di Dirac del §10, dove il confronto è più netto per la
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

Il meccanismo è classificato. In linguaggio di equazioni di confronto [17],
§4.2, la soluzione WKB non è valida presso un turning point classico, dove
$p^2(x)$ ha uno zero, "perché la mappa non è più biunivoca": il termine
trascurato nella costruzione diventa dominante e la condizione di validità è
violata. Al massimo di barriera i due zeri di $q$ si fondono in un **turning
point del secondo ordine**, il caso $\nu=2$ della classificazione di Langer, la
cui equazione di confronto è quella del cilindro parabolico — la stessa che [7]
usa per costruire i $\Lambda_j$. La degradazione misurata nella seconda riga
della tabella è dunque la firma numerica di un fenomeno noto e classificato: il
suo interesse non è la scoperta, ma il confronto quantitativo con la prima riga.

L'exact WKB dà a questo la formulazione più tagliente. In [19] la rottura non è
descritta come divergenza della serie, ma come **pinzatura del cammino di
normalizzazione**: gli sviluppi WKB ben normalizzati "sono singolari per quei
valori dell'energia per cui il cammino di normalizzazione è *pinched* dalla
confluenza di alcuni turning point". Localmente si può sempre scegliere una base
non pinzata; ciò che si perde è la regolarità **delle mappe di connessione** fra
regioni. La singolarità è dunque nel problema di connessione, non nel singolo
sviluppo.

La cura è una riscalatura: si pone $E=E_{\rm crit}+E_r\hbar$ presso un punto
critico quadratico del potenziale, e il problema di connessione riscalato dipende
regolarmente da $E_r$ su tutto il piano complesso. Vale la pena notare che la
condizione di quantizzazione di [7] ha esattamente quella forma —
$q_0/\sqrt{2q_0''}=i\varepsilon(n+\tfrac12)$, cioè $q_0=O(\varepsilon)$ — per
cui **i $\Lambda_j$ sono il contenuto troncato del problema di connessione
riscalato a un turning point doppio**. Questo chiude il §7 dall'altro capo: i
$Q_{2j}$ sono coefficienti dello sviluppo *non* riscalato, che è precisamente
quello che diventa singolare; i $\Lambda_j$ appartengono al problema riscalato,
che è regolare. Sono espansioni di oggetti diversi.

---

## 13. Il fluido di un QNM è aperto

Il passaggio dalla chiusura formale a un'interpretazione idrodinamica richiede
una decomposizione a densità reale e positiva. Per una frequenza complessa questo
è ostacolato in modo strutturale.

Poniamo $\psi=Ae^{iS}$ con $A,S$ reali e $\Omega^2=E+i\Gamma$. Separando parte
reale e immaginaria dell'equazione master si ottiene il sistema

$$A''-A\,S'^2+(E-V)A=0, \tag{13.1}$$
$$2A'S'+AS''+\Gamma A=0 \iff (\rho v)'=-\Gamma\rho,\qquad \rho=A^2,\;v=S'. \tag{13.2}$$

La (13.1) è l'equazione di Hamilton–Jacobi con il termine di Madelung; la (13.2)
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

In termini geometrici, la sorgente è precisamente ciò che impedisce di
collocare il problema su $T^*\mathrm{Dens}(M)$: lo spazio delle densità di [31]
è definito da $\varrho>0$ e $\int\varrho=1$, e la seconda condizione è
preservata dal flusso solo se la continuità è omogenea. Il termine
$-\Gamma\rho$ misura quindi l'ostruzione, e $\Gamma=\operatorname{Im}\omega^2$
la quantifica: si annulla esattamente per i modi normali, e mai per i
quasi-normali.

---

## 14. Discussione

### 14.1 Il bilancio, nei termini della domanda iniziale

La domanda era quanta informazione la riscrittura di Madelung contenga, e a
quale ordine. La risposta si separa nettamente lungo la divisione fra le Parti
II e III.

**Sui valori, nulla, e per una ragione strutturale.** Il potenziale quantistico
amplifica di $\sim10^2$ l'errore sulla frequenza (§5), quindi non è una presa
indipendente su $\omega$ ma una lente d'ingrandimento su una frequenza che si
deve già conoscere. La conseguenza si misura: il diagnostico d'ampiezza è
proporzionale a $|\Lambda_3|$ allo 0.15% (§6), e non esiste corrispondenza locale
fra i coefficienti della gerarchia e i $\Lambda_j$ (§7). Tre enunciati, un solo
contenuto: un'espansione della Riccati non contiene più della WKB che essa è.

**Sull'ordinamento, qualcosa, e in tre geometrie.** La potenza di $\varepsilon$
alla quale compare il primo termine subprincipale è informazione sulla
geometria, non sulla frequenza, ed è **robusta** proprio dove i valori sono mal
condizionati: le pendenze sopravvivono a errori su $\omega$ di ordine $10^{-3}$.
Su Schwarzschild vale 2 e spin e Madelung sono degeneri (§8); su Kerr la
rotazione la porta a 1 nel settore bosonico (§9); per Dirac la connessione di
spin la porta a 1 nel fermionico (§10). Su Vaidya la stessa tesi si presenta in
forma non asintotica: il termine misto è una derivata della famiglia congelata
(§11).

È in questa asimmetria — **ordinamento robusto, valori mal condizionati** — che
si colloca il potere predittivo effettivo del formalismo, e riteniamo sia questo
il dato che un lettore tentato da questa via deve avere prima di percorrerla.

### 14.2 Dove vive il contenuto invariante

Il verdetto complessivo è che l'ipotesi iniziale è **confermata in senso
formale-locale e complessificato**, e **non confermata in senso forte, canonico e
fisico**. Per il funzionale di chiusura una formulazione invariante esiste ed è
nota: è l'informazione di Fisher dell'ampiezza (§3). Per la **gerarchia** invece
il candidato robusto non è una collezione di potenziali locali. La formulazione
appropriata esiste ed è quella dell'exact WKB: il momento vive sul rivestimento
doppio $\dot{\mathbb C}_2$, superficie di Riemann di $p=(E-V)^{1/2}$, e gli
sviluppi vanno trattati "come oggetti impliciti, il cui interesse principale sta
nei piccoli esponenziali che generano attraverso il processo di risorgenza"
[19]. È un contenuto invariante di natura diversa da un coefficiente locale: non
sopravvive alla continuazione analitica *malgrado* le linee di Stokes, ma è
definito *attraverso* di esse. I singoli $Q_{2j}$, coefficienti di uno sviluppo
non riscalato e pinzato ai turning point doppi, non hanno quella proprietà.

Va aggiunto che il problema quasi-normale rientra nell'ambito di [19] senza
forzature: la ricerca di risonanze vi è definita come la ricerca dei valori
complessi dell'energia per cui la componente riflessa dell'onda si annulla, che
è la condizione QNM.

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

### 14.3 Direzioni aperte

**Direzioni aperte.** (a) Sostituire la serie locale con la forma normale
parabolico-cilindrica nella regione di coalescenza, e verificare se la chiusura
di Madelung ammetta lì un analogo uniforme. (b) Montare il funzionale $Q_M$ sul profilo radiale di Kerr del §9, dove la
ricorsione proiettiva del sistema di Chandrasekhar è già stata verificata
simbolicamente fino a $O(\varepsilon^2)$, e stabilire se la connessione di spin
mantenga l'ordine $\varepsilon^1$ in presenza di rotazione. (c) Chiarire se il
caso $s=1$, in cui $v_2$ si riduce al puro termine di Langer, ammetta una
caratterizzazione invariante — è l'unico dei tre spin bosonici in cui il termine
subprincipale non contiene curvatura. (d) Estendere il numeratore del §11.3 agli
overtoni: per $n\geq1$ la primitiva esterna cresce come
$e^{2|\operatorname{Im}\omega|r_*}$ e la cancellazione consuma quattordici delle
sedici cifre della doppia precisione, sicché servirebbe aritmetica
multiprecisione all'estremo esterno. È l'unico punto del lavoro in cui la
precisione di macchina, e non la struttura, sia il limite.

---


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
| Proposizione 5 e tabella §8.1 | `verification/scalar_eikonal_scaling.py` | pendenze 2.0000 / ≈1.87 |
| Legge con sorgente, §13 | `verification/open_continuity.py` | residuo $\sim10^{-6}$ |
| Eq. (10.1)–(10.2), pendenze §12 e §10 | `../core/dirac_madelung_profile.py --scaling` | pendenza spin 1.0000 |
| §6, ridondanza di $\mathcal E_M$ | `../../calculations/robust_indicator_test.py` | $\mathcal E_M/|\Lambda_3|$ = 1.0426 ± 0.15% |
| §9.1–9.3, termine $\varepsilon$ di Kerr | `../../calculations/kerr_eikonal_order_test.py` | $A_1\neq0$; Langer esatto per $a=0$ |
| §9.2, pendenze di Kerr | `../../calculations/kerr_radial_order_profile.py` | 1.951 statico, 0.963–1.014 con rotazione |
| Figura 1 | `figures/make_fig1.py` | pendenze 2.00 / 1.00 / 1.00 |
| §9.6, frequenza autoconsistente | `../../calculations/kerr_wkb3_selfconsistent.py` | WKB3 scalare a $3\times10^{-11}$ per $a=0$ |
| Frequenze esatte di riferimento | `../../core/leaver_qnm.py` | 10 cifre contro Leaver 1985 |
| §5, previsione e sensibilità | `../../calculations/madelung_wkb_prediction.py` | $1.7\times10^{-7}$ a $L=300$ |
| §11.1, residuo adiabatico di Vaidya | `../../calculations/vaidya_adiabatic_residual.py` | residuo $=2\,\partial_r\partial_M Z$, simbolico |
| §11.2, denominatore regolarizzato | `../../calculations/vaidya_solvability.py` | indipendente da $L_+$ allo 0.2% |
| §11.3, numeratore ai due bordi | `../../calculations/vaidya_numerator_factored.py` | $N=20.666542-40.326537\,i$; disp. $1.1\times10^{-7}$ |
| Suite completa | `core/` e `calculations/` | **56 test superati** |

I QNM di riferimento sono confrontati con Iyer–Will ($s=2$, $\ell=2$, $n=0$:
$0.3732-0.0892i$) e con Cho ($\kappa=1,2$).

Nota di convenzione: [7] costruisce i $\Lambda_j$ dalle derivate di
$Q=\Omega^2-V$, mentre la nostra implementazione usa quelle di $V$. Le due
coincidono perché $Q^{(k)}=-V^{(k)}$ per $k\geq1$ e le formule dei $\Lambda_j$
contengono solo rapporti pari in cui il segno si cancella, come
$V_4/V_2$ e $(V_3/V_2)^2$. Il test contro i valori di [8] lo conferma.


---

## Appendice B. Nota sulla scelta di Langer


L'uso di $L=\ell+\tfrac12$ anziché $L=\ell$ nella (8.1) non è cosmetico, ed è
anzi **unico**. Berry e Mount [17], §5.1, discutono la famiglia di sostituzioni
$\ell(\ell+1)\to\ell(\ell+1)+\alpha$ associate a una relazione lineare generale
$L=(\ell+\beta)\hbar$ fra numero quantico e momento angolare classico. Richiedere
che il termine centrifugo diventi esattamente $L^2$ impone
$\ell^2+\ell+\alpha=(\ell+\beta)^2$, cioè $\beta=\tfrac12$ e $\alpha=\tfrac14$:
non c'è altra scelta. È la modifica di Langer, che [17] indica come "la più
naturale per problemi che coinvolgono il limite semiclassico".

La sua origine non è di comodo: Langer (1937) mostrò che la WKB ingenua non è
applicabile presso $r=0$, e che la sostituzione $r=e^x$ ne ripristina la
validità. Il $-\tfrac14$ che compare nel nostro $v_2$ è quel termine.

Con $L=\ell$ si otterrebbe invece $\ell(\ell+1)=L^2+L$, che genera un termine di
ordine $\varepsilon^1$ nel potenziale riscalato. Quel termine è però un artefatto della
parametrizzazione, non una connessione di spin: scompare con la sostituzione di
Langer, mentre il termine $\tau\varepsilon h'$ del caso di Dirac non scompare per
nessuna riparametrizzazione di $\kappa$, essendo fissato dalla separazione dei
partner di Darboux. La distinzione fra i due è essenziale perché il criterio del
§10.1 sia ben posto, e va dichiarata esplicitamente in ogni applicazione.
Lo stesso test, applicato all'autovalore sferoidale, e' quello che nel §9.3
mostra che il termine di rotazione non e' rimovibile.

---


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
23. S. R. Dolan, *Quasinormal mode spectrum of a Kerr black hole in the eikonal limit*, Phys. Rev. D **82**, 104003 (2010). doi:10.1103/PhysRevD.82.104003 — Publisher's Note: doi:10.1103/PhysRevD.82.109906
24. H. Yang, D. A. Nichols, F. Zhang, A. Zimmerman, Z. Zhang, Y. Chen, *Quasinormal-mode spectrum of Kerr black holes and its geometric interpretation*, Phys. Rev. D **86**, 104006 (2012). doi:10.1103/PhysRevD.86.104006
25. R. A. Konoplya, A. Zhidenko, A. F. Zinhailo, *Higher order WKB formula for quasinormal modes and grey-body factors: recipes for quick and accurate calculations*, Class. Quantum Grav. **36**, 155002 (2019). doi:10.1088/1361-6382/ab2e25, arXiv:1904.10333
26. G. Aminov, A. Grassi, Y. Hatsuda, *Black hole quasinormal modes and Seiberg–Witten theory*, Ann. Henri Poincaré **23**, 1951–1977 (2022). doi:10.1007/s00023-021-01137-x, arXiv:2006.06111
27. T. Miyachi, R. Namba, H. Omiya, N. Oshita, *Path to an exact WKB analysis of black hole quasinormal modes*, Phys. Rev. D **111** (2025). doi:10.1103/1gmr-9f1g, arXiv:2503.17245
28. T. Miyachi, R. Namba, H. Omiya, N. Oshita, *Spectral instability of parametrized black hole quasinormal modes in the high-overtone limit via the exact WKB analysis*, Phys. Rev. D **113** (2026). doi:10.1103/gsht-7s9l, arXiv:2512.18631
29. J. Meza-Domínguez, T. Matos, *A covariant chiral-hydrodynamic formulation of the Dirac equation in curved spacetime*, arXiv:2605.28887 (2026). Precedente diretto per il ramo Dirac; v1, nessuna pubblicazione su rivista rilevata.
30. E. Abdalla, C. B. M. H. Chirenti, A. Saa, *Quasinormal modes for the Vaidya metric*, Phys. Rev. D **74**, 084029 (2006). doi:10.1103/PhysRevD.74.084029, arXiv:gr-qc/0609036
31. B. Khesin, G. Misiołek, K. Modin, *Geometric hydrodynamics via Madelung transform*, Proc. Natl. Acad. Sci. **115**, 6165–6170 (2018). doi:10.1073/pnas.1719346115, arXiv:1711.00321
32. B. Khesin, G. Misiołek, K. Modin, *Geometry of the Madelung transform*, Arch. Ration. Mech. Anal. **234**, 549–573 (2019). doi:10.1007/s00205-019-01397-2 — dimostrazioni del precedente
33. E. Seidel, S. Iyer, *Black-hole normal modes: a WKB approach. IV. Kerr black holes*, Phys. Rev. D **41**, 374–382 (1990). doi:10.1103/PhysRevD.41.374
34. L. Capuano, L. Santoni, E. Barausse, *Perturbations of the Vaidya metric in the frequency domain: quasinormal modes and tidal response*, Phys. Rev. D **110**, 084081 (2024). doi:10.1103/PhysRevD.110.084081, arXiv:2407.06009
35. P. T. Leung, Y. T. Liu, W.-M. Suen, C. Y. Tam, K. Young, *Perturbative approach to the quasinormal modes of dissipative systems*, J. Phys. A **31**, 3271–3290 (1998). arXiv:physics/9712037 — **DOI da verificare.** Nota di attribuzione: gli autori sono **cinque**; il lavoro è correntemente citato come “Leung–Liu–Young”.
