# Interpretazione idrodinamica dell’espansione WKB per i modi quasi-normali dei buchi neri tramite il formalismo di Madelung

## Rapporto di ricerca — stato della letteratura, formulazione matematica e programma di verifica

**Data di chiusura della ricerca:** 23 agosto 2026  
**Ambito:** equazioni master del tipo Schrödinger, WKB locale ed esatta, quantizzazione WKB al massimo della barriera, modi quasi-normali (QNM), trasformazione di Madelung.  
**Criterio metodologico:** l’esistenza di un “potenziale di Madelung a ogni ordine” è trattata come ipotesi, non come premessa.

---

## Sintesi esecutiva

Il risultato principale è deliberatamente più preciso — e più limitato — della congettura iniziale.

1. **Sì, esiste una gerarchia locale a tutti gli ordini.** Per un’equazione canonica
   \[
   \varepsilon^2\psi''+q(x,\omega)\psi=0,
   \]
   la serie WKB completa può essere riorganizzata, lontano da zeri e punti di svolta, mediante un unico momento di fase pari \(u\). Esso soddisfa formalmente
   \[
   q=u^2-\varepsilon^2\frac{(u^{-1/2})''}{u^{-1/2}}.
   \]
   Il secondo termine è esattamente il funzionale di potenziale quantistico di Madelung/Bohm valutato sull’ampiezza WKB completa \(A=u^{-1/2}\). Espandendolo in potenze pari di \(\varepsilon\) si ottiene una gerarchia ricorsiva compatta \(Q_{2},Q_{4},\ldots\).

2. **No, questo non dimostra una famiglia di potenziali fisici indipendenti.** La lettura più parsimoniosa è che vi sia **un solo** funzionale di Madelung, la cui valutazione auto-consistente genera tutti i coefficienti WKB. I singoli \(Q_{2j}\) dipendono dalla coordinata, dalla variabile master e dalla convenzione con cui si introduce il parametro formale \(\varepsilon\).

3. **Non esiste una corrispondenza canonica nota \(Q_{2j}(x)\leftrightarrow\Lambda_j\).** Le correzioni \(\Lambda_j\) di Iyer–Will, Konoplya e Matyjasek–Opala sono scalari della condizione di quantizzazione uniforme, costruiti dal getto del potenziale al massimo della barriera e dall’indice di overtone. Non sono potenziali locali. Inoltre, proprio nella regione di coalescenza dei turning point, la serie locale in inverse potenze di \(q\) diventa singolare e va sostituita dalla forma normale parabolico-cilindrica.

4. **Per i QNM il fluido reale è aperto, oppure il fluido è complesso.** Con frequenza complessa, una decomposizione con densità reale positiva porta a una legge di continuità con sorgente. Mantenendo invece la forma WKB esatta si ottengono ampiezza e momento complessi, ma si perde l’interpretazione ordinaria di densità e velocità reali.

5. **La ricerca bibliografica non ha identificato un precedente diretto** che reinterpreti la gerarchia WKB dei QNM gravitazionali come gerarchia di potenziali di Madelung. Esistono tre filoni adiacenti ma distinti: (i) WKB di ordine elevato e Padé per QNM; (ii) equazione di Hamilton–Jacobi quantistica/exact WKB; (iii) buchi neri acustici in condensati di Bose–Einstein, dove il termine di pressione quantistica è fisicamente microscopico. L’assenza riscontrata è un risultato di ricerca, non una prova matematica di inesistenza.

**Verdetto:** l’ipotesi è **confermata in senso formale-locale e complessificato**, ma **non confermata in senso forte, canonico e fisico**. Il candidato più robusto per una formulazione invariante non è una collezione di potenziali locali, bensì il momento quantistico come 1-forma e i suoi periodi/Voros symbols.

---

## 1. Domande di ricerca e livelli dell’ipotesi

La domanda iniziale contiene almeno tre ipotesi diverse, che vanno valutate separatamente.

| Ipotesi | Contenuto | Esito |
|---|---|---|
| H1 — chiusura formale | La serie WKB completa può essere espressa mediante un funzionale di tipo Madelung e una ricorsione compatta. | **Confermata localmente** lontano da turning point e zeri. |
| H2 — potenziali per ordine | A ogni correzione WKB/QNM corrisponde un potenziale locale autonomo e canonico. | **Non confermata; in generale non unica.** |
| H3 — fluido fisico | La gerarchia descrive correzioni dispersive di un fluido reale efficace associato ai QNM. | **Aperta ma fortemente vincolata.** Per QNM complessi serve un flusso con sorgente o una complessificazione. |

Questa separazione impedisce di trasformare un’identità algebrica utile in una dichiarazione ontologica non dimostrata.

---

## 2. Convenzioni e dominio di validità

Si considera una forma canonica senza derivata prima,
\[
\varepsilon^2\frac{d^2\psi}{dx^2}+q(x,\omega)\psi=0,
\qquad q=\omega^2-V(x)
\]
nel caso più semplice. Per perturbazioni più generali \(q\) può dipendere da \(\omega\) in modo non quadratico o essere complesso. Il tempo è assunto proporzionale a \(e^{-i\omega t}\); per un QNM stabile \(\operatorname{Im}\omega<0\).

Il parametro \(\varepsilon\) è un parametro formale di conteggio. Nella letteratura WKB dei QNM viene spesso posto uguale a uno al termine del calcolo. L’exact WKB recente sottolinea che la sua introduzione non è unica: una qualunque interpretazione “fisica” dei coefficienti deve quindi specificare coordinata tortoise, variabile master, normalizzazione e scala semiclassica.

I risultati locali sotto assumono un dominio semplicemente connesso in cui \(q\neq0\), una scelta di ramo di \(\sqrt q\) e l’assenza di zeri di \(\psi\). Turning point, linee di Stokes, poli e condizioni al contorno globali richiedono continuazione analitica e metodi uniformi.

---

## 3. La gerarchia WKB standard: cosa viene usato e cosa non viene ridimostrato

Schutz e Will hanno adattato la WKB di barriera ai QNM di buchi neri. Iyer e Will hanno costruito la sistematica di ordine superiore mediante matching uniforme presso il massimo. Konoplya l’ha estesa fino al sesto ordine; Matyjasek e Opala fino al tredicesimo con risommazione Padé; lavori successivi hanno sviluppato ordini molto elevati, ricette operative e codici efficienti.

La condizione di quantizzazione è usualmente scritta nella forma schematica
\[
i\frac{Q_0}{\sqrt{2Q_0''}}-\sum_{k=2}^{N}\Lambda_k=n+\frac12,
\]
dove \(Q=\omega^2-V\), il pedice zero indica il massimo della barriera, le derivate sono rispetto alla coordinata tortoise e ogni \(\Lambda_k\) è un polinomio razionale nel getto \(Q_0^{(3)},\ldots,Q_0^{(2k)}\), con dipendenza da \(\alpha=n+1/2\).

Questo rapporto non ripete le lunghe derivazioni dei \(\Lambda_k\). Ne usa invece tre proprietà strutturali:

- sono coefficienti **globali di matching/quantizzazione**, non campi locali;
- dipendono dal getto del potenziale al massimo della barriera e dall’overtone;
- la serie è asintotica e, a ordine finito, richiede selezione dell’ordine o risommazione.

---

## 4. Proposizione: chiusura locale di Madelung a tutti gli ordini

### 4.1 Equazione di Riccati esatta

Definiamo il momento quantistico
\[
p(x,\varepsilon)=-i\varepsilon\frac{\psi'}{\psi}.
\]
L’equazione d’onda è equivalente, dove \(\psi\neq0\), a
\[
p^2-i\varepsilon p'=q. \tag{4.1}
\]
Ponendo \(p=\sum_{n\ge0}\varepsilon^n p_n\), si ottiene la ricorsione standard
\[
p_0^2=q,\qquad
p_n=\frac{i p_{n-1}'-\sum_{k=1}^{n-1}p_kp_{n-k}}{2p_0}.
\tag{4.2}
\]

### 4.2 Separazione pari/dispari

Scriviamo \(p=u+v\), con \(u\) somma dei termini pari e \(v\) somma dei termini dispari in \(\varepsilon\). La parte dispari della (4.1) dà
\[
v=\frac{i\varepsilon}{2}\frac{u'}{u}. \tag{4.3}
\]
Segue immediatamente
\[
\psi=C\,u^{-1/2}\exp\!\left(\frac{i}{\varepsilon}\int^x u(s)\,ds\right). \tag{4.4}
\]
La parte pari fornisce l’identità formale esatta
\[
q=u^2+\frac{\varepsilon^2}{2}\frac{u''}{u}
-\frac{3\varepsilon^2}{4}\left(\frac{u'}{u}\right)^2. \tag{4.5}
\]
Se \(A=u^{-1/2}\), allora
\[
\boxed{\;q=u^2+Q_{\rm M}[A],\qquad
Q_{\rm M}[A]=-\varepsilon^2\frac{A''}{A}\;} \tag{4.6}
\]
che è precisamente la forma Hamilton–Jacobi del potenziale quantistico di Madelung, nelle unità e convenzioni adottate.

**Interpretazione rigorosa.** Tutta la gerarchia WKB locale è generata da un solo funzionale \(Q_{\rm M}\), valutato sull’ampiezza completa e dipendente da \(\varepsilon\). Non è necessario postulare un nuovo funzionale a ogni ordine.

### 4.3 Ricorsione compatta per la gerarchia efficace

Poniamo
\[
u=\sum_{j\ge0}\varepsilon^{2j}u_j,\qquad u_0=\sqrt q,
\]
e definiamo
\[
F[u]=\frac{(u^{-1/2})''}{u^{-1/2}}
=\frac34\left(\frac{u'}u\right)^2-\frac12\frac{u''}{u}.
\]
La (4.6) equivale a \(u^2=q+\varepsilon^2F[u]\). Per \(n\ge1\),
\[
\boxed{\;u_n=
\frac{[\varepsilon^{2n-2}]F[u]-\sum_{j=1}^{n-1}u_ju_{n-j}}
{2u_0}\;} . \tag{4.7}
\]
Qui \([\varepsilon^m]\) denota l’estrazione del coefficiente. Definendo
\[
Q_{\rm eff}\equiv Q_{\rm M}[u^{-1/2}]
=\sum_{n\ge1}\varepsilon^{2n}Q_{2n},
\]
si ha
\[
Q_{2n}=-[\varepsilon^{2n-2}]F[u]. \tag{4.8}
\]
Il primo livello è
\[
Q_2=\frac{q''}{4q}-\frac{5(q')^2}{16q^2},\qquad
u_1=\frac{5(q')^2}{32q^{5/2}}-\frac{q''}{8q^{3/2}}. \tag{4.9}
\]
La ricorsione (4.7) è la formulazione matematica generale cercata: compatta, algoritmica e controllabile simbolicamente.

### 4.4 Che cosa dimostra — e che cosa non dimostra

La (4.6) dimostra una **chiusura Madelung all-order** della WKB locale. Non dimostra:

- che i coefficienti \(Q_{2n}\) siano osservabili separati;
- che siano invarianti per trasformazioni di Liouville;
- che restino regolari ai turning point;
- che coincidano con i \(\Lambda_n\) della quantizzazione di barriera;
- che descrivano nuovi termini microscopici nella relazione di dispersione.

---

## 5. Frequenze complesse: fluido con sorgente o fluido complesso

Per mostrare l’ostacolo senza ambiguità, imponiamo una decomposizione con ampiezza e fase reali:
\[
\psi=\sqrt\rho\,e^{i\theta/\varepsilon},\qquad q=q_R+iq_I.
\]
Separando parti reale e immaginaria si ottiene esattamente
\[
(\theta')^2-\varepsilon^2\frac{(\sqrt\rho)''}{\sqrt\rho}=q_R, \tag{5.1}
\]
\[
(\rho\theta')'=-\frac{q_I}{\varepsilon}\rho. \tag{5.2}
\]
La prima è un’equazione di Hamilton–Jacobi con potenziale quantistico; la seconda è una legge di continuità con sorgente/pozzo. Per \(q=\omega^2-V\), con \(V\) reale, \(q_I=\operatorname{Im}(\omega^2)\).

Quindi un QNM non è in generale un flusso stazionario conservativo. Le condizioni outgoing/ingoing e \(\operatorname{Im}\omega<0\) rendono inoltre le autofunzioni divergenti alle estremità del dominio tortoise: non sono stati normalizzabili in \(L^2\). Una densità positiva globale va perciò maneggiata come oggetto locale o come variabile di un sistema aperto.

Alternativamente si mantiene la (4.4) con \(u\) complesso. L’identità exact-WKB sopravvive, ma \(A=u^{-1/2}\) e il “momento del fluido” sono complessi. Questa è una descrizione matematica coerente, non ancora una idrodinamica fisica ordinaria.

---

## 6. Perché \(Q_{2n}(x)\) non coincide canonicamente con \(\Lambda_n\)

### 6.1 Non uniformità presso il massimo

I coefficienti (4.9) contengono inverse potenze di \(q\). Per i QNM di barriera a basso overtone, i turning point coalescono vicino al massimo e \(q(x_0)\) è piccolo nell’opportuna scalatura semiclassica. La serie locale diventa singolare proprio dove il matching determina la frequenza.

Iyer–Will risolvono il problema espandendo il potenziale al massimo, riscalando \(x-x_0\sim\sqrt\varepsilon\) e riducendo l’equazione a una forma normale di cilindro parabolico perturbata. Le \(\Lambda_n\) sono invarianti della procedura di forma normale/matching, non valori regolari dei \(Q_{2n}(x)\) locali.

### 6.2 Non unicità inversa

A ordine \(N\), la condizione di quantizzazione usa un numero finito di derivate di \(Q\) in \(x_0\). Infinite funzioni \(\delta V_N(x)\) possono condividere lo stesso getto fino all’ordine richiesto e produrre la stessa correzione scalare. Pertanto il problema inverso “dato \(\Lambda_N\), trovare il potenziale efficace locale” è sottodeterminato senza una scelta di ansatz o di gauge.

### 6.3 Dipendenza da trasformazioni di Liouville

Sotto un cambio di coordinata \(x=x(z)\) e la ridefinizione canonica \(\phi(z)=\psi(x(z))/\sqrt{x'(z)}\),
\[
\varepsilon^2\phi''+\widetilde q(z)\phi=0,
\]
con
\[
\widetilde q=(x')^2q+\frac{\varepsilon^2}{2}\{x,z\}, \tag{6.1}
\]
dove \(\{x,z\}\) è la derivata schwarziana. Un termine di ordine \(\varepsilon^2\) può dunque essere spostato tra “potenziale di fondo” e “correzione quantistica” mediante una trasformazione ammessa che lascia equivalente l’equazione.

Ne segue che un singolo \(Q_{2n}(x)\) non è, da solo, un oggetto geometrico invariante. Oggetti più robusti sono la 1-forma del momento quantistico, le sue integrazioni su cicli e i dati di Stokes/Voros dell’exact WKB.

### 6.4 Proposizione di non unicità — status logico

La non unicità appena esposta è una deduzione matematica di questo rapporto, non un teorema attribuito agli autori WKB citati. Essa segue da tre fatti verificabili: singolarità della gerarchia locale ai turning point, carattere finito del getto usato da \(\Lambda_N\), termine schwarziano nelle trasformazioni di Liouville.

---

## 7. “Correzioni dispersive”: significato consentito e significato eccessivo

Nel formalismo di Madelung della Schrödinger standard, il termine di pressione/potenziale quantistico è già un termine di gradiente di ordine due. Iterarlo nella soluzione WKB genera coefficienti contenenti derivate sempre più alte di \(q\). Questa è una **gerarchia asintotica della soluzione**.

Non segue automaticamente che l’equazione microscopica contenga nuove potenze di \(k\) nella relazione di dispersione. Una vera gerarchia costitutiva di gradienti superiori — come in modelli non locali di fase o in teorie efficaci con operatori ad alte derivate — modifica l’equazione dinamica. I due concetti possono essere confrontati, ma non identificati senza un’ulteriore derivazione.

La formulazione prudente è quindi: i \(Q_{2n}\) **organizzano correzioni di gradiente/dispersive efficaci della soluzione**, ma non costituiscono di per sé nuovi coefficienti di dispersione del mezzo.

---

## 8. Revisione della letteratura e matrice delle lacune

### 8.1 Filone WKB/QNM

- Schutz–Will: applicazione semianalitica WKB ai QNM di barriera.
- Iyer–Will: sviluppo sistematico di ordine superiore e matching uniforme.
- Konoplya: sesto ordine e regime pratico di accuratezza, tipicamente migliore per bassi overtone e \(\ell>n\).
- Matyjasek–Opala: tredicesimo ordine e Padé; Matyjasek–Telecka: ordini molto elevati.
- Hatsuda: Borel summability in casi di Schwarzschild/Reissner–Nordström e connessione con oscillatori anarmonici.
- Exact WKB recente: geometria globale di Stokes, periodi quantistici e ambiguità dell’introduzione del parametro semiclassico.

### 8.2 Filone Madelung/exact WKB

Madelung e Takabayasi stabiliscono la trasformazione idrodinamica; la letteratura moderna su quantum Hamilton–Jacobi ed exact WKB usa la stessa equazione di Riccati e la decomposizione pari/dispari. Il collegamento al funzionale \(-\varepsilon^2A''/A\) è quindi strutturalmente naturale. I lavori matematici su Madelung–Gross–Pitaevskii–Korteweg chiariscono anche il ruolo dei gradienti e i problemi agli zeri della densità.

### 8.3 Filone analog-gravity

Nei buchi neri acustici in condensati di Bose–Einstein esistono QNM con dispersione di Bogoliubov e una pressione quantistica microscopica reale. Questi lavori mostrano che QNM e idrodinamica dispersiva possono coesistere fisicamente, ma studiano un mezzo analogico, non reinterpretano le correzioni WKB dei QNM gravitazionali come potenziali di Madelung.

### 8.4 Esito della ricerca d’intersezione

Sono state cercate combinazioni di termini fra “Madelung”, “Bohm quantum potential”, “quasinormal mode”, “black hole”, “higher-order WKB”, “Schutz–Will” e “Iyer–Will” in motori bibliografici e archivi preprint. Sono stati inoltre seguiti i riferimenti dei lavori exact-WKB e analog-gravity più vicini.

**Nessun lavoro diretto identificato al 23 agosto 2026** formula e verifica la gerarchia qui proposta per i QNM gravitazionali. Il risultato è qualificato: l’indicizzazione incompleta, terminologie diverse o lavori non digitalizzati possono produrre falsi negativi.

| Tema | Cosa esiste | Lacuna residua | Confidenza |
|---|---|---|---|
| WKB QNM ad alto ordine | Formule fino ad alto/altissimo ordine, Padé, Borel, codici | Interpretazione Madelung dei coefficienti | Alta |
| Madelung + WKB | Riccati, quantum Hamilton–Jacobi, exact WKB | Applicazione sistematica ai QNM di barriera | Alta |
| QNM + fluidi | Buchi neri acustici/BEC, dispersione Bogoliubov | Ponte con i \(\Lambda_n\) gravitazionali | Alta |
| Gerarchie di potenziali quantistici | Modelli non locali/gradienti superiori | Equivalenza con gerarchia WKB | Alta che l’equivalenza non sia dimostrata |
| Precedente diretto | Non identificato | Ricerca bibliografica mai logicamente esaustiva | Media |

---

## 9. Programma di ricerca falsificabile

### Fase A — fissare il gauge matematico

1. Scegliere metrica, coordinata tortoise e variabile master.
2. Dichiarare come viene introdotto \(\varepsilon\).
3. Specificare classe di trasformazioni di Liouville ammesse.
4. Distinguere esplicitamente quantità locali da periodi e dati di connessione.

### Fase B — costruire la gerarchia locale

1. Generare \(u_n\) e \(Q_{2n}\) con la (4.7).
2. Verificare simbolicamente l’equazione di Riccati a ogni ordine.
3. Studiare crescita fattoriale, singolarità e trasformazioni di gauge.
4. Costruire un indicatore locale \(|Q_{\rm M}|/|q|\) come diagnostica, senza presentarlo come bound d’errore.

### Fase C — effettuare l’uniformizzazione

1. Passare alla forma normale parabolico-cilindrica al massimo.
2. Tradurre i coefficienti della forma normale in combinazioni della gerarchia di Riccati.
3. Verificare se le combinazioni finite riproducono \(\Lambda_2,\Lambda_3,\ldots\).
4. Se possibile, definire non “potenziali locali” ma **invarianti di stress/forma normale**.

### Fase D — complessità e interpretazione idrodinamica

1. Confrontare il sistema reale con sorgente (5.1)–(5.2) e la chiusura complessa (4.6).
2. Identificare un flusso o bilineare che rimanga finito sotto regolarizzazione dei QNM.
3. Testare covarianza rispetto a trasformazioni di Darboux/isospectralità fra potenziali master.

### Fase E — benchmark e validazione

- Barriera di Pöschl–Teller, con spettro QNM esatto.
- Regge–Wheeler e Zerilli per Schwarzschild, come test di isospectralità.
- Reissner–Nordström e un caso con potenziale dipendente dalla frequenza.
- Doppia barriera o geometria “dirty” per separare effetti locali e globali.
- Confronto con continued fractions/integrazione diretta e con Padé/Borel–Padé.

### Criteri di confutazione

La versione forte dell’ipotesi va respinta se:

1. i presunti potenziali cambiano sotto trasformazioni ammesse mentre lo spettro resta invariato;
2. dopo aver fissato il gauge e uniformizzato, la costruzione non riproduce almeno \(\Lambda_2\) e \(\Lambda_3\);
3. non esiste una quantità di densità/flusso con significato coerente per frequenze complesse;
4. nel benchmark esatto la gerarchia non migliora né organizza le approssimazioni rispetto alla WKB ordinaria.

---

## 10. Applicazioni plausibili, ordinate per solidità

1. **Organizzazione simbolica e diagnostica della WKB** — alta plausibilità. La ricorsione compatta può generare coefficienti e mettere in evidenza dove i gradienti diventano grandi.
2. **Descrizione invariante tramite periodi quantistici** — alta plausibilità matematica. Può collegare la lettura di Madelung all’exact WKB senza attribuire realtà fisica ai singoli coefficienti locali.
3. **Sensibilità dello spettro alla geometria della barriera** — plausibilità media. I termini della gerarchia separano famiglie di derivate e possono aiutare a classificare correzioni EFT o ambientali.
4. **Controllo di isospectralità e scelte di variabile master** — plausibilità media. Il fallimento della covarianza locale diventa un test, non un difetto nascosto.
5. **Ponte con analog gravity** — esplorativo. Nei BEC il potenziale quantistico ha significato microscopico; il confronto può suggerire osservabili, ma non prova che il buco nero gravitazionale sia un fluido di Madelung.
6. **Nuova idrodinamica fisica dei QNM** — bassa plausibilità allo stato attuale. Richiede una teoria del mezzo, un’azione o una derivazione costitutiva aggiuntiva.

---

## 11. Conclusione

La congettura iniziale contiene un nucleo matematico solido: la serie WKB locale completa ammette una chiusura Madelung esatta, e la famiglia \(Q_{2n}\) è generata ricorsivamente dall’espansione di un unico funzionale di potenziale quantistico. Questa formulazione è nuova come sintesi applicata ai QNM, per quanto risulta dalla ricerca svolta.

Il passaggio da questa identità a una gerarchia fisica e canonica non è però giustificato. Le correzioni di barriera \(\Lambda_n\) sono oggetti di matching, la regione decisiva è non uniforme, la decomposizione locale è dipendente dal gauge di Liouville e le frequenze QNM rendono il flusso non conservativo o complesso.

La direzione più robusta è quindi duplice: usare la gerarchia locale come strumento di organizzazione e diagnostica; cercare la formulazione fisicamente significativa negli invarianti globali dell’exact WKB — momenti quantistici, periodi e dati di Stokes — oppure in un sistema idrodinamico aperto esplicitamente definito. Il programma proposto è falsificabile: può confermare una versione ristretta e utile dell’ipotesi oppure dimostrare che l’analogia idrodinamica non sopravvive oltre il livello formale.

---

## 12. Registro delle affermazioni

| Affermazione | Tipo | Base | Confidenza |
|---|---|---|---|
| La WKB locale all-order si chiude nel singolo funzionale \(Q_{\rm M}\) | Identità matematica | Eq. (4.1)–(4.8); exact-WKB pari/dispari | Alta |
| La mappa \(\Lambda_n\to Q_{2n}(x)\) non è canonica | Deduzione | Non uniformità, getto finito, Schwarziana | Alta |
| Un QNM complesso induce una continuità con sorgente per \(\rho\) reale | Identità matematica | Eq. (5.1)–(5.2) | Alta |
| I \(Q_{2n}\) non sono automaticamente nuovi termini microscopici dispersivi | Distinzione concettuale | Ordine dell’equazione vs ordine dell’espansione | Alta |
| Non esiste un precedente diretto | Esito bibliografico qualificato | Ricerche incrociate e catene di citazioni | Media |
| I periodi quantistici sono candidati più invarianti | Proposta motivata | Exact WKB e trasformazioni di Liouville | Medio-alta |
| Esiste una idrodinamica fisica nuova dei QNM | Ipotesi | Non dimostrata | Bassa/da testare |

---

## Bibliografia essenziale

[1] Schutz, B. F.; Will, C. M., “Black Hole Normal Modes: A Semianalytic Approach”, *Astrophysical Journal Letters* 291 (1985) L33. [DOI](https://doi.org/10.1086/184453) · [testo ADS](https://adsabs.harvard.edu/pdf/1985ApJ...291L..33S)

[2] Iyer, S.; Will, C. M., “Black-hole normal modes: A WKB approach. I. Foundations and application of a higher-order WKB analysis”, *Physical Review D* 35 (1987) 3621. [DOI](https://doi.org/10.1103/PhysRevD.35.3621)

[3] Iyer, S., “Black-hole normal modes: A WKB approach. II. Schwarzschild black holes”, *Physical Review D* 35 (1987) 3632. [DOI](https://doi.org/10.1103/PhysRevD.35.3632) · [testo APS](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.35.3632/fulltext)

[4] Konoplya, R. A., “Quasinormal behavior of the D-dimensional Schwarzschild black hole and higher order WKB approach”, *Physical Review D* 68 (2003) 024018. [arXiv](https://arxiv.org/abs/gr-qc/0303052) · [DOI](https://doi.org/10.1103/PhysRevD.68.024018)

[5] Matyjasek, J.; Opala, M., “Quasinormal modes of black holes: The improved semianalytic approach”, *Physical Review D* 96 (2017) 024011. [arXiv](https://arxiv.org/abs/1704.00361) · [DOI](https://doi.org/10.1103/PhysRevD.96.024011)

[6] Konoplya, R. A.; Zhidenko, A.; Zinhailo, A. F., “Higher order WKB formula for quasinormal modes and grey-body factors: recipes for quick and accurate calculations”, *Classical and Quantum Gravity* 36 (2019). [arXiv](https://arxiv.org/abs/1904.10333)

[7] Matyjasek, J.; Telecka, M., “Quasinormal modes of black holes. II. Padé summation of the higher-order WKB terms”, *Physical Review D* 100 (2019) 124006. [arXiv](https://arxiv.org/abs/1908.09389)

[8] Hatsuda, Y., “Quasinormal modes of black holes and Borel summation”, *Physical Review D* 101 (2020) 024008. [arXiv](https://arxiv.org/abs/1906.07232)

[9] Hatsuda, Y.; Kimura, M., “Spectral Problems for Quasinormal Modes of Black Holes”, *Universe* 7 (2021) 476. [arXiv](https://arxiv.org/abs/2111.15197)

[10] Miyachi, T.; Namba, R.; Omiya, H.; Oshita, N., “Path to an exact WKB analysis of black hole quasinormal modes”, (2025). [arXiv](https://arxiv.org/abs/2503.17245)

[11] Hatsuda, Y.; Shiga, Y., “Exact WKB and Quantum Periods for Extremal Black Hole Quasinormal Modes”, (2026). [arXiv](https://arxiv.org/abs/2605.01321)

[12] Konoplya, R. A.; Matyjasek, J.; Zhidenko, A., “An efficient higher-order WKB code for quasinormal modes and greybody factors”, (2026). [arXiv](https://arxiv.org/abs/2603.12466) · [DOI](https://doi.org/10.53941/ijgtp.2026.100005)

[13] Madelung, E., “Quantentheorie in hydrodynamischer Form”, *Zeitschrift für Physik* 40 (1927) 322–326. [DOI](https://doi.org/10.1007/BF01400372) · [traduzione inglese](https://www.neo-classical-physics.info/uploads/3/0/6/5/3065888/madelung_-_hydrodynamical_interp..pdf)

[14] Takabayasi, T., “On the Formulation of Quantum Mechanics associated with Classical Pictures”, *Progress of Theoretical Physics* 8 (1952) 143. [DOI](https://doi.org/10.1143/ptp/8.2.143)

[15] Carles, R.; Danchin, R.; Saut, J.-C., “Madelung, Gross–Pitaevskii and Korteweg”, *Nonlinearity* 25 (2012) 2843. [arXiv](https://arxiv.org/abs/1111.4670) · [DOI](https://doi.org/10.1088/0951-7715/25/10/2843)

[16] Türe, M.; Ünsal, M., “Quantum Hamilton–Jacobi Theory, Spectral Path Integrals and Exact-WKB”, *Physical Review D* 111 (2025) 105010. [arXiv](https://arxiv.org/abs/2406.07829) · [DOI](https://doi.org/10.1103/PhysRevD.111.105010)

[17] Silva, H. O. et al., “Quasinormal-mode filters: A new approach to analyze the black-hole ringdown”, *Physical Review D* 110 (2024) 024042. [arXiv](https://arxiv.org/abs/2404.11110)

[18] Barceló, C. et al., “Quasi-normal mode analysis in BEC acoustic black holes”, *Physical Review D* 75 (2007) 084024. [arXiv](https://arxiv.org/abs/gr-qc/0701173) · [DOI](https://doi.org/10.1103/PhysRevD.75.084024)

[19] Daghigh, R. G.; Green, M. D., “High Overtone Quasinormal Modes of Analog Black Holes and the Small Scale Structure of the Background Fluid”, (2015). [arXiv](https://arxiv.org/abs/1411.7066)

[20] Mauri, R., “A Non-local Phase Field Model of Bohm’s Quantum Potential”, *Foundations of Physics* 51 (2021). [DOI](https://doi.org/10.1007/s10701-021-00454-9)

[21] NIST Digital Library of Mathematical Functions, §1.13, “Differential Equations — Liouville Transformation”. [DLMF](https://dlmf.nist.gov/1.13)

---

## Nota metodologica sulla ricerca

La revisione ha privilegiato articoli originali, DOI, arXiv e documentazione matematica primaria. Le affermazioni negative sono formulate come “non identificato” e non come “non esiste”. Le deduzioni originali — in particolare la non unicità della ricostruzione locale e il sistema con sorgente per \(q\) complesso — sono esplicitamente marcate e non attribuite retroattivamente alle fonti. Non è stata assunta validità uniforme della WKB locale nei pressi dei turning point.
