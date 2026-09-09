# Fonti primarie lette, e cosa hanno cambiato

**Data:** 9 settembre 2026.
**Natura:** **letteratura**. Nessun risultato nuovo nostro; il valore è nel
passaggio di tre sezioni del manoscritto da argomentate a fondate, e in una
verifica di non precedenza.
**PDF:** in `Schw-QNM-WKB-Fluid/papers/` (cartella in `.gitignore`: i lavori sotto
paywall non vengono ridistribuiti dal repo pubblico).

---

## 1. Khesin, Misiołek, Modin — PNAS 115, 6165 (2018), arXiv:1711.00321

**Contenuto.** La trasformazione di Madelung è un **simplettomorfismo**
$T^*\mathrm{Dens}(M)\to PC^\infty(M,\mathbb C\setminus\{0\})$ (Thm 4.2) e
un'**isometria** fra metrica di Sasaki–Fisher-Rao e Fubini–Study (Thm 4.9).
Il termine di Bohm compare nel potenziale dell'equazione di Newton come
$4I(\varrho)$, con $I$ il funzionale di **informazione di Fisher** (Prop. 4.3).

**Ipotesi, e come il problema QNM le viola.**

| ipotesi di [KMM] | caso QNM |
|---|---|
| $\int_M\varrho=1$ | i QNM **non sono normalizzabili** |
| $M$ compatta | dominio radiale non compatto |
| flusso hamiltoniano, norma conservata | $\omega$ complessa: $(\rho v)'=-\Gamma\rho$, $\Gamma\neq0$ (§8) |
| $\mathbb C\setminus\{0\}$ | ✓ coincide con l'ipotesi $\psi\neq0$ del nostro Teorema 1 |

**Effetto sul manoscritto.** §1 riscritto: i §7 e §8 non sono limiti del nostro
trattamento, sono le **ostruzioni identificate** che separano il problema QNM dal
caso in cui la geometria è nota. §3: $Q_M$ è, variazionalmente, l'informazione di
Fisher — risposta a una domanda che il §12 poneva come aperta. §12 aggiornato di
conseguenza.

**Errore mio da registrare.** Avevo proposto di riposizionare il lavoro dicendo
che «non aggiunge informazione» diventa un teorema. **Non è lecito**: le ipotesi
del teorema le violiamo. Tentazione da evitare anche in futuro: usare l'isometria
per «spiegare» la ridondanza di $\mathcal E_M$.

**Appendice sulla momentum map (§4.6, da Fusca 2015, arXiv:1512.04611).** La
Madelung inversa è la momentum map dell'azione di
$S=\mathrm{Diff}(M)\ltimes C^\infty(M)$ su mezze densità:
$\mathbf M(\psi)=(\operatorname{Im}\bar\psi\nabla\psi,\ \bar\psi\psi)=(\rho\nabla\theta,\rho)$.
Cioè $\rho$ è la carica di Noether della **fase locale** e $\rho v$ quella dei
**diffeomorfismi**. Madelung è anche mappa di Poisson verso Lie–Poisson su
$\mathfrak s^*$.

*Non usato nel manoscritto, e la ragione va ricordata*: il quadro è per
Schrödinger, primo ordine in $t$; la nostra corrente è di tipo Klein–Gordon.
Si vede dal fattore: Schrödinger darebbe $(A^2S')'=-2\operatorname{Im}(\omega)A^2$,
noi abbiamo $-\operatorname{Im}(\omega^2)A^2$, che differiscono per
$\operatorname{Re}\omega$. L'idea che $\Gamma$ sia un'anomalia della momentum map
è **plausibile e non dimostrata**: richiede di costruire la versione KG.

---

## 2. Iyer & Will — Phys. Rev. D 35, 3621 (1987), Parte I

**Perché serviva.** Il §6 argomentava in tre punti che non esiste corrispondenza
$Q_{2j}\leftrightarrow\Lambda_j$, senza aver letto la fonte dei $\Lambda_j$.

**Cosa dice.** Il metodo è "carried to third order beyond the eikonal
approximation" — coerente con la riga (O3) del nostro §2.1. Nella regione fra i
turning point approssimano $-Q$ con uno sviluppo di Taylor **fino alla derivata
sesta inclusa**: sono quelle sei derivate al picco a produrre $\Lambda_2$ e
$\Lambda_3$. La modifica del WKB serve "per valori della frequenza tali che i
turning point $x_1$ e $x_2$ siano vicini fra loro, presso il picco", dove "il
matching standard non è più valido". La soluzione interna, senza i termini in
$\varepsilon$, è una **funzione del cilindro parabolico** $D_\nu(t)$, cercata poi
come $f(t)D_\nu[g(t)]$.

**Effetto.** I punti (i) e (iii) del §6 diventano citazioni. La frase chiave: le
due gerarchie vivono in regioni **disgiunte per costruzione**.

**Convenzione chiarita.** [7] costruisce i $\Lambda_j$ dalle derivate di
$Q=\Omega^2-V$, il nostro codice da quelle di $V$. Coincidono perché
$Q^{(k)}=-V^{(k)}$ per $k\geq1$ e le formule contengono solo rapporti in cui il
segno si cancella. Era un'ambiguità silenziosa nel codice; ora annotata in
Appendice A.

---

## 3. Seidel & Iyer — Phys. Rev. D 41, 374 (1990), Parte IV (Kerr)

**Perché serviva.** Verificare che il §9 non fosse preceduto.

**Esito: non lo è.** Applicano il WKB di barriera di [7] a Kerr per i modi
**bassi**, espandendo il potenziale **in potenze di $a\omega$**. Occorrenze in
tutto l'articolo: `spheroidal` 0, `eikonal` 0, `Langer` 0, `separation constant`
0, `large l` 0.

Il loro parametro piccolo è la rotazione; il nostro è $1/L$. È il terzo caso nel
progetto di espansioni distinte sullo stesso problema, dopo Capuano *et al.*
(tasso di massa) e la nostra.

**Contrasto utile, aggiunto al §9.5.** Le difficoltà che dichiarano ad $a$ alto
derivano proprio dall'espansione in $a\omega$. Le nostre pendenze ad $a=0.9$ sono
pulite perché il problema angolare è risolto esattamente, senza espandere in
$a\omega$.

---

## 4. Berry & Mount — Rep. Prog. Phys. 35, 315 (1972)

Il più produttivo dei quattro: ha promosso tre sezioni.

**(a) Unicità della sostituzione di Langer (§5.1) → nostra Appendice B.**
Considerano $\ell(\ell+1)\to\ell(\ell+1)+\alpha$ con $L=(\ell+\beta)\hbar$.
Richiedere che il centrifugo diventi esattamente $L^2$ impone
$\ell^2+\ell+\alpha=(\ell+\beta)^2$, cioè $\beta=\tfrac12$, $\alpha=\tfrac14$:
**unica**. La indicano come "la più naturale per problemi che coinvolgono il
limite semiclassico". Origine: Langer (1937) mostrò che la WKB ingenua non è
applicabile presso $r=0$; la sostituzione $r=e^x$ ne ripristina la validità. Il
$-\tfrac14$ del nostro $v_2$ è quel termine.

**(b) Principio di equivalenza dei turning point (§4.1) → nostro §6.**
*"In the semiclassical limit all problems are equivalent which have the same
classical turning-point structure."* Ne segue che una costruzione per equazione
di confronto codifica la **struttura dei turning point**, non il potenziale
locale; i $Q_{2j}$ distinguono potenziali che il metodo uniforme identifica.
È l'argomento che il §6 cercava: gli oggetti differiscono per **ciò che sono
costruiti per registrare**, non solo per tipo.

**(c) Meccanismo della rottura (§4.2) → nostro §7.**
La WKB "non è valida presso un turning point classico dove $p^2(x)$ ha uno zero,
perché la mappa non è più biunivoca": il termine trascurato diventa dominante.
Al massimo di barriera i due zeri si fondono in un **turning point del secondo
ordine**, il caso $\nu=2$ di Langer, con equazione di confronto del cilindro
parabolico — la stessa di [7].

**Ridimensionamento onesto, aggiunto al §7.** La degradazione che misuriamo è la
firma numerica di un fenomeno **noto e classificato**. L'interesse non è la
scoperta, ma il confronto quantitativo fra le due righe della tabella.

---

## 5. Delabaere, Dillinger & Pham — J. Math. Phys. 38, 6126 (1997)

Il più consequenziale dei cinque: tocca §1.1, §3, §4, §7 e §12.

**(a) La nostra chiusura è la loro eq. (1.2).** Scrivono
$w=P(q,\hbar^2)^{-1/2}\exp\big((i/\hbar)\int P\,dq\big)$ con $P$ "definita come
la **parte pari in $\hbar$ della soluzione dell'equazione di Riccati**". Il
nostro Teorema 1 è quella formula, il nostro Teorema 2 quella caratterizzazione.

L'attribuzione che avevamo ("noto nel formalismo di Hamilton–Jacobi quantistico")
era corretta ma vaga: ora è puntuale. Nessuna perdita — non lo rivendicavamo —
ma è la differenza fra citare un'area e citare un'equazione.

**(b) La rottura è una pinzatura, non una divergenza (§III.A) → nostro §7.**
Gli sviluppi ben normalizzati "sono singolari per quei valori dell'energia per
cui il cammino di normalizzazione è *pinched* dalla confluenza di alcuni turning
point". Localmente si può sempre scegliere una base non pinzata: **ciò che si
perde è la regolarità delle mappe di connessione**. La singolarità è nel
problema di connessione, non nel singolo sviluppo. È più preciso di "serie
asintotica valutata nel proprio punto singolare", che è ciò che scrivevamo.

**(c) La riscalatura, e un ponte fra §6 e §7 che non avevo visto.**
La cura è $E=E_{\rm crit}+E_r\hbar$ presso un punto critico quadratico; il
problema riscalato dipende regolarmente da $E_r$ su tutto il piano complesso.

Ma la condizione di quantizzazione di Iyer–Will ha **esattamente quella forma**:
$q_0/\sqrt{2q_0''}=i\varepsilon(n+\tfrac12)$, cioè $q_0=O(\varepsilon)$.
Quindi **i $\Lambda_j$ sono il contenuto troncato del problema di connessione
riscalato a un turning point doppio**. I $Q_{2j}$ sono invece coefficienti dello
sviluppo *non* riscalato, quello che si pinza. Sono espansioni di oggetti
diversi — argomento indipendente da quello di Berry–Mount e che chiude il §6
dall'altro capo.

**(d) §12: il contenuto invariante ha un nome.**
Il momento vive sul rivestimento doppio $\dot{\mathbb C}_2$, superficie di
Riemann di $p=(E-V)^{1/2}$; e nello spirito della risorgenza di Écalle gli
sviluppi vanno trattati "come oggetti impliciti, il cui interesse principale sta
nei **piccoli esponenziali** che generano attraverso il processo di risorgenza".
Il §12 diceva "periodi e simboli di Voros — oggetti che sopravvivono alla
continuazione analitica": impreciso. Non sopravvivono *malgrado* le linee di
Stokes, sono definiti *attraverso* di esse.

**(e) Le risonanze sono nell'ambito.** Vi sono definite come i valori complessi
dell'energia per cui la componente riflessa si annulla: è la condizione QNM. Non
c'è forzatura nel citarli.

**Non usato:** risommazione di Borel, condizione di Zinn-Justin, espansioni
multi-istantoniche, iperasintotica di Berry. Sono il contenuto principale
dell'articolo e non ci servono, ma vanno ricordati se un giorno si volesse
davvero estrarre numeri dalla struttura risorgente.

---

## 6. Leung, Liu, Suen, Tam & Young — J. Phys. A 31, 3271 (1998)

**Attribuzione da correggere.** Citato ovunque come "Leung–Liu–Young" o "Leung
*et al.*", inclusa la valutazione di novità dell'8 settembre. Gli autori sono
**cinque**: Leung, Liu, **Suen, Tam**, Young.

**Perché serviva.** L'integrale della condizione di solvibilità di Vaidya diverge,
come tutti gli integrali QNM. Serviva la norma bilineare regolarizzata.

**Cosa dà.** L'eq. (2.16): norma generalizzata con termini di superficie, con la
proprietà che numeratore e denominatore della (2.14) sono **separatamente
indipendenti** dai punti di raccordo \(L_\pm\). La (2.19), più semplice, vale
solo per potenziali **senza code** — e loro lo dichiarano; Schwarzschild ha code
di potenza e la richiede nella forma generale.

**Effetto.** Regolarizzazione implementata e validata: vedi
[vaidya_solvability](vaidya_solvability_2026-09-09.md) §4, dove la norma risulta
indipendente da \(L_+\) allo 0.2% fra punti di raccordo che differiscono di tre
ordini di grandezza nell'integrale nudo.

**Nota strutturale utile.** Sviluppano la teoria perturbativa logaritmica proprio
perché **i QNM non formano un insieme completo**, quindi non serve sommare su
stati imperturbati. È esattamente la nostra situazione, dove l'operatore
congelato è singolare sul modo. Discutono anche Pöschl–Teller, la nostra famiglia
di controllo.

---

## 7. Stato bibliografico

`references.bib`: 37 voci. Manoscritto: 34 riferimenti numerati, dodici citati
nel testo. In `papers/`: 24 PDF.

**Letto e integrato:** Leung *et al.* 1998, Voros 1983, KMM 2018, Fusca 2015 (parziale), Iyer & Will I,
Iyer II, Seidel & Iyer IV, Berry & Mount, Cho 2003, Konoplya 2003, Cardoso 2009,
Miyachi 2025, Capuano 2024, Yoo 2025, Abdalla 2006, Meza-Domínguez 2026.

**Ancora da leggere:** Olver, capitoli sui punti di transizione — **priorità
bassa**: Berry & Mount ha dato l'inquadramento uniforme e DDP la versione esatta.
Eventuale, se si volesse la parte risorgente: Écalle, e l'iperasintotica di
Berry–Howls per estrarre numeri dalle serie divergenti.
