# Originalità e pubblicabilità: valutazione dopo l'audit spettrale

**Data della ricerca:** 8 settembre 2026.
**Oggetto:** WKB-QNM-Fluid, base Git ec97d56 più modifiche locali dell'audit.
**Natura:** ricerca bibliografica mirata e giudizio scientifico provvisorio,
non peer review, prova di priorità o previsione di accettazione.

## 1. Verdetto operativo

**Non proporrei ancora il materiale attuale come articolo che introduce
un nuovo metodo predittivo per i QNM.** I risultati sono utili e riproducibili,
ma il diagnostico usa ancora informazioni dalla soluzione di riferimento,
mentre diverse strutture presentate come possibili sviluppi hanno precedenti
molto vicini.

Non ho identificato, nelle fonti consultate, la stessa identica prescrizione
numerica del nostro difetto uniforme di curvatura. Questo non dimostra
che sia inedita e, soprattutto, non basta a dimostrare che sia utile.
Originalità di una definizione e avanzamento scientifico sono due cose diverse.

La strada più difendibile è un **articolo metodologico circoscritto ai
fondi statici**, se si ottiene uno stimatore calcolabile dalla WKB che
mostri un vantaggio verificato rispetto ai controlli già disponibili.
Kerr e Vaidya restano estensioni successive, non condizioni necessarie per
completare il primo articolo.

## 2. Che cosa è stato cercato

Sono state confrontate le note locali con fonti degli autori su arXiv,
articoli editoriali e pagine ufficiali delle riviste. Sono stati usati
anche termini matematicamente equivalenti, perché una ricerca della sola
parola «Madelung» perde molti precedenti pertinenti.

| Gruppo | Esempi di query effettivamente usate | Scopo |
|---|---|---|
| sovrapposizione diretta | “Madelung” “quasinormal”; “Bohm” “quasinormal”; “Madelung” “Vaidya” | trovare il lessico esplicito |
| WKB e diagnostica | quasinormal modes WKB amplitude phase a posteriori error estimator residual | confronti operativi già proposti |
| equivalenze | “quasinormal” “logarithmic perturbation”; resonance eigenvalue Wronskian Jost semiclassical WKB | trovare la stessa struttura con altri nomi |
| analisi numerica | “resonances” “a posteriori error”; “quasinormal” “residual” “condition number” | residui, sensibilità e problemi non autoaggiunti |
| uniformità | “residual” “WKB” “error bounds” turning points Olver | confine fra resto locale e controllo spettrale |
| sede editoriale | pagine ufficiali CQG e PRD | distinguere pertinenza e requisiti dall'accettazione |

Limiti: ricerca non esaustiva del grafo delle citazioni; nessun accesso
sistematico a Scopus/Web of Science; possibili lavori non indicizzati;
lettura delle sezioni pertinenti, non verifica integrale di tutti i paper.
I risultati di motori generalisti sono stati usati per scoprire fonti:
le conclusioni sotto si appoggiano alle fonti primarie collegate.
Le ricerche letterali poco fruttuose **non sono prova di assenza**.

## 3. Precedenti che restringono concretamente la novità

### 3.1 Perturbazione logaritmica dei QNM: il precedente aggiuntivo decisivo

[Leung et al., preprint 1997, J. Phys. A 31 (1998) 3271](https://arxiv.org/abs/physics/9712037)
usano la Riccati della derivata logaritmica per calcolare correzioni QNM
con condizioni radiative. Le eq. (2.14)–(2.16), pagina 5 del preprint,
contengono una risposta spettrale pesata con norma bilineare e termini
di bordo; §4 tratta anche Pöschl–Teller. Pagina 3 discute già lo zero
centrale del settore dispari dei potenziali simmetrici.

**Impatto sul progetto:** il collegamento generale fra difetto locale,
correzione della frequenza e normalizzazione radiativa non è una nuova
idea del nostro audit. L'eventuale contributo deve essere nella costruzione
WKB approssimata, nel controllo dell'errore o nel costo, non nella sola
riscrittura della formula perturbativa.

Verifica: testo delle sezioni pertinenti e controllo visivo delle pagine
3 e 5 del PDF originale. La copia locale è solo materiale di consultazione,
non un allegato da redistribuire nel repository.

### 3.2 Ampiezza–fase applicata a Schwarzschild e Kerr

[Glampedakis–Andersson (2003)](https://arxiv.org/abs/gr-qc/0304030),
§II–III, impiegano funzioni di fase/Prüfer, Riccati e contorni complessi.
Le eq. (6)–(8) collegano ampiezza entrante e Wronskiano.

**Impatto:** né decomporre il QNM in ampiezza e fase, né annullare
un Wronskiano, né integrare una Riccati costituisce da solo una novità.
Il nostro passo di Newton su C_in è un controllo standard, non il
risultato da rivendicare.

### 3.3 Esistono già stimatori pratici dell'errore WKB

[Konoplya–Zhidenko–Zinhailo (2019)](https://arxiv.org/abs/1904.10333),
eq. (18), usano

\[
\Delta_k=\frac{|\omega_{k+1}-\omega_{k-1}|}{2}.
\]

Discutono anche dispersione fra approssimanti Padé e casi in cui
le stime falliscono, in particolare §VI.B.

**Impatto:** battere soltanto ε² o il criterio ℓ>n non è sufficiente.
Δ_k e i controlli Padé sono confronti obbligatori.
Δ_k è una stima assoluta: per confrontarla con il nostro errore relativo
va normalizzata coerentemente. Per k=3 servono almeno WKB2 e WKB4:
lo stato attuale WKB1/WKB3 non implementa ancora questo confronto.
Una differenza |ω3−ω1| non va etichettata come Δ3.

### 3.4 Residui pesati e risonanze non sono un territorio inesplorato

[Gopalakrishnan et al. (2024), §3](https://arxiv.org/html/2403.19485v2)
costruiscono stimatori dual-weighted residual per autovalori complessi
di modi ottici aperti, con problema aggiunto e PML.

[Araujo-Cabarcas–Engström (2016/2017)](https://arxiv.org/abs/1606.09635)
propongono un test basato su Lippmann–Schwinger e pseudospectro per
distinguere risonanze numeriche da soluzioni spurie.

**Impatto:** non possiamo rivendicare il primo stimatore a residuo per
risonanze aperte. L'applicazione WKB ai QNM gravitazionali resta distinta,
ma deve produrre un miglioramento specifico, non una semplice
ridenominazione del residuo aggiunto.

### 3.5 Sensibilità spettrale: un vincolo, non un dettaglio

[Jaramillo–Macedo–Al Sheikh (2021)](https://arxiv.org/abs/2004.06434)
studiano il pseudospectro QNM, usando Pöschl–Teller e Schwarzschild.
[Gasperin–Jaramillo (2021)](https://arxiv.org/abs/2107.12865)
discutono il ruolo del prodotto scalare nella misura delle perturbazioni.

**Impatto, per inferenza:** la sola piccolezza di una norma di residuo non
equivale a una garanzia sull'errore dell'autovalore. Occorre specificare
dominio dell'operatore, norma, sensibilità e classe di perturbazioni.
Un indicatore condizionato per D'(ω) è coerente con questa esigenza,
ma non diventa per questo una certificazione rigorosa.

### 3.6 Madelung, Weber e ramo Dirac

[Kumar (2026), §2–3](https://arxiv.org/html/2602.00507v1)
discute Bohm–Madelung, Ermakov–Pinney e funzioni di Weber in meccanica
quantistica stazionaria. È un precedente per il collegamento formale,
non una validazione del nostro problema a frequenza complessa.

[Meza-Domínguez–Matos (2026), §10](https://arxiv.org/html/2605.28887v1)
presentano una formulazione idrodinamica chirale con applicazione a QNM
fermionici in Schwarzschild. È una sovrapposizione diretta del ramo Dirac.
La pagina arXiv consultata mostra v1 e non una pubblicazione su rivista.
Le criticità annotate nel nostro progetto richiedono una verifica separata:
non autorizzano a ignorare l'esistenza del precedente.

### 3.7 Quantizzazione globale e Vaidya

[Miyachi et al. (2025), §II–IV](https://arxiv.org/html/2503.17245v2)
trattano condizioni QNM globali e strutture di Stokes nell'exact WKB.
I preprint [Hatsuda–Shiga (2026)](https://arxiv.org/abs/2605.01321)
e [Lo Chiatto et al. (2026)](https://arxiv.org/abs/2609.02816),
dei quali qui sono stati controllati abstract e metadati, confermano
l'attività su quantizzazione esatta in geometrie estreme/quasi estreme.
Questo non coincide automaticamente col nostro limite L→∞ a n fissato.

In Vaidya, [Abdalla–Chirenti–Saa (2006)](https://arxiv.org/abs/gr-qc/0609036)
studiano già il comportamento inerziale dei QNM.
[Capuano–Santoni–Barausse (2024), §II–III](https://arxiv.org/html/2407.06009v2)
forniscono il confronto a tasso di variazione della massa costante.
[Yoo et al., v3 del 2026, §4.2](https://arxiv.org/html/2510.25062v3)
mostrano l'importanza dello scattering durante la propagazione nel caso
dinamico.

**Impatto:** chiamare «memoria Madelung» una deviazione dalla frequenza
congelata non basta. Va separata da propagazione, scelta dell'osservatore,
mescolamento dei modi e variazione della massa.

## 4. Stato delle rivendicazioni

| Rivendicazione possibile | Giudizio attuale |
|---|---|
| nuova fisica introdotta dalla trasformazione Madelung | non sostenibile: è una riscrittura dell'equazione |
| prima ampiezza–fase o idrodinamica dei QNM | non sostenibile come rivendicazione ampia |
| primo legame tra residuo pesato e correzione QNM | precedenti sostanziali, §3.1 e §3.4 |
| indicatore E_M specifico e sottrazione uniforme | stessa prescrizione non identificata nelle fonti viste; utilità predittiva non dimostrata |
| formule chiuse e limiti del benchmark Pöschl–Teller | risultati locali dell'audit; da soli sono soprattutto analisi di consistenza |
| stimatore WKB economico con condizioni radiative e controllo dei nodi | candidato da costruire e confrontare, non risultato già ottenuto |
| diagnostica aggiuntiva della risposta dinamica in Vaidya | ipotesi aperta, senza simulazione di validazione nel progetto |

Non assegno percentuali di originalità o probabilità di accettazione:
non avrebbero una base empirica.

## 5. Che cosa manca per un articolo metodologico

Il seguente è un protocollo proposto per il progetto, non una regola imposta
da una rivista.

1. **Definizione operativa.** Lo stimatore deve usare potenziale,
   frequenza e profilo approssimati. La soluzione esatta/Leaver va riservata
   alla verifica, senza entrare nella calibrazione dei dati esclusi.
2. **Identificazione dell'errore.** Distinguere errore WKB, integrazione
   numerica, condizioni ai bordi e troncamento del dominio. Separare
   frequenza oscillatoria e smorzamento, oltre all'errore complesso totale.
3. **Matching e nodi.** Costruire entrambi i rami radiativi alla frequenza
   di prova; dichiarare come si attraversano gli zeri senza mascherarli
   con soglie arbitrarie. Dimostrare la stabilità rispetto alla regione
   di raccordo, finestra e precisione.
4. **Confronti competitivi.** Includere ε, n/ℓ, Δ_k, dispersione Padé,
   Newton/Jost quando disponibile e un controllo con residuo complesso
   non reinterpretato in variabili Madelung.
5. **Validazione indipendente.** Usare famiglie con parametri di forma
   indipendenti da ε. Escludere intere famiglie dalla calibrazione;
   mantenere i confronti a ε e overtone fissati. Una LOOCV punto per punto
   sulla stessa famiglia non risolve questo problema.
6. **Misure di utilità.** Riportare rapporto stima/errore vero,
   sottostime, scelta dell'ordine e costo complessivo, includendo
   costruzione del profilo e del problema ausiliario.
   La sola correlazione Pearson non basta.
7. **Valore aggiunto.** Mostrare almeno uno fra: migliore previsione
   a costo comparabile; stesso controllo a costo inferiore; dominio
   utile più ampio; limite teorico nuovo con conseguenza concreta.
8. **Riproducibilità.** Archiviare codice, parametri e tabelle del test
   finale, tenendo distinti risultati pilota e conferme successive.
   Una copia pubblica del codice non è una pubblicazione peer-reviewed.

Prova preliminare consigliata: deformare una barriera controllando
separatamente altezza, curvatura al massimo e forma più lontana.
Per perturbazioni distanti, seguire le radici per continuazione e dichiarare
quando cambia l'identificazione del modo. Il test è motivato anche dal
precedente [“Elephant and Flea”](https://arxiv.org/abs/2111.05415):
non va presentato come scoperta dell'esistenza di sensibilità a perturbazioni
lontane. Un diagnostico calcolato da un profilo globale esatto può già
contenere tale informazione: non va confuso con una funzione del solo
potenziale locale.

**Criterio di arresto:** se il vantaggio scompare togliendo la soluzione
esatta dagli input o confrontando a scala fissata, non rivendicare il
predittore. Se il risultato coincide con una formula perturbativa nota,
presentarlo come applicazione solo quando l'applicazione aggiunge utilità
misurabile.

## 6. Pubblicabilità e possibile formato

| Formato | Valutazione del materiale presente | Condizione per rivalutarlo |
|---|---|---|
| appunti tecnici/repository | utile già adesso, con limiti espliciti | nessuna pretesa di revisione scientifica esterna |
| breve nota metodologica critica | possibile ma ancora debole come contributo autonomo | risultato generale o esempio significativo per una difficoltà della letteratura, non solo correzione di una nostra ipotesi |
| articolo completo su uno stimatore WKB | non pronto | soddisfare il protocollo di §5 e chiarire il vantaggio sui precedenti |
| articolo dinamico Vaidya | prematuro | evoluzione temporale convergente e separazione delle spiegazioni alternative |

**Classical and Quantum Gravity** è una sede tematicamente pertinente.
Le sue [indicazioni ufficiali](https://publishingsupport.iopscience.iop.org/journals/classical-and-quantum-gravity/about-classical-quantum-gravity/)
richiedono avanzamenti originali per i Research Papers e prevedono Notes
per risultati brevi ma utili e nuovi. La presenza di questo formato
non rende automaticamente pubblicabile il nostro audit.

**Physical Review D** include metodi e risultati di gravitazione nel
proprio [ambito ufficiale](https://journals.aps.org/prd/about).
Lo considererei per un lavoro con validazione robusta su buchi neri e
un risultato di interesse fisico/metodologico chiaro. Questa è una valutazione
di pertinenza, non un pronostico editoriale. Non suggerisco ora una Letter.

Titolo di lavoro, solo se il metodo verrà validato:
“Boundary-aware error diagnostics for WKB black-hole quasinormal modes”.
Il ruolo di Madelung va spiegato nel contenuto; il titolo non deve
promettere una nuova teoria del fluido o una certificazione assente.

## 7. Ordine consigliato del lavoro

Prima chiarire l'equivalenza/differenza con perturbazione logaritmica e
residui aggiunti. Poi eseguire il test minimo di §5 su Pöschl–Teller
deformato e su Schwarzschild, con una seconda famiglia statica esclusa
dalla calibrazione. Solo se passa, investire nel benchmark esteso e
nella stesura. Il percorso statico → Kerr → Vaidya resta sensato,
ma non serve accumulare geometrie prima di validare il nucleo.

## 8. Documenti e tracciabilità

- [Audit matematico e numerico](spectral_audit_2026-09-08.md):
  formule, 18 nuovi controlli Wolfram, risultati numerici.
- [Consegna autosufficiente a Claude](../CLAUDE_HANDOFF.md):
  contesto, conclusioni da non perdere e domande per una revisione indipendente.
- [Mappa storica della letteratura](prior_art.md):
  da leggere insieme al presente aggiornamento.

In questa fase sono stati aggiunti documenti e riferimenti: **nessun nuovo
benchmark fisico è stato eseguito**. I test riportati nel file per Claude
sono quelli già completati e registrati dall'audit precedente.
Nessun invio a riviste, nessun messaggio a Claude, nessun commit o push.
