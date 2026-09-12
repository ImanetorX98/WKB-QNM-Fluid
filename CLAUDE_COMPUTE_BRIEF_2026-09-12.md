# Compito delimitato per Claude Code — proposta, non ancora assegnata

## Aggiornamento prioritario: protocollo delle prove analitiche

Leggere PRIMA `research/analytic_core_response_2026-09-12.md`. Questa sezione
precisa e prevale sui punti generici sotto. Non avviare grandi sweep prima
di chiudere A e B. Il compito è prescritto per una futura sessione Claude;
Codex non l'ha inviato automaticamente a un'altra applicazione.

### A. Secondo problema spettrale e norma generalizzata

- Riprodurre la risonanza circa3.423295488-0.607387592i, senza chiamarla
  fondamentale finché non si è classificato lo spettro.
- Metodo raccomandato: matching di soluzioni uscenti dai DUE lati, con
  Wronskiano normalizzato in un punto interno, oppure collocazione con
  Robin dipendenti da omega. Cambiare linguaggio/CAS allo stesso shooting
  non basta come indipendenza di metodo.
- Usare eta=0,+/-0.001,+/-0.01; tolleranze/gradi su tre livelli. Primo gate:
  accordo delle frequenze fra metodi <1e-8 assoluto e variazioni risolte
  almeno10 volte sopra l'errore stimato. È un criterio richiesto, non un
  risultato già ottenuto.
- Verificare d omega/d eta=B/N senza coniugazione, con termini superficiali.
  Traslare i bordi in regioni dove V=0, mantenendo la stessa normalizzazione
  fisica quando si confrontano B e N; il rapporto è indipendente dal fattore
  costante. Confrontare passi eta1e-2,3e-3,1e-3,3e-4,1e-4 finché domina
  il rumore. Target relativo1e-5 o spiegazione quantitativa del limite.

### B. Chiudere PRIMA il problema di quadratura della derivata

- Il nuovo script trova per dE/deta circa0.00189 nel continuo non convergente,
  contro0.001873624265 nel funzionale DISCRETO originale. Non forzare
  coincidenza: sono livelli di approssimazione diversi.
- Localizzare tutti gli zeri di Q_M in J mediante bracketing raffinato;
  verificare stabilità del numero di zeri. Spezzare le quadrature a tali
  punti, usando Gauss/adattiva o alta precisione in ciascun tratto.
- Verificare assenza di zeri di psi e registrare min|psi| nella normalizzazione
  adottata. Distinguere zeri di psi da zeri di Q_M.
- Calcolare k con due metodi: quadratura della formula integrale e ODE
  k'+2zk=-2omega, k(a)=-i. Controllare il residuo e la convergenza.
- Confrontare la derivata analitica continua con differenze centrali di E
  calcolata a sua volta con quadratura convergente. Tre raffinamenti;
  target relativo1e-4 iniziale,1e-5 se raggiungibile. Documentare ogni plateau.
- Ripetere il test sulla stessa trapezoidale2001 punti come regressione,
  senza scambiarlo per la verifica del limite continuo. Non aumentare solo
  WorkingPrecision lasciando invariata una quadratura sottorisolta.

### C. Test della fattorizzazione e del rango (solo dopo A e B)

- Famiglia b_{c,w}=B((x-c)/w), con c=1.6,2.0,2.4 e w=0.2,0.4;
  tutte esterne a J=[-0.8,0.8] e alla barriera centrale [-1,1]. Allargare
  il dominio per includere il supporto fino a2.8. Non imporre il bordo
  libero dentro una bump! Controllare massimo dominante invariato.
- Per ogni perturbazione stimare d=omega_eta e h=partial_eta z con
  differenze centrali e confrontare h(x)=d k(x). Usare profili normalizzati
  e una metrica relativa L2, più un errore assoluto dove la risposta è piccola.
- Costruire colonne di risposte REALI dot Q_M campionate nello stesso J.
  Predire ogni colonna con le due basi della §3, senza rifittare d.
  La terza singolare deve essere compatibile con l'errore numerico e
  di differenziazione stimato, non semplicemente sotto una soglia arbitraria.
  Nessuna richiesta che la seconda singolare sia non nulla.
- Controllo negativo: una perturbazione INTERNA, ad esempio centro0.3,
  larghezza0.1, cambia l'ODE interna: non pretendere h=d k. In questo test
  non rivendicare getto o massimo invarianti; è solo un controllo delle ipotesi.

### D. Utilità predittiva, distinta dalla prova strutturale

- Calcolare separatamente rapporto integrale e rapporto delle mediane pesate,
  dichiarando la definizione esatta di quantile, peso e denominatore.
  La formula analitica di dot E NON vale automaticamente per la mediana.
- Calcolare Lambda3 e frequenza WKB3 dal getto centrale. Valutare errore
  spettrale indipendente, non solo variazione di E. Nessuna calibrazione e
  valutazione sullo stesso campione; dichiarare training/test prima dei fit.
- Se E richiede omega esatta, chiamarla diagnosi a posteriori; una pretesa
  predittiva richiede un esperimento separato con soli input approssimati.

### Consegna e condizioni di arresto

Un unico `research/claude_compact_validation.md`: metodi, risultati dei gate
A–D, intervalli di stabilità, costi di esecuzione, fallimenti e comandi.
Salvare dati numerici separati. Se A/B falliscono, fermare C/D e riportare
il primo limite riproducibile; non selezionare soltanto i casi riusciti.
Non estendere automaticamente a Kerr, Vaidya o a pretese di originalità.

## Obiettivo

### Controllo F: residuo, norma e matching

Leggere `research/matching_invariance_proof_2026-09-12.md`.
Con normalizzazione psi(a)=1 e residuo non normalizzato
F=psi'(b)-D_plus psi(b), verificare N=-psi(b)F_omega al QNM.
Calcolare F_omega sia tramite equazione variazionale sia differenze centrali
in omega lungo direzioni reale/immaginaria. Usare almeno tre passi e
riportare la normalizzazione: un Wronskiano normalizzato ha fattori diversi.
Confrontare la stima locale delta omega=-delta F/F_omega con l'effetto
di un raffinamento numerico; non usare il solo residuo come errore spettrale.
Per lo spostamento dei bordi distinguere invarianza del rapporto B/N da
quella di B,N separati, che richiede la stessa normalizzazione globale.

### Estensione E, dopo i gate A–B: derivata spettrale Kerr

Leggere `research/kerr_nonlinear_spectral_proof_2026-09-12.md`.
Non avviare nuovi grandi sweep né modificare i solutori storici.

1. Per ell,m=(28,19),(61,41) e c/L=0.2,0.4,0.3+0.1i, confrontare A_c
   dalla formula bilineare con differenze centrali lungo direzione reale
   e immaginaria del piano c. Seguire il medesimo autovalore e normalizzare
   coerentemente; nei casi complessi NON sostituire S^2 con |S|^2.
2. Usare tre passi relativi, indicativamente1e-3,1e-4,1e-5 di max(1,|c|),
   e variazione della dimensione della base. Registrare il denominatore
   bilineare normalizzato, non soltanto l'autovalore. Target relativo1e-6
   dove la derivata non è piccola; altrimenti usare errore assoluto.
3. Per a=0,0.6 e frequenze complesse dei benchmark, controllare q_omega
   della nota a tre raggi esterni all'orizzonte mediante differenziazione
   dell'intero potenziale, ricalcolando A a ogni omega. Confrontare anche
   la variante con A congelato, etichettandola come controllo negativo.
4. Nel caso a=0 recuperare esattamente il limite2omega, entro errore
   numerico. Non chiamare validata la norma globale Kerr finché non sono
   inclusi i D_plus/minus corretti ai bordi finiti e il loro error budget.

Salvare questa estensione in `research/claude_kerr_derivative_validation.md`
e nuovi script `claude_kerr_derivative_*`, senza interferire con file esistenti.

Verificare indipendentemente il test on-shell in
`calculations/compact_barrier_core_test.py` e il resoconto
`research/core_proof_scope_2026-09-12.md`. Non ricominciare la ricerca
bibliografica e non riscrivere il paper. Non dichiarare originale la formula
perturbativa: proviene da Leung et al., §2, già presente nei papers locali.

## Lavoro richiesto

1. Secondo metodo per la stessa risonanza (collocazione spettrale, matching
   da entrambi i lati o Mathematica), evitando di copiare semplicemente lo
   shooting e chiamarlo indipendente. Confrontare eta=0,+/-0.001,+/-0.01.
2. Controllare convergenza in precisione, passi/grado, quadratura della norma
   generalizzata e griglia del diagnostico. Usare almeno due raffinamenti.
3. Aggiungere la mediana pesata, mantenendola separata dal rapporto integrale.
   Ripetere il confronto su tre posizioni e due larghezze di bump; lasciare
   invariato un intorno del massimo centrale e seguire il medesimo ramo.
4. Calcolare Lambda3 dal getto centrale e verificare che rimanga identico;
   calcolare l'errore della frequenza WKB3 rispetto a quella indipendente.
   Una risposta del diagnostico non equivale a buona predizione dell'errore.
5. Restituire un unico Markdown con formule/convenzioni, comando riproducibile,
   tabella compatta, error budget e casi falliti. Salvare i dati completi,
   ma non riversare log estesi nella conversazione.

## Vincoli di coordinamento

Non modificare i file numerici già in uso: scrivere nuovi file con prefisso
`claude_compact_` e un resoconto `research/claude_compact_validation.md`.
Non toccare manoscritto, handoff precedenti o file non pertinenti. Non fare
commit/push o installazioni senza richiesta. Se serve un pacchetto mancante,
segnalarlo. Controllare prima cosa è già disponibile, incluso WolframKernel.

## Perché questo può evitare sprechi

Un incarico circoscritto con dati e criteri di accettazione evita di trasmettere
tutta la cronologia o duplicare analisi. Non è una stima di token né una
garanzia economica: dipende da piano, modello e iterazioni. Moltiplicare i
campioni dentro uno script non implica chiedere un nuovo ragionamento al
modello per ogni campione. Il risultato utile per Codex è un breve resoconto
con evidenze verificabili, non l'intera trascrizione della sessione.
