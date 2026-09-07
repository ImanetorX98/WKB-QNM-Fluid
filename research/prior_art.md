# Mappa della letteratura e confine di originalità

**Ricerca in corso — aggiornata al 7 settembre 2026.** La mancata
identificazione di un precedente non dimostra che il precedente non esista.
Questa mappa serve a formulare rivendicazioni controllabili e a guidare una
successiva ricerca sistematica per citazioni.

## Risultato provvisorio

Non è emerso un lavoro che combini tutti gli elementi seguenti:

- perturbazioni QNM di buchi neri relativistici, non soltanto analoghi;
- decomposizione Madelung della perturbazione stessa;
- trattamento coerente della corrente PDE e del problema spettrale aperto;
- confronto progressivo statico → stazionario non statico → non stazionario;
- termine dinamico normalizzato come diagnostico del ritardo QNM in Vaidya.

Sono però presenti precedenti strettamente adiacenti. La novità non può essere
attribuita alla trasformazione ampiezza–fase o alla sola comparsa di
\(-A''/A\).

## Precedenti più vicini

| Filone | Cosa è già stato fatto | Differenza rispetto al progetto |
|---|---|---|
| WKB di barriera | QNM di Schwarzschild e altri fondi a ordini elevati, Padé e matching uniforme | nessuna interpretazione Madelung necessaria |
| exact WKB / Riccati | momento quantistico, Stokes, periodi complessi e QNM di Schwarzschild/Kerr/RN | la chiusura locale e i periodi non sono nuovi come algebra |
| phase-amplitude QNM | Prüfer e funzioni di fase applicate a Schwarzschild e Kerr | usate come metodo numerico, non come diagnostica di memoria o simmetria |
| analogue gravity | QNM di metriche acustiche e pressione quantistica microscopica | il fluido è il mezzo fisico che genera la metrica analoga |
| correnti QNM | bilinear forms, regolarizzazione e slicing iperboloidale | chiariscono l'apertura, senza Madelung |
| QNM di Vaidya | inerzia, frequenza dipendente dal tempo/raggio, scattering causale | nessun termine di Madelung identificato come memoria |

### Phase-amplitude e risonanze

K. Glampedakis e N. Andersson, *Quick and dirty methods for studying black-hole
resonances*, applicano trasformazioni di Prüfer e metodi ampiezza–fase a QNM di
Schwarzschild e Kerr. È il precedente matematicamente più vicino e deve essere
citato quando si discute la polarizzazione della ODE:
[arXiv:gr-qc/0304030](https://arxiv.org/abs/gr-qc/0304030).

Un uso moderno della stessa famiglia di tecniche per frequenze ed excitation
factors oltre GR è H. O. Silva et al., *Quasinormal modes and their excitation
beyond general relativity*:
[arXiv:2404.11110](https://arxiv.org/abs/2404.11110).

### Exact WKB: aggiornamento decisivo

La linea exact-WKB è arrivata direttamente ai QNM di buchi neri e restringe
molto qualunque rivendicazione basata soltanto su ricorsioni all-order, curve di
Stokes o periodi quantistici:

- Miyachi, Namba, Omiya e Oshita trattano con exact WKB le condizioni QNM da
  orizzonte a infinito, incluse spirali logaritmiche delle curve di Stokes e
  branch cut, verificando Schwarzschild:
  [arXiv:2503.17245](https://arxiv.org/abs/2503.17245).
- Hatsuda e Shiga calcolano periodi quantistici/simboli di Voros ad alto ordine
  e li risommano con Borel–Padé per Kerr e Reissner–Nordström estremi:
  [arXiv:2605.01321](https://arxiv.org/abs/2605.01321).
- Lo Chiatto, Schenk, Wagner e Yu costruiscono geometria di Stokes e condizione
  di quantizzazione esatta per i modi zero-damped di Reissner–Nordström quasi
  estremo:
  [arXiv:2609.02816](https://arxiv.org/abs/2609.02816).

Ne segue che «periodi di Voros dei QNM» non è più una direzione originale di
per sé. Il possibile contributo del progetto è collegare quegli oggetti globali
alla decomposizione Madelung e a diagnostici di flusso o non adiabaticità.

### Buchi neri analoghi

R. G. Daghigh e M. D. Green mostrano che, per QNM ad alto overtone di un buco
nero acustico in un condensato, il termine di quantum potential della
Gross–Pitaevskii linearizzata domina la struttura microscopica:
[arXiv:1411.7066](https://arxiv.org/abs/1411.7066).

L. C. N. Santos et al. costruiscono una metrica acustica tramite Madelung e
calcolano poi QNM con WKB. È una combinazione esplicita delle parole chiave, ma
la trasformazione descrive il condensato di fondo, non il profilo QNM di un
buco nero relativistico:
[EPJC 85, 1036 (2025)](https://link.springer.com/article/10.1140/epjc/s10052-025-14789-4).

### Stati risonanti e regolarità

I QNM divergono sulle usuali sezioni che raggiungono la biforcazione e
l'infinito spaziale, mentre risultano regolari verso l'orizzonte futuro e
l'infinito nullo in coordinate/slicing appropriati. Questo rende pericolosa
l'interpretazione globale di \(A^2\) come densità:
[Macedo e Zenginoglu, arXiv:2409.11478](https://arxiv.org/abs/2409.11478).

Per Kerr, correnti bilineari conservate e ortogonalità richiedono strutture e
regolarizzazioni specifiche:
[Green et al., arXiv:2210.15935](https://arxiv.org/abs/2210.15935).

### Vaidya

- Abdalla, Chirenti e Saa identificano un comportamento di «inerzia» dei QNM
  durante variazioni rapide della massa:
  [arXiv:gr-qc/0609036](https://arxiv.org/abs/gr-qc/0609036).
- Lin, Sun e Zhang trovano frequenze dipendenti da tempo e posizione per la
  propagazione causale fra orizzonte e infinito:
  [arXiv:2104.06631](https://arxiv.org/abs/2104.06631).
- Capuano, Santoni e Barausse trattano \(\dot M\) costante in frequenza; è il
  benchmark naturale prima del caso pienamente dinamico:
  [arXiv:2407.06009](https://arxiv.org/abs/2407.06009).
- Yoo, Kimura, Ishibashi e Ohashi confrontano photon sphere/Penrose limit e
  waveform temporale, mostrando il ruolo della propagazione e dello
  scattering:
  [arXiv:2510.25062](https://arxiv.org/abs/2510.25062).

Un precedente concettualmente vicino alla variazione di massa, ma non formulato
su Vaidya né in termini di Madelung, è l'«absorption-induced mode excitation»:
l'assorbimento cambia la massa del buco nero e proietta il modo iniziale sul
nuovo spettro. Nel regime studiato l'effetto è non adiabatico ed è approssimato
come salto improvviso:
[Sberna et al., arXiv:2112.11168](https://arxiv.org/abs/2112.11168).
Questo lavoro deve essere usato come confronto obbligatorio per distinguere
«memoria del profilo» da semplice mode mixing dovuto al cambiamento di massa.

Le ricerche testuali mirate per combinazioni di *Vaidya*, *Madelung*, *Bohm
potential* e *quasinormal modes* non hanno restituito un precedente diretto.

### Barriera esatta di Pöschl–Teller

Lo spettro e le autofunzioni QNM della barriera di Pöschl–Teller sono noti e
non costituiscono una rivendicazione del progetto. Cardona e Molina ottengono
la torre QNM di potenziali di Pöschl–Teller generalizzati con metodi algebrici:
[arXiv:1711.00479](https://arxiv.org/abs/1711.00479). Beyer studia completezza
e convergenza dell'espansione QNM:
[arXiv:gr-qc/9803034](https://arxiv.org/abs/gr-qc/9803034).

La ricerca mirata non ha identificato in questi lavori, né nei risultati
testuali adiacenti, l'uso di \(-A''/A\) come diagnostico quantitativo
dell'errore WKB o l'analisi del suo fallimento sugli zeri QNM. Questa resta una
conclusione provvisoria di ricerca, non una prova di priorità.

## Rivendicazioni da evitare

- «Prima decomposizione ampiezza–fase di un QNM».
- «Prima relazione fra QNM e fluidi».
- «Il potenziale di Madelung modifica lo spettro»: una trasformazione esatta di
  variabili non introduce nuova fisica.
- «\(A^2\) è una probabilità QNM» senza una costruzione funzionale e di bordo.
- «Il termine misto di Vaidya è covariante» se viene definito soltanto come
  \(A_{vr}/A\) in una coordinata particolare.

## Rivendicazioni plausibili da verificare

1. Un unico quadro Madelung distingue in modo controllato PDE conservativa e
   riduzione spettrale aperta dei QNM.
2. Un funzionale Madelung localizzato sulla barriera predice il dominio di
   affidabilità della WKB meglio del solo criterio \(\ell>n\).
3. La scala statico → stazionario → dinamico separa rispettivamente curvatura,
   shift rotazionale e risposta temporale.
4. Una curvatura logaritmica normalizzata rispetto al fondo congelato misura
   una componente del ritardo QNM di Vaidya non contenuta nella sola frequenza
   istantanea.
5. Un diagnostico globale deve trattare esplicitamente gli zeri del modo e la
   regione uniforme dei turning point; la serie locale di Madelung non basta.

## Prossime ricerche bibliografiche

- seguire tutte le citazioni del lavoro phase-amplitude del 2003;
- cercare *Milne/Ermakov equation* insieme a black-hole resonances;
- cercare letteratura su Gamow/Siegert states e Madelung complessa;
- seguire citazioni e sviluppi dei tre lavori exact-WKB 2025–2026;
- seguire citazioni e lavori correlati ai quattro articoli Vaidya sopra;
- distinguere sempre buchi neri astrofisici, analog gravity e modelli
  idrodinamici olografici.
