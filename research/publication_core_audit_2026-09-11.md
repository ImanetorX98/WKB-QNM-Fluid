# Nucleo pubblicabile: audit del nuovo handoff e priorità

11 settembre 2026. Letti CODEX_HANDOFF_CONDIZIONAMENTO.md (prioritario),
CODEX_HANDOFF.md, il programma di condizionamento e le sezioni rilevanti del
manoscritto. Branch analisi-kerr-vaidya-2026-09, worktree pulito all'inizio.
Il manoscritto non è stato riscritto in questo intervento.

## 1. Verdetto

Non proporrei l'invio del manoscritto attuale. Esiste materiale per un lavoro
metodologico, ma il nucleo va delimitato e verificato: non abbiamo ancora una
garanzia di originalità o pubblicabilità. La conclusione del §5 secondo cui
nessun funzionale di Q_M possa predire omega non segue dai suoi dati.

Il nuovo handoff ha ragione a contestarla, ma 'amplificazione universale 2'
è a sua volta troppo forte. Il codice misura una sola direzione di variazione
di una frequenza complessa. Un controllo aggiuntivo distingue chiaramente
sensibilità radiale, sensibilità alle singole componenti e condizionamento
rispetto alla norma complessa.

## 2. Risultato analitico: 2 e 10^2 possono coesistere

Fissati L, potenziale e convenzione per r*, sia omega=a+ib, b non nullo.
All'infinito il ramo uscente ha u -> omega/L, le derivate spaziali di u
svaniscono, e dalla formula di ln A del handoff segue esattamente nel limite:

~~~text
Q_infinity = -b^2/L^2,
dQ/Q = 2 db/b.
~~~

Q è reale e dipende da Im omega: non è olomorfa in omega. Una scrittura
indifferenziata 'partial ln Q / partial ln omega' non definisce quindi
un'unica derivata complessa. Per le sensibilità relative usare ln|Q|.

Tre domande diverse:

| Perturbazione e misura dell'errore | Sensibilità asintotica |
|---|---:|
| omega -> (1+delta)omega, delta reale | 2 |
| b -> (1+delta)b, errore relativo rispetto a b | 2 |
| errore complesso arbitrario, misurato come abs(domega)/abs(omega) | massimo 2 abs(omega)/abs(b) |

La terza quantità è il condizionamento relativo normwise della funzione
scalare puntuale, con norma euclidea sull'ingresso (a,b):

~~~text
kappa = |omega| sqrt[(partial_a Q)^2+(partial_b Q)^2] / |Q|.
~~~

Per ell=70,s=2: omega=13.564180...-0.096217...i, quindi kappa_infinity=281.9565,
pur avendo sensibilità radiale esattamente 2. Usando la formula eikonale
Schwarzschild, a ~ L/(3 sqrt(3)), |b| ~ (n+1/2)/(3 sqrt(3)), si ottiene
kappa_infinity ~ 4L/(2n+1) a overtone fissato e L grande.

È una conseguenza elementare della dipendenza quadratica dal damping:
**non va promossa da sola a scoperta originale**. Non conferma i numeri del
vecchio §5 e non recupera la sua conclusione universale. Mostra quale problema
matematico deve essere definito prima di confrontare i numeri.

Le correzioni a r finito non sono automaticamente O(epsilon): contano anche
il potenziale locale e la distanza dalla barriera. I limiti r -> infinito
e L -> infinito vanno distinti.

## 3. Riproduzione e nuovo controllo numerico

Comando originale eseguito:

~~~sh
python3.13 calculations/madelung_conditioning_schwarzschild.py
~~~

Riprodotti i pavimenti 1.735e-4,4.418e-5,1.457e-5,7.092e-6 per ell=20,40,70,100,
e i fattori circa 2 sopra il pavimento. Però (errore(delta)-errore(0))/delta
non è una derivata del funzionale: sottrae due mediane di errori assoluti.
Il programma produce anche valori negativi per piccole delta, per cancellazione
con l'errore WKB di base. Non usare questo estimatore come numero di condizione.

Nuovo comando:

~~~sh
python3.13 calculations/conditioning_directional_audit.py
~~~

Calcola differenze centrate della **stessa mappa WKB**, separatamente rispetto
ad a e b; non sottrae il pavimento del riferimento esatto. La griglia è uniforme
in r*, r è ottenuto invertendo la coordinata con Lambert W. Finestra 20–50,
s=2, 40001 punti, incremento assoluto h=1e-5 abs(omega).

| ell | Direzione radiale | Componente b, relativa a b | Mediana kappa puntuale normwise | Limite asintotico |
|---:|---:|---:|---:|---:|
| 20 | 1.9821 | 2.0070 | 82.1371 | 81.8502 |
| 40 | 1.9824 | 2.0070 | 162.4902 | 161.9243 |
| 70 | 1.9828 | 2.0070 | 282.9420 | 281.9565 |
| 100 | 1.9829 | 2.0070 | 403.3742 | 401.9695 |

Questa è la mediana dei massimi **puntuali**, non il numero di condizione
di un funzionale integrale o di una norma globale del profilo.
Il valore 282.942 resta stabile a circa 2e-6 relativi variando punti
20001–80001 e h/abs(omega)=1e-4–1e-5. La componente reale, più piccola,
è meno stabile per la differenziazione numerica ripetuta: non rivendicarne
molte cifre. Questi test misurano la mappa WKB troncata, non una famiglia
di soluzioni esatte fuori frequenza QNM con entrambe le condizioni al bordo.

Due nuovi test superati: limite a potenziale nullo, con Q esatto; confronto
Schwarzschild direzionale e di risoluzione. Nessuna nuova verifica Mathematica
in questo intervento e nessuna riesecuzione dell'intera suite del progetto.

## 4. Stato dei quattro compiti del handoff

1. **(a), fattore 2:** chiarito analiticamente all'infinito e confermato
   numericamente. Va specificata la direzione: non è universale.
2. **(b), picco:** un primo test formale 2.5<r<4, ell=70 dà sensibilità radiale
   circa 196, contro 1.98 fuori. Ma qui la WKB non uniforme è in difficoltà:
   questo numero NON certifica sensibilità fisica. Serve confronto con soluzione
   esatta fuori shell ben definita e/o approssimazione uniforme, convergenza,
   scelta dei rami e trattamento degli zeri di Q. Non salva il vecchio §5.
3. **(c), Kerr:** non ricontrollato in questo turno. Priorità successiva, usando
   variazioni separate di Re omega e Im omega, e ricalcolando coerentemente
   l'autovalore sferoidale per ciascuna frequenza perturbata.
4. **(d), Figura 2:** ancora sospesa. La figura utile deve distinguere direzioni,
   norma dell'errore e validità della WKB, non raffigurare un solo fattore.

## 5. Altri punti che un referee contesterebbe

- 'Proporzionalità esatta' in §6 non è dimostrata da una dispersione 0.15%
  su un campione finito. Scrivere 'proporzionalità approssimata nel campione'.
  Handoff e manoscritto riportano anche coefficienti/dispersioni diversi:
  fissare una versione dei dati. Distinguere correlazione, comune scaling
  epsilon^2 e assenza di informazione addizionale fuori campione.
- Un cambio esatto di variabili conserva l'informazione ma può migliorare un
  algoritmo. La dipendenza della soluzione dalla frequenza non impedisce di
  cercare gli zeri di un residuo spettrale. Anche un alto condizionamento diretto
  non è un teorema di impossibilità per tutti i problemi inversi/funzionali.
- Nei §§11.3–11.4 compare ancora N non nullo => Z1 non nullo. Una proiezione
  risonante non nulla non dimostra che la componente complementare del profilo
  sia non nulla: può essere assorbita in trasporto di ampiezza/fase e correzione
  della frequenza, secondo le convenzioni. Serve separare queste componenti.
- 'Autoaggiunto' va distinto da formalmente simmetrico bilineare. Il dominio,
  i bordi e la regolarizzazione fanno parte dell'argomento. Il prodotto
  sesquilineare non impedisce in generale la solvibilità: richiede il modo
  aggiunto appropriato, invece della proiezione semplificata.
- 7.8e24 contro un risultato di ordine 10 implica una sottrazione di circa
  **24 cifre**, non 14; la doppia precisione non può conservare il risultato.

## 6. Nucleo candidato e condizioni per chiamarlo pubblicabile

Proposta: **validità e sensibilità direzionale dei diagnostici di ampiezza
per QNM, con benchmark controllati e confronto con WKB standard**.
Il contributo utile dovrebbe essere un criterio quantitativo che dica quando
il diagnostico è affidabile e se offre vantaggi rispetto a controlli economici,
non l'affermazione che Madelung introduce nuova fisica.

Minimo necessario prima di scegliere questa tesi centrale:

1. Enunciato preciso della mappa, normalizzazione degli errori, direzioni,
   regime di validità e distinzione tra profilo, funzionale e spettro.
2. Convergenza del Jacobiano con derivate più stabili; verifica su un benchmark
   esatto e su Schwarzschild; test uniforme vicino al picco.
3. Per E_M, derivazione dell'eventuale coefficiente dominante rispetto a
   Lambda3, oppure test fuori campione su deformazioni controllate della
   barriera a L fissato. Non basta ripetere la stessa legge in L.
4. Un vantaggio concreto: soglia d'errore validata, criterio di esclusione
   degli artefatti o previsione che batta un controllo esistente.
5. Solo allora aggiornare abstract, §§5–6, discussione e figura, con ipotesi
   limitate ai risultati verificati. Non forzare Kerr/Vaidya dentro il nucleo
   se servono solo ad ampliare la casistica.

Vaidya è una seconda opzione, più ambiziosa: collegare la solvibilità a una
quantità indipendente dalla normalizzazione e verificarla con evoluzione nel
tempo. Un numeratore finito da solo non basta. Non è necessario risolvere
tutti gli overtoni prima di stabilire se il fondamentale predice qualcosa.

## 7. Controllo bibliografico mirato, non esaustivo

Fonti primarie consultate oggi (pagine abstract; nessuna pretesa di lettura
integrale o di ricerca di priorità completa):

- [Glampedakis–Andersson 2003](https://arxiv.org/abs/gr-qc/0304030): metodi
  fase-ampiezza usati concretamente per trovare risonanze, con test Schwarzschild
  e Kerr. Precedente che impedisce di liquidare l'utilità algoritmica di tutte
  le riformulazioni fase-ampiezza; non è necessariamente lo stesso Q_M.
- [Capuano–Santoni–Barausse 2024](https://arxiv.org/abs/2407.06009): QNM e
  risposta mareale su Vaidya con tasso di variazione della massa costante.
- [Yoo et al. 2025/26](https://arxiv.org/abs/2510.25062): confronto tra
  descrizione tramite sfera fotonica dinamica/limite di Penrose e forme d'onda
  numeriche su Vaidya. È un confronto da affrontare prima di rivendicare un
  nuovo effetto non adiabatico o una nuova previsione di ringdown.

Conclusione editoriale: il lavoro può essere rifocalizzato, ma nessuna delle
semplici identità qui derivate prova da sola originalità. La decisione di
investire sul nucleo va basata sull'esito dei test discriminanti sopra.
