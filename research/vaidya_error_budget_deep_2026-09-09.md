# Analisi approfondita dello scarto 1.65e-4 di Claude

Data: 9 settembre 2026. Seguito di vaidya_quadrature_audit_2026-09-09.md.

## Esito e limite della diagnosi

Abbiamo riprodotto **il meccanismo del paradosso**, non il preciso esperimento
storico: un integrando verificato, una primitiva verificata e un errore del
trapezio inferiore a 1e-7 possono convivere con uno scarto integrale di ordine
1e-4, se gli estremi effettivi della quadratura differiscono da quelli della
primitiva. Questo è ora il candidato prioritario, non una causa accertata del
numero 1.65e-4. Manca il driver che lo produsse (non compare fra i file di
solvibilità/asintotica tracciati consultati, né fra i programmi elencati nel
workspace). Non abbiamo cercato log privati di altre applicazioni.

Separatamente è stato diagnosticato e rimosso **nel driver di audit** il
plateau di circa 1e-7: nasce prevalentemente dai dati iniziali all'orizzonte
limitati all'onda piana. Con dati di Frobenius il confronto esterno arriva
a circa 1e-11. Non abbiamo sostituito automaticamente i dati iniziali del
modulo di produzione: l'esperimento resta isolato e confrontabile.

## 1. Esperimento sugli estremi: il trapezio non spiega tutto

Caso ell=2, n=0, s=0, M=1; omega=0.483643872210713-0.09675877597828786i.
Integrazione in r* da r=2.0001 a r=70; selezione 40 <= r <= 60.
Usiamo la medesima funzione asintotica p=W exp(2i omega r*) per tutti i
confronti di questa tabella: **nessun errore del profilo radiale è necessario
per produrre il difetto**. Gli errori sono |rapporto complesso - 1|.

| Campioni globali | Trapezio / primitiva, estremi identici | Simpson / primitiva, estremi identici | Simpson / primitiva, estremi nominali |
|---:|---:|---:|---:|
| 24000 | 1.267e-6 | 1.095e-11 | 2.217e-3 |
| 100000 | 7.299e-8 | 8.307e-15 | 7.941e-4 |
| 200000 | 1.825e-8 | 7.097e-15 | 2.739e-4 |

Con 100000 punti, la maschera integra realmente da 40.0002440801 a
59.9992260194. Confrontarla con F(60)-F(40) introduce

~~~text
J_grid / [F(60)-F(40)] - 1
 = -1.8038492925e-4 - 7.7334952436e-4 i.
~~~

Lo spostamento dei bordi predice questo risultato al primo ordine:

~~~text
errore ~ [p(b) delta_b - p(a) delta_a] / [F(b)-F(a)].
~~~

Il resto complesso della previsione lineare è 3.13e-7. Il residuo differenziale
della primitiva, max |(F'-p)/p| sui campioni, è circa 3.2e-15.
Quindi il grosso difetto non nasce dall'algebra dell'antiderivata.

Per r=60, F'/F = 0.2168472933 + 0.9972690554i. Quando domina il bordo superiore,
uno spostamento di soli **1.617e-4 in r** basta a produrre uno scarto relativo
1.65e-4 al primo ordine. È una stima di sensibilità, NON un fit o la
ricostruzione degli estremi usati da Claude.

L'errore dei bordi non necessariamente migliora monotonamente raffinando:
nella finestra 40–50 passa da 3.80e-5 a 3.60e-4 passando da 100000 a 200000
punti, perché cambia l'allineamento. Una traslazione sistematica del bordo può
anche dare un errore relativo quasi costante nella regione esponenziale;
la selezione su griglia, invece, può mostrare oscillazioni a dente di sega.

## 2. Un vincolo che restringe le cause possibili

Se su **tutta** la finestra |n/p-1| <= epsilon, allora

~~~text
|integrale(n-p)| / |integrale(p)| <= epsilon * kappa,
kappa = integrale(|p|) / |integrale(p)|.
~~~

È la disuguaglianza triangolare, non un'ipotesi fisica. Vale anche esattamente
per le somme di trapezi a pesi positivi sugli stessi nodi. Nel test 40–60:

- kappa = 4.6403;
- massimo errore puntuale misurato = 3.286e-7;
- limite discreto = 1.525e-6;
- errore integrato del solo profilo misurato = 9.804e-8.

Perciò un errore puntuale uniformemente 2e-7 NON può trasformarsi da solo in
1.65e-4 su questa finestra (richiederebbe kappa almeno 825).
Questa esclusione riguarda il confronto degli integrali esterni: **non** la
parte finita dopo sottrazione, che può essere molto peggio condizionata.
Il controllo di Claude '0.9999998' va ripetuto come massimo su tutti i nodi,
non su un singolo punto o soltanto sull'estremo esterno.

## 3. Il plateau 1e-7: dati all'orizzonte, non quadratura

La routine originaria impone R=exp(-i omega r*) e dR/dr*=-i omega R a
r=2+delta. Questo è solo il termine iniziale dello sviluppo entrante.
Scrivendo R=exp(-i omega r*)h(r), l'equazione dà

~~~text
f h'' + (f' - 2i omega) h' - U h = 0,
h = 1 + h1 delta + h2 delta^2 + ...,
h1 = (7/4)/(1/2 - 2i omega),
h2 = [(9/4)h1 - 15/8]/(2 - 4i omega)
~~~

per ell=2, s=0. Entrambi i coefficienti sono stati verificati simbolicamente
con il kernel Mathematica locale. Nel driver indipendente inizializziamo sia
R sia la derivata coerentemente, mantenendo la stessa equazione e tolleranza.

Errore del confronto integrale esterno, 24000 punti, Simpson, estremi identici:

| delta | Solo onda piana | Frobenius fino a delta | Frobenius fino a delta^2 |
|---:|---:|---:|---:|
| 1e-3 | 4.153e-7 | 7.007e-10 | 2.418e-11 |
| 1e-4 | 9.805e-8 | 1.873e-11 | 1.677e-11 |
| 1e-5 | 2.396e-8 | 1.186e-10 | 1.017e-10 |
| 1e-6 | 5.810e-9 | 5.204e-11 | 2.544e-10 |

Avvicinarsi indefinitamente all'orizzonte non è la strategia migliore: con
Frobenius il caso delta=1e-4 è già più accurato di delta=1e-6 in doppia precisione.
Con 100000 punti e ordine 2, variando rtol=2e-9,2e-11,2e-13 si ottiene
3.688e-9,2.281e-11,1.462e-11. Non attribuiamo a questi ultimi valori più cifre
affidabili di quanto consentano le cancellazioni numeriche.

Indizio analitico sulla convergenza lenta dell'onda piana: il ramo entrante
locale è delta^(-2i omega), quello uscente delta^(+2i omega). Trascurare il
termine lineare del ramo entrante induce, al primo ordine locale, una
contaminazione relativa dei coefficienti proporzionale a delta^(1-4i omega),
di modulo delta^(1+4 Im omega) = delta^0.612965. Ciò è compatibile con la
lenta riduzione misurata. Non è una stima rigorosa dell'errore globale.

## 4. Ipotesi escluse o distinguibili

- **Coordinate r* incoerenti:** con i parametri principali, max |r*_ODE-
  r*(r_ODE)| è 2.03e-11, insufficiente a spiegare 1e-4. Il controllo è nel codice.
- **Trapezio puro:** per una singola esponenziale su griglia uniforme,
  T/I=(kh/2)coth(kh/2)=1+(kh)^2/12+... . Per k=2i omega servirebbe h circa
  0.0451 per uno scarto 1.65e-4; è incompatibile con l'errore 1e-7 misurato
  sulla griglia fine di questo audit. È una stima asintotica, non la formula
  esatta per W variabile o una griglia non uniforme.
- **Normalizzazione incoerente:** usare C diversi nell'integrando e nella
  primitiva produce circa 2 delta_C/C. La convergenza di C a molte cifre non
  dimostra che lo stesso C, con la stessa fase, sia usato in entrambi.
- **Parte finita:** anche i residui 1e-11 qui ottenuti possono venire amplificati
  sottraendo quantità enormi. Non abbiamo certificato il numeratore globale,
  la norma globale o uno shift fisico.

Per distinguere le cause serve registrare **parte reale e immaginaria** dello
scarto: un errore da estremi ha firma proporzionale a (F'/F)delta, quello
del trapezio a k^2 h^2/12, quello di normalizzazione a delta_C/C.
Il solo modulo 1.65e-4 non permette di scegliere fra queste firme.

## Riproduzione e controlli

~~~sh
python3.13 calculations/audit_vaidya_error_budget.py
cd calculations
python3.13 -m unittest test_vaidya_quadrature.py test_outgoing_asymptotics.py
~~~

Il driver stampa JSON per ciascuna configurazione, inclusi estremi effettivi,
rapporti, condizionamento, residuo della primitiva e sweep dei dati al bordo.
**11 test superati**, fra cui tre nuovi controlli su estremi, limite dell'errore
puntuale e dati di Frobenius. git diff --check senza errori.

Dalla radice, verifica simbolica:

~~~sh
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_vaidya_error_budget.wl
~~~

Esito: **4/4 True**, codice di uscita 0. Usato il kernel Mathematica con
licenza locale; non il wrapper wolframscript. I test simbolici verificano i
due coefficienti di Frobenius e due identità di sensibilità/quadratura, non
certificano la soluzione numerica completa.

## Cosa chiedere al driver originale di Claude

Per ogni finestra stampare a,b nominali; r[mask][0],r[mask][-1]; C complesso;
integrale(n), integrale(p), F(b_eff)-F(a_eff), F(b)-F(a); massimo |n/p-1| e
kappa. Usare la stessa funzione primitiva e lo stesso C per entrambi gli
estremi. Questo separa gli errori senza dover indovinare il numero storico.
