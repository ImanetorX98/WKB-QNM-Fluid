# Audit spettrale: cosa dimostrano davvero i diagnostici di Madelung

**Data:** 8 settembre 2026. **Base analizzata:** commit ec97d56.

Questa nota aggiorna l'interpretazione dello snapshot del 7 settembre.
Non modifica i solver esistenti e non cancella i risultati precedenti:
separa i dati riproducibili dalle conclusioni che quei dati non bastano a
dimostrare. Le nuove formule sono nuove analisi nel repository, **non una
rivendicazione di priorità scientifica**.

## 1. Risposta in breve e mappa del progetto

Il progetto studia se la curvatura dell'ampiezza di una perturbazione QNM
possa diagnosticare l'errore WKB. Il «fluido» è una rappresentazione
ampiezza–fase dell'equazione d'onda, non un nuovo mezzo materiale.

| Ramo | Cosa fanno i programmi | Stato dopo questo audit |
|---|---|---|
| Schwarzschild bosonico | frequenze WKB1/WKB3, confronto Leaver, profili e scansione di 69 modi | benchmark numerico; non ancora stimatore utilizzabile senza soluzione di riferimento |
| Schwarzschild Dirac | partner radiali, gerarchia di spin e curvatura dell'ampiezza | contributi subprincipali da tenere distinti dal Madelung scalare |
| Pöschl–Teller | frequenze/autofunzioni esatte, residui locali e riferimento quadratico | banco di prova; qui sono stati ricavati i nuovi controlli analitici |
| Kerr scalare | riduzione simbolica esatta e limite di lenta rotazione | formalismo verificato; non ancora validazione numerica del diagnostico in Kerr |
| Vaidya | riduzione della PDE, split Madelung e bilancio di Kodama | formalismo verificato; non ancora simulazione che dimostri memoria aggiuntiva |

I tre risultati principali dell'audit sono:

1. Per il fondamentale Pöschl–Teller l'indicatore completo è esattamente
   una funzione del solo parametro eikonale. La correlazione non contiene
   informazione indipendente su questa famiglia.
2. Il riferimento uniforme attuale usa la frequenza esatta e dati di Cauchy
   locali; a quella frequenza non è un QNM della barriera quadratica.
   Il suo difetto misura differenze fra profili, non certifica ancora
   l'errore di una frequenza ignota.
3. Un controllo spettrale che incorpora le condizioni radiative è disponibile
   analiticamente in Pöschl–Teller. Fornisce un confronto severo per il futuro
   diagnostico Madelung, incluso il modo nodale n=1.

## 2. Convenzioni

Si considera

\[
\psi''+[\omega^2-L^2\operatorname{sech}^2y]\psi=0,\qquad
L>\tfrac12,\quad e^{-i\omega t}.
\]

Le condizioni sono uscenti verso entrambi gli estremi:
entrante nel buco nero sul lato sinistro nella corrispondente interpretazione
radiale. Poniamo

\[
\varepsilon=L^{-1},\quad N=n+\tfrac12,\quad B=N^2+\tfrac14,\quad
\lambda=\sqrt{L^2-\tfrac14},\quad
\omega_{\rm ex}=\lambda-iN,\quad \widehat\omega=\omega/L.
\]

Lo spettro esatto è noto: non è un risultato originale del progetto.
Si veda [Cardona–Molina](https://arxiv.org/abs/1711.00479).

## 3. Fondamentale: indicatore completo in forma chiusa

Per n=0, l'ipergeometrico della soluzione esatta si riduce a 1 e

\[
\psi_0\propto(\cosh y)^{i\omega_0},\qquad
A_0\propto\sqrt{\cosh y}.
\]

La forma dell'ampiezza è quindi indipendente da L; una normalizzazione
complessiva eventualmente dipendente da L si cancella nei rapporti.
Derivando, senza utilizzare la chiusura dell'ODE per definire Q, si trova

\[
\frac{A_0''}{A_0}=\frac{1+\operatorname{sech}^2y}{4},\qquad
Q_M=-\frac{\varepsilon^2}{4}(1+\operatorname{sech}^2y),
\]

\[
P=\varepsilon\operatorname{Im}(\psi_0'/\psi_0)
=\sqrt{1-\varepsilon^2/4}\tanh y,\qquad
|\widehat\omega_0|^2=1.
\]

La densità al denominatore del diagnostico implementato è dunque

\[
|\widehat\omega_0|^2+\operatorname{sech}^2y+P^2
=2-\frac{\varepsilon^2}{4}\tanh^2 y.
\]

Per qualsiasi finestra non negativa **fissata indipendentemente da L**,
definiamo

\[
c_w=\frac{\int w\tanh^2y\,dy}{\int w\,dy}.
\]

Allora l'indicatore completo del benchmark è esattamente

\[
\boxed{\mathcal E_M=
\frac{\varepsilon^2(2-c_w)}{8-\varepsilon^2 c_w}.}
\]

Con la griglia originale di 4001 punti,
y in ±2.25/√2 e w=exp(-y²), la quadratura dà
c_w=0.258076281399089. Il confronto con il programma originale per
L=1.5,2,4,8,16,32,64 ha discrepanza relativa massima 2.641×10⁻¹³.

**Conseguenza:** per questo modo e questa famiglia, conoscere E_M equivale
a conoscere una funzione monotona di ε². La correlazione prossima a 1 con
l'errore WKB non dimostra un vantaggio rispetto a un predittore basato su ε.
Questo argomento non annulla automaticamente i risultati multispin di
Schwarzschild, che richiedono controlli propri.

La formula vale per l'indicatore completo con finestra fissa. Per la finestra
ristretta dello script uniforme, c_w dipende a sua volta da ε; rimane però
determinato da ε e dalla finestra, non da un nuovo dato spettrale.

## 4. Errori WKB: coefficienti asintotici, non soltanto fit

La sostituzione delle derivate della barriera nelle formule implementate
dà esattamente

\[
\omega_1^2=L^2-2iNL,\qquad
\omega_3^2=L^2-B-2iNL\left(1-\frac{1}{8L^2}\right).
\]

Invece

\[
\omega_{\rm ex}^2=L^2-B-2iN\sqrt{L^2-\tfrac14}.
\]

Sul ramo a parte reale positiva, espandendo il rapporto delle frequenze
per L→∞ **a n fissato**, si ottiene

\[
\frac{\omega_1}{\omega_{\rm ex}}-1
=\frac{B}{2}\varepsilon^2+O(\varepsilon^3),
\qquad
\frac{\omega_3}{\omega_{\rm ex}}-1
=-\frac{iN}{128}\varepsilon^5+O(\varepsilon^6).
\]

Quindi, per l'errore relativo in modulo,

\[
\boxed{E_1\sim\frac{B}{2L^2},\qquad E_3\sim\frac{N}{128L^5}.}
\]

Queste formule spiegano le pendenze osservate. L'ordine L⁻⁵ discende da
cancellazioni dello spettro di questa barriera: non dimostra che una
correzione quartica/sestica della sola ampiezza produca un indicatore con
quell'ordine, né che L⁻⁵ sia universale.

Controllo a 80 cifre, L=256:

| n | L² E1 misurato | limite B/2 | L⁵ E3 misurato | limite N/128 |
|---:|---:|---:|---:|---:|
| 0 | 0.249999642376 | 0.25 | 0.00390625745060 | 0.00390625 |
| 1 | 1.24995014952 | 1.25 | 0.0117184147342 | 0.01171875 |
| 2 | 3.24962215812 | 3.25 | 0.0195294992739 | 0.01953125 |

Ripetendo il fit del programma uniforme originale sulla coda
L=16,24,32,48,64:

| n | difetto uniforme | Madelung completo nella finestra ristretta | errore WKB1 | errore WKB3 |
|---:|---:|---:|---:|---:|
| 0 | -1.98339 | -1.99220 | -1.99976 | -5.00031 |
| 2 | -1.94235 | -1.00471 | -1.98130 | -4.98539 |

Sono pendenze effettive log–log rispetto a L. La precedente pendenza
circa -1.76 per n=2 su L=4…16 non è un nuovo esponente asintotico
dimostrato. Il difetto uniforme tende numericamente verso una pendenza -2;
qui non ne viene fornita una dimostrazione analitica. La colonna WKB3 di
questo fit usa ancora il codice originale a precisione doppia; le formule
e la tabella precedente sono controllate separatamente ad alta precisione.

## 5. Il riferimento quadratico non soddisfa la quantizzazione QNM

Nel programma uniforme

\[
\varepsilon^2\psi_{\rm PC}''+(q_0+y^2)\psi_{\rm PC}=0,\qquad
q_0=\widehat\omega_{\rm ex}^2-1.
\]

Il codice calcola esplicitamente la frequenza esatta e impone al centro
ψ=1 e ψ'/ψ uguale al valore esatto. Per i modi pari, la derivata nulla
segue già dalla simmetria, ma la frequenza rimane un input esatto.

Con z=e^{-iπ/4}√(2/ε)y, l'indice di Weber è

\[
\nu=\frac{iq_0}{2\varepsilon}-\frac12.
\]

Il Wronskiano delle soluzioni uscenti dai due lati è proporzionale a

\[
\mathcal W_z\{D_\nu(z),D_\nu(-z)\}
=\frac{\sqrt{2\pi}}{\Gamma(-\nu)}.
\]

L'identità segue da Dν(z)=U(-ν-1/2,z) e dal
[Wronskiano DLMF 12.2.11](https://dlmf.nist.gov/12.2.E11).
La condizione QNM della parabola è quindi ν=0,1,2,… .
Valutando invece ν alla frequenza esatta Pöschl–Teller:

\[
\boxed{\nu-n
=N\left(\sqrt{1-\varepsilon^2/4}-1\right)
-\frac{iB\varepsilon}{2}.}
\]

Per ε>0 e n≥0 la parte immaginaria è strettamente negativa: la condizione
non è soddisfatta. Alla frequenza WKB1 si ottiene invece esattamente ν=n.

**Interpretazione corretta:** il programma costruisce una soluzione locale
esatta dell'equazione quadratica, adattata ai dati di Cauchy del modo esatto
della barriera completa. Questo è lecito come confronto di profili e non
era presentato dal codice come un solver globale. Non basta però a
identificarla con il QNM uscente della forma normale. Le condizioni
asintotiche della barriera completa richiedono il matching appropriato.

Il difetto attuale utilizza inoltre la curvatura della soluzione esatta
nel numeratore: è un diagnostico di riferimento, non ancora una quantità
disponibile quando la soluzione è ignota.

## 6. Controllo globale: coefficiente entrante e correzione di Newton

È possibile costruire un confronto spettrale che non richieda in input
la frequenza esatta. Per una frequenza di prova ω, definiamo

\[
a=\tfrac12+i\lambda-i\omega,\quad
b=\tfrac12-i\lambda-i\omega,\quad c=1-i\omega,\qquad
z=(1+\tanh y)/2.
\]

La soluzione normalizzata uscente a sinistra è

\[
\psi_L=[z(1-z)]^{-i\omega/2}\,{}_2F_1(a,b;c;z).
\]

Alla destra si decompone in

\[
\psi_L\sim C_{\rm out}e^{i\omega y}+C_{\rm in}e^{-i\omega y},
\quad
C_{\rm in}(\omega)=
\frac{\Gamma(1-i\omega)\Gamma(-i\omega)}
{\Gamma(a)\Gamma(b)}.
\]

La formula discende dalla
[connessione ipergeometrica DLMF 15.10.21](https://dlmf.nist.gov/15.10.E21).
Vicino ai QNM qui considerati i fattori al numeratore sono regolari e
non nulli; C_in=0 seleziona a=-n sul ramo di frequenza positiva.
Questa specificazione evita di applicare ingenuamente il quoziente in
corrispondenza di altre singolarità delle funzioni gamma.

Per una radice semplice e una prova sufficientemente vicina, poniamo

\[
\delta\omega=-\frac{C_{\rm in}}{\partial_\omega C_{\rm in}},
\qquad
\eta_{\rm spec}=\frac{|\delta\omega|}{|\omega|}.
\]

La correzione è calcolata senza usare ω_ex:

\[
\partial_\omega\log C_{\rm in}
=-i[\psi_D(1-i\omega)+\psi_D(-i\omega)-\psi_D(a)-\psi_D(b)],
\]

dove ψ_D è la digamma, non la funzione d'onda.
Il codice usa il reciproco negativo di questa derivata; il calcolo con
digamma va inteso fuori dalla radice, dove essa ha un polo.

### 6.1 Risultati numerici riproducibili

E è l'errore rispetto allo spettro esatto, usato **solo per validare**.
La terza colonna numerica è E dopo una sola correzione.

| n | L | ordine | E iniziale | η_spec/E | E dopo Newton |
|---:|---:|---:|---:|---:|---:|
| 0 | 4 | 1 | 1.553710573e-2 | 1.086251017 | 2.236021499e-3 |
| 0 | 4 | 3 | 3.844792014e-6 | 1.000022061 | 1.287317886e-10 |
| 0 | 16 | 1 | 9.762057942e-4 | 1.022872698 | 4.837483676e-5 |
| 0 | 16 | 3 | 3.727110398e-9 | 1.000000161 | 6.915225486e-16 |
| 0 | 64 | 1 | 6.103375949e-5 | 1.005943325 | 1.040759563e-6 |
| 0 | 64 | 3 | 3.638089834e-12 | 1.000000001 | 3.677785980e-21 |
| 1 | 4 | 1 | 6.767685962e-2 | 1.400973444 | 3.950602755e-2 |
| 1 | 4 | 3 | 1.025274883e-5 | 1.000041230 | 9.524441541e-10 |
| 1 | 16 | 1 | 4.833674627e-3 | 1.100030954 | 9.276817739e-4 |
| 1 | 16 | 3 | 1.109465420e-8 | 1.000000322 | 4.719362792e-15 |
| 1 | 64 | 1 | 3.049812195e-4 | 1.027681771 | 2.091392412e-5 |
| 1 | 64 | 3 | 1.090894287e-11 | 1.000000002 | 2.621541317e-20 |
| 2 | 4 | 1 | 1.423877008e-1 | 1.913536131 | 1.524684165e-1 |
| 2 | 4 | 3 | 1.398098038e-5 | 1.000073626 | 2.185210208e-9 |
| 2 | 16 | 1 | 1.233219044e-2 | 1.193292335 | 4.957542389e-3 |
| 2 | 16 | 3 | 1.820878358e-8 | 1.000000420 | 1.155704252e-14 |
| 2 | 64 | 1 | 7.919846333e-4 | 1.064198582 | 1.258311712e-4 |
| 2 | 64 | 3 | 1.816384198e-11 | 1.000000003 | 6.384204148e-20 |

Il modo n=1 è incluso senza dividere per l'ampiezza nel suo zero.
Questo elimina il problema nodale per **questo controllo spettrale**,
non regolarizza automaticamente Q_M.

Il caso n=2,L=4,WKB1 peggiora: da errore 0.14239 a 0.15247.
Newton è locale, non una garanzia globale né un limite superiore certificato.
Una normalizzazione analitica dipendente da ω modifica il passo lontano
dalla radice; la stima rimane equivalente al primo ordine vicino a una
radice semplice. Qui la normalizzazione è quella della soluzione di Jost.

Questa verifica sfrutta una barriera risolvibile: non dimostra che un
controllo altrettanto economico sia già disponibile per Schwarzschild/Kerr.
Newton, le funzioni di Jost e il matching QNM sono metodi standard.
I precedenti ampiezza–fase per Schwarzschild e Kerr includono
[Glampedakis–Andersson (2003)](https://arxiv.org/abs/gr-qc/0304030).

## 7. Ponte formale verso un diagnostico spettrale Madelung

Il seguente passaggio è una base esatta, non uno stimatore già validato.
Per una funzione di prova non nulla
\(\widetilde\psi=\widetilde A e^{i\widetilde S/\varepsilon}\), definiamo

\[
R=(\varepsilon^2\partial_y^2+q)\widetilde\psi,\quad
\widetilde P=\widetilde S',\quad
\widetilde Q=-\varepsilon^2\widetilde A''/\widetilde A.
\]

Si ha l'identità

\[
\boxed{\frac{R}{\widetilde\psi}
=\underbrace{q_R-\widetilde P^2-\widetilde Q}_{\text{difetto reale}}
+i\underbrace{\left[q_I+\varepsilon
\left(\widetilde P'+2\widetilde P\frac{\widetilde A'}{\widetilde A}\right)
\right]}_{\text{difetto di trasporto}}.}
\]

Quindi una differenza fra le sole curvature reali non è, in generale,
il residuo complesso dell'equazione.
In particolare il trasporto \(\widetilde A=\widetilde P^{-1/2}\) con
ampiezza e fase reali implica corrente radiale costante; per una soluzione
esatta non nulla richiederebbe q_I=0. Il momento ausiliario complesso della
chiusura WKB/Schwarziana non va identificato automaticamente con
P=ε Im(ψ'/ψ) del Madelung reale a frequenza QNM complessa.
Questo chiarisce il cambio di variabili implicito in §2.2 della nota covariante.

Se χ risolve esattamente la stessa equazione omogenea alla frequenza di prova,
il Wronskiano bilineare soddisfa

\[
\frac{d}{dy}\mathcal W(\chi,\widetilde\psi)
=\frac{\chi R}{\varepsilon^2},
\qquad
\mathcal W\big|_{y_a}^{y_b}
=\frac{1}{\varepsilon^2}\int_{y_a}^{y_b}\chi R\,dy.
\]

Non compare un complesso coniugato: è l'identità bilineare di Lagrange
per l'operatore differenziale, non un prodotto di probabilità.
Vale direttamente su un intervallo finito. Estenderla ai bordi QNM
richiede controllo dei termini asintotici e, se necessario, continuazione
complessa o regolarizzazione; le autofunzioni QNM divergono sulle sezioni
spaziali ordinarie.

Questo mostra quali ingredienti mancano a una norma locale di |ΔQ|:
fase del residuo, trasporto, matching ai bordi e sensibilità della
condizione spettrale D(ω) attraverso D'(ω).
Se χ è solo approssimata compare anche il residuo di χ: non si può
usare l'identità per promettere a costo nullo un controllo esatto.

## 8. Priorità aggiornata e criterio di originalità

La priorità non è imporre a un indicatore la pendenza L⁻⁵ osservata.
È definire prima quale errore spettrale il residuo debba controllare.

1. Costruire la forma uniforme alla **frequenza di prova**, con i rami
   radiativi e le condizioni di matching espliciti. Correggere
   ampiezza, fase e quantizzazione in modo coerente.
2. Confrontare almeno: solo |ΔQ|; residuo complesso completo; difetto
   spettrale normalizzato per D'(ω). Separare quantità ottenibili dalla
   soluzione approssimata da quelle che richiedono una soluzione esatta.
3. In Pöschl–Teller usare C_in e le formule di §4 come riferimento,
   includendo n=1 e punti fuori dal buon regime WKB.
4. Passare a famiglie di potenziali con parametri indipendenti:
   controllare la previsione a ε fissato, validare su intere famiglie
   escluse dalla calibrazione e confrontare anche costo/accuratezza
   con il calcolo spettrale diretto. La sola rimozione lineare di log ε²
   non elimina ogni dipendenza non lineare dallo stesso parametro.
5. Solo dopo questo test, estendere lo stesso campo scalare a Kerr
   e poi a Vaidya. In Vaidya la domanda resta se la curvatura del profilo
   aggiunga informazione rispetto a massa istantanea, derivate della massa,
   propagazione e mescolamento dei modi; qui non è ancora dimostrato.

**Possibile contributo originale da verificare:** uno stimatore spettrale
economico, risolto per ordine, che traduca difetti di ampiezza e trasporto
in un errore QNM e rimanga utilizzabile vicino ai turning point e ai nodi.
Le identità e il metodo di Newton di questa nota non costituiscono da soli
tale contributo. Non è stata svolta qui una nuova ricerca bibliografica
esaustiva né è dimostrata la pubblicabilità del risultato.

## 9. Verifiche e riproduzione

File aggiunti:

- [audit numerico](../calculations/audit_poschl_teller.py);
- [audit simbolico Wolfram](../calculations/verify_spectral_audit.wl).

Dalla radice del repository:

~~~bash
python3.13 calculations/audit_poschl_teller.py --uniform-tail
/Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noprompt -script calculations/verify_spectral_audit.wl
~~~

Eseguiti il giorno dell'audit:

- **18 nuovi controlli Wolfram superati**: ampiezza, normalizzazione,
  indicatore integrato, errori WKB, coefficienti Iyer–Will, indice di Weber,
  riduzione ipergeometrica, derivata di Jost, residuo polare e Wronskiano.
- **19 controlli Wolfram preesistenti superati**, compresi Kerr e Vaidya.
- **32 test Python preesistenti superati**: 20 core e 12 calculations.
- Audit numerico: tutte le asserzioni superate, mpmath 1.3.0 a 80 cifre.
  Connessione ipergeometrica verificata in 24 combinazioni (L=4,16;
  n=0,1,2; WKB1/WKB3; z=0.37,0.8), discrepanza scalata massima
  1.0368×10⁻⁸⁰. Derivata logaritmica confrontata con differenziazione
  indipendente: discrepanza relativa massima 1.4210×10⁻⁷⁹.

Le cifre di precisione aritmetica non sono una misura dell'accuratezza
fisica dell'approssimazione WKB. I controlli simbolici verificano algebra
locale; le condizioni al bordo sono controllate separatamente.

È stato utilizzato direttamente **Mathematica/WolframKernel 13.3.1 ARM64**,
non il launcher wolframscript. Non è stata modificata la licenza.
I risultati sono salvati localmente; questo audit non effettua commit o push.
