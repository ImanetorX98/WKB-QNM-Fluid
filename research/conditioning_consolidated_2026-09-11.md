# Condizionamento direzionale: risultato consolidato

11 settembre 2026. Seguito di publication_core_audit_2026-09-11.md.
Scopo: consolidare la mappa WKB nel settore esterno di Schwarzschild, non
certificare il comportamento presso i turning point, Kerr o la pubblicabilità.

## 1. Enunciato esatto al bordo esterno

Fissati L>0, coordinata r*, potenziale e ramo uscente, sia omega=a+ib con b<0.
Supponiamo u -> omega/L, u' -> 0 e u'' -> 0 al bordo esterno. La formula

~~~text
B = (ln A)'  = -Re(u'/u)/2 - L Im(u),
C = (ln A)'' = -Re[u''/u-(u'/u)^2]/2 - L Im(u'),
Q = -(C+B^2)/L^2
~~~

implica Q_infinity=-b^2/L^2. Per perturbazioni reali (da,db) della frequenza,
a L e potenziale fissati, il differenziale relativo è dQ/Q=2 db/b.

Ne conseguono:

- lungo omega -> (1+delta)omega, delta reale, sensibilità relativa 2;
- per la componente b, misurando l'errore rispetto a |b|, sensibilità 2;
- con errore d'ingresso ||(da,db)||_2/|omega|, il massimo condizionamento
  relativo puntuale è kappa_infinity=2|omega|/|b|;
- rispetto alla sola a la sensibilità asintotica è zero.

Non sono risultati contraddittori: cambiano la direzione e la scala di
normalizzazione dell'errore. Q non è olomorfa in omega. Usiamo ln|Q| per
il differenziale relativo; il caso Q=0 o b=0 richiede un condizionamento
assoluto e non rientra nell'enunciato relativo.

La stessa conclusione vale per l'onda uscente esatta al bordo asintotico,
se le correzioni al suo fattore esponenziale hanno derivate logaritmiche
che svaniscono. A r finito i numeri qui sotto riguardano invece la mappa
WKB troncata specificata: non una famiglia off-shell di QNM esatti.

Questa legge segue dalla regola della catena: non ne rivendichiamo originalità.
Non dimostra impossibilità di predire frequenze con ogni funzionale di Q.

## 2. Derivate indipendenti, senza griglia di differenziazione

Nuovo modulo: calculations/conditioning_analytic.py.
Costruisce simbolicamente il momento uscente agli ordini 0 e 2, poi applica
D_*=f d/dr per ottenere u',u'' e le loro derivate rispetto a omega. Non usa
np.gradient né differenze finite per il Jacobiano in frequenza.

Posti v=u'/u e z=u''/u-v^2, la derivata olomorfa del momento fornisce

~~~text
v_w = u'_w/u - u' u_w/u^2,
z_w = u''_w/u - u'' u_w/u^2 - 2 v v_w.
~~~

Per d=1 oppure i si calcolano

~~~text
dB = -Re(v_w d)/2 - L Im(u_w d),
dC = -Re(z_w d)/2 - L Im(u'_w d),
dQ = -(dC+2 B dB)/L^2.
~~~

Il ramo della radice resta quello uscente esterno. Non si applica
automaticamente questa scelta all'altro lato della barriera.

### Verifiche indipendenti

1. Differenziazione annidata mpmath a 45 cifre, ricostruendo il momento
   indipendentemente dalle espressioni SymPy. A r=22,35,60, per ell=70,s=2,
   confronto di Q, partial_a Q e partial_b Q: tutti entro 1e-9 relativo.
2. Differenze centrate in frequenza delle nuove formule: entro 1e-6 relativo
   per entrambe le componenti al punto di prova.
3. Confronto con il precedente algoritmo a griglia: kappa entro 1e-5 relativo,
   sensibilità radiale entro 1e-3 assoluto nella finestra 20–50.
4. Limite r=1e6: Q e partial_b Q concordano con le formule asintotiche entro
   1e-8 relativo; la componente a si sopprime.
5. Kernel Mathematica locale: quattro identità simboliche asintotiche verificate
   (4/4 True, uscita 0). Usato WolframKernel, non il wrapper wolframscript.

Le soglie sopra sono quelle dei test, non una stima universale d'errore del modello.

## 3. Quantità globali: non confondere una mediana con una norma

Per la finestra finita W usiamo ||Q||_2^2=integrale_W Q^2 dr*.
Il Jacobiano J porta (da,db) nel profilo da Q_a + db Q_b. Definiamo

~~~text
G_ij = integrale_W Q_i Q_j dr*, i,j in {a,b},
kappa_L2 = |omega| sqrt(lambda_max(G)) / ||Q||_2.
~~~

È il massimo su una **singola perturbazione complessa comune a tutta la
finestra**, non la mediana dei massimi puntuali. La sensibilità radiale
globale è ||a Q_a+b Q_b||_2/||Q||_2.
Le quadrature sono in r* e hanno controllo di convergenza 501–3001 campioni
(accordo di kappa_L2 entro 1e-5 relativo nel test). Il dominio è fisso.

## 4. Valori riproducibili

Fondamentali s=2, finestra 20<r<50, momento WKB fino a ordine 2.
Le mediane sono su 3001 campioni uniformi in r*.

| ell | Sensibilità radiale, mediana | kappa puntuale, mediana | kappa del profilo L2 | kappa_infinity |
|---:|---:|---:|---:|---:|
| 20 | 1.982813 | 82.137285 | 82.200930 | 81.850249 |
| 40 | 1.982867 | 162.490588 | 162.616187 | 161.924256 |
| 70 | 1.982879 | 282.942036 | 283.160617 | 281.956495 |
| 100 | 1.982882 | 403.374298 | 403.685873 | 401.969483 |

Per ell=70 la sensibilità radiale L2 è 1.981879. Cambiando spin a ell=70,
kappa puntuale mediana è 282.992137 (s=0), 282.979613 (s=1), 282.942036 (s=2):
stesso andamento, non identità esatta tra spin.

### Distanza dalla barriera e ordine della WKB

ell=70,s=2:

| Finestra | Radiale ordine 0 | Radiale ordine 2 | kappa mediana ordine 0 | kappa mediana ordine 2 |
|---|---:|---:|---:|---:|
| 10–25 | 1.963496 | 1.963493 | 284.690419 | 284.690178 |
| 30–70 | 1.988912 | 1.988912 | 282.369591 | 282.369586 |
| 200–300 | 1.999243 | 1.999243 | 281.961153 | 281.961153 |
| 2.5–4 | 92.363302 | 195.597390 | 292.386079 | 385.650939 |

La regione esterna dà risultati stabili rispetto ai due ordini e tende al
limite analitico. Presso il picco la sensibilità cambia fortemente con il
troncamento: non è un semplice rumore di differenziazione spaziale, perché
qui le derivate sono analitiche. Ma **non è ancora una misura controllata
del problema esatto**. Servono approssimazione uniforme, scelta del ramo e
validazione del Jacobiano rispetto a una famiglia di soluzioni ben definita.

## 5. Conseguenza pratica, limitata al regime verificato

Con ell=70 e norma L2 del profilo esterno, un errore relativo complesso
generico eta produce al primo ordine un errore relativo al più circa
283.16 eta, nel caso peggiore. È una stima locale del differenziale, non
un limite certificato per perturbazioni finite arbitrariamente grandi.
Un errore massimo sul profilo di 1e-3 richiede quindi, nella stima lineare
peggiore, eta circa 3.5e-6. Se l'incertezza riguarda invece b in termini
relativi alla stessa b, bastano circa 5e-4 nella regione asintotica.

Il grande fattore normwise esprime il fatto che |b| è piccolo rispetto a
|omega|: dire 'la frequenza è accurata' senza specificare le due componenti
è insufficiente. Non implica che un metodo fase-ampiezza sia inutilizzabile.

## 6. Test e comandi

~~~sh
python3.13 calculations/conditioning_analytic.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_conditioning.wl
~~~

Da calculations:

~~~sh
python3.13 -m unittest test_conditioning_analytic.py test_conditioning_directional_audit.py test_repro_numerator_discrepancy.py test_vaidya_quadrature.py test_outgoing_asymptotics.py
~~~

La suite selezionata copre 7 test del condizionamento, più 13 regressioni
dei precedenti audit. Non è l'intera suite di tutti i moduli Kerr/Vaidya.
Esito dell'esecuzione: **20/20 superati**, senza errori.

## 7. Testo candidato per sostituire la conclusione del §5

> La sensibilità del potenziale di ampiezza dipende dalla direzione della
> perturbazione complessa e dalla norma adottata per l'errore. Nel limite
> uscente asintotico, Q_M=-(Im omega)^2/L^2: la sensibilità a una variazione
> proporzionale di omega è 2, mentre il condizionamento relativo rispetto
> alla norma complessa è 2|omega|/|Im omega|. Le due quantità sono compatibili
> e differiscono fortemente nel regime eikonale. I benchmark esterni su
> Schwarzschild confermano questa distinzione mediante derivate analitiche
> e controlli multiprecisione. Il risultato quantifica la propagazione
> dell'errore nella mappa WKB d'ampiezza, ma non costituisce un teorema di
> impossibilità per stimatori spettrali costruiti a partire da essa.

Questo testo non è stato ancora inserito nel manoscritto: occorre sostituire
coerentemente anche abstract, elenco dei contributi, richiamo su Kerr e
discussione. La figura vicino al picco resta sospesa. Il risultato esterno è
consolidato nei limiti dichiarati; originalità e vantaggio rispetto ai
diagnostici esistenti restano da dimostrare, come indicato nell'audit editoriale.
