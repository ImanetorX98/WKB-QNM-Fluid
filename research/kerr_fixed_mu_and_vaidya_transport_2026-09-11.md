# Kerr corretto e trasporto covariante di Vaidya

11 settembre 2026. Analisi successive a kerr_vaidya_core_status_2026-09-11.md.
Questa nota distingue misure numeriche riprodotte, identità locali esatte e
ipotesi ancora necessarie per una condizione di solvibilità globale.

## 1. Kerr: la misura radiale corretta dà ordine 2

Nuovo programma: calculations/kerr_fixed_mu_radial.py.
Non usa interpolazione in m: seleziona solo coppie intere (ell,m) con
mu=m/(ell+1/2) **esattamente** pari a 2/5, 2/3 oppure 4/5.
Il caso mu=1/2 non ammette tali coppie intere: non viene simulato arrotondando.

Per ogni chat=a omega_hat fissato si stima A0 da autovalori a ell circa
240,300,360,420, poi si misura la radiale a ell circa 40,60,80,120,160.
I due insiemi sono disgiunti. Il fit di A/L^2 include 1,1/L,1/L^2,1/L^3:
**il coefficiente 1/L è libero**, non imposto nullo. Il termine geometrico
h''/h della radiale è derivato analiticamente con D_*=(Delta/H)d/dr.

### Sweep a omega_hat=0.225 reale, finestra 2.5<r<4

| mu | a=0 | a=0.3 | a=0.6 | a=0.9 |
|---|---:|---:|---:|---:|
| 2/5 | 2.000000 | 2.000000 | 1.999999 | 1.999997 |
| 2/3 | 2.000000 | 2.000000 | 2.000001 | 2.000004 |
| 4/5 | 2.000000 | 2.000001 | 2.000003 | 2.000009 |

Le celle riportano la pendenza del massimo spaziale di |q/L^2-Q0| contro 1/L.
I coefficienti A1 del fit di testa hanno modulo inferiore a 1e-9 in questo
sweep; non sono cifre significative di un A1 fisico diverso da zero.

Per a=0.6,mu=2/3, spostando i riferimenti a ell circa 180–360 oppure 300–480,
la pendenza passa da 2.0000014499 a 2.0000014578. Non è un artefatto dominante
dell'estrapolazione del riferimento. Questa verifica non prova una parità
esatta a tutti gli ordini o per ogni valore dei parametri.

### Controllo sul picco eikonale autoconsistente

Per non limitarsi a una frequenza reale scelta a priori, è stato implementato
un secondo controllo: risolvere Q0=partial_r Q0=0 con A0 estrapolato a mu
fisso per ogni frequenza tentata. Mu=2/3; finestra r_peak +/- 0.2, esterna
all'orizzonte in tutti i casi riportati.

| a | omega_hat reale | r_peak | pendenza |
|---:|---:|---:|---:|
| 0.0 | 0.192450090 | 3.000000000 | 2.000000000 |
| 0.3 | 0.209928011 | 2.734032593 | 2.000000296 |
| 0.6 | 0.236005163 | 2.370738789 | 2.000001574 |
| 0.9 | 0.287199475 | 1.770672242 | 2.000005148 |

Residui della radice doppia inferiori a 7e-15 nelle esecuzioni riportate.
Questa è una misura del potenziale a frequenza eikonale reale congelata,
**non** della famiglia completa di frequenze QNM complesse a overtone fissato.

### Conseguenza

La vecchia pendenza 1 non sopravvive al controllo che mantiene davvero fissi
i parametri dichiarati. Il §9 e la curva Kerr della Figura 1 non possono
essere conservati come prova che la rotazione scalare introduca un A1
geometrico analogo al termine dei partner di Dirac.

Resta distinta la possibilità di termini 1/L lungo una **traiettoria fisica**
omega_hat(L)=Omega0+Omega1/L+... . Anche con A1 intrinseco nullo si ottiene
un termine partial_chat A0 * a Omega1/L. Questo include, per esempio, la
dipendenza dovuta allo smorzamento; non è la stessa affermazione del vecchio
test a chat costante. Non confondere effetti della traiettoria, arrotondamento
di m e struttura del potenziale a parametri congelati.

La validità del solutore WKB3 complesso, il suo condizionamento e il caso
quasi estremale non sono stati rivalidati in questo intervento.

## 2. Vaidya: legge esatta della normalizzazione

Indichiamo con G(r,M) il modo congelato, omega0(M) la sua frequenza, e

~~~text
psi = A(v) G(r,M(v)) exp[-i integral(omega0+delta_omega) dv],
alpha = A'/A,
S = -2 partial_r partial_M G.
~~~

La PDE, dopo fattorizzazione della fase e di A, contiene all'ordine considerato

~~~text
L_omega0 G + 2(alpha-i delta_omega) G_r - Mdot S.
~~~

Questa identità viene dalla PDE, non da una scelta di prodotto scalare.
Per G_tilde=c(M)G e A_tilde=A/c(M), con c non nullo e differenziabile,

~~~text
S_tilde = c S - 2 c_M G_r,
alpha_tilde = alpha - Mdot c_M/c.
~~~

Il residuo fattorizzato si moltiplica per c; il campo completo e la PDE
rimangono invariati. Sono identità locali esatte, verificate con Mathematica.
Non sono una simmetria di gauge gravitazionale: è libertà di scelta della
normalizzazione complessa del modo e della sua ampiezza.

## 3. Perché N isolato non è un osservabile

Assumendo una proiezione bilineare regolarizzata linearmente e coerentemente,
definiamo schematicamente

~~~text
P = Reg integral mu G G_r dr,
N = Reg integral mu G S dr.
~~~

Sotto la normalizzazione precedente:

~~~text
P_tilde = c^2 P,
N_tilde = c^2 N - 2 c c_M P.
~~~

Queste leggi seguono dall'algebra degli integrandi e dalla linearità della
regolarizzazione. Non dimostrano ancora che questi due integrali siano tutti
i termini richiesti dalla proiezione Fredholm del problema con bordi QNM.
I contributi di bordo devono essere inclusi e trasformati coerentemente.

Se la proiezione completa si riduce alla forma sopra e P non è nullo,
la solvibilità dà

~~~text
alpha - i delta_omega = Mdot K,
K = N/(2P),
K_tilde = K - c_M/c.
~~~

In un intervallo regolare di M si può scegliere c_M/c=K e annullare il termine
di trasporto K in quella normalizzazione. Ciò non annulla il campo né i suoi
effetti fisici: li redistribuisce nel profilo G. Si può anche scegliere
c(M0)=1, lasciando il profilo invariato alla massa M0, ma cambiarne la derivata
rispetto a M e dunque il numeratore. Questo mostra perché N non nullo non
prova da solo uno shift osservabile o una componente complementare Z1.

## 4. Una quantità candidata indipendente dalla normalizzazione

La combinazione

~~~text
Xi(r,M) = K(M) + partial_M ln G(r,M)
~~~

è invariata sotto G -> c(M)G e K -> K-c_M/c, dove G non si annulla.
È stata verificata simbolicamente anche questa identità.
Per l'ansatz adiabatico con trasporto, a r fisso:

~~~text
partial_v ln psi = -i omega0 + Mdot Xi.
~~~

Nel caso complesso a singolo modo, -Im(partial_v ln psi) definisce una
frequenza di fase locale e -Re(partial_v ln psi) un tasso di decadimento.
Il candidato di correzione di fase è dunque -Mdot Im Xi, non N isolato.

**Limiti:** la formula per il campo con solo G è un'identità dell'ansatz.
Per applicarla alla soluzione perturbata completa servono una gerarchia
adiabatica (anche per Mddot), il trattamento di Z1, le condizioni iniziali e
gli eventuali transitori. Occorre inoltre fissare coordinata temporale,
posizione dell'osservatore e definizione della fase del segnale. Non è stata
calcolata numericamente Xi, né validata come correzione della forma d'onda.
Non ne rivendichiamo originalità: è una combinazione di trasporto covariante.

Questo è tuttavia un obiettivo più definito per il prossimo calcolo: verificare
prima la proiezione globale e il denominatore appropriato, poi valutare Xi,
e infine confrontarlo con evoluzioni temporali a piccolo tasso di variazione
della massa. Un termine locale in Mdot non costituisce già memoria non locale.

## 5. Riproducibilità e stato

~~~sh
python3.13 calculations/kerr_fixed_mu_radial.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_vaidya_normalization.wl
~~~

Da calculations:

~~~sh
python3.13 -m unittest test_kerr_fixed_mu_radial.py test_vaidya_numerator_factored.py test_conditioning_analytic.py test_conditioning_directional_audit.py
~~~

**23 test superati**: quattro nuovi test Kerr (con sottocasi), dodici test
Vaidya del numeratore, sette test di condizionamento. Mathematica restituisce
**4/4 True** sulle identità di normalizzazione. Non sono test di una soluzione
non adiabatica completa, che resta da costruire.

I moduli storici e i loro vecchi test non sono stati riscritti: riproducono
ancora il campionamento arrotondato e non devono essere citati come validazione
del nuovo enunciato a mu fissato. I nuovi calcoli sono separati per confrontabilità.
Il manoscritto non è pronto all'invio: la classificazione Kerr va corretta e
il numeratore Vaidya deve essere collegato a una quantità fisica verificata.
