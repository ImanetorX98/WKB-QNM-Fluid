# Stato della ricerca e documento di ripresa

**Snapshot:** 7 settembre 2026

**Repository pubblica:** https://github.com/ImanetorX98/WKB-QNM-Fluid

**Base scientifica consolidata fino al commit:** 09ef057

**Scopo del file:** permettere di riprendere il lavoro senza dipendere dalla
cronologia della conversazione.

Questo documento distingue sempre:

- identità esatte, cioè riscritture verificate delle equazioni;
- risultati numerici effettivamente misurati;
- risultati negativi, che escludono una strada;
- ipotesi e candidati ancora da validare;
- precedenti noti, per evitare rivendicazioni di novità eccessive.

I dettagli completi restano nelle note specialistiche:

- [formalismo covariante](covariant_formalism.md);
- [benchmark statico](static_benchmark.md);
- [barriera esatta e residuo uniforme](exact_barrier_benchmark.md);
- [mappa della letteratura](prior_art.md);
- [nota Dirac](../calculations/dirac_wkb_hydrodynamics_notes.md).

## 1. Idea e obiettivo scientifico

Il progetto studia una rappresentazione ampiezza–fase delle perturbazioni
quasi-normali:

\[
\Phi=Ae^{iS/\varepsilon}.
\]

La curvatura dell'ampiezza produce il termine

\[
Q_g=-\varepsilon^2\frac{\Box_gA}{A},
\]

o, dopo una riduzione radiale,

\[
Q_M=-\varepsilon^2\frac{A''}{A}.
\]

La trasformazione di Madelung, la WKB e i QNM non sono nuovi separatamente.
La domanda di ricerca controllabile è più stretta:

> Un residuo di curvatura dell'ampiezza, definito rispetto alla soluzione WKB
> uniforme appropriata, può stimare il dominio di affidabilità e l'errore di
> troncamento della WKB per risonanze QNM?

Il programma segue tre livelli geometrici:

1. spazio-tempi statici: Schwarzschild e barriere esattamente risolvibili;
2. spazio-tempi stazionari non statici: slow Kerr e Kerr;
3. spazio-tempi non stazionari: Vaidya e risposta ritardata rispetto al fondo
   congelato.

“Fluido QNM” significa rappresentazione ampiezza–fase. Non significa che la
perturbazione gravitazionale sia un fluido materiale, né che \(A^2\) sia una
probabilità globale normalizzabile.

## 2. Convenzioni

- Unità: \(G=c=\hbar=1\).
- Schwarzschild: \(x=r/M\), \(x_*=r_*/M\).
- Dipendenza temporale: \(e^{-i\omega t}\).
- Un QNM smorzato soddisfa \(\operatorname{Im}\omega<0\).
- Condizioni QNM: entrante all'orizzonte e uscente all'infinito.
- Parametro eikonale bosonico:
  \[
  L=\ell+\frac12,\qquad \varepsilon=L^{-1}.
  \]
- Parametro eikonale Dirac:
  \[
  K=|\kappa|,\qquad \varepsilon=K^{-1},\qquad \omega=K\Omega.
  \]

Per la ODE

\[
\varepsilon^2\psi''+q\psi=0,\qquad q=q_R+iq_I,
\]

con \(\psi=Ae^{iS/\varepsilon}\), \(A,S\) reali e \(P=S'\), valgono
esattamente

\[
P^2+Q_M=q_R,\qquad Q_M=-\varepsilon^2\frac{A''}{A},
\tag{2.1}
\]

\[
\varepsilon(A^2P)'=-q_IA^2.
\tag{2.2}
\]

La seconda equazione mostra perché il profilo radiale separato di un QNM è un
sistema aperto: la frequenza complessa genera una sorgente o un pozzo radiale.
La corrente della PDE covariante completa resta invece conservata.

## 3. Risultati statici: Schwarzschild bosonico

### 3.1 Implementazione WKB e controllo campione

Il potenziale massless di spin \(s=0,1,2\) è

\[
V_s=f\left[\frac{\ell(\ell+1)}{r^2}
+\frac{2M(1-s^2)}{r^3}\right],
\qquad f=1-\frac{2M}{r}.
\]

Il codice calcola simbolicamente le derivate rispetto a \(r_*\) fino al sesto
ordine e applica WKB1 e WKB3 di Iyer–Will.

Per il modo gravitazionale assiale \((s,\ell,n)=(2,2,0)\):

\[
\frac{r_0}{M}=3.280776406166,
\]

\[
V_0=0.151286699570,\qquad
V_0''=-0.009919411728,
\]

\[
\Lambda_2=-0.141970003076,\qquad
\Lambda_3=-0.054526599347,
\]

\[
\Omega_{\rm WKB1}=0.398849605924-0.088285381437i,
\]

\[
\boxed{\Omega_{\rm WKB3}=0.373162064857-0.089217447231i}.
\]

Il valore WKB3 riproduce il benchmark Iyer–Will entro la tolleranza dei test.

### 3.2 Benchmark indipendente Madelung–Leaver

È stata eseguita una scansione di 69 modi:

- \(s=0,1\), \(1\leq\ell\leq8\);
- \(s=2\), \(2\leq\ell\leq8\);
- \(n=0,1,2\);
- 23 modi per overtone.

La frequenza è ottenuta con la frazione continua di Leaver. Il profilo è
integrato dall'ODE di Regge–Wheeler e non è costruito con la WKB di barriera.
Il diagnostico fissato prima della scansione è

\[
\mathcal E_M=
\frac{\int w|Q_M|\,dx_*}
{\int w\left(|\varepsilon\omega|^2+
\varepsilon^2|V|+|P|^2\right)\,dx_*}.
\tag{3.1}
\]

La finestra gaussiana è centrata al massimo della barriera, con
\(\sigma=\sqrt{V_0/(-V_0'')}\), e troncata a \(2.25\sigma\).

| \(n\) | errore | Pearson | Spearman | Pearson parziale dato \(\varepsilon^2\) | RMSE LOOCV con \(\mathcal E_M\) [dex] | RMSE LOOCV con \(\varepsilon^2\) [dex] |
|---:|---|---:|---:|---:|---:|---:|
| 0 | WKB1 | 1.0000 | 1.0000 | 0.9997 | 0.0014 | 0.0387 |
| 0 | WKB3 | 0.9998 | 1.0000 | 0.9855 | 0.0234 | 0.0608 |
| 1 | WKB1 | 0.9978 | 0.9872 | 0.9641 | 0.0305 | 0.0517 |
| 1 | WKB3 | 0.9987 | 0.9872 | 0.9214 | 0.0527 | 0.0953 |
| 2 | WKB1 | 0.9999 | 1.0000 | 0.9959 | 0.0059 | 0.0670 |
| 2 | WKB3 | 0.9989 | 1.0000 | 0.9925 | 0.0421 | 0.1129 |

Controlli numerici:

- residuo massimo della ricostruzione a differenze finite di \(A''/A\):
  \(4.99\times10^{-6}\) per \(n=0\),
  \(3.13\times10^{-4}\) per \(n=1\),
  \(4.74\times10^{-5}\) per \(n=2\);
- dimezzare la profondità della frazione continua cambia le frequenze di meno
  di \(2.2\times10^{-16}\) in relativo;
- per \(n=0\), cambiare la larghezza della finestra fra
  \(0.75\sigma\) e \(1.25\sigma\) conserva Spearman \(=1\) e porta la RMSE
  WKB1 soltanto da 0.0013 a 0.0015 dex.

Interpretazione supportata:

> A overtone fissato, \(\mathcal E_M\) ordina gli errori WKB in Schwarzschild
> meglio del solo parametro \(1/L^2\), includendo informazione sullo spin del
> potenziale.

Limite osservato: mescolando gli overtone, Pearson scende a 0.921 per WKB1 e
0.673 per WKB3. Il \(Q_M\) completo misura soprattutto la correzione totale
all'eikonale e non è automaticamente il residuo del metodo già troncato a
ordine superiore.

## 4. Barriera esatta di Pöschl–Teller

### 4.1 Soluzione esatta

Per

\[
V(y)=L^2\operatorname{sech}^2y,
\]

le frequenze QNM sono

\[
\omega_n=\sqrt{L^2-\frac14}-i\left(n+\frac12\right).
\tag{4.1}
\]

Con \(z=(1+\tanh y)/2\) e
\(\lambda=\sqrt{L^2-1/4}\), il modo esatto è

\[
\psi_n=[z(1-z)]^{-i\omega_n/2}
{}_2F_1\left(-n,-n-2i\lambda;
\frac12-n-i\lambda;z\right).
\tag{4.2}
\]

L'ipergeometrico termina. Frequenza e autofunzione non dipendono quindi da
Leaver o da un'integrazione approssimata.

### 4.2 Diagnostico completo

Per \(L=1.5,2,3,4,6,8,12,16\):

| overtone | dominio | Pearson con WKB1 | Pearson con WKB3 | Spearman |
|---:|---|---:|---:|---:|
| 0 | regolare | 1.0000 | 1.0000 | 1.0000 |
| 1 | nodale | non definito | non definito | non definito |
| 2 | regolare | 0.9994 | 0.9989 | 1.0000 |

Il controllo a differenze finite di \(-L^{-2}A''/A\) ha residuo relativo
massimo inferiore a \(9\times10^{-9}\) sui modi regolari.

Questo dimostra che la correlazione osservata in Schwarzschild non è un
artefatto specifico del metodo di Leaver o del potenziale di Regge–Wheeler.
Non dimostra ancora potere predittivo indipendente dallo scaling in \(1/L\),
perché la famiglia Pöschl–Teller usata ha un solo parametro.

### 4.3 Risultato nodale

Per \(n=1\):

\[
{}_2F_1(-1,-1-2i\lambda;-1/2-i\lambda;z)
=1-2z=-\tanh y,
\]

quindi

\[
\psi_1(0)=0.
\]

Poiché \(A=|\psi|\) si annulla, \(Q_M=-\varepsilon^2A''/A\) non è definito
nel centro della barriera. Il codice marca il caso come nodale e restituisce
NaN invece di un valore spurio.

Possibili trattamenti futuri, non equivalenti:

1. carte locali prive di zeri;
2. deflazione esplicita dei fattori nodali;
3. fase e ampiezza complesse invece di \(A=|\psi|\).

## 5. Risultato negativo: la gerarchia locale non è uniforme

La chiusura formale

\[
q=u^2-\varepsilon^2
\frac{(u^{-1/2})''}{u^{-1/2}}
\tag{5.1}
\]

è stata sviluppata con

\[
u^{(0)}=u_0,\qquad
u^{(1)}=u_0+\varepsilon^2u_2,\qquad
u^{(2)}=u_0+\varepsilon^2u_2+\varepsilon^4u_4.
\]

Mathematica verifica che i residui
\(\mathfrak R_0,\mathfrak R_1,\mathfrak R_2\) iniziano rispettivamente a
\(\varepsilon^2,\varepsilon^4,\varepsilon^6\).

Sui fondamentali Schwarzschild, però, gli integrali normalizzati sulla
barriera risultano tipicamente

\[
\mathfrak R_0\sim0.06\text{--}0.17,\qquad
\mathfrak R_1\sim0.18\text{--}0.50,\qquad
\mathfrak R_2\sim0.80\text{--}0.98.
\tag{5.2}
\]

La successione peggiora aumentando l'ordine. Non è un errore dell'algebra:
l'espansione locale è singolare quando \(q\to0\) e i turning point
coalescono. L'ordine formale a \(q\) fissato non è uniforme nella regione che
determina la quantizzazione QNM di barriera.

Conclusione operativa:

> Non usare \(Q_M-Q_M^{(N)}\) costruito dalla serie WKB locale come
> diagnostico dell'errore di barriera.

## 6. Primo risultato positivo del residuo uniforme

### 6.1 Forma parabolico-cilindrica

Per

\[
\varepsilon^2\psi''+
[\widehat\omega^2-\operatorname{sech}^2y]\psi=0,
\qquad \varepsilon=L^{-1},
\]

la forma quadratica al massimo è \(q_{\rm PC}=q_0+y^2\). Per la forma
generale \(q=q_0+q_2x^2/2\), la trasformazione

\[
\zeta=e^{-i\pi/4}\frac{(2q_2)^{1/4}}{\sqrt\varepsilon}x,
\qquad
\nu=\frac{iq_0}{\varepsilon\sqrt{2q_2}}-\frac12
\tag{6.1}
\]

porta esattamente all'equazione della funzione parabolico-cilindrica
\(D_\nu(\zeta)\). Questa mappa è verificata simbolicamente in Mathematica.

Il calcolo numerico integra la forma quadratica con gli stessi dati di Cauchy
del modo esatto al massimo. Ciò equivale a scegliere la corretta combinazione
dei rami \(D_\nu(\zeta)\) e \(D_\nu(-\zeta)\), evitando errori di parità.

Il residuo proposto è

\[
\boxed{
\Delta Q_M^{\rm unif}
=-\varepsilon^2\left(
\frac{A''}{A}-\frac{A_{\rm PC}''}{A_{\rm PC}}
\right)
}.
\tag{6.2}
\]

La regione uniforme è lo strato \(y=O(\sqrt\varepsilon)\). Infatti

\[
\frac{1-\operatorname{sech}^2(\sqrt\varepsilon\eta)}
{\varepsilon}
=\eta^2-\frac23\varepsilon\eta^4
+\frac{17}{45}\varepsilon^2\eta^6+O(\varepsilon^3),
\tag{6.3}
\]

anch'esso verificato con Mathematica.

### 6.2 Dati numerici

| \(L\) | \(n\) | \(\Delta Q_M^{\rm unif}\) | \(Q_M\) completo | rapporto | errore WKB1 | errore WKB3 |
|---:|---:|---:|---:|---:|---:|---:|
| 1.5 | 0 | 4.341e-02 | 1.007e-01 | 4.310e-01 | 1.076e-01 | 5.450e-04 |
| 2 | 0 | 2.587e-02 | 5.742e-02 | 4.505e-01 | 6.125e-02 | 1.260e-04 |
| 3 | 0 | 1.228e-02 | 2.605e-02 | 4.716e-01 | 2.751e-02 | 1.630e-05 |
| 4 | 0 | 7.141e-03 | 1.484e-02 | 4.811e-01 | 1.554e-02 | 3.845e-06 |
| 6 | 0 | 3.287e-03 | 6.696e-03 | 4.909e-01 | 6.927e-03 | 5.041e-07 |
| 8 | 0 | 1.884e-03 | 3.798e-03 | 4.960e-01 | 3.901e-03 | 1.194e-07 |
| 12 | 0 | 8.535e-04 | 1.703e-03 | 5.013e-01 | 1.735e-03 | 1.571e-08 |
| 16 | 0 | 4.849e-04 | 9.622e-04 | 5.040e-01 | 9.762e-04 | 3.727e-09 |
| 1.5 | 2 | 1.926e-02 | 6.568e-01 | 2.932e-02 | 4.088e-01 | 7.431e-04 |
| 2 | 2 | 1.608e-02 | 5.320e-01 | 3.023e-02 | 3.234e-01 | 2.521e-04 |
| 3 | 2 | 1.018e-02 | 3.743e-01 | 2.720e-02 | 2.094e-01 | 4.891e-05 |
| 4 | 2 | 6.724e-03 | 2.849e-01 | 2.360e-02 | 1.424e-01 | 1.398e-05 |
| 6 | 2 | 3.486e-03 | 1.909e-01 | 1.826e-02 | 7.529e-02 | 2.160e-06 |
| 8 | 2 | 2.113e-03 | 1.431e-01 | 1.477e-02 | 4.556e-02 | 5.460e-07 |
| 12 | 2 | 1.011e-03 | 9.507e-02 | 1.064e-02 | 2.146e-02 | 7.542e-08 |
| 16 | 2 | 5.900e-04 | 7.114e-02 | 8.293e-03 | 1.233e-02 | 1.821e-08 |

Correlazioni logaritmiche:

| overtone | confronto | Pearson | Spearman |
|---:|---|---:|---:|
| 0 | \(\Delta Q_M^{\rm unif}\) vs WKB1 | 0.9999 | 1.0000 |
| 0 | \(\Delta Q_M^{\rm unif}\) vs WKB3 | 0.9999 | 1.0000 |
| 2 | \(\Delta Q_M^{\rm unif}\) vs WKB1 | 0.9998 | 1.0000 |
| 2 | \(\Delta Q_M^{\rm unif}\) vs WKB3 | 0.9970 | 1.0000 |

Fit di potenza sui punti \(L=4,6,8,12,16\):

| \(n\) | \(\Delta Q_M^{\rm unif}\) | \(Q_M\) completo | errore WKB1 | errore WKB3 |
|---:|---:|---:|---:|---:|
| 0 | \(L^{-1.9408}\) | \(L^{-1.9739}\) | \(L^{-1.9964}\) | \(L^{-5.0051}\) |
| 2 | \(L^{-1.7599}\) | \(L^{-1.0018}\) | \(L^{-1.7719}\) | \(L^{-4.7999}\) |

Variando la larghezza gaussiana uniforme fra 0.75 e 1.5 volte quella
nominale:

- \(n=0\): esponente fra -1.949 e -1.941;
- \(n=2\): esponente fra -1.780 e -1.759.

Conclusioni:

1. Il residuo uniforme recupera lo scaling dell'errore WKB1 entro circa 0.06
   nell'esponente per entrambi gli overtone regolari.
2. L'elevata correlazione con WKB3 non implica che il suo ordine asintotico sia
   corretto.
3. La forma quadratica non riproduce il residuo WKB3, che scala circa come
   \(L^{-5}\).
4. Per WKB3 occorre calcolare le correzioni dell'ampiezza uniforme dovute ai
   termini quartico, sestico e successivi della forma normale.

Questo è il risultato più promettente, al momento, per un futuro articolo.

## 7. Risultati Dirac su Schwarzschild

Per Dirac massless:

\[
W=K\frac{\sqrt f}{r},\qquad
V_\pm=W^2\pm W',
\]

\[
\boxed{
V_\pm=
K^2\frac f{r^2}
\pm K\frac{\sqrt f(3M-r)}{r^3}
}.
\tag{7.1}
\]

I partner sono di Darboux e sono esattamente isospettrali per
\(\omega\neq0\), anche se una WKB troncata può rompere leggermente
l'equivalenza.

Inserendo il parametro semiclassico prima del disaccoppiamento:

\[
\varepsilon^2Z_\sigma''+
\left[\Omega^2-h^2-\sigma s\varepsilon h'\right]Z_\sigma=0,
\qquad h=\frac{\sqrt f}{r}.
\]

La chiusura è

\[
\boxed{
\operatorname{Re}\Omega^2
=P^2+h^2+
\underbrace{\sigma s\varepsilon h'}_{O(\varepsilon)}
+\underbrace{Q_M}_{O(\varepsilon^2)}
}.
\tag{7.2}
\]

Risultato strutturale:

> La gerarchia Dirac non è una successione di soli potenziali di Madelung.
> Il termine di spin connection compare già a \(O(\varepsilon)\); il
> Madelung scalare comincia a \(O(\varepsilon^2)\), seguito da termini misti
> spin–gradiente.

Il primo coefficiente di Madelung è

\[
Q_2=\frac{q_0''}{4q_0}
-\frac{5(q_0')^2}{16q_0^2}.
\tag{7.3}
\]

Il primo termine misto è

\[
Q_3=
\frac{q_1''}{4q_0}
-\frac{5q_0'q_1'}{8q_0^2}
-\frac{q_1q_0''}{4q_0^2}
+\frac{5q_1(q_0')^2}{8q_0^3}.
\tag{7.4}
\]

Scaling misurato per \(K=2,4,8,16\), \(\tau=+1\), \(n=0\):

| regione | pendenza spin | pendenza \(Q_M\) | \(\min|P|\) |
|---|---:|---:|---:|
| lontano dai turning point, \(20<x<50\) | 1.0000 | 1.84 | circa 0.187 |
| massimo della barriera | 1.14 | 1.40 | circa \(10^{-5}\) |

La pendenza di Madelung tende verso 2 spostando il bordo iniziale più vicino
all'orizzonte:

\[
x_{\min}=2.02\to1.557,\quad
2.001\to1.749,\quad
2.0001\to1.839,\quad
2.00005\to1.865
\]

La contaminazione residua è dovuta al ramo riflesso introdotto da una
condizione ingoing applicata a distanza finita. Al massimo della barriera la
gerarchia si rompe perché \(P\to0\), confermando indipendentemente la necessità
della forma uniforme.

Frequenze WKB3 Dirac, \(M=1,n=0\):

| \(K\) | \(r_{\rm peak}/M\) | \(M\omega_{\rm WKB3}\) |
|---:|---:|---:|
| 1 | 2.420483 | \(0.176452-0.100109i\) |
| 2 | 2.617914 | \(0.378627-0.096542i\) |
| 3 | 2.727147 | \(0.573685-0.096324i\) |
| 4 | 2.790186 | \(0.767194-0.096276i\) |
| 5 | 2.830203 | \(0.960215-0.096256i\) |

Questi valori riproducono Cho (2003). Per \(K=1\), il continued fraction di
Jing dà circa \(0.182963-0.0969825i\), mostrando il limite quantitativo della
WKB3 al modo più basso.

L'idrodinamica fisica del bispinore richiede i bilineari completi. Nella base
\(H=-i\sigma_3D_*+W\sigma_1\):

\[
N^2=J^2+\Sigma_r^2+\Pi_r^2,
\]

\[
\begin{aligned}
N'&=2W\Pi_r-2\omega_IJ,\\
J'&=-2\omega_I N,\\
\Sigma_r'&=2\omega_R\Pi_r,\\
\Pi_r'&=2WN-2\omega_R\Sigma_r .
\end{aligned}
\tag{7.5}
\]

La Madelung di una singola componente non sostituisce la corrente spinoriale.

## 8. Kerr: risultati già derivati

### 8.1 Campo scalare

Con

\[
\Delta=r^2-2Mr+a^2,\qquad
H=r^2+a^2,\qquad
\frac{dr_*}{dr}=\frac H\Delta,\qquad
\Psi=\sqrt H\,R,
\]

l'equazione radiale scalare diventa esattamente

\[
\Psi''+\mathcal Q_K\Psi=0,
\]

\[
\boxed{
\mathcal Q_K=
\left(\omega-\frac{am}{H}\right)^2
-\frac{\Delta\lambda}{H^2}
-\frac1{\sqrt H}\frac{d^2\sqrt H}{dr_*^2}
}.
\tag{8.1}
\]

All'orizzonte:

\[
\mathcal Q_K\longrightarrow(\omega-m\Omega_H)^2.
\]

Nel limite slow Kerr:

\[
\boxed{
\mathcal Q_K=
\omega^2
-f\left[\frac{\ell(\ell+1)}{r^2}+\frac{2M}{r^3}\right]
-\frac{4amM\omega}{r^3}
+O(a^2)
}.
\tag{8.2}
\]

Per \(\omega=\omega_R+i\omega_I\):

\[
\operatorname{Im}\mathcal Q_K=
2\omega_I\left(\omega_R-\frac{2amM}{r^3}\right)+O(a^2).
\tag{8.3}
\]

La rotazione entra dunque sia nel simbolo principale co-rotante sia nella
sorgente del flusso radiale aperto. Il termine di frame dragging non deve
essere rinominato “potenziale di Madelung”.

Le riduzioni (8.1)–(8.3), il limite Schwarzschild e il simbolo all'orizzonte
sono verificati indipendentemente con SymPy e Mathematica.

### 8.2 Dirac proiettivo su Kerr

Per il rapporto \(z=P_+/P_-\) è stata ottenuta l'equazione esatta

\[
\varepsilon z'=w_0+2i\kappa_rz-w_0z^2.
\tag{8.4}
\]

Con \(z=\sum_n\varepsilon^nz_n\):

\[
k=\sqrt{\kappa_r^2-w_0^2},\qquad
z_0=\frac{i(\kappa_r+\tau k)}{w_0},
\]

\[
\boxed{
z_n=
\frac{z_{n-1}'+w_0\sum_{j=1}^{n-1}z_jz_{n-j}}
{-2i\tau k}
}.
\tag{8.5}
\]

Questa è una gerarchia WKB formale all-order che trasporta insieme fase,
ampiezza e polarizzazione. Resta da includere in modo geometrico il trasporto
di Berry e l'accoppiamento con l'autovalore angolare dipendente da
\(a\omega\).

## 9. Vaidya e non conservazione dell'energia

Per la metrica entrante

\[
ds^2=-f(v,r)dv^2+2\,dv\,dr+r^2d\Omega^2,
\qquad f=1-\frac{2M(v)}r,
\]

e \(\Phi=\psi(v,r)Y_{\ell m}/r\), la Klein–Gordon massless si riduce a

\[
2\psi_{vr}+\partial_r(f\psi_r)-U_\ell\psi=0,
\qquad
U_\ell=\frac{f_r}{r}+\frac{\ell(\ell+1)}{r^2}.
\tag{9.1}
\]

Con \(\psi=Ae^{iS/\varepsilon}\):

\[
2S_vS_r+fS_r^2+\varepsilon^2U_\ell+Q_V=0,
\tag{9.2}
\]

\[
\boxed{
Q_V=-\varepsilon^2
\frac{2A_{vr}+\partial_r(fA_r)}A
=Q_r+Q_\times
},
\tag{9.3}
\]

\[
Q_r=-\varepsilon^2\frac{\partial_r(fA_r)}A,\qquad
Q_\times=-2\varepsilon^2\frac{A_{vr}}A,
\]

\[
\partial_v(A^2S_r)+
\partial_r[A^2(S_v+fS_r)]=0.
\tag{9.4}
\]

Le equazioni (9.1)–(9.4) sono verificate simbolicamente.

Il solo \(Q_\times\) non è ancora una misura pulita della memoria, perché
contiene l'inviluppo di decadimento e dipende dalla normalizzazione temporale.
Il candidato normalizzato è

\[
\boxed{
\mathcal M_{vr}=
-2\varepsilon^2\partial_v\partial_r
\ln\left(\frac A{A_{\rm fr}}\right)
},
\tag{9.5}
\]

dove \(A_{\rm fr}(v,r)=A_{\rm QNM}(r;M(v))\) è il profilo a fondo congelato.
La (9.5) è invariante sotto normalizzazioni dipendenti soltanto da \(v\), come
verificato simbolicamente.

Una candidata versione geometrica con il vettore di Kodama è

\[
\mathcal M_K=
-2\varepsilon^2K^a(\nabla^br)\nabla_a\nabla_b
\ln\left(\frac A{A_{\rm fr}}\right).
\tag{9.6}
\]

La sua normalizzazione e la sua reale indipendenza dalla variabile master
restano aperte.

### 9.1 Equazione del “selettore di rotaia”

Per un osservatore \(u^\mu\),

\[
\widehat E[u]=-u^\mu W_\mu,\qquad W_\mu=\nabla_\mu S,
\]

è l'energia localmente misurata, non una quantità automaticamente conservata.
Per un campo vettoriale \(X^\mu\):

\[
W^\nu\nabla_\nu(-X\cdot W)
=-\frac12W^\mu W^\nu\mathcal L_Xg_{\mu\nu}
+\frac12X^\mu\nabla_\mu(\mathcal U+Q_g).
\tag{9.7}
\]

La conservazione deriva da una simmetria di Killing, non dalla scelta di un
osservatore. Imporre \(-u\cdot W=\widehat E\) costante selezionerebbe una
congruenza di osservatori accelerati; non conserverebbe l'energia dell'intero
fluido QNM.

In Vaidya, con \(K=\partial_v\):

\[
\boxed{
\frac{dE_K}{d\lambda}
=-\frac{\dot M(v)}r(W^v)^2
}.
\tag{9.8}
\]

Questo drift è la firma fisica dell'assenza di una Killing temporale e non va
cancellato artificialmente con un vincolo o con una propulsione immaginaria.

## 10. Stato della novità

### 10.1 Cose già note e da non rivendicare

- trasformazione di Madelung e termine \(-A''/A\);
- Riccati, Schwarziana, Prüfer e metodi phase-amplitude;
- WKB di barriera, WKB ad alto ordine, Padé ed exact WKB per QNM;
- spettro e autofunzioni QNM di Pöschl–Teller;
- legame Bohm–Madelung–Ermakov–Pinney–base di Weber;
- QNM e quantum pressure nei buchi neri analoghi;
- idrodinamica covariante del campo di Dirac;
- QNM di Vaidya e loro risposta ritardata.

Precedenti particolarmente vicini:

- Glampedakis–Andersson, phase-amplitude per risonanze di buchi neri:
  https://arxiv.org/abs/gr-qc/0304030
- Kumar, Bohm–Madelung, Ermakov e base di Weber:
  https://arxiv.org/abs/2602.00507
- Cardona–Molina, QNM Pöschl–Teller:
  https://arxiv.org/abs/1711.00479
- Beyer, completezza dei QNM Pöschl–Teller:
  https://arxiv.org/abs/gr-qc/9803034
- Daghigh–Green, quantum potential nei QNM di buchi neri analoghi:
  https://arxiv.org/abs/1411.7066
- Capuano–Santoni–Barausse, Vaidya a \(\dot M\) costante:
  https://arxiv.org/abs/2407.06009

Per il ramo fermionico esiste inoltre una sovrapposizione diretta con
Meza-Domínguez–Matos:
https://arxiv.org/abs/2605.28887. La versione v1 discute idrodinamica chirale,
Schwarzschild, QNM e greybody factors, ma presenta problemi interni già
registrati nella nota Dirac: limite massless formulato con velocità contenenti
\(1/m\), conservazione separata delle correnti chirali nel caso massive,
potenziale radiale non coincidente con i partner standard e condizioni
asintotiche incoerenti con il segno di \(\operatorname{Im}\omega\). Va citata
come precedente diretto, ma ogni formula utile deve essere verificata
indipendentemente e confrontata con Chandrasekhar, Cho e Jing.

La ricerca di anteriorità resta aperta: l'assenza di un risultato nelle
ricerche testuali non dimostra priorità.

### 10.2 Nucleo potenzialmente originale

Non è stato finora identificato un lavoro che combini:

1. il profilo della perturbazione QNM di un buco nero relativistico, non il
   fluido materiale di un modello analogo;
2. un residuo fra curvatura di ampiezza esatta e curvatura della soluzione
   uniforme;
3. l'uso quantitativo di tale residuo per stimare l'errore di troncamento WKB;
4. il trattamento esplicito degli zeri QNM e della non uniformità ai turning
   point;
5. l'estensione coerente statico → Kerr → Vaidya.

Il risultato uniforme di §6 è quindi una candidata novità concreta, ma non è
ancora una rivendicazione di priorità definitiva e non è ancora una
pubblicazione scientifica sottoposta a peer review.

## 11. Verifiche automatiche eseguite

Ambiente usato:

- Python 3.13.2;
- NumPy 2.2.2;
- SciPy 1.17.1;
- SymPy 1.13.3;
- Mathematica kernel 13.3.1 ARM64.

Stato al momento di questo snapshot:

- 20 test Python nel ramo core: superati;
- 12 test Python nei benchmark calculations: superati;
- 19 controlli Wolfram Language in verify_formalism.wl: superati;
- 5 controlli SymPy scalari Kerr: superati;
- 6 controlli SymPy Vaidya: superati.

I controlli Mathematica comprendono:

1. chiusura Madelung 1D esatta;
2. ricorsione WKB pari fino a \(\varepsilon^4\);
3. ordini formali dei residui
   \(\varepsilon^2,\varepsilon^4,\varepsilon^6\);
4. split di Langer del potenziale Regge–Wheeler;
5. modi Pöschl–Teller \(n=0,1\) e zero nodale;
6. mappa parabolico-cilindrica;
7. sviluppo nello strato \(y=\sqrt\varepsilon\eta\);
8. riduzione radiale scalare Kerr;
9. limite Schwarzschild e termine slow Kerr;
10. split Madelung di Vaidya;
11. drift dell'energia di Kodama.

## 12. Riproduzione

Dalla radice della repository:

~~~bash
python3.13 -m unittest \
  calculations/test_static_madelung_benchmark.py \
  calculations/test_poschl_teller_madelung.py \
  calculations/test_order_resolved_madelung.py \
  calculations/test_uniform_madelung_defect.py

python3.13 calculations/static_madelung_benchmark.py --n 0 --n 1 --n 2
python3.13 calculations/poschl_teller_madelung_benchmark.py --n 0 --n 1 --n 2
python3.13 calculations/order_resolved_madelung_benchmark.py --n 0
python3.13 calculations/uniform_madelung_defect.py --n 0 --n 2

python3.13 calculations/kerr_scalar_madelung_symbolic.py
python3.13 calculations/vaidya_madelung_symbolic.py
~~~

I test core importano i moduli relativamente alla cartella core:

~~~bash
cd core
python3.13 -m unittest test_wkb.py test_leaver.py test_dirac_madelung.py
~~~

Verifica Mathematica:

~~~bash
wolframscript -file calculations/verify_formalism.wl
~~~

Nell'ambiente automatizzato Codex il launcher wolframscript si è bloccato
prima di produrre output, anche con TTY e opzione local. Non è un problema di
licenza o del file: lo stesso script passa usando direttamente il kernel:

~~~bash
/Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noprompt -script calculations/verify_formalism.wl
~~~

Sul terminale interattivo dell'utente wolframscript funziona regolarmente.

## 13. Mappa dei file principali

| File | Ruolo |
|---|---|
| core/schwarzschild_wkb.py | WKB1/WKB3 Schwarzschild |
| core/leaver_qnm.py | frequenze indipendenti con frazione continua |
| core/madelung_profile.py | profilo bosonico ampiezza–fase |
| core/dirac_madelung_profile.py | profilo e scaling Dirac |
| calculations/static_madelung_benchmark.py | scansione dei 69 modi |
| calculations/poschl_teller_madelung_benchmark.py | barriera esatta |
| calculations/order_resolved_madelung_benchmark.py | fallimento dei residui locali |
| calculations/uniform_madelung_defect.py | residuo uniforme e fit asintotici |
| calculations/verify_formalism.wl | controllo indipendente Mathematica |
| calculations/kerr_scalar_madelung_symbolic.py | Kerr scalare e slow Kerr |
| calculations/vaidya_madelung_symbolic.py | riduzione e Madelung su Vaidya |
| calculations/dirac_schwarzschild_wkb.py | gerarchia Dirac locale |
| calculations/dirac_kerr_projective_wkb.py | ricorsione proiettiva Kerr |
| research/prior_art.md | letteratura e confine di novità |
| research/references.bib | bibliografia |

## 14. Prossimi conti, in ordine

### Priorità 1: forma uniforme corretta

Partire dalla (6.3) e costruire perturbativamente

\[
A_{\rm unif}
=A_{\rm PC}
+\varepsilon A_1
+\varepsilon^2A_2+\cdots
\]

includendo almeno il termine quartico e poi il sestico. Definire

\[
\Delta Q_M^{[N]}
=-\varepsilon^2
\left(\frac{A''}{A}
-\frac{(A_{\rm unif}^{[N]})''}{A_{\rm unif}^{[N]}}\right).
\]

Criterio di successo:

- il residuo corretto deve mostrare un esponente compatibile con
  \(L^{-5}\) quando viene confrontato con WKB3;
- il risultato deve restare stabile sotto variazioni moderate della finestra;
- il confronto deve essere ripetuto almeno per \(n=0,2\);
- gli zeri devono essere trattati separatamente.

### Priorità 2: controllo fuori Pöschl–Teller

Ripetere il protocollo senza riadattarlo su:

1. Reissner–Nordström oppure Schwarzschild–de Sitter;
2. più variabili master, se disponibili;
3. un secondo metodo indipendente per le frequenze.

Questo è necessario per separare un vero diagnostico da una coincidenza
monoparametrica.

### Priorità 3: slow Kerr scalare numerico

Usare il termine

\[
-\frac{4amM\omega}{r^3}
\]

come perturbazione controllata, calcolare frequenze e profili per \(m\) di
segno opposto e verificare:

- la separazione dispari in \(am\);
- la risposta di \(Q_M\);
- il comportamento in funzione di
  \(\omega_R-m\Omega_H\).

### Priorità 4: Vaidya

1. Implementare \(M(v)=M_0+\dot Mv\) e riprodurre il limite in frequenza
   noto.
2. Costruire \(A_{\rm fr}(v,r)\) con le stesse coordinate e variabile master.
3. Misurare \(\mathcal M_{vr}\) e confrontarlo con
   \(\Delta\omega=\omega_{\rm dyn}-\omega_{\rm fr}[M(v)]\).
4. Confrontare il potere predittivo con \(|\dot M|\) e \(|\ddot M|\).
5. Se \(\mathcal M\) non aggiunge informazione, respingere la congettura di
   memoria Madelung.

### Priorità bibliografica

- seguire le citazioni di Glampedakis–Andersson;
- cercare Milne/Ermakov insieme a Gamow, Siegert e black-hole resonances;
- seguire Kumar 2026;
- cercare esplicitamente curvature di ampiezza come stimatori dell'errore
  semiclassico;
- verificare versioni successive dei lavori exact-WKB 2025–2026 e del
  preprint Dirac direttamente sovrapposto.

## 15. Errori concettuali da non reintrodurre

1. Non chiamare nuova la trasformazione di Madelung.
2. Non interpretare \(A^2\) come probabilità QNM globale senza una
   regolarizzazione e un prodotto bilineare adeguati.
3. Non confondere il potenziale geometrico con il potenziale di Madelung.
4. Non confondere la corrente conservata della PDE con il flusso radiale della
   ODE a frequenza complessa.
5. Non usare la serie WKB locale sui turning point coalescenti.
6. Non inferire l'ordine asintotico da una sola correlazione elevata.
7. Non nascondere gli zeri di \(A\) con una regolarizzazione numerica
   arbitraria.
8. Non imporre \(-u\cdot W=\mathrm{costante}\) come se fosse una legge di
   conservazione: è una scelta di osservatori accelerati.
9. In Kerr, non chiamare Madelung il termine di frame dragging.
10. In Vaidya, non dichiarare covariante \(A_{vr}/A\) senza specificare
    coordinate e variabile master.

## 16. Stato Git al momento della ripresa

Sequenza dei commit scientifici:

| commit | contenuto |
|---|---|
| dc63ed4 | prima versione pubblica |
| dab10a6 | fondamenti covarianti Kerr–Vaidya |
| d16fabf | energia di osservatore e bilancio di Kodama |
| 179d07c | benchmark statico e controlli Mathematica |
| 6103c8f | ponte scalare Kerr |
| bdf198f | barriera esatta, nodi e residui locali |
| 09ef057 | diagnostico uniforme di Madelung |

La repository è pubblica, il branch predefinito è main e, prima
dell'aggiunta di questo documento, il branch locale era sincronizzato con
origin/main.
