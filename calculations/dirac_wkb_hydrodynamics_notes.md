# Dirac WKB idrodinamica su Schwarzschild, con estensione a Kerr

## Nota di calcolo verificabile

**Data:** 26 agosto 2026  
**Convenzioni:** unità \(G=c=\hbar=1\), firma \((-+++ )\), dipendenza temporale \(e^{-i\omega t}\).  
**Campo di prova:** Dirac massless. Il caso massive viene rinviato perché la riduzione a potenziali partner introduce già una coordinata e un potenziale dipendenti da \(\omega\).

## 1. Risultato in breve

Su Schwarzschild il problema radiale massless è esattamente un sistema di Dirac in una dimensione con massa spaziale

\[
W(r)=\kappa\frac{\sqrt{f(r)}}{r},\qquad f(r)=1-\frac{2M}{r}.
\]

Le due componenti soddisfano equazioni di Schrödinger con potenziali di Darboux

\[
V_\pm=W^2\pm\frac{dW}{dr_*}.
\]

Il punto strutturalmente nuovo emerge quando il parametro semiclassico viene inserito **prima** di disaccoppiare il sistema. Con

\[
K=|\kappa|,\qquad \varepsilon=K^{-1},\qquad \omega=K\Omega,
\]

le componenti scalari obbediscono a

\[
\varepsilon^2 Z_\sigma''+
\left[\Omega^2-h^2-\sigma s\,\varepsilon h'\right]Z_\sigma=0,
\qquad
h=\frac{\sqrt f}{r},\quad s=\operatorname{sgn}\kappa.
\]

Quindi la gerarchia efficace non comincia dal potenziale di Madelung:

\[
\boxed{
\Omega^2=P^2+h^2
+\underbrace{\sigma s\,\varepsilon h'}_{\text{spin connection}}
+\underbrace{Q_{\rm M}}_{O(\varepsilon^2)}
}
\]

La correzione di spin è subprincipale, di ordine \(\varepsilon\); il funzionale scalare di Madelung comincia a ordine \(\varepsilon^2\). Ai successivi ordini compaiono termini misti spin-gradiente. Questo fornisce una risposta concreta alla domanda del progetto: per Dirac la famiglia efficace esiste localmente, ma **non è una successione di soli potenziali di Madelung**.

## 2. Schwarzschild: sistema radiale e potenziali partner

Definiamo

\[
ds^2=-fdt^2+f^{-1}dr^2+r^2d\Omega_2^2,
\qquad
r_*=r+2M\log\left(\frac{r}{2M}-1\right),
\qquad
\frac{d}{dr_*}=f\frac{d}{dr}.
\]

Per \(K=|\kappa|=j+\tfrac12=1,2,\ldots\), una scelta di fasi coerente dà

\[
(D_*-W)F=i\omega G,
\qquad
(D_*+W)G=i\omega F,
\qquad
W=K\frac{\sqrt f}{r}.
\tag{2.1}
\]

Applicando gli operatori nell'ordine opposto,

\[
F''+(\omega^2-V_+)F=0,
\qquad
G''+(\omega^2-V_-)G=0,
\tag{2.2}
\]

\[
V_\pm=W^2\pm W',
\qquad
W'=K\frac{\sqrt f(3M-r)}{r^3}.
\]

Pertanto

\[
\boxed{
V_\pm(r)=K^2\frac{f}{r^2}
\pm K\frac{\sqrt f(3M-r)}{r^3}
}
\tag{2.3}
\]

e, ponendo \(\Delta=r(r-2M)\),

\[
V_+=\frac{K\sqrt\Delta}{r^4}
\left[K\sqrt\Delta-(r-3M)\right].
\]

Questa è precisamente la forma usata da Cho. Il cambio \(\kappa\mapsto-\kappa\), oppure una diversa fase delle componenti, scambia le etichette \(+\) e \(-\) senza cambiare la fisica.

Poiché \(W\to0\) sia per \(r_*\to-\infty\) sia per \(r_*\to+\infty\), l'intertwiner di Darboux preserva, per \(\omega\ne0\), le condizioni QNM. I due spettri esatti sono isospettrali; una WKB troncata può tuttavia rompere leggermente questa equivalenza.

Con la convenzione temporale adottata,

\[
Z\sim e^{-i\omega r_*}\quad(r_*\to-\infty),
\qquad
Z\sim e^{+i\omega r_*}\quad(r_*\to+\infty).
\tag{2.4}
\]

Per \(\operatorname{Im}\omega<0\), le funzioni QNM crescono spazialmente agli estremi e non sono stati \(L^2\). Una loro "densità" idrodinamica è quindi locale oppure appartiene a un sistema aperto.

## 3. Due conteggi semiclassici non equivalenti

### 3.1 Parametro inserito dopo il disaccoppiamento

Si può formalmente scrivere

\[
\varepsilon^2 Z''+[\omega^2-V_\pm]Z=0.
\]

In questo conteggio l'intero \(V_\pm\), incluso \(\pm W'\), è trattato come potenziale di ordine principale. È una WKB scalare legittima a modo angolare fissato, ma nasconde l'origine spinoriale del termine derivativo.

### 3.2 Parametro inserito nel sistema di Dirac

Scriviamo \(\kappa=Ks\), \(K\gg1\), \(\omega=K\Omega\) e

\[
h=\frac{\sqrt f}{r}.
\]

Il sistema (2.1) equivale alla Hamiltoniana radiale

\[
\left[-i\varepsilon\sigma_xD_*+s h\,\sigma_y\right]\Psi
=\Omega\Psi,
\qquad \Psi=(F,G)^T.
\tag{3.1}
\]

Quadrando si ottiene

\[
\varepsilon^2Z_\sigma''+
\left[q_0+\varepsilon q_1\right]Z_\sigma=0,
\tag{3.2}
\]

\[
q_0=\Omega^2-h^2,
\qquad
q_1=-\sigma s\,h',
\qquad \sigma=\pm1.
\tag{3.3}
\]

Questo è il conteggio fisicamente naturale per un limite eikonale di Dirac. Mostra che \(h^2\) è il simbolo principale e \(\sigma s h'\) è il simbolo subprincipale di spin.

## 4. Chiusura WKB-Madelung della componente scalare

In una regione oscillatoria reale, scegliamo

\[
Z=P^{-1/2}\exp\left(\frac{i}{\varepsilon}\int^xP(y)\,dy\right).
\tag{4.1}
\]

La sostituzione in (3.2) è esatta e produce

\[
q_0+\varepsilon q_1=P^2+Q_{\rm M},
\qquad
Q_{\rm M}=-\varepsilon^2
\frac{(P^{-1/2})''}{P^{-1/2}}.
\tag{4.2}
\]

Equivalentemente,

\[
P^2=q_0+\varepsilon q_1+\varepsilon^2\mathcal F[P],
\qquad
\mathcal F[P]=\frac34\left(\frac{P'}P\right)^2-
\frac12\frac{P''}{P}.
\tag{4.3}
\]

Ponendo \(P=\sum_{n\ge0}\varepsilon^nP_n\), i primi coefficienti sono

\[
P_0=\sqrt{q_0},
\qquad
P_1=\frac{q_1}{2P_0}
=-\frac{\sigma s\,h'}{2\sqrt{q_0}},
\tag{4.4}
\]

\[
P_2=
\frac{5(q_0')^2}{32q_0^{5/2}}
-\frac{q_0''}{8q_0^{3/2}}
-\frac{q_1^2}{8q_0^{3/2}}.
\tag{4.5}
\]

Il coefficiente iniziale del potenziale di Madelung è

\[
Q_{\rm M}=\varepsilon^2Q_2+O(\varepsilon^3),
\qquad
\boxed{
Q_2=\frac{q_0''}{4q_0}-\frac{5(q_0')^2}{16q_0^2}
}.
\tag{4.6}
\]

Poiché \(q_1\neq0\), il funzionale di Madelung completo non resta una serie
puramente pari. Il primo termine misto spin-gradiente è

\[
Q_{\rm M}=\varepsilon^2Q_2+\varepsilon^3Q_3+\cdots,
\]

\[
\boxed{
Q_3=
\frac{q_1''}{4q_0}
-\frac{5q_0'q_1'}{8q_0^2}
-\frac{q_1q_0''}{4q_0^2}
+\frac{5q_1(q_0')^2}{8q_0^3}
}.
\tag{4.7}
\]

La ricorsione generale è compatta. Se \(q=q_0+\varepsilon q_1\) e \(q_n=0\) per \(n\ge2\), allora

\[
2P_0P_n=q_n+[\varepsilon^{n-2}]\mathcal F[P]
-\sum_{j=1}^{n-1}P_jP_{n-j},
\qquad n\ge1,
\tag{4.8}
\]

dove il termine con \(\mathcal F\) è assente per \(n<2\). Questa relazione genera sia le correzioni puramente Madelung sia quelle miste con il simbolo di spin \(q_1\).

Nel caso Schwarzschild,

\[
h'=\frac{\sqrt f(3M-r)}{r^3},
\]

\[
q_0=\frac{2M+\Omega^2r^3-r}{r^3},
\qquad
q_0'=\frac{2(r-3M)(r-2M)}{r^5},
\]

\[
q_0''=-\frac{2(r-2M)(30M^2-20Mr+3r^2)}{r^7},
\tag{4.9}
\]

dove tutti i primi indicano \(D_*=fD_r\). Il primo contributo di spin alla fase si integra localmente:

\[
\int P_1dr_*=-\frac{\sigma s}{2}
\arcsin\left(\frac{h}{\Omega}\right),
\tag{4.10}
\]

per un ramo fissato e \(q_0>0\). È una fase di polarizzazione/spin, non una pressione quantistica scalare.

## 5. Idrodinamica del bispinore: la corrente fisica

La Madelung di una singola componente decoupled non coincide automaticamente con l'idrodinamica della corrente di Dirac. Dalla versione dipendente dal tempo di (3.1),

\[
i\varepsilon\partial_t\Psi=
\left[-i\varepsilon\sigma_xD_*+s h\sigma_y\right]\Psi,
\]

seguono esattamente

\[
\rho=\Psi^\dagger\Psi,
\qquad
j=\Psi^\dagger\sigma_x\Psi,
\qquad
\partial_t\rho+\partial_{r_*}j=0.
\tag{5.1}
\]

All'ordine eikonale, con \(\Psi=\sqrt\rho\,\chi\,e^{iS/\varepsilon}\) e \(p=S'\),

\[
(p\sigma_x+s h\sigma_y)\chi=\Omega\chi,
\qquad
\Omega^2=p^2+h^2.
\tag{5.2}
\]

Per il ramo positivo, il vettore di Bloch locale è

\[
\mathbf n=\left(\frac p\Omega,\frac{s h}{\Omega},0\right),
\]

e quindi

\[
j=\rho\frac p\Omega,
\qquad
j'=0\quad\Longrightarrow\quad \rho\propto\frac{\Omega}{p}
\tag{5.3}
\]

per scattering stazionario reale. Questo fattore di trasporto è diverso dall'ampiezza scalare \(|Z|^2\propto P^{-1}\), perché include la polarizzazione del bispinore.

Per un modo QNM \(e^{-i\Omega t/\varepsilon}\),

\[
j_0'=-\frac{2\operatorname{Im}\Omega}{\varepsilon}\rho_0
=-2\operatorname{Im}\omega\,\rho_0.
\tag{5.4}
\]

Il profilo radiale QNM si comporta dunque come un flusso aperto con sorgente/pozzo nella descrizione stazionaria; non come un fluido conservativo normalizzabile.

Una rotazione unitaria costante porta la stessa Hamiltoniana non scalata alla forma

\[
H=-i\sigma_3D_*+W\sigma_1.
\]

Questa base rende visibile un sistema idrodinamico radiale reale ed esatto. Definiamo

\[
N=\Psi^\dagger\Psi,\qquad
J=\Psi^\dagger\sigma_3\Psi,\qquad
\Sigma_r=\Psi^\dagger\sigma_1\Psi,\qquad
\Pi_r=\Psi^\dagger\sigma_2\Psi.
\]

L'identità pura di Bloch e le equazioni radiali sono

\[
N^2=J^2+\Sigma_r^2+\Pi_r^2,
\tag{5.5}
\]

\[
\boxed{
\begin{aligned}
N'&=2W\Pi_r-2\omega_IJ,\\
J'&=-2\omega_I N,\\
\Sigma_r'&=2\omega_R\Pi_r,\\
\Pi_r'&=2WN-2\omega_R\Sigma_r .
\end{aligned}}
\tag{5.6}
\]

Per \(\omega\) reale, \(J=J_0\) è costante. Se \(J_0>0\), la parametrizzazione

\[
N=J_0\cosh\eta,
\quad
\Sigma_r=J_0\sinh\eta\cos\beta,
\quad
\Pi_r=J_0\sinh\eta\sin\beta
\]

riduce (5.6) a

\[
\boxed{
\eta'=2W\sin\beta,
\qquad
\beta'=-2\omega+2W\coth\eta\cos\beta .
}
\tag{5.7}
\]

Questo sistema è un candidato più fedele per l'"idrodinamica radiale di Dirac": \(N,J\) descrivono densità e flusso, mentre \(\Sigma_r,\Pi_r\), oppure \(\eta,\beta\), conservano la coerenza interna e la polarizzazione che una singola Madelung scalare perde. Nella base che diagonalizza il quadrato,

\[
Z_\sigma=\frac{F-\sigma iG}{\sqrt2},
\qquad
|Z_\sigma|^2=\frac{N+\sigma\Pi_r}{2}.
\tag{5.8}
\]

Occorrono dunque entrambi i partner e il vincolo del primo ordine per ricostruire il bispinore.

| Oggetto | Equazione | Significato |
|---|---|---|
| \(P^{-1}\) | componente scalare disaccoppiata | ampiezza WKB locale |
| \(\rho=\Psi^\dagger\Psi\) | sistema di Dirac 2x2 | densità spinoriale |
| \(j=\Psi^\dagger\sigma_x\Psi\) | continuità esatta | flusso fisico radiale |
| \(Q_{\rm M}\) | chiusura scalare | correzione di ampiezza, non intera dinamica di spin |

## 6. Barriera, photon sphere e controllo numerico WKB

Per

\[
V_\sigma=K^2h^2+\sigma sK h',
\]

il massimo della barriera ha l'espansione

\[
r_{\rm peak}=3M-\frac{\sigma s\sqrt3\,M}{2K}+O(K^{-2}),
\tag{6.1}
\]

\[
V_{\rm peak}=\frac{K^2}{27M^2}+\frac{1}{108M^2}+O(K^{-1}).
\tag{6.2}
\]

Inoltre

\[
\sqrt{-2V_0''}=\frac{2K}{27M^2}+O(K^0).
\]

La formula di Schutz-Will al primo ordine dà, per \(\alpha=n+\tfrac12\),

\[
\boxed{
M\omega_{Kn}=\frac{K-i\alpha}{3\sqrt3}+O(K^{-1})
}.
\tag{6.3}
\]

Il coefficiente completo di ordine \(K^{-1}\) non si può ricavare dal solo massimo: riceve anche le correzioni uniformi WKB.

Come controllo indipendente, lo script `dirac_schwarzschild_wkb3.py` applica la formula Iyer-Will/Cho di terzo ordine, iterando correttamente \(D_*=fD_r\). Per \(M=1\), \(n=0\):

| \(K\) | \(r_{\rm peak}/M\) | \(M\omega_{\rm WKB3}\) |
|---:|---:|---:|
| 1 | 2.420483 | \(0.176452-0.100109i\) |
| 2 | 2.617914 | \(0.378627-0.096542i\) |
| 3 | 2.727147 | \(0.573685-0.096324i\) |
| 4 | 2.790186 | \(0.767194-0.096276i\) |
| 5 | 2.830203 | \(0.960215-0.096256i\) |

I valori riproducono la tabella WKB3 di Cho. Il continued fraction di Jing dà per \(K=1\), nelle stesse unità, \(M\omega\simeq0.182963-0.0969825i\): il modo più basso mostra chiaramente il limite quantitativo della WKB3.

**Avvertenza decisiva.** I coefficienti (4.5)-(4.6) divergono quando \(q_0\to0\). Non vanno valutati ingenuamente sui turning point coalescenti presso la barriera QNM. Per collegarli alle \(\Lambda_n\) di Iyer-Will serve una forma normale uniforme.

## 7. Kerr: il passo successivo non è banale

Per Kerr,

\[
\Delta=r^2-2Mr+a^2,
\qquad
\frac{dr_*}{dr}=\frac{r^2+a^2}{\Delta},
\]

e definiamo

\[
\mathcal K(r)=(r^2+a^2)\omega-am,
\qquad
v(r)=\frac{am}{r^2+a^2},
\qquad
w(r)=\frac{\lambda\sqrt\Delta}{r^2+a^2}.
\tag{7.1}
\]

Con \(\Psi\propto e^{-i\omega t+im\phi}\) e
\(P_-=R_{-1/2}\), \(P_+=\sqrt\Delta\,R_{+1/2}\), una convenzione
di Chandrasekhar coerente è

\[
(D_*+i[\omega-v])P_-=wP_+,
\qquad
(D_*-i[\omega-v])P_+=wP_-.
\tag{7.2}
\]

Con \(P=(P_-,P_+)^T\), essa equivale a

\[
\boxed{
\left[vI+i\sigma_zD_*+w\sigma_y\right]P=\omega P
}.
\tag{7.3}
\]

Nel limite eikonale poniamo

\[
m=K\mu,
\quad \lambda=K\Lambda(a\Omega,\mu),
\quad \omega=K\Omega,
\quad \varepsilon=K^{-1},
\]

\[
v_0=\frac{a\mu}{r^2+a^2},
\qquad
w_0=\frac{\Lambda\sqrt\Delta}{r^2+a^2}.
\]

Il simbolo principale è

\[
H_0=v_0I-p\sigma_z+w_0\sigma_y,
\]

perciò

\[
\boxed{
(\Omega-v_0)^2=p^2+w_0^2,
\qquad
j=-P^\dagger\sigma_zP
=\rho\frac{p}{\Omega-v_0}
}.
\tag{7.4}
\]

Il termine \(v_0\) ha una lettura idrodinamica pulita: è il trascinamento azimutale. All'orizzonte,

\[
\omega-v(r_+)\longrightarrow\omega-m\Omega_H,
\]

mentre \(w\to0\); all'infinito \(v,w\to0\). La dispersione interpola quindi tra il momento co-rotante vicino all'orizzonte e quello asintotico.

La struttura matriciale ammette già una ricorsione all-order senza introdurre
una coordinata dipendente dall'autovalore. Posto

\[
\kappa_r(r)=\Omega-\frac{a\mu}{r^2+a^2},
\qquad
z=\frac{P_+}{P_-},
\]

il rapporto proiettivo soddisfa esattamente

\[
\boxed{
\varepsilon z'=w_0+2i\kappa_r z-w_0z^2
}.
\tag{7.5}
\]

Con \(z=\sum_{n\ge0}\varepsilon^nz_n\),

\[
k=\sqrt{\kappa_r^2-w_0^2},
\qquad
z_0=\frac{i(\kappa_r+\tau k)}{w_0},
\qquad \tau=\pm1,
\tag{7.6}
\]

e, per \(n\ge1\),

\[
\boxed{
z_n=
\frac{z_{n-1}'+w_0\sum_{j=1}^{n-1}z_jz_{n-j}}
{-2i\tau k}
}.
\tag{7.7}
\]

Infine

\[
(\log P_-)'=\frac{i\tau k}{\varepsilon}
+w_0z_1+\varepsilon w_0z_2+\cdots .
\tag{7.8}
\]

Questa è una gerarchia WKB concreta, compatta e polarizzata. I \(z_n\)
trasportano insieme fase, ampiezza e polarizzazione; ridurli a un singolo
\(Q_{\rm eff}\) richiede la scelta di componente, gauge spinoriale e ramo.
La domanda corretta su Kerr è quindi se esista una chiusura **matriciale o
proiettiva** geometricamente naturale, non se ogni ordine sia un nuovo
potenziale scalare.

Esiste una seconda via, utile ma meno geometrica. Definendo una coordinata spettrale

\[
\frac{d\widehat r_*}{dr}=\frac{\mathcal K(r)}{\omega\Delta},
\]

si ottiene un superpotenziale

\[
\widehat W(r;\omega)=
\frac{\lambda\omega\sqrt\Delta}{\mathcal K(r)},
\qquad
\widehat V_\pm=\widehat W^2
\pm\frac{d\widehat W}{d\widehat r_*}.
\tag{7.9}
\]

Questa forma presenta però tre ostacoli reali:

1. \(\widehat r_*\) e \(\widehat W\) dipendono dall'autovalore complesso \(\omega\);
2. \(\lambda\) dipende dal problema angolare e da \(a\omega\);
3. \(\mathcal K=0\) produce una singolarità di coordinata/superpotenziale.

Per costruire una gerarchia idrodinamica Kerr conviene quindi partire dalla Hamiltoniana matriciale (7.3), diagonalizzarne il simbolo principale e trattare spin connection e Berry transport come termini subprincipali, senza forzare prematuramente una singola equazione scalare.

L'autovalore \(\lambda\) proviene inoltre dall'equazione angolare
sferoidale e dipende da \(a\omega\); per una QNM è in generale complesso.
Il problema spettrale è dunque una coppia non lineare radiale-angolare.
La WKB scalare su Kerr è già stata studiata fino al sesto ordine: la novità
plausibile è proprio la gerarchia matriciale/proiettiva. Infine, per il campo
di Dirac neutro non c'è superradianza classica, pur restando essenziali
frame dragging e dinamica quasi-estremale.

## 8. Letteratura direttamente rilevante e controllo critico

### 8.1 Precedenti solidi

- Takabayasi ha sviluppato la rappresentazione idrodinamica del campo di Dirac in termini di bilineari, densità, momento, spin e pseudoscalare. Questa è la cornice concettuale più robusta per non confondere la fase di una componente con la corrente del bispinore.
- Matos, Gallegos e Chavanis (2022) hanno pubblicato una trasformazione generalizzata di Madelung per fermioni di Dirac e Weyl in spazio-tempo curvo, con equazioni di continuità e bilancio energetico. Il progetto non può quindi rivendicare come nuova l'idea generale di una idrodinamica di Dirac curva.
- Cho (2003) e Jing (2005) forniscono rispettivamente la barriera WKB e un benchmark continued-fraction per i QNM di Dirac su Schwarzschild.
- Chandrasekhar e Page forniscono la separazione di Dirac su Kerr/Kerr-Newman.

Lo spazio originale rimasto è più preciso: **ordinamento semiclassico del sistema di Dirac, ricorsione WKB-Madelung con simbolo subprincipale di spin, compatibilità con la corrente spinoriale e uniformizzazione QNM; quindi estensione matriciale a Kerr.**

Anche il lavoro del 2022 va usato come precedente concettuale, non come
scorciatoia tecnica per il caso massless: la sua velocità di fase
\(mv_\mu=\nabla_\mu S+qA_\mu\) è singolare per \(m=0\), e la formula generale
ottenuta quadrando Dirac non mostra esplicitamente il termine \(R/4\) della
formula standard di Lichnerowicz. Questi punti non annullano il valore della
rappresentazione proposta, ma impongono di controllare le identità covarianti
con le convenzioni standard prima di trasferirle ai QNM.

### 8.2 Preprint 2026 direttamente sovrapposto

Meza-Domínguez e Matos, arXiv:2605.28887v1, dichiarano una formulazione chiral-idrodinamica covariante e la applicano a Schwarzschild, QNM e greybody factors. È una sovrapposizione diretta e deve essere citata. Tuttavia la versione v1, così come scritta, non può essere usata come base senza verifiche indipendenti:

1. le velocità (3.10)-(3.11) sono definite con fattori \(1/m\), ma varie applicazioni e la tabella QNM includono il limite massless;
2. l'eq. (4.20) postula conservazione separata delle correnti sinistra/destra anche per \(m\ne0\). Definendo \(\Pi_{\rm cov}=i\bar\psi\gamma^5\psi\), per il Dirac massive e \(J_{R,L}=(J\pm J_5)/2\),
   \[
   \nabla_\mu J_R^\mu=+\frac{m}{\hbar}\Pi_{\rm cov},
   \qquad
   \nabla_\mu J_L^\mu=-\frac{m}{\hbar}\Pi_{\rm cov}
   \]
   a livello classico, salvo convenzioni di segno; si conserva la somma, non ciascuna chiralità;
3. il potenziale massless stampato nell'eq. (7.36),
   \[
   V_1=f\left(\frac{\kappa^2}{r^2}
   -\frac{\kappa M}{r^3\sqrt f}\right),
   \]
   non coincide con nessuno dei partner standard (2.3), che contengono \(\pm\kappa\sqrt f(3M-r)/r^3\);
   un controllo numerico per \(M=\kappa=1\) dà
   \(r_{\max}=3.85313\), \(V_{\max}=0.0202712\) per il potenziale
   stampato, contro \(r_{\max}=2.42048\), \(V_{\max}=0.0466839\)
   per \(V_+\);
4. l'eq. (10.41) contiene letteralmente un punto interrogativo, mentre la (10.43) divide per \(\sqrt{f(r_s)}=0\);
5. il testo afferma che il ramo outgoing \(e^{+i\omega r_*}\) decade all'infinito per \(\operatorname{Im}\omega<0\), ma il suo modulo cresce;
6. l'eq. (11.47) uguaglia \((r_s/2M)^2(\omega^2/\omega^2)\) a \(16\pi M^2\omega^2\): i due membri non sono uguali e il fattore di matching non è derivato nella formula mostrata.

Questi sono controlli interni sulle equazioni della v1, non un giudizio definitivo sul programma degli autori. La versione può essere corretta in futuro. La scelta prudente è citare il preprint come precedente diretto, riprodurre indipendentemente ogni risultato e usare come baseline i potenziali standard (2.3).

## 9. Cosa è dimostrato e cosa resta aperto

**Dimostrato in questa nota**

- disaccoppiamento Schwarzschild e forma esplicita dei partner \(V_\pm\);
- ordinamento \(O(1)\), \(O(\varepsilon)\), \(O(\varepsilon^2)\) di fondo, spin e Madelung;
- ricorsione locale per \(P_n\) con \(q=q_0+\varepsilon q_1\);
- continuità esatta della Hamiltoniana radiale e sua modifica per \(\omega\) complessa;
- espansione del massimo a grande \(K\) e riproduzione numerica della WKB3 di Cho.
- ricorsione proiettiva radiale Kerr (7.5)-(7.8) a tutti gli ordini formali.

**Formulato ma non ancora dimostrato all-order**

- una gerarchia covariante in termini di bilineari completi di Dirac;
- l'equivalenza, ordine per ordine, tra i due partner dopo troncamento e risommazione;
- il ponte uniforme tra i coefficienti locali \(P_n,Q_n\) e i \(\Lambda_n\) di barriera;
- la chiusura covariante Kerr comprendente connessione di Berry, equazione angolare e \(a\omega\) complesso.

## 10. Prossimi conti consigliati

1. Derivare la forma normale uniforme di (3.2) vicino al massimo e verificare esplicitamente \(\Lambda_2\) e \(\Lambda_3\) mantenendo \(q_1\).
2. Ricostruire \(J^\mu=\bar\Psi\gamma^\mu\Psi\) dalle soluzioni WKB \(F,G\) fino a \(O(\varepsilon^2)\).
3. Testare l'isospectralità \(V_+\leftrightarrow V_-\) a ogni ordine troncato e dopo Padé.
4. Su Kerr, sviluppare il trasporto adiabatico delle autofibre di \(H_0\), calcolando connessione di Berry e primo termine subprincipale.
5. Confrontare Schwarzschild \(K=1,2\) con continued fractions prima di attribuire significato fisico alle correzioni.

## Riferimenti primari essenziali

1. S. Chandrasekhar, *The solution of Dirac's equation in Kerr geometry*, Proc. R. Soc. A **349** (1976), DOI: 10.1098/rspa.1976.0090.
2. D. N. Page, *Dirac equation around a charged, rotating black hole*, Phys. Rev. D **14** (1976), DOI: 10.1103/PhysRevD.14.1509.
3. H. T. Cho, *Dirac quasinormal modes in Schwarzschild black hole spacetimes*, Phys. Rev. D **68** (2003), arXiv:gr-qc/0303078.
4. J. Jing, *Dirac quasinormal modes of Schwarzschild black hole*, Phys. Rev. D **71** (2005), arXiv:gr-qc/0502023.
5. B. F. Schutz and C. M. Will, *Black Hole Normal Modes: A Semianalytic Approach*, ApJL **291** (1985), DOI: 10.1086/184453.
6. S. Iyer and C. M. Will, *Black-hole normal modes: A WKB approach. I*, Phys. Rev. D **35** (1987), DOI: 10.1103/PhysRevD.35.3621.
7. T. Takabayasi, *Relativistic Hydrodynamics of the Dirac Matter*, Prog. Theor. Phys. Suppl. **4** (1957), DOI: 10.1143/PTPS.4.2.
8. T. Matos, O. Gallegos and P.-H. Chavanis, *Hydrodynamic representation and energy balance for Dirac and Weyl fermions in curved space-times*, Eur. Phys. J. C **82**, 898 (2022), arXiv:2106.15803.
9. J. Meza-Domínguez and T. Matos, *A Covariant Chiral-Hydrodynamic Formulation of the Dirac Equation in Curved Spacetime*, arXiv:2605.28887v1 (2026).
10. M. A. Oancea and A. Kumar, *Semiclassical analysis of Dirac fields on curved spacetime*, Phys. Rev. D **107**, 044029 (2023), DOI: 10.1103/PhysRevD.107.044029.
11. N. Carlson, A. S. Cornell and B. Jordan, *Fermion quasinormal modes of the Kerr black hole*, arXiv:1201.3267.
12. W. G. Unruh, *Separability of the Neutrino Equations in a Kerr Background*, Phys. Rev. Lett. **31**, 1265 (1973), DOI: 10.1103/PhysRevLett.31.1265.

## File riproducibili

- `dirac_schwarzschild_wkb.py`: identità simboliche, gerarchia \(P_0,P_1,P_2,Q_2\), massimo eikonale.
- `dirac_schwarzschild_wkb3.py`: massimo esatto del partner \(V_+\) e frequenze WKB di terzo ordine.
- `dirac_kerr_projective_wkb.py`: verifica simbolica della ricorsione proiettiva Kerr fino a \(O(\varepsilon^2)\).
