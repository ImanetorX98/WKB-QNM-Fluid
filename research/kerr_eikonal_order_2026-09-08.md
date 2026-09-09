# La rotazione produce un termine di ordine ε nel settore bosonico

**Data:** 8 settembre 2026.
**Natura:** risultato **misurato**, con la parte di classificazione indicata come
derivata. L'espansione eikonale dell'autovalore sferoidale è nota in letteratura
(Dolan 2010; Yang *et al.* 2012): **non è questa la rivendicazione**.
**Riproduzione:** `python3.13 calculations/kerr_eikonal_order_test.py`;
test in `calculations/test_kerr_eikonal_order.py` (7 test).

## 1. La domanda

Il criterio di ordinamento del manoscritto afferma: pendenza 2 in \(\varepsilon\)
identifica il regime puramente di Madelung, pendenza 1 segnala un contributo
geometrico che lo precede. Nel settore bosonico statico la pendenza è 2 perché
la sostituzione di Langer dà \(\ell(\ell+1)=L^2-\tfrac14\), senza termine di
ordine \(L\).

In Kerr l'autovalore sferoidale \(A_{\ell m}(c)\), con \(c=a\omega\), sostituisce
\(\ell(\ell+1)\). Se la sua espansione eikonale

\[
\frac{A}{L^2}=A_0+\frac{A_1}{L}+\frac{A_2}{L^2}+\dots,
\qquad \mu=\frac mL,\quad \hat c=\frac cL \ \text{fissati},
\]

ha \(A_1\neq0\), la rotazione produce un termine di ordine \(\varepsilon\) e la
pendenza 1 **non è una firma fermionica**.

## 2. Metodo

Problema agli autovalori sferoidale risolto come matrice nella base delle
armoniche sferiche: \(M=\operatorname{diag}[\ell(\ell+1)]-c^2\langle\cos^2\theta\rangle\).

Le convenzioni alternative dell'equazione angolare differiscono per termini
\(a^2\omega^2=L^2\hat c^2\) e \(2am\omega=2L^2\mu\hat c\), entrambi puramente di
ordine \(L^2\): **spostano \(A_0\), non \(A_1\)**. Il risultato è quindi
indipendente dalla convenzione.

### 2.1 Un errore di metodo e la sua correzione

Il primo tentativo tracciava il modo per prossimità di autovalore lungo una
continuazione da \(c=0\). Con \(c\) grande e complessa il tracciamento **salta
modo**: per \(\mu=0.3,\ \hat c=0.9\) la successione di \(A/L^2\) risultava
0.646, 0.711, 0.956, 1.033, 0.986, 0.995 — erratica, e il fit dava
\(A_1=-12.1\) con residuo \(10^{-1}\).

La troncatura della matrice non c'entrava: a `pad` crescente l'autovalore è
stabile a \(10^{-8}\). La correzione è in due fasi: per \(c\) reale la matrice è
reale simmetrica, gli autovalori non si incrociano e l'indice \(\ell-|m|\) della
lista ordinata identifica il modo senza ambiguità; la parte immaginaria si
continua poi seguendo l'**autovettore** anziché l'autovalore. Le successioni
diventano lisce e monotone, e i residui scendono a \(10^{-9}\)–\(10^{-7}\).

## 3. Risultati

**Controllo, \(\hat c=0\):**
\(A_1=-1.4\times10^{-14}\) e \(A_2=-0.25000000\), cioè esattamente il termine di
Langer. Il metodo riproduce il caso statico.

**Con rotazione, \(\hat c\) complessa:**

| \(\mu\) | \(\hat c\) | \(A_1\) | \(\delta_{\rm Langer}=A_1/(2A_0)\) | residuo |
|---|---|---|---|---|
| 0.5 | \(0.4-0.04i\) | \(-0.01913+0.00372i\) | \(-0.01014+0.00211i\) | \(5\times10^{-9}\) |
| 0.8 | \(0.6-0.06i\) | \(-0.11408+0.02280i\) | \(-0.06095+0.01316i\) | \(4\times10^{-8}\) |
| 0.3 | \(0.9-0.09i\) | \(-0.02878+0.00465i\) | \(-0.02162+0.00587i\) | \(2\times10^{-7}\) |
| 0.9 | \(0.75-0.05i\) | \(-0.24143+0.03379i\) | \(-0.12821+0.01921i\) | \(2\times10^{-7}\) |
| 0.6 | \(0.5-0.15i\) | \(-0.04015+0.02499i\) | \(-0.02089+0.01458i\) | \(1\times10^{-8}\) |

**Non rimovibilità.** Ridefinendo \(L\to L+\delta\) si ha
\(A_1\to A_1-2\delta A_0\). Una costante unica assorbirebbe \(A_1\) solo se
\(\delta=A_1/(2A_0)\) fosse comune ai modi. I valori richiesti vanno da
\(-0.010\) a \(-0.128\), con dispersione del 162% attorno alla media: **nessuna
sostituzione di Langer elimina il termine.** È lo stesso test che nell'Appendice
B del manoscritto separa l'artefatto di parametrizzazione dal termine di spin
genuino di Dirac.

**Limite fisico.** A overtone fissato \(\operatorname{Im}\hat c=O(\varepsilon)\).
Facendo tendere \(\operatorname{Im}\hat c\to0\), \(\operatorname{Re}A_1\) resta
stabile a \(-0.0193\) e \(\operatorname{Im}A_1\to0\): il termine sopravvive nel
limite fisico ed è reale all'ordine dominante.

## 4. Conseguenza: la formulazione del criterio va corretta

Il manoscritto (§9) attribuisce la pendenza 1 al **settore fermionico**. Kerr
scalare è bosonico e ha un termine \(\varepsilon^1\): l'attribuzione è troppo
stretta. La formulazione corretta è più generale:

> Il criterio non distingue spin intero da semintero. Distingue **presenza o
> assenza di struttura geometrica che precede il termine di ampiezza**. Se ne
> conoscono almeno due sorgenti indipendenti: la connessione di spin dei partner
> di Darboux (Dirac, statico) e la separazione sferoidale indotta dalla
> rotazione (Kerr, bosonico). Entrambe danno pendenza 1.

Predizione verificabile con la macchina già disponibile: la pendenza 1 deve
comparire in Kerr scalare come in Schwarzschild Dirac.

## 5. Il potenziale radiale: misura non circolare

**Riproduzione:** `python3.13 calculations/kerr_radial_order_profile.py`;
test in `calculations/test_kerr_radial_order.py` (6 test).

Con \(\Psi=\sqrt{r^2+a^2}R\), \(H=r^2+a^2\), \(\Delta=r^2-2r+a^2\),
\(\lambda=A+a^2\omega^2-2am\omega\) e la scalatura \(\omega=L\hat\Omega\),
\(m=\mu L\), \(\hat c=a\hat\Omega\), il potenziale efficace si decompone
**esattamente** come

\[
\frac{q}{L^2}=Q_0(r)
-\varepsilon\,\frac{\Delta A_1}{H^2}
-\varepsilon^2\left[\frac{\Delta A_2}{H^2}+\frac{h''}{h}\right],
\qquad
Q_0=\left(\hat\Omega-\frac{\mu a}{H}\right)^2-\frac{\Delta\bar A_0}{H^2},
\tag{5.1}
\]

con \(\bar A_0=A_0+\hat c^2-2\mu\hat c\). È il parallelo esatto della (5.1) di
Schwarzschild e della (9.1) di Dirac nel manoscritto. I due termini di
convenzione entrano in \(\bar A_0\), non in \(A_1\).

**Misura non circolare.** Definire \(A_1\) come coefficiente di \(1/L\) e poi
verificarne lo scaling sarebbe tautologico. Si calcola invece a ogni \(\ell\) il
potenziale con l'autovalore sferoidale **esatto** a quel \(\ell\), lo si
confronta con la forma di testa \(Q_0\), e si misura la pendenza in
\(\varepsilon\) del residuo \(D(r)=q/L^2-Q_0(r)\) su una finestra centrata sul
picco.

| \(a\) | \(\mu\) | \(\hat\Omega\) | \(r_{\rm picco}\) | \(A_0\) | **pendenza** |
|---|---|---|---|---|---|
| 0.0 | 0.5 | 0.192450 | 3.00000 | 1.00000 | **1.951** |
| 0.0 | 0.9 | 0.192450 | 3.00000 | 1.00000 | **1.951** |
| 0.3 | 0.5 | 0.205614 | 2.78835 | 0.99857 | **0.963** |
| 0.6 | 0.5 | 0.225120 | 2.47487 | 0.99315 | **0.996** |
| 0.9 | 0.5 | 0.261204 | 1.93474 | 0.97924 | **1.006** |
| 0.6 | 0.9 | 0.252767 | 2.23949 | 0.99778 | **1.014** |
| 0.9 | 0.9 | 0.330631 | 1.60599 | 0.99135 | **0.997** |

Il controllo statico dà 2, la rotazione dà 1, per ogni spin e ogni \(\mu\)
provati. **Il criterio del manoscritto, applicato a Kerr, misura pendenza 1 nel
settore bosonico.**

### 5.1 Validazioni

- Per \(a=0\): \(\hat\Omega=1/(3\sqrt3)\) e \(r_{\rm picco}=3\), e nessuna
  dipendenza da \(\mu\), come deve essere.
- Condizione di radice doppia \(Q_0(r_0)=Q_0'(r_0)=0\) soddisfatta a
  \(10^{-16}\)–\(10^{-12}\).
- I raggi trovati stanno sopra l'orbita fotonica equatoriale prograda
  \(r_{\rm ph}=2[1+\cos(\tfrac23\arccos(-a))]\) e vi si avvicinano al crescere
  di \(\mu\): 2.788 e 2.239 contro 2.630 e 2.189, 1.935 e 1.606 contro 1.558.
  È il comportamento corretto delle orbite fotoniche sferiche.

### 5.2 Due errori di metodo incontrati e corretti

1. **Tracciamento del modo** (§2.1): risolto con ordinamento a \(c\) reale più
   continuazione sull'autovettore.
2. **Solutore eikonale**: partendo da un ansatz analitico grezzo, `fsolve` non
   convergeva per \(a\gtrsim0.9\) o \(\mu\) grande, restituendo
   \(r_{\rm picco}\sim2500\). Risolto con continuazione in \(a\) a passi di
   0.05 dal caso di Schwarzschild, più una penalizzazione fuori dall'intervallo
   ammissibile. Entrambi gli errori producevano numeri plausibili ma sbagliati.

## 6. Limiti dichiarati

1. Che \(A_1\neq0\) per l'autovalore sferoidale **è noto** ed è contenuto nella
   letteratura eikonale su Kerr (Dolan 2010; Yang *et al.* 2012). Il contenuto
   eventualmente nuovo è la classificazione del §4 e la misura operativa del
   §5, non l'espansione.
2. \(A_0\) è stimato a \(\ell_{\rm rif}=400\), dove vale \(1-1/(4L^2)\): l'errore
   residuo su \(\hat\Omega\) è \(\sim1.5\times10^{-7}\). Per un uso di precisione
   andrebbe estrapolato.
3. La pendenza del controllo statico è 1.951, non 2.000: è l'effetto del
   troncamento a \(\ell\le160\), non un difetto strutturale.
4. Il termine di Madelung \(Q_M\) **non è ancora stato montato** sul profilo
   radiale di Kerr: servirebbe la frequenza QNM completa, cioè Leaver per Kerr
   con il sistema angolare-radiale accoppiato. Qui si misura l'ordinamento del
   potenziale, che è ciò che il criterio richiede, non il profilo dell'ampiezza.
5. Nessun modo di overtone \(n>0\) è stato considerato.
