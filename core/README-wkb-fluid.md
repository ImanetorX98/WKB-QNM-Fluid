# QNM di Schwarzschild, WKB di ordine superiore e fluido di Madelung

Questo piccolo progetto calcola frequenze quasi-normali (QNM) della metrica di
Schwarzschild con il metodo WKB al primo ordine e al terzo ordine di Iyer--Will.
Inoltre riscrive la funzione radiale in variabili di Madelung e rende visibile il
termine

\[
Q_M=-\frac{1}{2}\frac{A''}{A},
\qquad \Psi=Ae^{iS},
\]

dove gli apici indicano derivate rispetto alla coordinata tortoise adimensionale
\(x_*=r_*/M\). Il progetto usa \(G=c=1\), \(x=r/M\) e restituisce
\(\Omega=M\omega\).

## 1. Equazione radiale e condizioni QNM

Per un campo massless di spin \(s=0,1,2\), il settore assiale per \(s=2\) si
scrive

\[
\frac{d^2\Psi}{dr_*^2}+[\omega^2-V_s(r)]\Psi=0,
\qquad
r_*=r+2M\log\left(\frac{r}{2M}-1\right),
\]

\[
V_s(r)=f(r)\left[\frac{\ell(\ell+1)}{r^2}
+\frac{2M(1-s^2)}{r^3}\right],
\qquad f(r)=1-\frac{2M}{r}.
\]

La convenzione temporale è \(e^{-i\omega t}\). Le condizioni QNM sono onda
entrante all'orizzonte, \(\Psi\sim e^{-i\omega r_*}\), e onda uscente
all'infinito, \(\Psi\sim e^{+i\omega r_*}\). Di conseguenza il modo smorzato
ha \(\operatorname{Im}\omega<0\).

Valori ammessi nel codice:

- `spin=0`: campo scalare, \(\ell\geq0\);
- `spin=1`: campo elettromagnetico, \(\ell\geq1\);
- `spin=2`: perturbazione gravitazionale assiale Regge--Wheeler, \(\ell\geq2\).

## 2. Il conto WKB

Si trova il massimo \(x_0\) della barriera e si valutano lì le derivate
\(V_0^{(k)}=d^k(M^2V)/dx_*^k\). Il codice non differenzia numericamente:
costruisce in modo simbolico l'operatore

\[
\frac{d}{dx_*}=f(x)\frac{d}{dx}
\]

fino alla sesta derivata. La frequenza al terzo ordine è

\[
\Omega^2=V_0+\sqrt{-2V_0''}\,\Lambda_2
-i\alpha\sqrt{-2V_0''}\,(1+\Lambda_3),
\qquad \alpha=n+\frac12,
\]

con

\[
\Lambda_2=\frac{1}{\sqrt{-2V_0''}}\left[
\frac18\frac{V_0^{(4)}}{V_0''}\left(\frac14+\alpha^2\right)
-\frac1{288}\left(\frac{V_0^{(3)}}{V_0''}\right)^2
(7+60\alpha^2)\right],
\]

\[
\begin{aligned}
\Lambda_3=\frac{1}{-2V_0''}\bigg[&
\frac5{6912}\left(\frac{V_0^{(3)}}{V_0''}\right)^4(77+188\alpha^2)
-\frac1{384}\frac{(V_0^{(3)})^2V_0^{(4)}}{(V_0'')^3}(51+100\alpha^2)\\
&+\frac1{2304}\left(\frac{V_0^{(4)}}{V_0''}\right)^2(67+68\alpha^2)
+\frac1{288}\frac{V_0^{(3)}V_0^{(5)}}{(V_0'')^2}(19+28\alpha^2)\\
&-\frac1{288}\frac{V_0^{(6)}}{V_0''}(5+4\alpha^2)\bigg].
\end{aligned}
\]

Nella notazione diffusa in alcuni articoli, le nostre \(\Lambda_2\) e
\(\Lambda_3\) sono indicate rispettivamente con \(\Lambda\) e \(\Omega\).
Si sceglie la radice quadrata con parte reale positiva e parte immaginaria negativa.
Ponendo \(\Lambda_2=\Lambda_3=0\) si ottiene il primo ordine.

### Conto campione: modo gravitazionale fondamentale

Per \((s,\ell,n)=(2,2,0)\), il potenziale Regge--Wheeler ha il massimo in

\[
\frac{r_0}{M}=\frac{9+\sqrt{17}}{4}=3.280776406166.
\]

Lo script trova

\[
V_0=0.151286699570,\qquad
V_0''=-0.009919411728,
\]

e, con \(\alpha=1/2\),

\[
\Lambda_2=-0.141970003076,\qquad
\Lambda_3=-0.054526599347.
\]

Le due approssimazioni sono quindi

\[
\Omega_{\rm WKB1}=0.398849605924-0.088285381437i,
\]

\[
\boxed{\Omega_{\rm WKB3}=0.373162064857-0.089217447231i}.
\]

La differenza tra gli ordini è \(2.5704\times10^{-2}\). È utile come allarme
di convergenza, ma non va interpretata come incertezza statistica.

### Limiti da non nascondere

Il WKB per una singola barriera funziona meglio quando \(\ell>n\) e migliora
nel limite eikonale \(\ell\gg1\). La serie è asintotica: un ordine più alto non
garantisce automaticamente un risultato migliore. Il confronto WKB1--WKB3 è
un diagnostico pratico, non una barra d'errore rigorosa. Per bassa \(\ell\), alti
overtone o precisione spettrale conviene validare con frazioni continue di
Leaver; per andare al sesto/tredicesimo ordine conviene aggiungere le correzioni
pubblicate e usare approssimanti di Padé, mantenendo questo stesso calcolo delle
derivate al picco.

## 3. Dove compare il potenziale di Madelung

Dividendo l'equazione radiale per due si ottiene una forma Schrödinger-like:

\[
-\frac12\Psi''+\frac{V_s}{2}\Psi=\frac{\Omega^2}{2}\Psi.
\]

Questa scelta equivale a porre \(\hbar=m_{\rm eff}=1\). Ripristinando le
costanti, il termine standard è
\(Q_M=-(\hbar^2/2m_{\rm eff})A''/A\).

Con \(\Psi=Ae^{iS}\), la parte reale dà

\[
\frac{(S')^2}{2}+\frac{V_s}{2}+Q_M
=\frac{\operatorname{Re}(\Omega^2)}{2},
\qquad
Q_M=-\frac12\frac{A''}{A}.
\]

Questa è l'identificazione cercata. Tre cautele sono essenziali:

1. \(V_s/2\) è il potenziale esterno geometrico; \(Q_M\) è invece costruito
   dalla curvatura dell'ampiezza della soluzione.
2. Il fattore `1/2` dipende dalla normalizzazione Schrödinger-like. Nell'equazione
   d'onda non divisa per due il termine corrispondente è \(-A''/A\).
3. Poiché un QNM ha energia complessa, la continuità stazionaria contiene una
   sorgente/pozzo:
   \((A^2S')'=-\operatorname{Im}(\Omega^2)A^2\). Non è quindi un autostato
   normalizzabile di un Hamiltoniano hermitiano ordinario.

Lo script evita seconde derivate numeriche rumorose. Se
\(y=\Psi'/\Psi=A'/A+iS'\), usa direttamente l'ODE:

\[
\frac{A''}{A}=\operatorname{Re}(V_s-\Omega^2)+(S')^2.
\]

## 4. Installazione ed esempi

Da questa cartella:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Modo gravitazionale fondamentale \((s,\ell,n)=(2,2,0)\):

```bash
python schwarzschild_wkb.py --spin 2 --ell 2 --n 0 --order both
```

Conversione per un buco nero di 10 masse solari:

```bash
python schwarzschild_wkb.py --spin 2 --ell 2 --n 0 --order 3 --mass-solar 10
```

Scansione di multipoli e overtones:

```bash
python scan_qnms.py --spin 2 --ell-max 6 --n-max 2
```

Profilo di Madelung, CSV e figura:

```bash
python madelung_profile.py --spin 2 --ell 2 --n 0 --order 3
```

La figura confronta \(V_s/2\), il termine cinetico di fase \((S')^2/2\),
\(Q_M\) e la loro somma. Il CSV contiene anche il residuo ricostruibile della
relazione Hamilton--Jacobi. L'integrazione parte vicino all'orizzonte con la
condizione ingoing; siccome una frequenza WKB approssimata non soddisfa
esattamente entrambe le condizioni globali, il profilo lontano può contenere una
piccola componente indesiderata. L'identità locale di Madelung resta comunque
esatta per la soluzione integrata.

Test rapidi:

```bash
python -m unittest -v test_wkb.py
```

## 5. Riferimenti

- B. F. Schutz e C. M. Will, [*Black hole normal modes: a semianalytic
  approach*](https://doi.org/10.1086/184453), Astrophysical Journal **291**,
  L33 (1985).
- S. Iyer e C. M. Will, [*Black-hole normal modes: a WKB approach.
  I*](https://doi.org/10.1103/PhysRevD.35.3621), Physical Review D **35**,
  3621 (1987).
- S. Iyer, [*Black-hole normal modes: a WKB approach. II. Schwarzschild black
  holes*](https://doi.org/10.1103/PhysRevD.35.3632), Physical Review D **35**,
  3632 (1987).
- R. A. Konoplya, A. Zhidenko e A. F. Zinhailo, [*Higher order WKB formula for
  quasinormal modes and grey-body factors*](https://arxiv.org/abs/1904.10333),
  Classical and Quantum Gravity **36**, 155002 (2019).
