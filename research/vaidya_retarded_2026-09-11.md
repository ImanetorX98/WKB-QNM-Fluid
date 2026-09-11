# Vaidya: il secolare di Codex è il tempo ritardato, e si toglie

**Data:** 11 settembre 2026.
**Origine:** `vaidya_forced_boundary_check_2026-09-11.md` §5 (Codex).
**Codice:** `calculations/vaidya_retarded_uniformity.py`, 3 test.

## 1. Il problema che Codex ha trovato

Scrivendo la correzione forzata come \(H=Gy\), l'equazione esterna dà
\(y\sim k r^2\) con \(k=2i\omega\). Quindi \(\dot M\,H/G\) cresce come
\(\dot M|\omega|r^2\): **lo sviluppo adiabatico non è uniformemente piccolo a
grande \(r\)**, e i limiti adiabatico e di grande distanza non commutano.

Verificato indipendentemente con sympy: il coefficiente quadratico è esattamente
\(k\).

È un difetto reale del §11 come scritto, e il manoscritto non lo dichiara.

## 2. Che cos'è

**Un termine secolare, cioè lo sviluppo di Taylor di un argomento ritardato.**

Nella regione esterna la riduzione caratteristica di Vaidya diventa
\(2\psi_{vr}+\psi_{rr}=0\), e si verifica per sostituzione diretta che

> **ogni** funzione \(F(u)\) della coordinata uscente \(u=v-2r_*\) la risolve
> **esattamente**, qualunque sia \(M(v)\).

Non «approssimativamente per \(\dot M\) piccolo»: esattamente, senza ipotesi
sulla massa. Nella regione esterna l'adiabaticità organizzata attorno a \(u\)
non è un'approssimazione.

Sviluppare \(F(v-2r)\) attorno a \(v\) a \(r\) fissato dà

\[
F(v-2r)=F(v)-2rF'(v)+2r^2F''(v)-\dots
\]

e il quadratico è il secolare osservato.

## 3. Il coefficiente coincide, non solo l'ordine

Questo è il controllo che separa una diagnosi da una coincidenza. Con
\(F=\exp(-i\!\int\!\omega\,dv)\) e \(\omega=\omega_1/M\), la parte lineare in
\(\dot M\) del termine \(2r^2F''(v)/F\) vale

\[
-2i\,\omega'(M)\Big|_{M=1} = +2i\omega = k .
\]

Due derivazioni indipendenti — l'equazione forzata di Codex e lo sviluppo del
ritardo — danno **lo stesso coefficiente**, differenza simbolica zero.

## 4. La cura

Il difetto è nella scelta della variabile, non nell'espansione.

> L'ansatz adiabatico va congelato al **tempo di emissione** \(v-2r_*\), non al
> tempo locale \(v\).

Il contenuto fisico è immediato: l'onda uscente osservata a raggio \(r\) è
partita dalla fotosfera quando la massa era diversa, e congelarla a \(v\)
accumula un errore di fase quadratico nella distanza percorsa. Con la variabile
ritardata il secolare non c'è.

## 5. Che cosa cambia per il manoscritto

Il §11 va corretto in due punti, e migliora.

- **Va dichiarato** che l'espansione in \(\dot M\) attorno al tempo locale è non
  uniforme, con il raggio di validità \(r\lesssim(\dot M|\omega|)^{-1/2}\).
- **Va aggiunto** che la non uniformità è risolta, non solo diagnosticata, e che
  la variabile corretta è \(u\). È un enunciato più forte di quello che il §11
  fa oggi: non «l'espansione vale in una regione», ma «la variabile sbagliata
  produce un secolare calcolabile, e quella giusta non lo produce».

La (11.2) — il residuo \(=2\partial_r\partial_M Z\) — **non cambia**: è un
enunciato a \(r\) fissato e ordine \(\dot M\), e il secolare riguarda il
comportamento a grande \(r\) della correzione, non la sorgente.

## 6. Che cosa resta aperto

Riscrivere \(Z_1\) nella variabile ritardata e verificare numericamente che il
secolare sparisca — qui è dimostrato analiticamente nella regione esterna
asintotica, non su tutto il dominio. Il raccordo con la regione della barriera,
dove \(u\) non è più caratteristica, è il punto delicato.

---

## 7. Seguito: nella zona esterna la correzione è **interamente** ritardo

Il conto è più forte di quanto scritto sopra. Con la sottrazione generale
\(\chi=a\,r^2+b\,r\), l'annullamento simultaneo del termine lineare e della
costante nella sorgente forzata dà

\[
a = 2i\omega,\qquad b = -2K,\qquad \text{residuo} = \mathbf{0}\ \text{esatto}.
\]

Non «più piccolo»: **zero simbolico**. A ordine \(\dot M\), nella zona esterna,
la correzione forzata non contiene nulla oltre la retrodatazione:

- il **quadratico** è la deriva di massa, \(\theta''=\omega_M\dot M\);
- il **lineare** è lo shift di frequenza di ordine \(\dot M\), anch'esso
  valutato al tempo di emissione.

### La struttura asintotica, in forma chiusa

**Correzione al primo tentativo.** La verifica sul dominio vero era fatta
adattando \(H/G\) a un polinomio in \(r_*\). Un adattamento su intervallo
finito non distingue \(r^2\) da \(r_*^2\) e assorbe un logaritmo nei
coefficienti polinomiali: ne usciva un quadratico basso dello 0.5% e la
conclusione, **sbagliata**, che \(H/G\) fosse un polinomio.

Rifatto analiticamente (`calculations/vaidya_secular_series.py`). L'equazione
per \(p=Dy\) è del **primo ordine**, quindi si inverte per serie senza
integrare nulla:

\[
D^2y+\Big(2i\omega+2f\tfrac{u'}{u}\Big)Dy=\frac{f\,(S-2KG_r)}{G},
\qquad D=\partial_{r_*} .
\]

I coefficienti sono serie pure in \(1/r\), quindi la forma asintotica è
**forzata**:

\[
\frac{H}{G}=2i\omega\,r^2+c_1 r+c_{\log}\ln r+\sum_k c_k r^{-k}.
\]

| | risultato |
|---|---|
| quadratico | \(2i\omega\) con scarto **zero macchina** (0.0, 2.2e-16, 1.1e-16) |
| in \(r\), non \(r_*\) | un \(r_*^2\) produrrebbe \(r\ln r\), assente |
| \(c_{\log}\neq0\) | \(H/G\) **non** è un polinomio |

### Universalità, e forma chiusa

I rapporti a \(2i\omega\) valgono **8 e 24 esatti** per tutti e sei i modi
provati (\(\ell=2,3,4\); \(s=0,1,2\)); da \(1/r\) in poi dipendono dal modo.
I tre termini di testa sono dunque **geometrici**, non del modo, e si
riconoscono:

\[
\frac{H}{G}=4i\omega\!\int\!\frac{r}{f}\,dr_*+O(1/r),
\]

\[
\int\frac{r\,dr}{f^2}=\frac{16M^3+(2M-r)\big(24M^2\ln(r-2M)+8Mr+r^2\big)}{2(2M-r)}
=\frac{r^2}{2}+4Mr+12M^2\ln r+O(1/r).
\]

Moltiplicando per \(4i\omega\) si ottengono esattamente \(2i\omega(r^2+8r+24\ln r)\)
a \(M=1\).

> Il secolare di testa è il **trasporto di fase lungo il raggio uscente**, pesato
> da \((dr_*/dr)^2\), e non conosce né il multipolo né lo spin.

### La conseguenza che conta

Nella zona esterna la correzione forzata è, ai tre ordini di testa, una quantità
puramente geometrica: nessuna informazione sul modo, quindi nessuna memoria.
Ciò che dipende dal modo comincia a \(O(1/r)\) — **decadente**.

Qualunque effetto genuino deve vivere nella regione della **barriera**. È un
restringimento utile: elimina il settore dove i numeri sono enormi e le
cancellazioni peggiori.
