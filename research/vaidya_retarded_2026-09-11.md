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
