# Vaidya: il termine misto è il difetto di adiabaticità

**Data:** 9 settembre 2026.
**Natura:** risultato **derivato** (simbolico), non ancora corroborato da
evoluzione numerica.
**Riproduzione:** `python3.13 calculations/vaidya_adiabatic_residual.py`.

## 1. La domanda, e perché si risponde senza PDE

L'ipotesi nulla del programma Vaidya era: se l'ampiezza dipende da \(v\) solo
attraverso \(M(v)\), allora
\(Q_\times=-2\varepsilon^2\dot M\,(\partial_M\partial_rA_0)/A_0+O(\dot M^2)\),
cioè il termine misto è determinato dal profilo congelato ed è ridondante come
lo era \(\mathcal E_M\).

È una domanda **all'ordine \(\dot M\)**, e a quell'ordine si risponde
simbolicamente: non serve un'evoluzione caratteristica, che è un progetto di
relatività numerica.

## 2. Impostazione

Vaidya entrante, \(ds^2=-f\,dv^2+2\,dv\,dr+r^2d\Omega^2\), \(f=1-2M(v)/r\). Per
\(\Phi=\psi(v,r)Y_{\ell m}/r\) la Klein–Gordon massless si riduce a

\[
2\psi_{vr}+\partial_r(f\psi_r)-U_\ell\psi=0,
\qquad U_\ell=\frac{f_r}{r}+\frac{\ell(\ell+1)}{r^2}.
\]

Il modo congelato è la soluzione statica a \(M\) fissato riscritta in coordinate
entranti: con \(v=t+r_*\), il modo \(e^{-i\omega t}R(r)\) diventa

\[
\psi_{\rm ad}(v,r)=Z(r;M)\exp\!\left(-i\!\int^v\!\omega(M(v'))\,dv'\right),
\qquad Z=e^{i\omega r_*}R .
\]

## 3. Risultato

Sostituendo e raccogliendo in potenze di \(\dot M\) (a tasso costante,
\(\ddot M=0\)):

**Ordine \(\dot M^0\).** Si ritrova l'equazione del modo congelato — validazione
dello schema.

**Ordine \(\dot M^1\).** Il residuo è

\[
\boxed{\;\mathcal R_1=2\,\partial_r\partial_M Z\;}
\]

senza altri termini: i coefficienti di \(Z\), \(\partial_rZ\) e
\(\partial_MZ\) sono tutti nulli.

## 4. Lettura

Il residuo ha **esattamente la struttura di \(Q_\times\)**. Ricordando
\(Q_\times=-2\varepsilon^2A_{vr}/A\) e che sotto l'ansatz adiabatico
\(A_{vr}\simeq\dot M\,\partial_M\partial_rA_0\), i due oggetti coincidono nella
forma a meno del passaggio da \(Z\) complessa ad \(A=|\psi|\).

Ne segue l'enunciato preciso:

> Il termine misto di Vaidya non è informazione dinamica indipendente: **è la
> misura del fallimento dell'approssimazione adiabatica**, ed è calcolabile
> dalla sola famiglia di modi congelati \(Z(r;M)\), senza risolvere il problema
> dinamico.

È lo stesso schema già incontrato due volte in questo progetto: l'oggetto di
Madelung riproduce una quantità già determinata da dati più semplici, invece di
aggiungerne. Con \(\mathcal E_M\) era \(|\Lambda_3|\); qui è la derivata mista
della famiglia statica.

## 5. Cosa NON è stato stabilito

1. Il passaggio da \(Z\) complessa ad \(A=|\psi|\) coinvolge la fase
   \(\exp(-i\int\omega)\) con \(\omega\) complessa, che contribuisce a \(|\psi|\)
   attraverso \(\operatorname{Im}\omega\). L'identificazione del §4 è
   strutturale; la costante di proporzionalità esatta va derivata.
2. La correzione \(\psi_1\) che ripara l'ansatz risolve l'operatore congelato con
   sorgente \(-2\partial_r\partial_MZ\). Il suo effetto retroagisce su
   \(A_{vr}\) a \(O(\dot M^2)\): non è stato calcolato.
3. Nessuna evoluzione numerica è stata eseguita. Il risultato è perturbativo in
   \(\dot M\) e non dice nulla sul regime di massa rapidamente variabile, che è
   proprio quello in cui Abdalla–Chirenti–Saa osservano l'effetto inerziale.
4. Il condizionamento misurato nell'Appendice C del manoscritto
   (\(Q_M\) amplifica di \(\sim10^2\) l'errore sulla frequenza) suggerisce che
   una verifica numerica di \(Q_\times\), che coinvolge una derivata **mista**
   su un campo evoluto, sarebbe ancora peggio condizionata. Va tenuto presente
   prima di investire in un codice caratteristico.

## 6. Conseguenza per il programma

La domanda decidibile del passo C ha ricevuto risposta, e la risposta è
negativa nel senso utile: **Vaidya non aggiunge un oggetto nuovo all'ordine
\(\dot M\)**. Resta aperta solo la parte che avevo già dichiarato probabilmente
fuori portata — il regime non adiabatico e il legame con l'effetto inerziale.

Il manoscritto può quindi riportare Vaidya come **terzo caso della stessa
tesi**, in forma analitica e breve, senza evoluzione numerica: la struttura
(6.1)–(6.5) del formalismo covariante, il residuo \(2\partial_r\partial_MZ\), e
la conclusione che il termine misto è il difetto di adiabaticità.
