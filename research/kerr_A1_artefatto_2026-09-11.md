# Il termine \(\varepsilon^1\) di Kerr non esiste: era l'arrotondamento di \(m\)

**Data:** 11 settembre 2026. **Gravità: massima.** Invalida il §9 del manoscritto
e la «correzione» che il §10.1 dichiarava di aver fatto al proprio criterio.

## 1. Il difetto

`kerr_eikonal_order_test.eikonal_fit` tiene fisso \(\mu\) e varia \(\ell\), ma
calcola \(m=\texttt{int(round}(\mu L))\). Poiché \(m\) deve essere intero,
\(\mu_{\rm eff}=m/L\) **non** è \(\mu\): differisce di \(\delta\mu=O(0.5/L)\).

E \(A_0\) dipende da \(\mu\). Quindi

\[
\frac{A}{L^2}\Big|_{\rm misurato}
= A_0(\mu) + \frac{\partial A_0}{\partial\mu}\,\delta\mu + \dots
= A_0(\mu) + \underbrace{\frac{\partial A_0}{\partial\mu}\,(\pm0.25)}_{\text{falso }A_1}\frac{1}{L}+\dots
\]

Per \(\mu=0.5\) e \(L=\ell+\tfrac12\): \(\mu L=\tfrac{\ell}{2}+\tfrac14\), che per
\(\ell\) **pari** arrotonda in giù di 0.25 e per \(\ell\) **dispari** in su di
0.25. Il modulo è sempre \(0.25/L\); cambia solo il segno.

## 2. La firma

| insieme di \(\ell\) | \(A_1\) misurato |
|---|---|
| 40, 60, 80, 120, 160, 240 (pari) | **−0.041441** |
| 41, 61, 81, 121, 161, 241 (dispari) | **+0.041443** |
| 42, 62, 82, … (pari) | −0.041443 |

**Cambia segno con la parità di \(\ell\)**, a modulo identico. Nessuna quantità
fisica si comporta così.

Ed è anche la ragione per cui il test di parità sulla pendenza *radiale* non
rivelava nulla: quella misura usa \(|\Delta q|\), e il modulo dell'errore di
arrotondamento è \(0.25/L\) in entrambe le parità — quindi pendenza 1 sempre.

## 3. La misura corretta

Si scelgono \(\mu\) e \(\ell\) tali che \(m=\mu(2\ell+1)/2\) sia intero
**esatto**: serve \(\mu=2a/b\) con \(b\mid(2\ell+1)\). Per \(\mu=2/3\) bastano gli
\(\ell\equiv1\pmod 3\); per \(\mu=2/5\) gli \(\ell\equiv2\pmod 5\).

Allora \(L^2\!\left(A/L^2-A_0\right)\) è **costante**:

| \(\mu\) | \(\hat c\) | \(A_2\) su una decade di \(L\) | dispersione |
|---|---|---|---|
| 2/3 | 0.6 | −0.2311 | 1.4e−5 |
| 2/3 | 0.3+0.2i | −0.2437 | 3.2e−6 |
| 2/5 | 0.6 | −0.2699 | 1.2e−5 |
| 4/5 | 0.8 | −0.1734 | 2.6e−4 |
| 2/3 | 0.9 | −0.2431 | 5.5e−5 |

> \(A_1=0\). Il primo termine è \(A_2/L^2\), anche per \(\hat c\) **complessa**.

Fit diretto con \(m\) esatto: \(A_1=-3\times10^{-7}\), residuo \(5\times10^{-11}\).

## 4. Conseguenza analitica per il potenziale radiale

Non serve rimisurare. Con \(\omega=L\hat\omega\), \(m=\mu L\) e
\(\Sigma=A+a^2\omega^2-2am\omega\):

\[
\frac{\hat q}{L^2}=\Big(\hat\omega-\frac{\mu a}{H^2}\Big)^2
-\frac{\Delta}{H^4}\Big(\frac{A}{L^2}+a^2\hat\omega^2-2a\mu\hat\omega\Big),
\]

e l'**unica** dipendenza da \(L\) è attraverso \(A/L^2\). Se \(A_1=0\) allora

\[
\frac{\hat q}{L^2}=Q_0+O(\varepsilon^2)\quad\text{identicamente}.
\]

**Kerr si comporta come il caso statico bosonico: pendenza 2.**

## 5. Che cosa cade

- §9.1–9.3, il termine di ordine \(\varepsilon\) da rotazione: **falso**
- §9.5, la discussione di attribuzione: priva di oggetto
- la serie di Kerr nella Figura 1 (pendenza 0.9964): artefatto
- l'abstract, dove Kerr è uno dei due casi che «rompono la degenerazione»
- il §10.1 nella forma attuale

## 6. Che cosa sopravvive — ed è più pulito di prima

Le altre due decomposizioni sono **identità algebriche esatte**, non misure:

\[
\frac{V_s}{L^2}=h^2+\varepsilon^2v_2 \quad(\S8),\qquad
\frac{V_\tau}{K^2}=h^2+\tau\varepsilon h' \quad(\S10).
\]

Il criterio torna quindi alla forma che il manoscritto dichiarava **superata**:

> pendenza 2 nel settore bosonico, pendenza 1 in quello fermionico.

E in questa forma è più forte, perché entrambi i membri sono identità esatte e
non pendenze adattate: la distinzione è fra un sistema del **primo** ordine, dove
la connessione di spin di Darboux compare a \(\varepsilon^1\), e uno del
**secondo**, dove non c'è nulla a quell'ordine — con o senza rotazione.

La «correzione» che il §10.1 rivendicava era essa stessa l'artefatto.

## 7. Lezione

Terzo caso nella stessa sessione in cui un **fit** produce un coefficiente
inesistente: il condizionamento \(10^2\), il polinomio secolare, e ora \(A_1\).
In tutti e tre il difetto era una quantità discreta o una variabile sbagliata
nascosta dentro la regressione.

Il controllo che li avrebbe trovati tutti: **variare un parametro che non
dovrebbe contare** — qui la parità di \(\ell\).
