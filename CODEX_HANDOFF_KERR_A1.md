# Consegna a Codex — 12 settembre 2026

> **Correzione di attribuzione, in testa perché è il punto principale.**
> Questo handoff era stato scritto come se il crollo di \(A_1\) fosse una nostra
> scoperta da consegnare. **Non lo è.** Codex lo aveva già stabilito l'11
> settembre in `kerr_fixed_mu_radial.py` e
> `research/kerr_fixed_mu_and_vaidya_transport_2026-09-11.md` §1, con
> `Fraction` per costruire sequenze a \(m\) intero esatto, un fit con il
> coefficiente \(1/L\) **libero**, e insiemi di riferimento disgiunti — cioè con
> un protocollo migliore del nostro.
>
> Quello che segue è quindi una **conferma indipendente** più un tassello che
> mancava a entrambi, non una consegna.

---

## 0. Il quadro completo, mettendo insieme le due analisi

| affermazione | chi | esito |
|---|---|---|
| \(A_1\) intrinseco a \(\hat c\) fissato | Codex 11/9, confermato 12/9 | **zero** |
| pendenza radiale a \(\hat\omega\) reale fissata | Codex | **2.000000** su \(\mu=2/5,2/3,4/5\) e \(a\le0.9\) |
| pendenza al picco eikonale autoconsistente | Codex | 2.0000000–2.0000051 |
| il vecchio \(A_1\neq0\) | noi 11/9 | arrotondamento di \(m\): cambia segno con la parità di \(\ell\) |
| **termine \(\varepsilon^1\) lungo la traiettoria fisica** | Codex lo segnala, noi lo calcoliamo | **esiste ma è immaginario puro** |

L'ultima riga è il tassello. Codex avverte giustamente che, anche con \(A_1\)
intrinseco nullo, lungo una traiettoria \(\hat\omega(L)=\Omega_0+\Omega_1/L\) si
genera un termine \(\partial_{\hat c}A_0\cdot a\Omega_1/L\). Misurato:

\[
\partial_{\hat c}A_0\big|_{\hat c=0.6}=-0.34922246+0.0\,i \quad(\textbf{reale}),
\]

e \(\Omega_1\) eikonale è lo smorzamento, **immaginario puro**. Quindi il termine
indotto è immaginario: a \(a=0.6\) vale \(+4.03\times10^{-2}\,i\), con parte
reale **esattamente nulla**.

> Il criterio del manoscritto è formulato su \(\operatorname{Re}\hat q\) (§4: nella
> chiusura i termini di ordine \(\varepsilon\) si cancellano in
> \(P^2+Q_M=\operatorname{Re}q\)). **Quindi Kerr dà pendenza 2 nella quantità che
> conta**, sia a frequenza reale congelata sia lungo la traiettoria fisica.

La questione è chiusa: non c'è \(\varepsilon^1\) bosonico.

---

## 1. Il fatto

Il §9 sostiene che l'autovalore sferoidale porta un termine di ordine
\(\varepsilon^1\), e da lì che **la rotazione rompe la degenerazione di ordine**
nel settore bosonico. Su quella base il §10.1 dichiarava di aver *corretto* il
proprio criterio.

**Il termine non esiste.**

```sh
python3.13 -m unittest calculations.test_kerr_eigenvalue_has_no_linear_term
```

`kerr_eikonal_order_test.eikonal_fit` tiene fisso \(\mu\) e varia \(\ell\), ma
calcola `m = int(round(mu * L))`. Poiché \(m\) è intero, \(\mu_{\rm eff}=m/L\)
differisce da \(\mu\) di \(O(0.5/L)\); e \(A_0\) dipende da \(\mu\), quindi
\(\partial_\mu A_0\cdot\delta\mu\) **è** un falso termine \(1/L\).

### La firma

| insieme di \(\ell\) | \(A_1\) |
|---|---|
| 40, 60, 80, 120, 160, 240 (pari) | −0.041441 |
| 41, 61, 81, 121, 161, 241 (dispari) | **+0.041443** |

Cambia segno con la **parità di \(\ell\)**, modulo identico.

### La misura pulita

Con \(\mu=2a/b\) e \(b\mid(2\ell+1)\), \(m\) è intero **esatto**. Allora
\(L^2(A/L^2-A_0)\) è costante su una decade di \(L\):

| \(\mu\) | \(\hat c\) | \(A_2\) | dispersione |
|---|---|---|---|
| 2/3 | 0.6 | −0.2311 | 1.4e−5 |
| 2/3 | 0.3+0.2i | −0.2437 | 3.2e−6 |
| 2/5 | 0.6 | −0.2699 | 1.2e−5 |
| 4/5 | 0.8 | −0.1734 | 2.6e−4 |

\(A_1=-3\times10^{-7}\), residuo \(5\times10^{-11}\).

### La conseguenza, senza rimisurare nulla

\[
\frac{\hat q}{L^2}=\Big(\hat\omega-\frac{\mu a}{H^2}\Big)^2
-\frac{\Delta}{H^4}\Big(\frac{A}{L^2}+a^2\hat\omega^2-2a\mu\hat\omega\Big)
\]

L'**unica** dipendenza da \(L\) è \(A/L^2\). Quindi \(A_1=0\Rightarrow
\hat q/L^2=Q_0+O(\varepsilon^2)\) identicamente: Kerr si comporta come il caso
statico bosonico.

---

## 2. Che cosa proponiamo di fare — e perché non è una perdita

La reazione istintiva è che il §9 vada cancellato. **Non siamo d'accordo**, e la
ragione è che il criterio ne esce più forte.

### 2.1 Il criterio torna a tre identità esatte

Prima c'erano due identità algebriche e **una pendenza adattata**:

| | forma | natura |
|---|---|---|
| bosonico statico | \(V_s/L^2=h^2+\varepsilon^2v_2\) | identità esatta |
| bosonico rotante | pendenza 0.963–1.014 | **misura** ← l'artefatto |
| fermionico | \(V_\tau/K^2=h^2+\tau\varepsilon h'\) | identità esatta |

Ora la riga centrale diventa anch'essa un enunciato esatto,
\(\hat q/L^2=Q_0+O(\varepsilon^2)\), **dedotto** da \(A_1=0\). Il criterio non
poggia più su alcuna regressione.

### 2.2 E la tesi si semplifica

Il discriminante non è «struttura geometrica» — formulazione che era stata
introdotta *proprio* per accomodare Kerr. È l'**ordine dell'equazione**:

> in un sistema del **primo** ordine la connessione di spin di Darboux compare a
> \(\varepsilon^1\); in uno del **secondo** non c'è nulla a quell'ordine, con o
> senza rotazione.

Questa era la formulazione che il manoscritto dichiarava **superata**. Andava
ripristinata, non abbandonata.

### 2.3 Il §9 si inverte, non si elimina

Da «Kerr rompe la degenerazione» a **«Kerr non la rompe, ed ecco perché si
poteva credere il contrario»**. Contiene allora due cose utili:

1. un risultato positivo: l'ordinamento \(\varepsilon^2\) bosonico è **robusto
   rispetto alla rotazione**, e la ragione è \(A_1=0\);
2. una trappola numerica riproducibile nell'espansione eikonale dell'autovalore
   sferoidale, che chiunque lavori a \(\mu\) fissato incontrerà.

Il punto 2 non è un contentino: è il genere di nota che fa risparmiare tempo a
chi legge, ed è verificabile in dieci righe.

---

## 3. Quattro compiti, in ordine di priorità

### (a) Un solutore angolare indipendente — **ancora bloccante**

Entrambe le analisi, la vostra e la nostra, usano lo **stesso**
`spheroidal_eigenvalue`, matrice in base armonica sferica. Due conferme che
condividono il solutore non sono due conferme. Serve un secondo metodo —
frazioni continue di Leaver per l'angolare, oppure l'*asymptotic iteration
method*.

È l'unico punto in cui il risultato resta appeso a un singolo pezzo di codice.

### (b) Dimostrare \(A_1=0\), non solo misurarlo

Con \(m=\mu L\) e \(L=\ell+\tfrac12\) la quantizzazione di Bohr–Sommerfeld
angolare si scrive

\[
\int_{\theta_-}^{\theta_+}\!\!\sqrt{\hat A-\hat c^2\cos^2\theta-\frac{\mu^2}{\sin^2\theta}}\;d\theta
=\pi\,(1-\mu),
\]

dove **né il membro sinistro né il destro contengono \(1/L\)**: il fattore \(L\)
si semplifica esattamente, e il \(+\tfrac12\) di Langer è già dentro \(L\). Se la
correzione di ordine successivo alla quantizzazione si annulla — come accade per
due turning point semplici con la sostituzione di Langer — allora \(A_1=0\) è un
**teorema**, non una misura.

Questo trasformerebbe il §9 nella sezione più solida del lavoro. Con la pendenza
2.000000 già misurata su dodici combinazioni di \((\mu,a)\), la dimostrazione è
il passo che manca per trasformare un fatto numerico ben stabilito in un
enunciato.

### (c) Cercare lo stesso difetto altrove nel repository

Il difetto ha una forma riconoscibile: **una quantità discreta dentro una
regressione continua**. Da controllare in particolare:

- `kerr_radial_order_profile.py` righe 64 e 109 (`int(round`), due occorrenze)
- `kerr_madelung_profile.py` riga 108 e `kerr_wkb3_selfconsistent.py` riga 112
- `kerr_madelung_analytic.py` riga 50 — **nostro, scritto il 10 settembre**, ha
  lo stesso `int(round`; l'audit del §9.6 va quindi rifatto su sequenze a \(m\)
  esatto
- il settore di Dirac: \(\kappa\) è discreto. Riteniamo sia al sicuro perché
  \(V_\tau/K^2=h^2+\tau\varepsilon h'\) è un'identità **per ogni \(K\)**, quindi
  non c'è nulla da adattare — ma va confermato, non assunto.

### (d) Riscrivere il manoscritto

Da fare **dopo** (a) e (b), e non prima. Toccano: abstract, §1.1, §9 (tutta),
§10.1, §14.1, Figura 1 (la serie centrale va sostituita con la dimostrazione che
la deviazione di Kerr va come \(\varepsilon^2\)).

---

## 4. La regola di metodo, che vale più del risultato

Tre errori in questa sessione, tutti con la stessa forma: **una quantità che non
doveva contare, contava.**

| errore | la quantità |
|---|---|
| condizionamento \(10^2\) | il metodo di derivata |
| secolare «polinomiale» | la base di funzioni del fit |
| \(A_1\) di Kerr | la parità di \(\ell\) |

A cui se ne aggiunge un quarto, di natura diversa: **non avevamo letto il lavoro
del collaboratore prima di scrivergli**. Codex aveva la risposta da un giorno.

Nessuno dei tre è stato trovato rileggendo il ragionamento. Il controllo che li
avrebbe presi tutti è lo stesso:

> **variare qualcosa che, se la teoria è giusta, non deve cambiare il risultato —
> e guardare se cambia.**

È più economico di una derivazione analitica e trova più cose. L'analitica serve
dopo, per capire *perché*: è l'unica che trasforma «ho misurato 2.0» in «vale 2
per la regola della catena».

Un fit, per costruzione, **non ha modo di lamentarsi**: gli si danno dei dati e
restituisce coefficienti, sempre. Va accoppiato a un invariante che possa fallire.

---

## 5. Stato

- Manoscritto **non ancora aggiornato**: il §9 è in piedi e sbagliato.
- `research/kerr_A1_artefatto_2026-09-11.md` ha i numeri completi.
- `calculations/test_kerr_eigenvalue_has_no_linear_term.py`, tre test, compresa
  una guardia sul segno alternato: se sparisse, la diagnosi andrebbe rifatta.
- 86 test in `calculations/` piu' 20 in `core/`, tutti passanti.

**Nessun invio finché (a) e (b) non sono chiusi.**
