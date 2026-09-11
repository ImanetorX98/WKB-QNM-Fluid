# Consegna a Codex — 12 settembre 2026

> **Precede tutti gli altri handoff.** Questo invalida il §9 del manoscritto, che
> era uno dei due pilastri della tesi. `CODEX_HANDOFF_CONDIZIONAMENTO.md` resta
> valido su Kerr §9.6 e sulla sensibilità; `CODEX_HANDOFF.md` su Vaidya.

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

### (a) Verificare \(A_1=0\) con un solutore indipendente — **bloccante**

Il nostro `spheroidal_eigenvalue` usa la matrice in base sferica armonica. Se
avesse un difetto suo, tutta la diagnosi cade. Serve un secondo metodo —
frazioni continue di Leaver per l'angolare, oppure l'*asymptotic iteration
method* — e la stessa tabella di \(A_2\).

**Non fidatevi del nostro numero perché è nostro.**

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

Questo trasformerebbe il §9 nella sezione più solida del lavoro. È il compito con
il rapporto valore/sforzo migliore.

### (c) Cercare lo stesso difetto altrove nel repository

Il difetto ha una forma riconoscibile: **una quantità discreta dentro una
regressione continua**. Da controllare in particolare:

- `kerr_radial_order_profile.py` righe 64 e 109 (`int(round`), due occorrenze)
- `kerr_madelung_profile.py` riga 108 e `kerr_wkb3_selfconsistent.py` riga 112
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
