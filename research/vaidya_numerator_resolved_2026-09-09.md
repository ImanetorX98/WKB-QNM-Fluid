# Il numeratore di Vaidya, risolto — e perché mpmath non serviva

**Data:** 9 settembre 2026.
**Codice:** `calculations/vaidya_numerator_factored.py`, test in
`calculations/test_vaidya_numerator_factored.py`.
**Precede:** [`vaidya_original_bug_resolved_2026-09-09.md`](vaidya_original_bug_resolved_2026-09-09.md).

## Risultato

\[
N \;=\; 40.294405 + 4.154483\,i ,\qquad \ell=2,\ s=0,\ n=0,\ a=4 .
\]

Indipendenza dal taglio esterno: **7.6e-6** di dispersione fra \(L_+=40\) e 90
(1.9e-7 escludendo \(L_+=90\), che tocca il bordo \(r_{\rm far}=95\)). Era 3.5e-1.
**Quattro ordini e mezzo**, senza multiprecisione.

Regge su tutti i modi provati:

| \(\ell\) | \(s\) | \(N\) | dispersione |
|---|---|---|---|
| 2 | 0 | +40.294405 + 4.154483i | 1.9e-7 |
| 3 | 0 | +33.469793 + 46.968735i | 5.0e-7 |
| 4 | 0 | −15.533330 + 73.646500i | 7.6e-7 |
| 2 | 2 | +3.494440 − 25.557505i | 1.6e-7 |
| 3 | 2 | +47.808055 − 1.028573i | 2.4e-7 |
| 4 | 2 | +33.379221 + 59.121712i | 7.4e-7 |
| 1 | 1 | +3.478993 − 18.328993i | 1.8e-5 |
| 2 | 1 | +36.305268 − 8.102899i | 3.1e-7 |

**\(N\neq0\) ovunque, di ordine 10.** La forzatura risonante ha componente sul
modo. Non è più un'affermazione al 10%: è a sei cifre.

## Le quattro cause, in ordine di scoperta

Nessuna era quella sospettata. Ogni volta la successiva era invisibile finché la
precedente dominava — è questa la ragione per cui le "esclusioni" documentate
erano valide *al loro livello di errore* e inutili al successivo.

| # | causa | ampiezza | come si è vista |
|---|---|---|---|
| 1 | quadratura disallineata di una cella | 1.0e-3 | Codex: `integrand[i:j]` esclude `j`, `F[j]-F[i]` no |
| 2 | `np.gradient` per \(S=2(rZ')'\) | 2.2e-7 | costante in \(r\), \(=(2\omega h)^2/6\) |
| 3 | trapezio invece di Simpson | 1.1e-7 | grid refinement |
| 4 | **condizione al bordo dell'orizzonte** | dominante | vedi sotto |

Il punto 2 non richiedeva codice nuovo: `vaidya_solvability.source_integrand`
ricavava già \(S\) dall'ODE. Era `setup()` del riproduttore a usare differenze
finite.

## Il punto 4, che era il limite vero

\(\psi=e^{-i\omega r_*}\) imposta a \(r=2+\delta\) sbaglia di \(O(\delta)\), e
l'errore **semina il modo entrante all'infinito**. Due misure lo identificano:

- **il tasso.** \(|g/C-1|\) decade esattamente come \(|e^{-2i\omega r_*}|\):
  3.3e-7 → 5.7e-9 → 1.0e-10 a \(r=40,60,80\), rapporto 3300 contro l'atteso
  \(e^{0.1935\,\Delta r_*}\).
- **l'insensibilità a \(\omega\).** Perturbando \(\omega\) di 1e-9 (Leaver ne dà
  dieci cifre) il residuo non si muove di una cifra. Non è la frequenza.

Ridurre \(\delta\) migliora la BC ma peggiora la cancellazione vicino
all'orizzonte, dove \(f=1-2/r\) diventa \(\delta/2\). Ottimo a
\(\delta=10^{-10}\), e valeva 2.4e-4:

| \(\delta\) | 1e-4 | 1e-6 | 1e-8 | 1e-10 | 1e-12 | 1e-14 |
|---|---|---|---|---|---|---|
| dispersione | 3.5e-1 | 3.2e-2 | 1.8e-3 | **2.4e-4** | 9.2e-3 | 7.5e-1 |

Una U: è il segno che si stanno scambiando due errori invece di eliminarne uno.

## La cura: togliere l'esponenziale, non comprarlo con le cifre

Due cambi di variabile, uno per regione:

\[
\text{dentro}\quad \psi=e^{-i\omega r_*}h(r)\ \Longrightarrow\ \boxed{Z=h},
\qquad
\text{fuori}\quad \psi=e^{+i\omega r_*}u(r)\,g(r),\quad g\to C .
\]

Dentro, \(Z=e^{i\omega r_*}\psi=h\): l'esponenziale sparisce **identicamente**,
non approssimativamente. L'equazione è \(f h''+(f'-2i\omega)h'-Uh=0\).
Fuori, \(g''+\text{drift}\,g'+\text{defect}\,g=0\) con i coefficienti costruiti
esattamente dai coefficienti della serie (`Asymptotic.at`). Nessuna delle due
variabili cresce: il range dinamico \(10^9\) del vecchio \(\psi\) è sparito.

La BC al bordo diventa **Frobenius a due termini**. Il punto \(r=2\) è singolare
regolare, e l'equazione indiciale dà
\(c_1=U(2)/(\tfrac12-2i\omega)\),
\(c_2=[c_1(\tfrac12+U(2))+U'(2)]/(2-4i\omega)\), errore \(O(\delta^3)\).
Con \(\delta=10^{-5}\): errore \(10^{-15}\) e nessuna cancellazione. La U sparisce
— \(\delta\) da \(10^{-3}\) a \(10^{-8}\) muove \(N\) di 2.8e-6.

Costo: **5 secondi**, doppia precisione. `mpmath` non è mai entrato.

## Robustezza

| parametro | intervallo | variazione di \(N\) |
|---|---|---|
| \(\delta\) | 1e-3 … 1e-8 | 2.8e-6 |
| \(r_{\rm match}\) | 15 … 35 | 5.5e-6 |
| punti interni | 4e4 … 3.2e5 | **0** (vedi sotto) |
| termini della serie | 12 … 24 | 3.9e-7 |

Lo zero esatto sulla risoluzione non è un test superato: è **quantizzazione**.
\(N\approx40\) è la differenza fra \(\int I\) e \(F\), entrambi \(\sim10^{10}\);
l'ulp di un double lì vale 1.9e-6, e i valori stampati cadono esattamente su
quella griglia. \(N\) ha circa **sette cifre significative**, e raffinare la
griglia non ne aggiunge. *Questo* sì sarebbe il caso per mpmath — ma sette cifre
sono già tre più del necessario.

## Che cosa questo non dice

Il limite inferiore \(a=4\) resta un taglio arbitrario: il bordo **interno**
della regolarizzazione di Leung non è trattato. \(N\neq0\) stabilisce che la
proiezione risonante non si annulla; non è ancora un coefficiente fisico, e non
dimostra memoria osservabile.

## La lezione, di nuovo la stessa

Quattro cause in fila, ciascuna nascosta dalla precedente. Le esclusioni della
sessione precedente — "quadratura esclusa: trapezio = Simpson a 7 cifre" — erano
**corrette e inutili**: vere a errore 1e-3, false a errore 1e-7. Un'esclusione
va datata con il livello di errore a cui è stata fatta, altrimenti scade in
silenzio.

E la diagnosi giusta non è arrivata guardando meglio i numeri, ma cambiando la
domanda: non «quanta precisione serve» ma «perché serve precisione». La risposta
— nove ordini di range dinamico — si toglieva con l'algebra.
