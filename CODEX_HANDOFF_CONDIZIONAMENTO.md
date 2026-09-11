# Consegna a Codex — 11 settembre 2026

> **Verifica successiva del bordo forzato:** integrata direttamente L H=S-2K G_r,
> con H analitica all'orizzonte. Il Wronskiano pesato meno la primitiva uscente
> si annulla al K candidato entro 2.5e-6 assoluti per r=30–40 (rtol=1e-12);
> due perturbazioni di K producono il residuo complesso previsto -2 deltaK P.
> Tre test e un'identità Mathematica superati. Limite emerso: H/G~2i omega r^2,
> quindi l'espansione adiabatico-radiale non è uniforme fino all'infinito.
> La prescrizione del problema congelato è verificata; resta il raccordo al
> problema temporale uscente. [Nota completa](research/vaidya_forced_boundary_check_2026-09-11.md).

> **Ultimo progresso Vaidya:** calcolato separatamente P=Reg integral mu G G_r
> e la norma generalizzata, verificando P=i omega Norm a circa 1e-10 relativo
> per cutoff 40–60. K=N/(2P)=1.1642938694+4.3489002268i è ancora un candidato:
> manca la verifica dei bordi della sorgente temporale per identificarlo con
> il trasporto fisico completo. Sono tabulati anche i primi valori candidati
> di Xi=K+partial_M ln G. 15 test e due identità Mathematica superati.
> [Nota e comandi](research/vaidya_denominator_progress_2026-09-11.md).

> **Aggiornamento successivo su Kerr/Vaidya:** non è più valida la rassicurazione
> storica del §5 di questo handoff sull'ordinamento Kerr. Il test arrotondava m,
> variando mu=m/L di ordine 1/L. Il nuovo confronto a mu esattamente fisso dà
> pendenza 2 su dodici casi scalari, confermata anche sul picco eikonale
> autoconsistente. Su Vaidya il numeratore supera i test, ma cambia sotto
> normalizzazioni dipendenti da M: il nuovo obiettivo è il trasporto combinato
> con la derivata del profilo, non N da solo. Dettagli e condizioni ancora
> aperte: [Kerr corretto e trasporto Vaidya](research/kerr_fixed_mu_and_vaidya_transport_2026-09-11.md).
> 23 test numerici e quattro identità Mathematica superati in questa analisi.

> **Aggiornamento Codex, risultato consolidato:** il fattore 2 è direzionale,
> non universale. Per omega=a+ib, Q_infinity=-b^2/L^2: lungo delta omega
> proporzionale a omega la sensibilità è 2; rispetto a perturbazioni complesse
> arbitrarie normalizzate a |omega| è 2|omega|/|b|. A ell=70,s=2 il limite è
> 281.9565; il calcolo esterno con derivate analitiche dà mediana puntuale
> 282.9420 e condizionamento L2 del profilo 283.1606. Non recupera il teorema
> di impossibilità del vecchio §5. **20 test numerici e 4 identità Mathematica
> superati.** Presso il picco i risultati cambiano fortemente fra ordini WKB
> e rimangono non validati; Kerr non ricontrollato.
> Stato e testo sostitutivo proposto:
> [conditioning_consolidated_2026-09-11.md](research/conditioning_consolidated_2026-09-11.md).
> Le affermazioni storiche sotto di universalità del fattore 2 sono superate.

> **Precedenza su `CODEX_HANDOFF.md`** per quanto riguarda il §5 del manoscritto.
> Quel documento e' ancora valido su Kerr, Vaidya e sul numeratore.

---

## 1. Il fatto, in tre righe

Il §5 del manoscritto — **il risultato che era stato promosso a enunciato di
testa** — afferma che il potenziale quantistico amplifica di ~10² l'errore sulla
frequenza, e ne deduce che nessun funzionale di \(Q_M\) possa predire \(\omega\).

**L'amplificazione misurata e' 2.** Universale in \(\ell\), spin e finestra.

```
python3.13 calculations/madelung_conditioning_schwarzschild.py
```

| \(\ell\) | 20 | 40 | 70 | 100 |
|---|---|---|---|---|
| pavimento | 1.7e-4 | 4.4e-5 | 1.5e-5 | 7.1e-6 |
| amplificazione | 2.0 | 2.0 | 2.0 | 2.0 |

Controllo incrociato: \(s=0\) e \(s=2\), finestre 10–25, 20–50, 30–70 — sempre
1.94–2.00.

**Ragione analitica.** Nella regione esterna \(Q_M\) e' dominato da
\(-\varepsilon^2(\operatorname{Im}u/\varepsilon)^2\) con
\(\operatorname{Im}u\to\operatorname{Im}(\omega)/L\), quindi
\(Q_M\propto(\operatorname{Im}\omega)^2\) e un errore relativo \(\delta\) su
\(\omega\) ne da' \(2\delta\) su \(Q_M\). Regola della catena. Nessuna ostruzione.

---

## 2. Perche' non se n'era accorto nessuno

**Il codice che genera la tabella del §5 non e' nel repository.** La misura
«Schwarzschild \(\ell=70\), finestra \(20<r<50\)» esiste solo come tabella in
`research/kerr_madelung_sensitivity_2026-09-09.md`.
`madelung_wkb_prediction.py` contiene solo la validazione su Poschl-Teller.

L'ha fatto emergere il tentativo di **disegnarne la figura**.

---

## 3. I tre difetti della ricostruzione, per non rifarli

Si sono nascosti l'uno dietro l'altro, come nel §11. Ogni volta il successivo era
invisibile finche' il precedente dominava.

| # | difetto | sintomo | correzione |
|---|---|---|---|
| 1 | BC troncata all'orizzonte | pavimento ×160 fra \(\delta=10^{-4}\) e \(10^{-10}\) | Frobenius completa, `horizon_values` |
| 2 | griglia uniforme in \(r\) | `predict` assume passo costante, ma \(dx_*=dr/f\) varia ×5 | integrazione **in \(x_*\)** con \(r\) di stato |
| 3 | troncamento della serie di Frobenius | \(Q_M\) **oscilla** di 5.9e-5 contro 1.9e-6, cambia segno | due termini per multipolo |

Il (3) e' il piu' insidioso: i coefficienti crescono con \(\ell(\ell+1)\), e a
\(\ell=70\) sessanta termini a \(x=0.5\) non bastano. Con 120 l'ampiezza
dell'oscillazione scende da 5.9e-5 a **1.7e-8**.

Diagnostico che li distingue: **il pavimento a frequenza esatta deve CALARE con
\(\ell\).** Se cresce, il difetto e' nel riferimento, non nel modello.

---

## 4. Che cosa va fatto — quattro compiti indipendenti

### (a) Verificare il fattore 2 analiticamente — priorita' alta

Derivare a mano \(\partial\ln Q_M/\partial\ln\omega\) nella regione esterna,
partendo da \(u=\sqrt{q}+\varepsilon^2u_2\) e
\(\ln A=-\tfrac12\ln|u|-\varepsilon^{-1}\!\int\operatorname{Im}u\). Attesa: 2
esatto al termine dominante, con correzioni \(O(\varepsilon)\). Se il calcolo
simbolico conferma, il §5 va **riscritto**, non corretto.

### (b) Cercare un regime dove l'amplificazione non sia 2

Il fattore 2 e' un artefatto della regione **asintotica**, dove \(Q_M\) e'
dominato dal decadimento esponenziale. Vicino al **picco della barriera** \(Q_M\)
ha struttura vera, e il condizionamento potrebbe essere diverso. Finestre da
provare: \(2.5<r<4\) (attorno a \(r=3M\)), e la regione di coalescenza dei
turning point del §12.

Se in quella regione l'amplificazione fosse grande, il §5 sarebbe **vero ma
enunciato nel posto sbagliato**. Se resta \(O(1)\), va ritirato.

### (c) Ricontrollare lo stesso schema su Kerr

`research/kerr_madelung_sensitivity_2026-09-09.md` §5 diagnostica un
comportamento erratico di \(Q_M\) su Kerr attribuendolo all'errore sulla
frequenza eikonale. **Quella diagnosi ora e' sospetta**: i difetti 1–3 sono gli
stessi, e il codice di Kerr li eredita.

### (d) Ricostruire la Figura 2

Non va prodotta finche' (a) e (b) non sono chiusi. Lo scheletro e'
`Schw-QNM-WKB-Fluid/figures/make_fig2.py`, **da non usare com'e'**: misura su
Poschl-Teller, dove \(Q_M=-\tfrac{\varepsilon^2}{4}(1+\operatorname{sech}^2y)\)
non dipende affatto dalla frequenza, quindi il test e' inconcludente per
costruzione (da' 1.41 costante in \(L\)).

---

## 5. Che cosa NON e' in discussione

- **§6, la ridondanza.** \(\mathcal E_M=c(n)|\Lambda_3|\), \(c=1.0437\) con
  dispersione **0.10%** su dodici modi. Misurato indipendentemente, riprodotto
  ieri nella Figura 3. Non dipendeva dal §5.
- **La ragione strutturale vera.** Non e' il condizionamento: e' l'argomento di
  Riccati. Un'espansione della Riccati non puo' contenere piu' della WKB che essa
  e'. Era gia' nel manoscritto (§1.1) e non ha bisogno del §5.
- **§§8–11**, l'ordinamento su tre geometrie, e il numeratore di Vaidya.

---

## 6. Conseguenza editoriale

Il manoscritto e' stato riorganizzato il 10 settembre attorno a **due** risultati
di testa, e il §5 era il primo. Se (a) e (b) confermano, la Parte II perde il suo
capofila e va ricostruita attorno al §6 da solo.

Non e' fatale — il §6 e' misurato e riproducibile — ma l'inquadramento «il
funzionale non puo' predire la frequenza, ed ecco perche'» va sostituito con «il
funzionale riproduce la correzione nota, ed ecco la misura».

**Nessun invio prima che questo sia chiuso.**
