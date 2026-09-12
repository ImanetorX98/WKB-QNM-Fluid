# Sincronizzazione Codex ↔ Claude — 12 settembre, sera

Due tracce, con le dipendenze dichiarate. **Analitico a Codex, numerico a me**,
e dove una traccia sblocca l'altra è scritto.

Stato dei controlli A–H del brief: **tutti eseguiti**. Quanto segue è ciò che
resta, più i difetti che l'esecuzione ha fatto emergere.

---

## Decisione che precede tutto — non tecnica

**A quale accuratezza serve la norma spettrale di Kerr, e per sostenere cosa?**

Il budget dell'ambiguità di ramo è già quantificato, \(|2\omega|e^{2\operatorname{Im}(\omega)r}\):

| \(r\) | 40 | 80 | 120 |
|---|---|---|---|
| fondamentale (\(\operatorname{Im}\omega\simeq-0.085\)) | 1.1e−3 | 1.3e−6 | **1.4e−9** |

* bersaglio \(10^{-6}\) → **spostare il bordo a \(r=120\) e basta**, zero lavoro;
* bersaglio \(10^{-15}\) → serve la traccia A1 qui sotto.

Va detto che la norma Kerr **non è sul percorso critico del manoscritto**: il §9
si regge su \(A_1=0\) e sulla (9.4), il §11 usa la norma di Schwarzschild, già
validata a \(5\times10^{-11}\). È un'estensione del nucleo analitico, non un
blocco della pubblicazione.

---

## Traccia analitica — Codex

### A1. Rappresentazione convergente per il ramo uscente di Kerr

**Solo se il bersaglio di accuratezza lo richiede.**

Il problema non è «selezionare» un ramo: è che una serie **asintotica** in
\(1/r\) non lo determina, perché la contaminazione entrante è oltre ogni ordine.
Con una rappresentazione **convergente** l'ambiguità non esiste per costruzione.

Cosa serve, concretamente:

1. la ricorrenza a tre termini di Jaffé/Leaver per la radiale di Kerr **nelle
   nostre convenzioni**, per il **solo bordo esterno** — quello interno è chiuso
   da `HeunC`, vedi N1 — variabile master \(\Psi=\sqrt{H}R\), \(e^{-i\omega t}\),
   \(\lambda=A+a^2\omega^2-2am\omega\) — con i prefattori espliciti;
2. il dominio di convergenza dichiarato, e la relazione fra convergenza della
   serie e condizione uscente;
3. \(D_+ = v\,\partial_r\log\psi_{\rm out}\) e \(D_{+,\omega}\) da quella serie,
   con le derivate in \(\omega\) ottenute **differenziando la ricorrenza**, non
   per differenze finite (è la stessa scelta che ha funzionato per \(D_{\pm}\)
   ai bordi e per \(A_c\)).

*Prima di derivare*, però, vale la pena guardare il punto N1 della mia traccia:
Mathematica ha `HeunC`/`HeunG` built-in e la radiale di Kerr è una Heun
confluente. Se le connessioni bastano, questo compito non si fa. L'ultima volta
proposi di riscrivere un solutore angolare che esisteva già.

### A2. \(I_2\) e il resto: fin dove arriva la struttura pari

Aperto e dichiarato tale nella vostra nota. Due domande distinte:

* la parità della quantizzazione vale **a tutti gli ordini**, o solo formalmente
  fino a quello calcolato? Serve per dire \(A_{2k+1}=0\) e non solo \(A_1=A_3=0\);
* esiste un controllo del **resto**, cioè un bound su \(\hat A-\sum_{k\le K}A_{2k}\varepsilon^{2k}\)?
  Senza, «serie asintotica» resta la descrizione corretta.

Non chiedo un teorema uniforme nei limiti di turning point coalescenti o
\(\mu\to0,1\): quelli sono dichiaratamente fuori.

### A3. Vaidya, la variabile ritardata dentro la barriera

Il §11.4 chiude la zona esterna: la correzione forzata è lì interamente
geometrica, e la separazione esatta \(H/G=4i\omega\!\int\!(r/f)dr_*+\tilde y\)
isola il contenuto non geometrico in \(\tilde y\), sorgentato solo da
\(\rho=G'/G-2i\omega/f\).

Resta da capire se \(\tilde y\) ammetta una forma chiusa o una caratterizzazione
sulla barriera, dove \(\rho=O(1)\). È il punto in cui l'eventuale memoria vive.

---

## Traccia numerica — mia

### N1. `HeunC` per i bordi di Kerr — **ESEGUITO: riuscito a metà**

> **Esito.** Bordo **interno risolto**: la radiale è confluente Heun, i parametri
> escono con residuo di matching nullo, e \(D_-\) da `HeunC` concorda con la serie
> di Frobenius a **undici cifre** dove quella vale, con residuo della Riccati
> \(10^{-13}\) fino a \(\rho=1\) — contro un raggio di \(1.7\times10^{-3}\) della
> serie, cioè un fattore \(\sim600\).
>
> Bordo **esterno non risolto**: `HeunC` è la soluzione al punto singolare
> **regolare**, e per l'infinito irregolare servirebbero i coefficienti di
> connessione, che *Mathematica* non espone.
>
> **A1 resta necessario ma dimezzato**: solo il bordo esterno, e solo se il
> bersaglio di accuratezza supera il budget. Dettagli in
> [`research/claude_kerr_heun_2026-09-12.md`](research/claude_kerr_heun_2026-09-12.md).

### N1 (testo originale)

Mathematica ha le Heun built-in e nessun pacchetto Teukolsky installato. La
radiale di Kerr è confluente Heun, con punti singolari regolari in \(r_\pm\) e
irregolare all'infinito.

* `HeunC` è la soluzione locale al punto regolare: dovrebbe dare \(D_-\)
  **convergente**, risolvendo il raggio \(1.7\times10^{-3}\) che la mia serie di
  Frobenius ha a \(a=0.6\), \((\ell,m)=(28,19)\);
* per \(D_+\) serve la connessione all'infinito irregolare, ed è lì che potrebbe
  non bastare. Da provare, non da assumere.

Esito che sblocca o annulla A1.

### N2. Rifare l'audit del §9.6 su sequenze a \(m\) intero esatto

**Difetto noto e non ancora corretto.** Sei file contengono ancora
`m = int(round(mu * large_l))`:

| file | riga | alimenta |
|---|---|---|
| `kerr_madelung_analytic.py` | 63 | la nota di precisione del §9.6 — **mia, del 10 settembre** |
| `kerr_madelung_profile.py` | 108 | la tabella del §9.6 |
| `kerr_wkb3_selfconsistent.py` | 112 | la frequenza autoconsistente del §9.6 |
| `kerr_radial_order_profile.py` | 64, 109 | non più citato dal manoscritto |
| `kerr_eikonal_order_test.py` | 100 | non più citato dal manoscritto |

La tabella del §9.6 usa \(\mu=0.5\), per cui \(m\) va **sempre** arrotondato di
\(0.25\). Va rifatta su \(\mu=2/3\) con \(\ell\equiv1\pmod 3\), o su \(\mu=2/5\)
con \(\ell\equiv2\pmod 5\).

Attenzione: la conclusione qualitativa del §9.6 — nulla di patologico
nell'ampiezza — non dipende da quelle cifre, ma le cifre sì.

### N3. Aggiornare il §9.2 con \(A_2\) in forma chiusa

Fatto in parte: la (9.4) è nel manoscritto e i controlli sono in
`research/claude_second_period_check.md`. Resta da sostituire ovunque i valori
misurati con quelli calcolati, ora che esistono.

### N4. Il rango cinque al secondo ordine

Il controllo G ha verificato il rango **due** al primo ordine. La vostra nota
prevede al più **cinque** al secondo — due componenti di \(e\) più i tre monomi
in \(d\). Non è stato testato, e il brief non lo chiedeva. È un test netto e lo
posso fare quando serve.

---

## Dipendenze, in una riga

```
N1 (HeunC)  ──┬─► se basta ──► A1 NON serve
              └─► se non basta ──► A1 (Codex) ──► N1bis (verifica numerica)

N2, N3, N4   indipendenti, posso procedere subito
A2, A3       indipendenti da tutto il resto
```

---

## Regola che vorrei tenere

Tre errori di questa sessione avevano la stessa forma — una quantità che non
doveva contare, contava: il metodo di derivata, la base del fit, la parità di
\(\ell\). E un quarto, mio, ieri: l'espressione analitica al posto della
definizione a tratti.

Nessuno trovato rileggendo il ragionamento. Il controllo che li prende è sempre
lo stesso, ed è a costo nullo:

> **variare qualcosa che, se la teoria è giusta, non deve cambiare il risultato.**

E una seconda, imparata la scorsa ora: **prima di reimplementare, guardare se
esiste già**. `SpheroidalEigenvalue` ha chiuso in dieci minuti un compito che
avevo classificato come bloccante.
