# Consegna a Codex: da dove ripartire

**Data:** 9 settembre 2026.
**Scopo:** riprendere l'analisi senza la cronologia della chat.
**Simmetrico a** `CLAUDE_HANDOFF.md`, che è dell'8 settembre e **non copre il
lavoro del 9**: frequenza autoconsistente di Kerr, figura, regolarizzazione di
Leung, algebra di serie. Leggere questo, non quello, per lo stato corrente.

> **Nota di stato Git.** Il lavoro dell'8 e 9 settembre è sul branch
> `analisi-kerr-vaidya-2026-09`, non ancora unito a `main` (base `ec97d56`).
> Verificare `git status` e il branch attivo prima di toccare qualcosa, e non
> sovrascrivere modifiche dell'utente. I PDF sotto paywall stanno in
> `Schw-QNM-WKB-Fluid/papers/`, che è in `.gitignore`: **non vanno committati**.

---

## 1. Che cos'è il progetto, in un paragrafo

Si studia che cosa la decomposizione di Madelung dica sui modi quasi-normali dei
buchi neri, e **a quale ordine**. La tesi che ha resistito a tutti i controlli:
la riscrittura è una trasformazione esatta, il suo contenuto informativo è
esattamente quello della WKB da cui proviene, e l'unico punto non ridondante è
l'**ordinamento** — a quale potenza di \(\varepsilon\) compare il primo termine
subprincipale, e se il termine di ampiezza sia solo o preceduto da struttura
geometrica.

**Il macchinario formale è l'equazione di Riccati in \(\varepsilon\)**, e non è
nostro: la chiusura è l'eq. (1.2) di Delabaere–Dillinger–Pham. Il contributo sta
in ciò che la geometria *immette* in \(q\), non nel risolverlo.

Il documento da leggere per primo è
[`research/RESULTS_FOR_MANUSCRIPT.md`](research/RESULTS_FOR_MANUSCRIPT.md):
indice consolidato, 23 sezioni, ogni voce dichiara se è derivata, misurata, nota
in letteratura o corretta rispetto a un'affermazione precedente.

---

## 2. Stato: il manoscritto è pronto, l'analisi ha un solo buco

`Schw-QNM-WKB-Fluid/manuscript.md` — 12 sezioni, 3 appendici, Figura 1,
34 riferimenti con DOI verificati su Crossref, ~8200 parole (≈14 pagine su due
colonne). **Non contiene sezioni che dichiarino "non stabilito".**

Test: **57**, tutti passanti.

```bash
# dalla RADICE del progetto (questi importano calculations.*)
python3.13 -m unittest calculations/test_static_madelung_benchmark.py \
  calculations/test_poschl_teller_madelung.py \
  calculations/test_order_resolved_madelung.py \
  calculations/test_uniform_madelung_defect.py     # 12

cd calculations && python3.13 -m unittest \
  test_kerr_eikonal_order.py test_kerr_radial_order.py \
  test_kerr_wkb3_selfconsistent.py test_outgoing_asymptotics.py   # 25

cd core && python3.13 -m unittest \
  test_wkb.py test_leaver.py test_dirac_madelung.py               # 20
```

**Interprete: `python3.13`.** Il `python3` di sistema non ha numpy. I primi
quattro test falliscono con `ModuleNotFoundError: calculations` se lanciati da
dentro `calculations/`: è la directory, non una regressione.

---

## 3. Il problema aperto, enunciato per essere attaccato a freddo

**Dove:** [`research/vaidya_solvability_2026-09-09.md`](research/vaidya_solvability_2026-09-09.md) §5.2.
**Codice:** `calculations/asymptotic_series.py`, `calculations/vaidya_solvability.py`,
`calculations/outgoing_asymptotics.py`.

### Contesto

Su Vaidya, a ordine \(\dot M\) la correzione \(Z_1\) risolve
\(L_M Z_1=-2\partial_r\partial_M Z\), dove \(L_M\) è l'operatore congelato,
**singolare** alla frequenza QNM perché \(Z\) sta nel suo nucleo. La componente
della sorgente parallela al modo dà uno spostamento di frequenza, quella
ortogonale dà \(Z_1\) — ed è lì che potrebbe entrare la dipendenza dalla storia,
cioè la "memoria". Serve quindi la proiezione, regolarizzata alla Leung *et al.*
(1998), eq. (2.14)–(2.16).

### Che cosa è già fatto e verificato

| pezzo | esito |
|---|---|
| residuo adiabatico a ordine \(\dot M\) | \(2\partial_r\partial_M Z\), simbolico |
| **denominatore** (norma generalizzata) | indipendente da \(L_+\) allo 0.2%: 10.21 a \(L_+=40\), 10.19 a \(L_+=70\) |
| serie asintotica uscente | residuo \(10^{-12}\); \(D_+\) sul ramo corretto |
| algebra di serie con derivate esatte | \(f\cdot(1/f)=1\); ricorsione dell'antiderivata a \(1.8\times10^{-15}\) |
| asintotica della sorgente | \(W=2u[(2i\omega/f)P+P']\), \(P=r(2i\omega u/f+u')\), grado principale \(r^{+1}\) |
| costante di normalizzazione \(C\) | convergente a 12 cifre (l'integrando è **quadratico** in \(Z\)) |

### La discrepanza da spiegare

Con tutti gli ingredienti verificati a \(10^{-7}\) o meglio:

- integrando numerico contro serie, rapporto **complesso**: \(0.9999998\)
- antiderivata contro integrale della serie su griglia fine: \(1.000000000\)

eppure **integrale definito numerico e differenza dell'antiderivata differiscono
di \(1.65\times10^{-4}\), costante sull'intervallo**. Un errore *relativo*
costante su un integrale che cresce esponenzialmente produce un residuo che
cresce esponenzialmente: è questo a impedire l'indipendenza da \(L_+\) del
numeratore.

Sospetto residuo: la quadratura sulla griglia dell'ODE. Ma la stima dell'errore
di trapezio a quella risoluzione è \(10^{-7}\), tre ordini sotto. **Non c'è
un'ipotesi**, ed è per questo che ci si è fermati.

### Primo passo suggerito

Confrontare \(\int\text{pred}\,dr\) calcolato **sulla griglia dell'ODE** con
\(F(b)-F(a)\), invece che su una griglia fine indipendente. Isola se la
discrepanza sia nella quadratura o altrove. È un test da dieci righe e non è
stato fatto.

---

## 4. Che cosa NON rifare

Verificato e chiuso. Rifarlo è tempo perso.

| ipotesi | esito |
|---|---|
| troncamento della matrice sferoidale | convergente a \(10^{-8}\); la causa era il tracciamento del modo |
| risoluzione della griglia insufficiente | converge al 2% raddoppiando i punti |
| quasi-zeri di \(\psi\) nella finestra Kerr | zero minimi in \(20<r<50\) |
| autovalore sferoidale troncato al reale | correggerlo non cambia nulla |
| degrado monotono con l'overtone | **falso**, sweep sistematica con frequenze di Leaver |
| \(\mathcal E_M\) come predittore | **ridondante**: \(\mathcal E_M=c(n)|\Lambda_3|\), dispersione 0.15% |
| geometria di Stokes come capitolo | **preceduta**: Miyachi *et al.* 2025, PRD 111 |
| Madelung come teorema di non-informazione | **non lecito**: le ipotesi di Khesin–Misiołek–Modin sono violate dai QNM |

---

## 5. Disciplina di metodo — la parte che conta di più

In questo progetto **ogni risultato intermedio sbagliato è stato plausibile**, e
nessuno è stato scoperto rileggendo il ragionamento. Nove casi documentati in
`research/RESULTS_FOR_MANUSCRIPT.md` §8.

Il controllo efficace è sempre stato guardare **due cose invece di una**:

| coppia | perché una sola non basta |
|---|---|
| valore **e** stabilità | una serie divergente dà valori concordi fra parametri diversi: entrambi dominati dallo stesso termine |
| modulo **e** fase | due complessi di uguale modulo integrano in modo diverso |
| segno **e** ordine di grandezza | un ramo sbagliato dà il modulo giusto col segno rovesciato |
| statistica aggregata **e** dati grezzi | un residuo di fit piccolo non esclude che la successione salti |
| convergenza numerica **e** fisica | raffinare la griglia non corregge una frequenza sbagliata |

Caso emblematico: l'antiderivata asintotica mostrava scarto fra punti di raccordo
in calo fino a \(10^{-5}\) **mentre il valore divergeva di undici ordini**.
Guardando solo lo scarto si sarebbe riportato \(3\times10^{9}\).

**Corollario operativo:** l'indipendenza dal parametro di regolarizzazione non è
una verifica sufficiente. Mai.

---

## 6. Trappole specifiche di questo codice

1. **Frequenza.** \(Q_M\) amplifica di \(\sim10^2\) l'errore su \(\omega\)
   (Appendice C del manoscritto). Per Kerr usare
   `kerr_wkb3_selfconsistent.selfconsistent_frequency`, **non** quella eikonale:
   l'errore di quest'ultima è \(10^{-3}\) e basta a rovinare tutto.
2. **Bordo all'orizzonte.** La condizione entrante è esatta solo per \(V\to0\).
   La pendenza di \(Q_M\) passa da 0.569 a 2.026 spostando l'offset da
   \(10^{-3}\) a \(10^{-6}\); sotto \(10^{-7}\) l'integratore cede.
3. **Finestra di misura.** Al massimo di barriera \(|P|\to0\) e la gerarchia si
   rompe (§7 del manoscritto). Misurare **lontano dai turning point**,
   convenzionalmente \(20<r<50\).
4. **Ramo di \(\sqrt q\).** Va scelto coerentemente con la condizione al
   contorno su **ciascun lato**. Usarne uno solo dà un errore del 10% che *non*
   svanisce al crescere di \(L\), quindi imita un difetto di teoria.
5. **Ampiezza con \(\omega\) complessa.** Non è \(|u|^{-1/2}\): il termine
   \(-\varepsilon^{-1}\!\int\operatorname{Im}u\) in \(\ln A\) è \(O(1)\).
   Ometterlo sbaglia \(Q_M\) di un fattore ~50. Usare
   `madelung_wkb_prediction.py`, che è validato su Pöschl–Teller a
   \(1.7\times10^{-7}\).
6. **Cancellazione.** Le quantità QNM regolarizzate richiedono di cancellare
   integrali che crescono esponenzialmente. Oltre \(L_+\approx70\) la doppia
   precisione non basta: servirebbe `mpmath`.

---

## 7. Se invece si vuole lavorare al manoscritto

Aperto, in ordine di utilità:

1. **Figure 2 e 3.** La 2 schematizzerebbe dove vivono \(Q_{2j}\) e
   \(\Lambda_j\) (§6, oggi tre paragrafi di prosa); la 3 la validazione
   dell'Appendice C, oggi solo tabella. Fig. 1 è fatta:
   `Schw-QNM-WKB-Fluid/figures/make_fig1.py`.
2. **Rivista.** Non JMP — il teorema centrale è l'eq. (1.2) di DDP, *pubblicata
   su JMP*, e il nostro contributo è una classificazione misurata in contesto
   relativistico. CQG o PRD. PRD non ha categoria "Note": Letters (5 pagine) o
   Articles (senza limite).
3. **Rapporto contenuto/lunghezza.** ~14 pagine di cui ~4 originali. O si
   comprime (§3 e §4 a mezza pagina con rimando a [19], §7 in appendice), o si
   completa con lavoro *nostro*. Allungare con altra letteratura peggiora.
4. **Citazioni.** 34 riferimenti, 13 citati nel testo. In stesura finale gli
   altri vanno agganciati o tolti.
