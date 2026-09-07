# WKB, QNM e potenziale di Madelung

Progetto unico dopo la fusione di `WKB-Fluid/WKB-QNM-Schw-Fluid` (numerica
bosonica) in questa cartella (teoria + Dirac). Unità \(G=c=\hbar=1\),
\(x=r/M\), \(x_*=r_*/M\), dipendenza temporale \(e^{-i\omega t}\), quindi un
QNM smorzato ha \(\operatorname{Im}\omega<0\).

## Struttura

| Percorso | Contenuto |
|---|---|
| `core/schwarzschild_wkb.py` | QNM Schwarzschild WKB1/WKB3 Iyer–Will, spin 0/1/2 |
| `core/madelung_profile.py` | profilo di Madelung numerico, caso scalare |
| `core/dirac_madelung_profile.py` | profilo di Madelung numerico, Dirac massless |
| `core/scan_qnms.py` | tabella WKB1 vs WKB3 |
| `core/test_wkb.py`, `core/test_dirac_madelung.py` | test |
| `core/README-wkb-fluid.md` | derivazione del ramo bosonico |
| `calculations/` | verifiche simboliche Dirac (Schwarzschild e Kerr proiettivo) |
| `calculations/dirac_wkb_hydrodynamics_notes.md` | nota di calcolo Dirac |
| `research/covariant_formalism.md` | base covariante statico → Kerr → Vaidya |
| `research/RESEARCH_STATE.md` | documento autosufficiente per riprendere risultati e prossimi conti |
| `research/prior_art.md` | mappa aggiornata dei precedenti e confine di originalità |
| `research/references.bib` | bibliografia verificata del nuovo programma |
| `research/static_benchmark.md` | benchmark di \(\mathcal E_M\) contro Leaver per 69 modi statici |
| `research/exact_barrier_benchmark.md` | controllo Pöschl–Teller, residui per ordine e limite nodale |
| `calculations/vaidya_madelung_symbolic.py` | verifica KG e split Madelung su Vaidya |
| `calculations/static_madelung_benchmark.py` | scansione Schwarzschild indipendente WKB–Leaver–Madelung |
| `calculations/order_resolved_madelung_benchmark.py` | residui locali \(\mathfrak R_0,\mathfrak R_1,\mathfrak R_2\) |
| `calculations/poschl_teller_madelung_benchmark.py` | benchmark QNM esattamente risolvibile fuori Schwarzschild |
| `calculations/uniform_madelung_defect.py` | residuo di Madelung uniforme nello strato parabolico-cilindrico |
| `calculations/kerr_scalar_madelung_symbolic.py` | riduzione scalare Kerr esatta e ponte slow-rotation |
| `calculations/verify_formalism.wl` | controlli simbolici con kernel Mathematica/Wolfram |
| `report-source.md` | rapporto di ricerca, sorgente del PDF |
| `build_report.py`, `build_dirac_note.py` | generatori PDF (`ROOT = parent` del file: tenerli qui) |
| `output/pdf/` | PDF prodotti |
| `tmp/` | render intermedi e paper di riferimento |

**Interprete:** usare `python3.13` (o `/usr/local/bin/python3`). Il `python3`
di default sul sistema non ha numpy.

```
cd core
python3.13 -m unittest test_wkb.py test_dirac_madelung.py
```

## Il ponte fra i due rami

Entrambi i rami calcolano lo stesso oggetto per vie diverse. Il ramo scalare
integra \(\psi\) e decompone \(\psi=Ae^{iS}\); il ramo Dirac fa lo stesso sui
partner di Darboux, ma **dopo** la scalatura eikonale, che è il punto in cui
la struttura cambia.

Scalare: \(\tfrac12(S')^2+\tfrac12 V+Q_M=\tfrac12\operatorname{Re}\Omega^2\),
con \(Q_M=-A''/(2A)\).

Dirac, con \(K=|\kappa|\), \(\varepsilon=1/K\), \(\omega=K\Omega\),
\(h=\sqrt f/x\), \(\tau=\sigma\operatorname{sgn}\kappa\):

\[
\operatorname{Re}\Omega^2=P^2+h^2+\underbrace{\tau\varepsilon h'}_{O(\varepsilon)}+\underbrace{Q_M}_{O(\varepsilon^2)},
\qquad Q_M=-\varepsilon^2\frac{A''}{A}.
\]

Entrambe le chiusure sono **identità algebriche**: valgono a precisione
macchina per costruzione (residuo misurato \(\sim10^{-18}\)) e non sono di per
sé una verifica. Il controllo indipendente è ricostruire \(A''/A\) per
differenze finite dall'ampiezza integrata: `fd_residual` nel modulo Dirac, che
si assesta a \(\sim10^{-7}\).

## Il criterio di riconoscimento

La domanda "come riconosciamo il regime di Madelung" ha una risposta
operativa: **guardare le pendenze in \(\varepsilon\)**, non i valori assoluti.

```
cd core
python3.13 dirac_madelung_profile.py --kappa 4 --scaling
```

Risultato misurato (\(\tau=+1\), \(n=0\), \(K=2,4,8,16\)):

| Regione | pendenza spin | pendenza \(Q_M\) | \(\min|P|\) |
|---|---|---|---|
| lontano dai turning point (\(20<x<50\)) | 1.0000 | 1.84 | ≈ 0.187 |
| al massimo di barriera | 1.14 | 1.40 | ≈ 10⁻⁵ |

Tre letture:

1. **Il termine di spin connection è genuinamente di ordine \(\varepsilon\)**,
   pendenza 1 a precisione macchina — è analitico, \(\tau\varepsilon h'\), e
   coincide con metà della separazione fra i partner: \(V_+-V_-=2Kh'\).
   Nel caso scalare questo termine **non esiste**. È la firma che distingue
   Dirac, e conferma numericamente la conclusione della nota: la gerarchia
   efficace fermionica non è una successione di soli potenziali di Madelung.

2. **Lontano dai turning point \(Q_M\) segue \(\varepsilon^2\)** (1.84 misurato
   contro 2 atteso). Il difetto residuo non è fisico: è contaminazione dal ramo
   riflesso, perché la condizione ingoing \(Z\sim e^{-i\Omega x_*}\) è esatta
   solo per \(V\to0\). Spingendo `--x-min` verso l'orizzonte la pendenza sale
   monotonamente: 2.02 → 1.557; 2.001 → 1.749; 2.0001 → 1.839; 2.00005 → 1.865.

3. **Al massimo di barriera la gerarchia si rompe.** Lì \(|P|\to0\): i due
   turning point coalescono, e \(Q_M\) smette di seguire \(\varepsilon^2\)
   (pendenza 1.40). È esattamente il punto 3 del rapporto — la serie locale in
   potenze inverse di \(q\) è singolare nella zona di coalescenza e va
   sostituita dalla forma normale parabolico-cilindrica — qui però *misurato*,
   non solo argomentato. Ed è la regione in cui la WKB di barriera costruisce
   i \(\Lambda_j\): un'altra ragione per cui non esiste corrispondenza
   \(Q_{2j}\leftrightarrow\Lambda_j\).

## Validazione

- `core/test_wkb.py`: WKB3 gravitazionale \(\ell=2,n=0\) contro Iyer–Will
  \(0.3732-0.0892i\).
- `core/test_dirac_madelung.py`: WKB3 Dirac contro Cho (2003) per \(K=1,2\);
  isospettralità dei partner; chiusura esatta; accordo differenze finite;
  le tre pendenze sopra.

Tutti e 10 + 4 i test passano.
