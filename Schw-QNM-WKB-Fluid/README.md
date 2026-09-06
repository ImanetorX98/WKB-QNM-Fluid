# Schw-QNM-WKB-Fluid

Preparazione di un manoscritto per rivista di fisica matematica sul potenziale
di Madelung come chiusura della gerarchia WKB per l'equazione di tipo
Schrödinger di un buco nero di Schwarzschild.

- `manuscript.md` — la bozza.
- `verification/` — gli script che sostengono ogni enunciato quantitativo.
- `references.bib` — bibliografia con DOI verificati via Crossref.
- `papers/` — gli otto PDF ad accesso libero, con testo estratto.
- `codex-download-manifest.tsv` — i 14 articoli da recuperare altrove.

## Tesi in una riga

La serie WKB si chiude su **un solo** funzionale di Bohm–Madelung, valutato
sull'ampiezza WKB completa; il funzionale compare dal **terzo termine della
serie di fase** in poi e da lì genera tutti gli ordini successivi.

## Riproduzione

Interprete: `python3.13` (il `python3` di sistema non ha numpy).

```
cd verification
python3.13 madelung_recursion.py      # chiusura, parita', ricorsione u2 u4
python3.13 scalar_eikonal_scaling.py  # pendenze in eps per s=0,1,2
python3.13 open_continuity.py         # legge di continuita' con sorgente
cd ../../core
python3.13 dirac_madelung_profile.py --kappa 4 --scaling   # contrappunto Dirac
```

## Mappa enunciato → verifica

| § | Enunciato | Script | Numero chiave |
|---|---|---|---|
| 3 | Teorema 1, chiusura esatta | `madelung_recursion.py` (a) | identità simbolica |
| 4 | Teorema 2, niente ordini dispari | `madelung_recursion.py` (b) | $u_1\equiv0$ |
| 4 | Eq. (4.1)–(4.2), $u_2$ e $u_4$ | `madelung_recursion.py` (c,d) | accordo con WKB classica |
| 4 | Corollario 3, singolo funzionale | `madelung_recursion.py` (e) | $u_4$ rigenerato |
| 5 | Proposizione 5, spin a ordine $\varepsilon^2$ | `scalar_eikonal_scaling.py` | pendenze 2.0000 |
| 7 | Rottura ai turning point | `dirac_madelung_profile.py --scaling` | 1.84 → 1.40 |
| 8 | Fluido aperto | `open_continuity.py` | residuo $\sim10^{-6}$ |
| 9 | Connessione di spin a ordine $\varepsilon^1$ | `dirac_madelung_profile.py --scaling` | pendenza 1.0000 |

## Nota sui conteggi di ordine

Nel manoscritto convivono tre nozioni distinte di "ordine" (§2.1): il numero del
termine nella serie di fase, la potenza di $\varepsilon$, e l'ordine WKB di
barriera di Iyer–Will usato altrove in questo repo. "Terzo ordine" nel titolo
della tesi è il **primo** dei tre. Confonderli è il modo più rapido di
sbagliare un enunciato in questo ambito.

## Stato della bibliografia

Verificati contro Crossref: 20 DOI su 22 voci. Le due eccezioni sono
Chandrasekhar 1983 (monografia, solo ISBN) e Voros 1983 (NUMDAM legacy, nessun
DOI mai assegnato; identificativo stabile `AIHPA_1983__39_3_211_0`).

Tre correzioni emerse dalla verifica, tutte già applicate al manoscritto:

1. "Iyer & Will 1987" sono **due** articoli distinti. La sistematica di ordine
   superiore è la Parte I (Iyer & Will, PRD 35, 3621); la Parte II (Iyer da
   solo, PRD 35, 3632) applica il metodo a Schwarzschild ed è la fonte dei
   valori di riferimento dei test.
2. *Analogue gravity* ha un'edizione **2026** (Living Rev. Rel. 29) che
   sostituisce quella 2011 citata nella bozza iniziale.
3. L'universalità del termine eikonale di testa è un risultato noto
   (Cardoso et al. 2009), non nostro. La Proposizione 5 è stata riscritta per
   rivendicare solo ciò che è effettivamente nuovo: la coincidenza di ordine
   fra spin e Madelung.

## Cosa il lavoro NON afferma

Nessuna previsione osservativa. La trasformazione di Madelung è esatta e non
sposta lo spettro dei QNM, quindi non tocca alcuna forma d'onda. L'unico
contesto in cui il termine è fisicamente microscopico è quello degli analoghi
acustici in condensati di Bose–Einstein, e la corrispondenza fra i due contesti
non è stabilita qui.
