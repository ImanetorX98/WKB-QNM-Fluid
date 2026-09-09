# Risultati consolidati per la stesura del manoscritto

**Aggiornato:** 9 settembre 2026.
**Scopo:** documento unico da consultare in stesura. Ogni voce dichiara se è
**derivata**, **misurata**, **nota in letteratura** o **corretta rispetto a una
nostra affermazione precedente**.

**Interprete:** `python3.13`. Il `python3` di sistema non ha numpy.

---

## 1. La tesi, nella forma che ha resistito ai controlli

> La riscrittura di Madelung è una trasformazione esatta e non può aggiungere
> informazione allo spettro. Il suo contenuto informativo è **esattamente**
> quello della WKB da cui proviene. L'unico punto in cui dice qualcosa di non
> ridondante è l'**ordinamento**: a quale potenza di \(\varepsilon\) compare il
> primo termine subprincipale, e se il termine di ampiezza sia solo o preceduto
> da struttura geometrica.

Due risultati la sostengono: uno positivo (criterio di ordinamento, §3) e uno
negativo (ridondanza del diagnostico, §5).

---

## 2. Che cosa NON rivendichiamo

Da mettere esplicitamente in §1.1 del manoscritto. In quest'area la letteratura
è densa e tre direzioni promettenti sono già state chiuse da precedenti.

| Elemento | Stato |
|---|---|
| Chiusura \(q=u^2-\varepsilon^2(u^{-1/2})''/u^{-1/2}\) | **Nota**: formalismo di Hamilton–Jacobi quantistico, momento quantistico su Riccati |
| Coefficienti \(u_2\), \(u_4\) | **Noti**: WKB classica (Bender–Orszag cap. 10) |
| Universalità del termine eikonale di testa | **Nota**: Cardoso *et al.* 2009 |
| Ampiezza–fase per QNM, Wronskiano, Riccati | **Noti**: Glampedakis–Andersson 2003 |
| WKB esatta dei QNM, geometria di Stokes | **Nota e attiva**: Aminov–Grassi–Hatsuda 2022; Miyachi *et al.* 2025 |
| Idrodinamica di Dirac con QNM su Schwarzschild | **Precedente diretto**: Meza-Domínguez–Matos 2026 (arXiv:2605.28887) |
| \(A_1\neq0\) per l'autovalore sferoidale | **Noto**: letteratura eikonale su Kerr (Dolan 2010; Yang *et al.* 2012) |
| Stimatori d'errore per QNM | **Noti**: \(\Delta_k\) di Konoplya *et al.* 2019; residui pesati |

---

## 3. Risultato positivo: il criterio di ordinamento

### 3.1 Enunciato (derivato, e **corretto** rispetto alla prima versione)

> Il criterio non distingue spin intero da semintero. Distingue **presenza o
> assenza di struttura geometrica che precede il termine di ampiezza**. Se ne
> conoscono due sorgenti indipendenti: la connessione di spin dei partner di
> Darboux (Dirac, statico) e la separazione sferoidale indotta dalla rotazione
> (Kerr, bosonico). Entrambe danno pendenza 1.

La prima formulazione attribuiva la pendenza 1 al **settore fermionico**. Kerr
scalare è bosonico e ha un termine \(\varepsilon^1\): l'attribuzione era troppo
stretta. La correzione va riportata nel manoscritto (§9).

### 3.2 Le tre decomposizioni esatte

Tutte e tre hanno la stessa forma, ed è questo il punto del capitolo.

**Schwarzschild bosonico** (derivato), con Langer \(L=\ell+\tfrac12\):
\[
\frac{V_s}{L^2}=h^2+\varepsilon^2v_2,\qquad
h^2=\frac{f}{x^2},\quad v_2=f\left[\frac{2(1-s^2)}{x^3}-\frac1{4x^2}\right].
\]
Nessun termine \(\varepsilon^1\). Lo spin entra solo in \(v_2\), **degenere in
ordine** con il Madelung.

**Dirac su Schwarzschild** (derivato), \(K=|\kappa|\), \(\tau=\sigma\,\mathrm{sgn}\,\kappa\):
\[
\frac{V_\tau}{K^2}=h^2+\tau\varepsilon h',\qquad h=\frac{\sqrt f}{x}.
\]
Stesso \(h^2\) di testa; correzione a \(\varepsilon^1\), non \(\varepsilon^2\).

**Kerr scalare** (derivato), \(H=r^2+a^2\), \(\hat c=a\hat\Omega\):
\[
\frac{q}{L^2}=Q_0(r)-\varepsilon\frac{\Delta A_1}{H^2}
-\varepsilon^2\left[\frac{\Delta A_2}{H^2}+\frac{h''}{h}\right],\quad
Q_0=\left(\hat\Omega-\frac{\mu a}{H}\right)^2-\frac{\Delta\bar A_0}{H^2}.
\]

### 3.3 Pendenze misurate

| Caso | Termine \(\varepsilon^1\) | Pendenza misurata | Script |
|---|---|---|---|
| Schwarzschild \(s=0,1,2\) | assente | \(\varepsilon^2v_2\): **2.0000**; \(Q_M\): 1.85–1.91 | `Schw-QNM-WKB-Fluid/verification/scalar_eikonal_scaling.py` |
| Dirac, lontano dai turning point | \(\tau\varepsilon h'\) | **1.0000** (precisione macchina) | `core/dirac_madelung_profile.py --scaling` |
| Kerr, \(a=0\) (controllo) | assente | **1.951** | `calculations/kerr_radial_order_profile.py` |
| Kerr, \(a=0.3\ldots0.9\), \(\mu=0.5,0.9\) | \(-\varepsilon\Delta A_1/H^2\) | **0.963 – 1.014** | idem |

Tabella Kerr completa in [kerr_eikonal_order](kerr_eikonal_order_2026-09-08.md) §5.

### 3.4 Il rapporto spin/Madelung nel settore bosonico (misurato)

A \(n=0\), finestra \(20<x<50\), \(\ell=2\ldots32\):

| \(s\) | \(v_2\) | pendenza \(\varepsilon^2v_2\) | pendenza \(Q_M\) | \(\varepsilon^2v_2/Q_M\) a \(\ell=2\) |
|---|---|---|---|---|
| 0 | \(f[+2/x^3-1/(4x^2)]\) | 2.0000 | 1.907 | 0.035 |
| 1 | \(f[\ \ 0\ \ -1/(4x^2)]\) | 2.0000 | 1.894 | 0.061 |
| 2 | \(f[-6/x^3-1/(4x^2)]\) | 2.0000 | 1.854 | **0.150** |

Il caso \(s=1\) è degenere in modo notevole: \(v_2\) si riduce al **puro termine
di Langer**, senza curvatura. Il gravitazionale è quello in cui spin e Madelung
competono più da vicino.

Il difetto di \(Q_M\) rispetto a 2 **non è fisico**: è contaminazione dal ramo
riflesso alla condizione ingoing. Spingendo il bordo verso l'orizzonte la
pendenza sale monotonamente: \(x_{\min}=\) 2.02 → 1.557; 2.001 → 1.749;
2.0001 → 1.839; 2.00005 → 1.865. Conferma indipendente al §3.5.

### 3.5 Conferma incrociata dall'overtone (misurato)

A overtone crescente \(Q_M\) cresce di due ordini e la contaminazione diventa
irrilevante: la pendenza va a **2.0001** a \(n=3\). Questo conferma
indipendentemente che il difetto del §3.4 era numerico, non strutturale.

---

## 4. Limiti strutturali (misurati, attesi ma qui quantificati)

### 4.1 Rottura alla coalescenza dei turning point

| Regione | \(\min|P|\) | pendenza \(\varepsilon^1\) | pendenza \(Q_M\) |
|---|---|---|---|
| lontano dai turning point | \(\approx0.187\) | **1.0000** | 1.84 |
| al massimo di barriera | \(2\times10^{-5}\)–\(2\times10^{-4}\) | 1.14 | **1.40** |

La serie locale è singolare dove i due turning point coalescono — che è dove la
WKB di barriera costruisce i \(\Lambda_j\). Ragione ulteriore per cui non esiste
corrispondenza \(Q_{2j}\leftrightarrow\Lambda_j\).

### 4.2 Il fluido di un QNM è aperto

\[
A''-AS'^2+(E-V)A=0,\qquad (\rho v)'=-\Gamma\rho,\quad \Gamma=\operatorname{Im}\Omega^2 .
\]

| \(\ell\) | \(\Gamma\) | residuo della legge con sorgente | violazione della legge conservata |
|---|---|---|---|
| 2 | \(-0.0666\) | \(8.5\times10^{-7}\) | 1.000 |
| 16 | \(-0.6073\) | \(2.2\times10^{-6}\) | 1.000 |

Nessuna terza via conserva insieme densità reale positiva e continuità.
**Nota**: in Vaidya la legge caratteristica è invece **conservativa** — vedi §6.

### 4.3 Zeri di \(\psi\) e overtone: risultato NEGATIVO

Con frequenze esatte di Leaver, il numero di quasi-zeri **non è monotono** in
\(n\): \(\ell=16\) dà 1, 6, 6, 3, 2, 2, 2 per \(n=0\ldots6\), e `fd_l1` salta in
modo sporadico. **L'ipotesi di un degrado monotono con l'overtone è falsa.**

Sopravvive solo che \(n=0\) è sistematicamente diverso: un minimo poco profondo
(0.87) e `fd_l1` \(\sim10^{-4}\), contro 3–6 minimi e fino a \(10^0\) per
\(n\ge1\). La sporadicità è coerente con la geometria di Stokes, che è
letteratura consolidata e che **non tocchiamo**.

---

## 5. Risultato negativo: il diagnostico di Madelung è ridondante

Dettaglio completo in [indicator_verdict](indicator_verdict_2026-09-08.md).

**Il risultato decisivo (misurato):**

| overtone | \(\mathcal E_M/|\Lambda_3|\), \(s=0,2\), \(\ell=3\ldots8\) | dispersione |
|---|---|---|
| \(n=0\) | 1.0426 | **0.15%** |
| \(n=1\) | 3.95 | 1.6% |

A overtone fissato l'indicatore è una **riscalatura costante di \(|\Lambda_3|\)**,
che la WKB3 calcola gratis dal getto di \(V\) al massimo. \(\mathcal E_M\)
richiede invece frequenza esatta, dati di Frobenius e integrazione dell'ODE.

Supporto:
- **Pöschl–Teller, \(n=0\)**: identità esatte \(|\Omega/L|^2=1\) e
  \(P^2=(1-\varepsilon^2/4)\tanh^2y\), da cui
  \(\mathcal E_M=\varepsilon^2(2-c_w)/(8-\varepsilon^2c_w)\): funzione del solo
  \(\varepsilon\), verificata contro il codice a \(4\times10^{-17}\).
- **Il segnale apparente a \(n=1\) è il nodo**: l'integrale è dominato all'89%
  dal picco, e con la mediana pesata lo spread crolla da 85.3% a **1.9%**.
- **A \(\varepsilon\) fissato** (spin variato a \(\ell\) fissato, 18 coppie):
  \(\mathcal E_M\) 18/18 con errore 0.0179; \(\Delta_2\) 18/18 con 0.0185;
  \(|\Lambda_3|\) 18/18 con 0.0187. **Indistinguibili.**

**Nota utile**: dato \(\Delta_k=|\omega_{k+1}-\omega_{k-1}|/2\), la quantità
\(|\omega_3-\omega_1|/2\) **è** \(\Delta_2\), non un surrogato. Il confronto è
già disponibile senza implementare nulla; solo \(\Delta_3\) richiede WKB2 e WKB4.

---

## 6. Vaidya: stato e ipotesi nulla

**Letteratura (misurata su arXiv):** Vaidya+quasinormal 15 lavori;
Vaidya+WKB 2, **nessuno sui QNM**; Vaidya+Madelung/Bohm **0**;
quasinormal+Madelung **0**.

**I quattro lavori pertinenti:**

| | contenuto | rapporto con noi |
|---|---|---|
| Abdalla–Chirenti–Saa 2006 | primi QNM su Vaidya, numerico caratteristico; **rottura del regime adiabatico** con «effetto inerziale» (overshoot di \(\omega_R\)) | fenomeno bersaglio |
| Lin–Sun–Zhang 2021 | metodo matriciale, frequenze diverse fra orizzonte e infinito | contesto |
| Capuano–Santoni–Barausse 2024 | \(\dot M\) costante, cambio di coordinate → dominio delle frequenze, eikonale **solo all'ordine di testa**, Leaver | il loro cambio di coordinate **elimina** la struttura mista \(\partial_v\partial_r\) che sarebbe il nostro oggetto; il loro \(\epsilon\) è il tasso di massa, non \(1/L\) |
| Yoo *et al.* 2025/26 | limite di Penrose, fit di fase \(\omega=dS/dV\) sulla forma d'onda | fit diagnostico, non decomposizione dell'equazione |

**Struttura già verificata simbolicamente** (`calculations/vaidya_madelung_symbolic.py`):
\[
2S_vS_r+fS_r^2+\varepsilon^2U_\ell+Q_V=0,\qquad
Q_V=Q_r+Q_\times,\quad Q_\times=-2\varepsilon^2\frac{A_{vr}}{A},
\]
con legge caratteristica **conservativa**
\(\partial_v(\rho S_r)+\partial_r[\rho(S_v+fS_r)]=0\).

Il termine \(Q_\times\) **non ha analogo statico** e non compare in nessuno dei
quattro lavori.

**Ipotesi nulla da battere (derivata).** Se \(A(v,r)=A_0(r;M(v))\),
\[
Q_\times=-2\varepsilon^2\dot M\,\frac{\partial_M\partial_rA_0}{A_0}+O(\dot M^2,\ddot M),
\]
cioè nel regime adiabatico \(Q_\times\) è determinato da \(\dot M\) e dal profilo
congelato: **ridondante come \(\mathcal E_M\)**. Il contenuto eventualmente nuovo
è tutto nel resto non adiabatico
\(Q_\times^{\rm na}=Q_\times+2\varepsilon^2\dot M\,\partial_M\partial_rA_0/A_0\),
e la domanda è se tracci l'effetto inerziale di Abdalla–Chirenti–Saa.

**Non ancora fatto:** l'evoluzione Vaidya a \(\dot M\) costante nel dominio
caratteristico, con estrazione di \(A\), \(S\) e del resto non adiabatico.

---

## 7. Infrastruttura validata

| Modulo | Contenuto | Validazione |
|---|---|---|
| `core/leaver_qnm.py` | frazioni continue di Leaver, ricorsione **derivata simbolicamente** | 10 cifre contro i valori tabulati; residuo CF \(10^{-14}\) |
| `core/schwarzschild_wkb.py` | WKB1/WKB3 Iyer–Will, \(s=0,1,2\) | Iyer–Will \(0.3732-0.0892i\) |
| `core/dirac_madelung_profile.py` | partner di Darboux, gerarchia spin/Madelung | Cho 2003 per \(\kappa=1,2\) |
| `calculations/kerr_eikonal_order_test.py` | autovalore sferoidale, \(c\) complessa | Langer esatto nel limite sferico |
| `calculations/kerr_radial_order_profile.py` | decomposizione radiale di Kerr | \(1/(3\sqrt3)\), \(r=3\), orbite fotoniche |

**Ricorrenza di Leaver, derivata e non trascritta:**
\[
\alpha_n=(n+1)(n+1-4i\omega),\quad
\gamma_n=(n-4i\omega)^2-s^2,
\]
\[
\beta_n=-\ell(\ell+1)-2n^2-2n+16in\omega+32\omega^2+8i\omega+s^2-1.
\]

**Test:** 12 (calculations base) + 7 (Kerr angolare) + 6 (Kerr radiale)
+ 20 (core) = **45**, tutti superati.

---

## 8. Errori di metodo, e il controllo che li ha presi

Non è un elenco di scuse: è il catalogo dei controlli che funzionano. In questo
progetto **ogni singolo risultato intermedio sbagliato è stato plausibile**, e
nessuno è stato scoperto rileggendo il ragionamento. Sono stati scoperti tutti
misurando qualcosa in più.

### 8.1 Il pattern

Il controllo efficace è sempre stato guardare **due cose invece di una**. Un
solo indicatore, per quanto buono, si lascia soddisfare da un risultato
sbagliato.

| coppia da controllare | perché uno solo non basta |
|---|---|
| **valore** e **stabilità** | una serie divergente può dare valori concordi fra parametri diversi, perché entrambi sono dominati dallo stesso termine divergente |
| **modulo** e **fase** | due complessi con lo stesso modulo integrano in modo diverso |
| **segno** e **ordine di grandezza** | un ramo sbagliato dà il modulo giusto col segno rovesciato |
| **statistica aggregata** e **dati grezzi** | un residuo di fit piccolo non esclude che la successione salti |
| **convergenza numerica** e **convergenza fisica** | raffinare la griglia non corregge una frequenza sbagliata |

### 8.2 Catalogo

| # | affermazione sbagliata | come è stata presa |
|---|---|---|
| 1 | degrado monotono con l'overtone | sweep sistematica contro quattro casi scelti *dopo* aver visto i dati |
| 2 | errore WKB3 stimato dallo scarto WKB1↔WKB3 | Leaver: 500× di sovrastima (5.4% contro \(1.1\times10^{-4}\)) |
| 3 | tracciamento autovalori per prossimità | successioni grezze erratiche: \(A_1=-12.1\) invece di \(-0.029\) |
| 4 | solutore eikonale convergente | \(r_{\rm picco}\sim2500\), fuori dall'orizzonte di tre ordini |
| 5 | \(Q_M\) dalla sola ampiezza \(|u|^{-1/2}\) | fattore ~50, corretto includendo l'esponenziale da \(\operatorname{Im}u\) |
| 6 | stesso ramo di \(\sqrt q\) sui due lati | errore del 10% **che non svaniva** al crescere di \(L\) |
| 7 | \(Q_M\) erratico su Kerr | tre ipotesi testate e scartate prima di trovare la frequenza |
| 8 | antiderivata asintotica convergente | scarto \(\to10^{-5}\) mentre il **valore** divergeva di 11 ordini |
| 9 | integrando che coincide con la serie | verificato sul modulo; il rapporto complesso era l'informazione utile |

### 8.3 Casi in cui la diagnosi stessa era sbagliata

Vale la pena registrarli separatamente, perché sono i più insidiosi: l'errore era
reale, ma la spiegazione che ne davo no.

- **Troncamento della matrice** (§ Kerr): ipotizzato come causa dell'erraticità
  degli autovalori. Era convergente a \(10^{-8}\); la causa era il
  tracciamento.
- **Profondità dei quasi-zeri mal condizionata**: diagnosticata come ostacolo che
  richiedeva Leaver. L'errore vero di WKB3 era 500 volte più piccolo di come
  l'avevo stimato, e le profondità erano stabili.
- **Quasi-zeri di \(\psi\) nella finestra Kerr**: ipotesi plausibile,
  smentita — zero minimi in \(20<r<50\).

### 8.4 Ciò che resta aperto con la stessa disciplina

La discrepanza di \(1.65\times10^{-4}\) del §5.2 di
[vaidya_solvability](vaidya_solvability_2026-09-09.md) **non è stata risolta**.
Ogni ingrediente è verificato a \(10^{-7}\) o meglio, e i pezzi non tornano
insieme. Non si prosegue senza un'ipotesi precisa: è esattamente la situazione in
cui, sopra, sono nati i nove casi.

---

## 9. Struttura proposta del manoscritto

1. Impostazione, convenzioni, i tre conteggi di «ordine» (O1 serie di fase,
   O2 potenza di \(\varepsilon\), O3 barriera di Iyer–Will).
2. Chiusura esatta e ricorsione — **richiamate, non rivendicate**.
3. Schwarzschild: decomposizione di Langer, degenerazione spin–Madelung,
   pendenza 2.
4. **Kerr: la rotazione dà pendenza 1 nel settore bosonico.** Corregge il §3.
5. Dirac: la connessione di spin dà pendenza 1 nel settore statico.
6. Il criterio generale, che unifica 4 e 5.
7. Limiti: coalescenza dei turning point; fluido aperto.
8. Ridondanza del diagnostico: \(\mathcal E_M\propto|\Lambda_3|\).
9. Vaidya: struttura, \(Q_\times\), ipotesi nulla adiabatica — **come programma
   dichiarato**, non come risultato, finché non c'è l'evoluzione.

I capitoli 3–5 non ripetono la stessa cosa: il 4 falsifica una lettura ingenua
del 3, e il 9 introduce una struttura che nei precedenti non esiste.

---

## 9bis. Aggiornamenti del 9 settembre

**Fondazione bibliografica** (dettaglio in [sources_read](sources_read_2026-09-09.md)):
Iyer & Will I fonda i punti (i) e (iii) del §6; Berry & Mount fonda l'unicità di
Langer (Appendice B), il principio di equivalenza dei turning point (§6) e il
meccanismo della rottura (§7); Seidel & Iyer IV **non** precede il §9, perché
espande in $a\omega$ e non in $1/L$; Khesin–Misiołek–Modin dà la
caratterizzazione geometrica e identifica $Q_M$ con l'informazione di Fisher.

**Previsione analitica e condizionamento**
(dettaglio in [kerr_madelung_sensitivity](kerr_madelung_sensitivity_2026-09-09.md)):
il Teorema 2 vale solo per $q$ indipendente da $\varepsilon$, e Kerr lo viola
($u_1=q_1/2\sqrt{q_0}\neq0$). Formula chiusa per $Q_M$ validata su
Pöschl–Teller a $1.7\times10^{-7}$. **$Q_M$ amplifica di $\sim10^2$ l'errore
sulla frequenza**: nessun funzionale costruito su di esso può prevederla.

**Vaidya** (dettaglio in [vaidya_adiabatic](vaidya_adiabatic_2026-09-09.md)):
il residuo dell'ansatz adiabatico a ordine $\dot M$ è esattamente
$2\,\partial_r\partial_M Z$. Il termine misto **è il difetto di adiabaticità**,
calcolabile dalla famiglia congelata. Terzo caso della tesi «riproduce, non
estende», in forma analitica e senza evoluzione numerica.

**Passo B chiuso** (dettaglio in [kerr_selfconsistent](kerr_selfconsistent_2026-09-09.md)):
frequenza autoconsistente col potenziale a quel \(\ell\), risolta per
\((\omega,r_0)\) complessi. A \(a=0\) riproduce il WKB3 scalare a
\(3\times10^{-11}\). L'errore previsione/numerica su Kerr scende da
\(2.3\times10^{-1}\) a \(1.5\times10^{-3}\): **Kerr è ora trattato
completamente**, su potenziale (§9.2) e su ampiezza (§9.6). Nuovo §9.6.

**Figura 1** in `Schw-QNM-WKB-Fluid/figures/`, generata da `make_fig1.py`:
log-log dei tre termini subprincipali contro \(\varepsilon\), ancorati a
\(\varepsilon=0.1\) perché resti visibile solo la pendenza. Palette Okabe-Ito
validata per daltonismo (peggior coppia \(\Delta E=11.0\) in deuteranopia),
identità codificata anche da marcatore e etichetta diretta.

**Regolarizzazione di Leung implementata e validata**
(dettaglio in [vaidya_solvability](vaidya_solvability_2026-09-09.md)): serie
asintotica uscente con residuo \(10^{-12}\), \(D_+\) sul ramo corretto, norma
generalizzata indipendente da \(L_+\) allo 0.2%. Il denominatore della (2.14) è
fatto; manca il numeratore.

**Ancora aperto:** il numeratore della (2.14), cioè i termini \(\Delta_\pm\)
dall'asintotica della sorgente \(-2\partial_r\partial_M Z\) (il confronto previsione/numerica su Kerr
richiede una frequenza autoconsistente al potenziale, non quella eikonale);
DDP 1997 da leggere per il §12.

## 10. Documenti collegati

- [indicator_verdict_2026-09-08.md](indicator_verdict_2026-09-08.md) — verdetto su \(\mathcal E_M\)
- [kerr_eikonal_order_2026-09-08.md](kerr_eikonal_order_2026-09-08.md) — Kerr, angolare e radiale
- [novelty_publication_assessment_2026-09-08.md](novelty_publication_assessment_2026-09-08.md) — precedenti e protocollo
- [spectral_audit_2026-09-08.md](spectral_audit_2026-09-08.md) — audit e formule chiuse
- [sources_read_2026-09-09.md](sources_read_2026-09-09.md) — fonti primarie lette e loro effetto
- [kerr_madelung_sensitivity_2026-09-09.md](kerr_madelung_sensitivity_2026-09-09.md) — previsione analitica e condizionamento
- [vaidya_adiabatic_2026-09-09.md](vaidya_adiabatic_2026-09-09.md) — residuo adiabatico
- [covariant_formalism.md](covariant_formalism.md) — formalismo statico → Kerr → Vaidya
- `../Schw-QNM-WKB-Fluid/manuscript.md` — bozza corrente
