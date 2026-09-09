# WKB-QNM-Fluid: aggiornamento per Claude

**Data:** 8 settembre 2026.
**Scopo:** consentire una revisione indipendente senza la cronologia della chat.
Questo file è preparato per essere condiviso dall'utente; non è già stato
trasmesso a Claude e non costituisce un parere ottenuto da Claude.

## 1. Obiettivo dell'utente e stato effettivo

L'utente vuole risultati originali e pubblicabili nel collegamento
fra QNM, WKB e potenziale di Madelung, con progressione:
spazi-tempi statici → stazionari non statici/Kerr → non stazionari/Vaidya.
Vuole conservare i risultati in Markdown e verificare l'algebra con
Mathematica, disponibile con licenza sul suo Mac.

La repository pubblica è [WKB-QNM-Fluid](https://github.com/ImanetorX98/WKB-QNM-Fluid).
Il progetto locale è:

~~~text
/Users/iman.rosignoli/Desktop/ISTRUZIONE/PROGRAMMAZIONE/WKBQNMMadelungFluid
~~~

Base Git dell'audit: ec97d56, branch main. L'audit e questi documenti sono
modifiche locali non ancora committate/pubblicate. Controllare lo stato
effettivo prima di lavorare; non sovrascrivere modifiche dell'utente.

**Verdetto corrente:** formalismo e benchmark utili, ma nessuno stimatore
economico nuovo è ancora validato. Originalità e pubblicabilità non sono
dimostrate. Non riprendere le formulazioni più ottimistiche delle prime note
senza le correzioni qui elencate.

## 2. Letture in ordine

1. [Audit spettrale dell'8 settembre](research/spectral_audit_2026-09-08.md).
2. [Valutazione di originalità e pubblicabilità](research/novelty_publication_assessment_2026-09-08.md).
3. [Stato esteso della ricerca](research/RESEARCH_STATE.md), snapshot storico.
4. [Formalismo covariante](research/covariant_formalism.md),
   [benchmark statico](research/static_benchmark.md),
   [barriera esatta](research/exact_barrier_benchmark.md).

Il presente file contiene il minimo necessario anche se viene condiviso
da solo; le note collegate contengono le tabelle complete e le derivazioni.

## 3. Che cosa fanno i programmi

| Percorso | Contenuto |
|---|---|
| core/schwarzschild_wkb.py | frequenze Schwarzschild WKB1/WKB3, spin 0/1/2 |
| core/leaver_qnm.py | frequenze di riferimento con frazioni continue |
| core/madelung_profile.py | profilo scalare ampiezza–fase |
| core/dirac_madelung_profile.py | partner Dirac, gerarchia di spin/Madelung |
| calculations/static_madelung_benchmark.py | confronto su 69 modi Schwarzschild |
| calculations/poschl_teller_madelung_benchmark.py | spettro e profili esatti Pöschl–Teller |
| calculations/order_resolved_madelung_benchmark.py | residui locali e loro fallimento vicino ai turning point |
| calculations/uniform_madelung_defect.py | riferimento quadratico locale alla frequenza esatta |
| calculations/audit_poschl_teller.py | nuovi controlli analitico-numerici e confronto Jost/Newton |
| calculations/verify_spectral_audit.wl | 18 controlli simbolici dell'audit |
| calculations/verify_formalism.wl | 19 controlli simbolici precedenti |
| calculations/kerr_scalar_madelung_symbolic.py | riduzione scalare Kerr e lenta rotazione |
| calculations/vaidya_madelung_symbolic.py | PDE e split Madelung su Vaidya |

Non esiste ancora, in questi risultati, una validazione numerica del nuovo
diagnostico in Kerr o una simulazione Vaidya che dimostri memoria aggiuntiva.

## 4. Convenzioni e distinzioni da preservare

Unità geometriche, tempo e^{-iωt}, QNM smorzato con Im ω<0.
Per ε²ψ''+qψ=0 e ψ=A exp(iS/ε), A,S reali, P=S':

\[
P^2+Q_M=q_R,\qquad Q_M=-\varepsilon^2A''/A,\qquad
\varepsilon(A^2P)'=-q_I A^2.
\]

La corrente radiale non è costante quando q_I≠0. Questo non contraddice
la corrente covariante conservata della PDE completa.

Il «fluido» è una rappresentazione del campo: non un fluido materiale,
né una densità di probabilità QNM globalmente normalizzabile.
L'identità locale non introduce una nuova interazione.
Il momento complesso ausiliario della WKB/Schwarziana non coincide
automaticamente con P reale del Madelung QNM.

L'idea discussa con l'utente di imporre -u·W=E misurata costante non è
una legge di conservazione: seleziona osservatori. In Vaidya non si deve
cancellare artificialmente il drift di energia per imporre tale vincolo.
Il frame dragging di Kerr è geometrico, non un termine da rinominare Madelung.

## 5. Risultati precedenti e correzioni dell'audit

### 5.1 Schwarzschild

69 modi: s=0,1 con ℓ=1…8; s=2 con ℓ=2…8; n=0,1,2.
Frequenze Leaver e profili indipendenti dalla WKB mostrano forti
correlazioni fra E_M ed errore WKB a overtone fissato.
Esempio n=0,WKB3: Pearson 0.9998, correlazione parziale 0.9855;
RMSE LOOCV 0.0234 dex con E_M contro 0.0608 con ε².

Limiti: E_M usa il modo di riferimento; rimozione lineare di log ε²;
LOOCV punto per punto, non validazione su una famiglia esclusa.
Questi dati non dimostrano ancora un predittore utilizzabile senza Leaver.

### 5.2 Pöschl–Teller: risultati chiusi

\[
V=L^2\operatorname{sech}^2y,\quad \varepsilon=L^{-1},\quad
N=n+\tfrac12,\quad B=N^2+\tfrac14,\quad
\omega_{\rm ex}=\sqrt{L^2-\tfrac14}-iN,\quad L>\tfrac12.
\]

Per n=0:

\[
A_0\propto\sqrt{\cosh y},\qquad
Q_M=-\frac{\varepsilon^2}{4}(1+\operatorname{sech}^2y).
\]

Per la finestra fissa del benchmark:

\[
c_w=\frac{\int w\tanh^2y\,dy}{\int w\,dy},\qquad
\mathcal E_M=\frac{\varepsilon^2(2-c_w)}{8-\varepsilon^2c_w}.
\]

Pertanto E_M non aggiunge informazione indipendente da ε in questa famiglia
fondamentale. Discrepanza formula/codice ≤2.641×10⁻¹³ sulla griglia testata.

A n fissato e L→∞, per gli errori relativi:

\[
E_1\sim\frac{B}{2L^2},\qquad E_3\sim\frac{N}{128L^5}.
\]

Queste sono espansioni ricavate nell'audit da spettri noti, non rivendicazioni
di nuova fisica. Non implicano una pendenza L⁻⁵ universale per un difetto
della sola ampiezza.

### 5.3 Problema del riferimento uniforme attuale

Lo script usa q_PC=q0+y², con q0=(ω_ex/L)²−1, e dati di Cauchy esatti
al massimo. L'indice ν=iq0/(2ε)−1/2 soddisfa

\[
\nu-n=N(\sqrt{1-\varepsilon^2/4}-1)-iB\varepsilon/2\ne0.
\]

La quantizzazione QNM della parabola richiederebbe ν=n.
È un valido riferimento locale di profilo, non un QNM globale della parabola.
Il difetto usa inoltre la curvatura esatta nel numeratore.
Il nodo n=1 rende A=|ψ| nullo al centro e Q_M lì non definito:
non nascondere il problema sostituendo A con una piccola costante.

### 5.4 Confronto spettrale standard aggiunto

Alla frequenza di prova, con λ=√(L²−1/4):

\[
a=\tfrac12+i\lambda-i\omega,\quad b=\tfrac12-i\lambda-i\omega,\qquad
C_{\rm in}=\frac{\Gamma(1-i\omega)\Gamma(-i\omega)}{\Gamma(a)\Gamma(b)}.
\]

Si calcola δω=−C_in/C'_in, senza usare ω_ex nel passo; ω_ex serve
solo alla validazione. È Newton/Jost standard, non una nuova invenzione.

- n=0,L=4,WKB3: errore 3.844792×10⁻⁶ → 1.287318×10⁻¹⁰ in un passo.
- n=2,L=4,WKB1: errore 0.1423877 → 0.1524684, quindi peggiora.
- n=1 è gestibile dal controllo spettrale anche se Q_M è indefinito al nodo.

Il passo è locale, dipende dalla normalizzazione fuori dal limite vicino
alla radice e non è un bound globale.

## 6. Ponte formale da confrontare con i precedenti

Per una soluzione approssimata non nulla, R=(ε²∂²+q)ψ̃:

\[
\frac{R}{\widetilde\psi}
=q_R-\widetilde P^2-\widetilde Q_M
+i\left[q_I+\varepsilon\left(\widetilde P'
+2\widetilde P\frac{\widetilde A'}{\widetilde A}\right)\right].
\]

Se χ è soluzione omogenea esatta alla stessa frequenza:

\[
\mathcal W(\chi,\widetilde\psi)\big|_{y_a}^{y_b}
=\varepsilon^{-2}\int_{y_a}^{y_b}\chi R\,dy.
\]

Questa identità bilineare, senza coniugazione, vale su un intervallo finito.
L'estensione ai bordi richiede trattamento delle divergenze QNM.
Se χ è approssimata compare anche il suo residuo.
Non è di per sé uno stimatore economico, né una formula nuova.

## 7. Aggiornamento bibliografico decisivo

La ricerca dell'8 settembre ha ristretto ulteriormente il margine di novità:

- [Leung et al. (1997/1998)](https://arxiv.org/abs/physics/9712037),
  eq. (2.14)–(2.16): perturbazione logaritmica QNM e norma radiativa.
- [Glampedakis–Andersson (2003)](https://arxiv.org/abs/gr-qc/0304030):
  ampiezza–fase/Prüfer per risonanze di buchi neri.
- [Konoplya et al. (2019)](https://arxiv.org/abs/1904.10333), eq. (18):
  Δ_k=|ω_{k+1}−ω_{k−1}|/2 e confronti Padé.
- [Gopalakrishnan et al. (2024), §3](https://arxiv.org/html/2403.19485v2):
  stimatori a residui pesati con soluzione aggiunta per modi ottici aperti.
- [Jaramillo et al. (2021)](https://arxiv.org/abs/2004.06434):
  sensibilità spettrale QNM.
- [Kumar (2026)](https://arxiv.org/abs/2602.00507):
  Bohm–Madelung/Ermakov/Weber; [Meza-Domínguez–Matos (2026)](https://arxiv.org/abs/2605.28887):
  precedente idrodinamico diretto nel ramo Dirac, da verificare tecnicamente.
- [Yoo et al., v3 (2026)](https://arxiv.org/html/2510.25062v3):
  propagazione/scattering in Vaidya, controllo necessario per qualsiasi
  rivendicazione di memoria.

La prescrizione esatta del nostro indicatore non è stata identificata nelle
fonti consultate, ma questo non prova priorità. Anche combinare strumenti
noti può essere utile: deve emergere un vantaggio specifico e verificato.

## 8. Giudizio editoriale e richiesta di revisione a Claude

Il materiale è oggi una buona base di lavoro, non un articolo pronto su
un nuovo stimatore. Una breve nota critica richiederebbe un risultato utile
oltre la correzione della nostra ipotesi iniziale. Un articolo metodologico
statico potrebbe diventare proponibile se passa i confronti sotto.
CQG e PRD sono sedi pertinenti, non garanzie di accettazione.

Richiesta proposta per Claude:

1. Controlla indipendentemente le formule di §5–6; cerca errori di segno,
   normalizzazione, ipotesi di realtà e condizioni radiative.
2. Confronta il ponte del §6 con LPT e residui aggiunti: indica con precisione
   che cosa, se qualcosa, resterebbe non equivalente.
3. Proponi il test minimo che possa smentire o sostenere l'utilità del
   diagnostico senza usare soluzione/frequenza esatta negli input.
4. Includi Δ_k e controlli Padé: per Δ3 occorrono WKB2 e WKB4,
   non soltanto i due ordini attualmente implementati.
5. Richiedi validazione a ε e n fissati, su deformazioni indipendenti e
   intere famiglie escluse dalla calibrazione; misura anche il costo.
6. Formula una rivendicazione circoscritta e falsificabile, oppure
   raccomanda di abbandonare lo stimatore se non aggiunge informazione.

Non è richiesto confermare le conclusioni di Codex. Questa consegna serve
a una critica indipendente, non a ottenere un secondo consenso automatico.
La richiesta non autorizza push, invii a riviste o comunicazioni esterne.

## 9. Verifiche già eseguite e riproduzione

Audit precedente dell'8 settembre: 18 nuovi +19 precedenti controlli
Wolfram superati; 32 test Python preesistenti superati (20 core, 12
calculations); verifiche numeriche mpmath 1.3.0 a 80 cifre.
Questi test non sono stati rilanciati durante la sola ricerca bibliografica.

Usare Python 3.13, perché il python3 predefinito sul Mac può non avere NumPy.
Mathematica 13.3.1 ARM64 funziona tramite kernel diretto. Il launcher
wolframscript era bloccato nell'ambiente automatizzato, pur funzionando
nel terminale dell'utente: non dedurre un problema di licenza.

~~~bash
python3.13 calculations/audit_poschl_teller.py --uniform-tail
/Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noprompt -script calculations/verify_spectral_audit.wl
/Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noprompt -script calculations/verify_formalism.wl

python3.13 -m unittest \
  calculations/test_static_madelung_benchmark.py \
  calculations/test_poschl_teller_madelung.py \
  calculations/test_order_resolved_madelung.py \
  calculations/test_uniform_madelung_defect.py

# In un terminale separato, dalla cartella core:
python3.13 -m unittest test_wkb.py test_leaver.py test_dirac_madelung.py
~~~

Le cartelle output e tmp sono ignorate da Git. Salvare ogni nuovo risultato
scientifico rilevante in Markdown, specificando se sia derivato, misurato,
ipotizzato o tratto dalla letteratura. Verificare i precedenti sul testo
originale prima di usare la parola «nuovo».
