# Fondamenti covarianti per QNM, WKB e Madelung

**Nota di lavoro — 7 settembre 2026.** Questo documento estende il ramo
Schwarzschild già verificato nel repository. Le formule etichettate come
*identità* seguono direttamente dall'equazione d'onda; gli oggetti etichettati
come *candidati* sono ipotesi da sottoporre a test numerico e non risultati.

## 1. Obiettivo e confine della rivendicazione

La trasformazione di Madelung, la WKB per i QNM e i metodi ampiezza–fase non
sono nuovi separatamente. Il programma di ricerca è più specifico:

1. formulare la decomposizione senza confondere la corrente conservata della
   PDE con il problema spettrale aperto ottenuto dopo separazione;
2. verificare se il funzionale di Madelung fornisca un diagnostico quantitativo
   del fallimento WKB, non soltanto una riscrittura dell'equazione;
3. seguire lo stesso oggetto passando da spazio-tempi statici a stazionari
   non statici e infine non stazionari;
4. nel caso Vaidya, isolare una misura normalizzata della risposta ritardata del
   profilo QNM rispetto alla soluzione a fondo congelato.

La frase «fluido QNM» indica quindi una rappresentazione ampiezza–fase. Non
implica che la perturbazione gravitazionale sia un fluido materiale o che
\(A^2\) sia una densità di probabilità normalizzabile.

## 2. Identità covarianti per un campo scalare

Si prenda come modello iniziale un campo complesso scalare su un fondo reale:

\[
\left(\varepsilon^2\Box_g-\mathcal U\right)\Phi=0,
\qquad \mathcal U\in\mathbb R,
\]

e, in un dominio privo di zeri, si ponga

\[
\Phi=A\exp(iS/\varepsilon),\qquad A>0,\quad S\in\mathbb R.
\]

Separando parte reale e immaginaria si ottengono le identità

\[
g^{\mu\nu}\nabla_\mu S\nabla_\nu S+\mathcal U+Q_g=0,
\qquad
Q_g=-\varepsilon^2\frac{\Box_g A}{A},
\tag{2.1}
\]

\[
\nabla_\mu J^\mu=0,
\qquad
J^\mu=A^2\nabla^\mu S.
\tag{2.2}
\]

La (2.1) è una Hamilton–Jacobi esatta con un funzionale di curvatura
dell'ampiezza; la (2.2) è la corrente locale della PDE completa. Entrambe
dipendono dalla variabile di campo scelta: per perturbazioni gravitazionali
una variabile master non è automaticamente uno scalare geometrico.

### 2.1 PDE conservativa contro ODE QNM aperta

Dopo separazione temporale, una frequenza QNM complessa rende complesso il
coefficiente della ODE radiale. Per

\[
\varepsilon^2\psi''+q\psi=0,
\qquad q=q_R+iq_I,
\qquad \psi=Ae^{iS/\varepsilon},
\]

con \(A,S\) reali e \(P=S'\), si ha

\[
P^2+Q_{1\mathrm D}=q_R,
\qquad Q_{1\mathrm D}=-\varepsilon^2\frac{A''}{A},
\tag{2.3}
\]

\[
\varepsilon(\rho P)'=-q_I\rho,
\qquad \rho=A^2.
\tag{2.4}
\]

Non c'è contraddizione fra (2.2) e (2.4): la derivata temporale della densità è
stata eliminata nella riduzione spettrale e ricompare come sorgente radiale.
L'apertura fisica è fissata dalle condizioni entrante all'orizzonte e uscente a
infinito; non richiede un potenziale fondamentale non hermitiano.

### 2.2 Dizionario con i formalismi equivalenti

Se si impone il trasporto \(A=P^{-1/2}\), con \(P=S'\), allora

\[
\frac{A''}{A}=rac34\left(\frac{P'}P\right)^2
-\frac12\frac{P''}P=-\frac12\{S,x\},
\]

dove \(\{S,x\}\) è la derivata Schwarziana. La chiusura statica può quindi
essere scritta anche come

\[
q=P^2+\frac{\varepsilon^2}{2}\{S,x\}.
\tag{2.5}
\]

Lo stesso contenuto compare, con variabili diverse, nella Riccati logaritmica,
nella quantum Hamilton–Jacobi, nei metodi Prüfer/phase-amplitude e nell'exact
WKB. Il termine «Madelung» seleziona l'interpretazione in ampiezza e corrente;
non rende nuova l'identità (2.5). Il progetto deve produrre un diagnostico o un
risultato dinamico che non segua dalla sola rinominazione dello Schwarziano.

## 3. La scala delle simmetrie

| Classe | Simmetria temporale | Nuovo ingrediente | Oggetto da misurare |
|---|---|---|---|
| statica | Killing ipersuperficie-ortogonale | barriera e turning point | \(Q_{1\mathrm D}\), errore WKB |
| stazionaria non statica | Killing con shift | frame dragging e superradianza | frequenza co-rotante e flussi |
| non stazionaria | nessuna Killing temporale | ritardo, mixing e propagazione | eccesso dinamico rispetto al fondo congelato |

Questa progressione deve essere implementata sullo stesso campo scalare prima
di cambiare spin: altrimenti gli effetti della geometria si confondono con i
termini subprincipali di spin già trovati nel ramo Dirac.

## 4. Livello I: fondo statico

Il caso Schwarzschild è già sviluppato in
`Schw-QNM-WKB-Fluid/manuscript.md`. I prossimi controlli statici sono:

- applicare la stessa convenzione a Reissner–Nordström o
  Schwarzschild–de Sitter;
- sostituire la serie locale singolare con una forma uniforme nella regione dei
  turning point coalescenti;
- cercare un funzionale integrato nella barriera che ordini correttamente
  l'errore della frequenza WKB rispetto a Leaver.

Un candidato prudente, da calibrare e non ancora canonico, è

\[
\mathcal E_M=
\frac{\int_{\mathcal B}w(r_*)|Q_{1\mathrm D}|\,dr_*}
{\int_{\mathcal B}w(r_*)
\left(|\omega|^2+|V|+|P|^2\right)\,dr_*},
\tag{4.1}
\]

dove \(w\) è una finestra regolare centrata sulla barriera. La scelta di
\(w\), della variabile master e della coordinata deve essere dichiarata: il
primo test è stabilire se la correlazione con l'errore WKB sopravviva a scelte
ragionevoli.

## 5. Livello II: fondo stazionario non statico

In decomposizione ADM il termine principale di Hamilton–Jacobi è

\[
g^{\mu\nu}S_\mu S_\nu=
-\frac{(\partial_tS-\beta^i\partial_iS)^2}{\alpha^2}
+\gamma^{ij}\partial_iS\partial_jS.
\tag{5.1}
\]

Per Kerr, \(S=-\omega_Rt+m\phi+S_r+S_\theta\), e lo shift azimutale
produce localmente la combinazione co-rotante \(\omega_R-m\Omega\). Il piano
minimo è:

1. campo scalare, per evitare inizialmente la connessione di spin;
2. slow Kerr come ponte perturbativo dal caso statico;
3. Kerr completo separato, distinguendo \(Q_r\) e \(Q_\theta\);
4. confronto dei flussi ai due lati della soglia
   \(\omega_R=m\Omega_H\).

Il termine di frame dragging appartiene al simbolo geometrico principale, non
va rinominato «potenziale di Madelung». La domanda è come \(Q_g\) corregga la
propagazione su quel simbolo e se presenti una firma robusta della regione
superradiante.

## 6. Livello III: Vaidya

Per la metrica entrante

\[
ds^2=-f(v,r)dv^2+2\,dv\,dr+r^2d\Omega^2,
\qquad f(v,r)=1-\frac{2M(v)}r,
\]

e \(\Phi=\psi(v,r)Y_{\ell m}/r\), la Klein–Gordon massless diventa

\[
2\psi_{vr}+\partial_r(f\psi_r)-U_\ell\psi=0,
\qquad
U_\ell=\frac{f_r}{r}+\frac{\ell(\ell+1)}{r^2}.
\tag{6.1}
\]

Nel conteggio eikonale si intende
\(L=\ell+\tfrac12\), \(\varepsilon=L^{-1}\): in questo modo il contributo
angolare in \(\varepsilon^2U_\ell\) rimane di ordine principale. Senza questa
scalatura, \(\varepsilon\) nella derivazione seguente è soltanto un parametro
formale di fase.

Ponendo \(\psi=Ae^{iS/\varepsilon}\), la parte reale è

\[
2S_vS_r+fS_r^2+\varepsilon^2U_\ell+Q_V=0,
\tag{6.2}
\]

\[
Q_V=-\varepsilon^2
\frac{2A_{vr}+\partial_r(fA_r)}{A}
=Q_r+Q_\times,
\tag{6.3}
\]

\[
Q_r=-\varepsilon^2\frac{\partial_r(fA_r)}A,
\qquad
Q_\times=-2\varepsilon^2\frac{A_{vr}}A.
\tag{6.4}
\]

Qui \(A\) è l'ampiezza della variabile ridotta \(\psi=r\Phi/Y_{\ell m}\):
\(Q_V\) non coincide automaticamente con lo scalare covariante \(Q_g\) della
(2.1), perché la riduzione armonica e il fattore \(r\) cambiano la variabile
master. Ogni confronto fra coordinate o fra spin diversi deve fissare questa
convenzione, oppure tornare all'ampiezza del campo covariante \(\Phi\).

La parte immaginaria è la legge conservativa caratteristica

\[
\partial_v(\rho S_r)
+\partial_r\!\left[\rho(S_v+fS_r)\right]=0,
\qquad \rho=A^2.
\tag{6.5}
\]

Le (6.1)–(6.5) sono verificate simbolicamente da
`calculations/vaidya_madelung_symbolic.py`.

### 6.1 Dal termine misto a una memoria misurabile

\(Q_\times\) è un termine esatto, ma da solo non è ancora una misura pulita di
memoria: contiene anche l'inviluppo di decadimento del QNM e dipende dalla
normalizzazione temporale scelta durante l'estrazione del modo.

Si introduca perciò un profilo congelato

\[
A_{\rm fr}(v,r)=A_{\rm QNM}\bigl(r;M(v)\bigr),
\]

costruito con le stesse coordinate e la stessa variabile master. Un primo
*candidato computazionale*, invariante sotto normalizzazioni puramente
temporali di ciascun profilo, è

\[
\boxed{\;
\mathcal M_{vr}=-2\varepsilon^2\,
\partial_v\partial_r\ln\!\left(\frac{A}{A_{\rm fr}}\right)
\;}. \tag{6.6}
\]

La (6.6) non è identica a \(Q_\times\): è la curvatura logaritmica usata per
isolarne l'eccesso non adiabatico. È ancora legata alla gauge di
Eddington–Finkelstein. In simmetria sferica si può testare la versione scalare

\[
\mathcal M_K=-2\varepsilon^2
K^a(\nabla^b r)\nabla_a\nabla_b
\ln\!\left(\frac{A}{A_{\rm fr}}\right),
\tag{6.7}
\]

dove \(K^a\) è il vettore di Kodama sullo spazio orbitale bidimensionale. La
normalizzazione e il rapporto fra (6.6) e (6.7) sono problemi aperti, non
risolti in questa nota. Nella (6.7), per avere davvero uno scalare, \(A\) deve
essere l'ampiezza del campo covariante o di una variabile master geometricamente
definita, non una riscalatura radiale lasciata implicita.

## 7. Ipotesi falsificabili

1. **Limite statico.** Per \(M(v)=M_0\), dopo la rimozione dell'inviluppo
   globale, \(\mathcal M_{vr}=0\).
2. **Accrescimento costante.** Le frequenze devono riprodurre il problema in
   frequenza di Capuano–Santoni–Barausse nel limite corrispondente.
3. **Profilo generico.** Per una transizione regolare di massa, il massimo o un
   integrale di \(\mathcal M\) nella barriera deve correlare con
   \(\Delta\omega=\omega_{\rm dyn}-\omega_{\rm fr}[M(v)]\).
4. **Causalità.** Il ritardo deve dipendere dalla posizione dell'osservatore e
   non può essere descritto soltanto dalla photon sphere istantanea.
5. **Esito negativo informativo.** Se \(\mathcal M\) non predice il ritardo
   meglio dei parametri adiabatici standard \(|\dot M|\) e \(|\ddot M|\), la
   congettura di memoria Madelung è respinta.

## 8. Ordine operativo

1. chiudere la ricerca di anteriorità sui formalismi equivalenti;
2. implementare un benchmark statico con soluzione indipendente;
3. costruire la forma uniforme vicino al massimo della barriera;
4. implementare slow Kerr scalare;
5. implementare Vaidya a \(\dot M\) costante;
6. passare a un profilo \(M(v)\) generico e misurare (6.6)–(6.7).
