# Risposta al handoff Kerr A1 del 12 settembre

## Procedura adottata e stato

Letto integralmente CODEX_HANDOFF_KERR_A1.md; worktree pulito all'avvio sul
branch analisi-kerr-vaidya-2026-09. Il lavoro precedente era già salvato in Git.

1. Secondo metodo angolare, indipendente dalla base armonica: implementato e
   verificato per casi reali e complessi, con controllo di precisione e grado.
2. Derivazione di A1=0: completata **nel senso dell'espansione formale WKB**
   scalare con due turning point semplici, nelle ipotesi specificate sotto.
   Non è un teorema uniforme per ogni settore bosonico o ogni regime Kerr.
3. Audit delle quantità discrete: classificate sei occorrenze del rounding;
   individuato e corretto anche un fattore mancante nel termine geometrico
   del modulo analitico Kerr. Il rifacimento quantitativo del §9.6 resta aperto.
4. Manoscritto: non riscritto automaticamente. La proposta del handoff di una
   regola universale basata sull'ordine differenziale non è giustificata.

## 1. Secondo solutore: equazione differenziale su nodi Chebyshev

La convenzione del progetto è

~~~text
(1-x^2) S'' - 2x S' + [A+c^2 x^2-m^2/(1-x^2)] S = 0.
~~~

La sostituzione S=(1-x^2)^(|m|/2)y elimina la potenza singolare agli estremi:

~~~text
(1-x^2)y'' - 2(|m|+1)x y'
 + [A-|m|(|m|+1)+c^2 x^2]y = 0.
~~~

Il nuovo angular_chebyshev_independent.py discretizza direttamente questa
equazione con matrici di derivazione Chebyshev-Lobatto. Le righe agli estremi
impongono la relazione regolare dell'ODE; non imponiamo y=0 artificialmente.
Per c reale il modo è individuato dall'indice ell-|m| nell'ordinamento degli
autovalori; per c complesso si continua il modo tramite autovettori.
Il metodo non costruisce elementi di matrice in armoniche sferiche.

### Un limite trovato, non nascosto

La matrice di collocazione può essere mal condizionata in doppia precisione.
Per ell=61,m=41,chat=0.6, il disaccordo con la base armonica è 4.3e-6 a grado
60 e peggiora a 3.0e-5 a grado 120. Aumentare soltanto il grado non basta.

La versione multiprecisione ricostruisce **tutte** le matrici a 55 cifre e
risolve l'autocoppia con iterazione inversa. La sola inizializzazione del ramo
proviene dalla collocazione double, mai dal solutore armonico. La verifica
non si ferma al residuo della matrice: varia anche il grado.

| ell,m | chat | grado MP | differenza relativa dalla base armonica |
|---|---|---:|---:|
| 40,27 | 0.6 | 40 | 6.28e-10 |
| 40,27 | 0.6 | 52 | 7.8e-16 |
| 61,41 | 0.6 | 52 | 4.87e-8 |
| 61,41 | 0.6 | 64 | 1.1e-16 |
| 40,27 | 0.3+0.2i | 52 e 64 | 1.14e-14 |

Il residuo dell'autocoppia discreta può essere sotto 1e-28 anche quando il
grado è insufficiente: non va scambiato per errore dell'ODE.

### Fit indipendente del coefficiente lineare

Solo collocazione MP, grado 64, 55 cifre; chat=0.6, mu=2/3 esatto;
ell=28,34,40,46,52,61. Fit libero con potenze da 1 a 1/L^4:

~~~text
A0 = 0.8974781207948589
A1 = 8.8335292e-8
A2 = -0.2310722605873
~~~

Gli autovalori MP sono convertiti in double per questo fit; le cifre piccole
di A1 non vanno interpretate come un valore fisico non nullo. Il test verifica
|A1|<1e-6. Non abbiamo ripetuto in MP l'intero sweep fino a ell=480.
Questo fornisce un controllo indipendente della precedente misura, non
una prova numerica di nullità esatta universale.

## 2. Derivazione formale del termine lineare nullo

Partiamo dall'equazione scalare in theta nella convenzione appena dichiarata.
Con S=psi/sqrt(sin theta), L=ell+1/2, epsilon=1/L,
m=mu L, c=chat L, A=L^2 Ahat, si ottiene **esattamente**

~~~text
epsilon^2 psi'' + [Q0(theta,Ahat)+epsilon^2 Q2(theta)] psi = 0,
Q0 = Ahat + chat^2 cos(theta)^2 - mu^2/sin(theta)^2,
Q2 = (1+csc(theta)^2)/4.
~~~

Le sostituzioni sono verificate con Mathematica. Nella convenzione del
codice il segno del termine chat^2 è **positivo**: il segno negativo scritto
nella formula Bohr-Sommerfeld del handoff non è coerente con questo A.

### Ipotesi

Settore scalare, mu reale con 0<|mu|<1 fissato esattamente; chat reale fissato;
un unico intervallo permesso limitato da due turning point semplici, separati
dagli estremi singolari; ramo dell'autovalore seguito senza degenerazioni.
Usiamo la quantizzazione formale WKB a due turning point con indice di Maslov
incluso. Non includiamo coalescenze, cambi di topologia, limiti mu=0 o |mu|=1,
né termini esponenzialmente piccoli oltre tutti gli ordini.

Il numero angolare è n_theta=ell-|m|. La condizione di quantizzazione è

~~~text
I(Ahat) + epsilon^2 I2(Ahat) + ...
 = pi epsilon (n_theta+1/2) = pi(1-|mu|),
I(Ahat) = integral(theta-,theta+) sqrt(Q0) dtheta.
~~~

Non si deduce la parità semplicemente dalla mancanza di epsilon nel membro
destro: occorre la struttura WKB del membro sinistro. Il primo termine di
trasporto determina l'ampiezza e l'indice di Maslov; dopo averlo incluso,
la prima correzione di fase è di ordine epsilon^2. Anche Q2 entra a tale
ordine. Non stiamo affermando che quella correzione quadratica si annulli.

Espandendo Ahat=A0+epsilon A1+... si trova all'ordine epsilon:

~~~text
A1 I_A(A0)=0,
I_A(A0) = 1/2 integral(theta-,theta+) 1/sqrt(Q0) dtheta > 0.
~~~

L'integrale è finito per turning point semplici; i termini dovuti allo
spostamento degli estremi nell'azione si annullano perché sqrt(Q0)=0 lì.
Pertanto **A1=0 nell'espansione formale, sotto queste ipotesi**.
Per chat complesso è possibile proseguire formalmente finché il ramo e il
ciclo restano regolari e I_A non si annulla: non è stata dimostrata qui una
continuazione globale o una stima uniforme del resto complesso.

### Traiettoria fisica: che cosa si può affermare

Se omega_hat=Omega0+epsilon Omega1+... con Omega0 reale e Omega1 immaginario
puro, allora il termine indotto è Omega1 partial_Omega Q0, immaginario puro,
poiché Q0 e la sua derivata sono reali sul ramo reale. Vale quindi
Re(qhat)=Q0+O(epsilon^2) in quelle ipotesi. La derivata include **tutta** la
dipendenza di Q0 da Omega, non soltanto quella dell'autovalore angolare.
Non abbiamo provato che ogni traiettoria fisica Kerr soddisfi quelle ipotesi;
non estendiamo il risultato ai campi di spin diverso o al regime quasi estremale.

## 3. Il criterio 'primo contro secondo ordine' non segue

Il risultato scalare precedente non dimostra una classificazione universale
per ordine differenziale. Nel nostro stesso settore Dirac, eliminando una
componente del sistema del primo ordine si ottiene un'equazione del **secondo**
ordine con V/K^2=h^2+tau epsilon h'. Il termine lineare resta presente.
Quindi l'ordine differenziale della rappresentazione non è il discriminante.

L'identità Dirac è per ogni K: non usa un fit di un rapporto discreto mantenuto
solo nominalmente fisso. L'audit dei file Dirac selezionati non ha trovato
lo stesso arrotondamento di m. Questo non sostituisce una rivalidazione completa
del settore Dirac su Kerr.

## 4. Audit delle sei occorrenze di rounding

| File / funzione | Conseguenza |
|---|---|
| kerr_eikonal_order_test.py / eikonal_fit | Fit dichiarato a mu fisso: contaminato; usare sequenze esatte |
| kerr_radial_order_profile.py / angular_leading | Anche il riferimento ha mu nominale e L finito |
| stesso file / exact_scaled_potential | Il potenziale finito usa mu_eff; confronto incoerente con Q0 nominale |
| kerr_wkb3_selfconsistent.py / selfconsistent_frequency | Per un singolo ell seleziona un modo intero legittimo; non è da solo un bug del root finder. Mu_eff va dichiarato e mantenuto nella sequenza |
| kerr_madelung_profile.py / kerr_profile | Modo intero, ma riferimenti e confronti a mu nominale possono contaminare la diagnosi |
| kerr_madelung_analytic.py / analytic_profile | Stesso problema per confronti/sequenze; inoltre corretto il difetto geometrico sotto |

Non è stato eliminato indiscriminatamente round: un singolo autovalore fisico
richiede m intero. È la pretesa di tenere fisso mu durante la regressione che
richiede sequenze compatibili o un diverso parametro di riferimento.

### Difetto geometrico aggiuntivo corretto

In kerr_madelung_analytic.py il vecchio d2h era partial_r(D_*h), ma la radiale
richiede D_*^2h. Mancava il fattore esterno Delta/H. Ora geometric_term calcola

~~~text
D_*^2 h/h = v [v_r r/H + v a^2/H^2],
v=Delta/H, H=r^2+a^2, h=sqrt(H).
~~~

A a=0 restituisce 2f/r^3; il vecchio codice restituiva 2/r^3.
Due test verificano il limite Schwarzschild e la derivata simbolica per
a=0.3,0.6,0.9. È stata corretta anche la premessa del docstring che chiamava
universale il fattore di condizionamento 2.
Le tabelle precedenti prodotte da questo modulo non sono state rigenerate:
il §9.6 richiede ancora l'audit con frequenze coerenti, dati al bordo accurati
e sequenze a mu esatto. Non va dichiarato validato per effetto di questa patch.

## 5. Comandi, test e limiti editoriali

~~~sh
python3.13 calculations/angular_chebyshev_independent.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_angular_liouville.wl
~~~

Da calculations:

~~~sh
python3.13 -m unittest test_angular_chebyshev_independent.py test_kerr_fixed_mu_radial.py test_kerr_analytic_geometry.py
~~~

Esiti: **8 test** angolari/radiali e **2 test** geometrici superati; **4/4**
controlli Mathematica superati. La prova formale della quantizzazione non
è 'certificata da Mathematica': il kernel controlla le identità algebriche.
Non è stata rieseguita tutta la suite del repository.

Il §9 può essere preparato come risultato scalare condizionato, sostituendo
il vecchio termine geometrico con l'assenza formale di A1 a parametri fissati
e la diagnosi della contaminazione discreta. Non va sostituito con un nuovo
teorema universale su tutti i bosoni o sull'ordine delle equazioni. La revisione
integrata di abstract, §§9–10 e figure resta da fare con questi limiti espliciti.

## Fonti di convenzione

La [DLMF §30.2](https://dlmf.nist.gov/30.2) riporta l'equazione con
lambda+gamma^2(1-x^2). Il confronto con il progetto richiede gamma^2=-c^2
e lambda=A+c^2; confondere A con lambda cambia il segno apparente del termine
sferoidale. Il [quadro eikonale di Yang et al.](https://arxiv.org/abs/1207.4253)
resta un precedente da confrontare prima di attribuire originalità a questa
derivazione. Il presente audit non è una nuova ricerca esaustiva di priorità.
