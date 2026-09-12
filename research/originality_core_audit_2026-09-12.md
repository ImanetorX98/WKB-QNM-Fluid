# Originalità del core e programma verificabile — 12 settembre 2026

## Verdetto operativo

Non abbiamo ancora dimostrato un nucleo originale sufficiente per un articolo
di ricerca. Abbiamo risultati corretti e riproducibili, ma correttezza,
originalità e utilità editoriale sono tre condizioni diverse. La tesi forte
del manoscritto attuale non regge: il termine Kerr lineare era contaminazione
discreta; la proporzionalità esatta del diagnostico con Lambda3 non segue dai
fit; l'invarianza sotto cambio di variabile master non segue dalla geometria
di Fisher. Non basta riunire Schwarzschild, Kerr e Vaidya sotto il nome Madelung.

La strada più concreta è una **nota metodologica su ciò che un diagnostico
d'ampiezza misura realmente**, con controlli di parametrizzazione, frequenza,
nodi e bordo. La strada fisicamente più ambiziosa resta Vaidya, ma richiede
un confronto quantitativo con una risposta temporale e con precedenti vicini.
Queste sono valutazioni motivate, non garanzie di accettazione editoriale.

## Ricerca bibliografica e limiti

Ricerca web mirata per le combinazioni `quasinormal Madelung`, `Madelung WKB`,
`Dolan Ottewill inverse multipole`, `Vaidya adiabatic quasinormal`,
`logarithmic perturbation quasinormal` e `phase integral Schwarzian`.
Usate come evidenza fonti primarie arXiv e pagine degli editori. Non è una
revisione sistematica esaustiva né una prova di assenza di precedenti.
Per Yang, Capuano e Yoo consultati anche passaggi del testo HTML; per le altre
fonti sotto il confronto si limita agli abstract/metadati, salvo indicazione.

| Affermazione candidata | Precedente e sovrapposizione | Valutazione |
|---|---|---|
| Espansione QNM in 1/L | [Dolan–Ottewill 2009](https://arxiv.org/abs/0908.0329), frequenze e funzioni d'onda | Non nuova come metodo |
| Ordini subleading e spin Kerr | [Dolan 2010](https://arxiv.org/abs/1007.5097), correzioni oltre l'eikonale | Non rivendicare una nuova gerarchia generale |
| Quantizzazione angolare a L=ell+1/2 | [Yang et al. 2012](https://arxiv.org/html/1207.4253), §II e App. A, eq.116–118 | Precedente diretto; A1=0 è coerente con la struttura standard |
| Formula algoritmica oltre eikonale | [Konoplya–Zhidenko 2023](https://arxiv.org/abs/2309.02560) | Serve una differenza concreta, non una notazione diversa |
| Riccati e perturbazione dei QNM | [Leung et al. 1998](https://arxiv.org/abs/physics/9712037) | Non rivendicare derivata logaritmica o regolarizzazione in sé |
| Madelung e geometria di Fisher | [Khesin et al. 2018](https://arxiv.org/abs/1711.00321) | Precedente matematico; non certifica l'invarianza del nostro Q_M |
| QNM in fluidi BEC | [Barceló et al. 2007](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.75.084024) | Contesto distinto: dispersione Bogoliubov e buchi neri acustici, non identico al nostro problema |
| Risposta non adiabatica Vaidya | [Abdalla et al. 2006](https://arxiv.org/abs/gr-qc/0609036), [Lin et al. 2021](https://arxiv.org/abs/2104.06631) | Il fenomeno generale non è nuovo |
| QNM Vaidya a tasso costante | [Capuano et al. 2024](https://arxiv.org/html/2407.06009v2) | Benchmark obbligatorio per una correzione in dot M |
| Ringdown dinamico e frequenza osservata | [Yoo et al., v3 febbraio 2026](https://arxiv.org/html/2510.25062v3), §IV–V | Già confronto adiabatico/time-domain, propagazione e redshift |

Yang usa la stessa quantizzazione Bohr–Sommerfeld e il rapporto m/L.
La sua espansione in piccolo a omega/L non è una prova globale della nostra
parità a chat arbitrario: non confondere espansioni differenti. Resta però
troppo vicino per presentare la nostra derivazione formale come nuova senza
un confronto integrale di formule e ipotesi. Non ho trovato nelle ricerche
mirate la precisa diagnosi A1 spurio=delta_m partial_mu A0; ciò non dimostra
che sia inedita. È comunque una Taylor expansion elementare: valore possibile
nella dimostrazione numerica e nella prevenzione di errori, non come nuova fisica.

Yoo distingue già previsione locale e frequenza ricavata dalla forma d'onda
a raggio finito, inclusi ritardo e redshift. Perciò non basta chiamare
"memoria" una correzione locale della famiglia congelata. Capuano tratta
il caso speciale a tasso costante con cambio di coordinate e dominio in
frequenza: occorre allineare coordinata temporale, raggio e normalizzazione
prima di confrontare coefficienti.

## Nuovi controlli logici sul core

### 1. Identità di variabili non significa equivalenza con WKB troncata

La trasformazione Riccati/Madelung esatta conserva l'informazione dell'ODE.
Non segue che un funzionale della soluzione esatta contenga soltanto
l'informazione del getto di ordine sei al picco usato da WKB3. La dispersione
0.10% del rapporto E_M/|Lambda3| in un campione ristretto è una correlazione
approssimata, non una proporzionalità esatta o un no-go universale.

Controllo analitico: con z=psi'/psi, z'+z^2+q=0. Variando q di eta b(x)
a frequenza e dato iniziale fissati, h=partial_eta z soddisfa

    h'+2 z h=-b,
    h(x)=-psi(x)^(-2) integral(x0,x) b(t) psi(t)^2 dt.

Una deformazione liscia supportata fuori da un intorno del picco lascia
invariato tutto il getto locale e dunque Lambda3, ma può cambiare z e Q_M
nella regione raggiunta dall'integrazione. Questo argomento riguarda il
problema a dati iniziali, NON costituisce da solo un controesempio on-shell
al rapporto osservato nella famiglia Schwarzschild. Per quello serve
ricalcolare la frequenza QNM con entrambi i bordi sulla famiglia deformata.

Esperimento discriminante: V_eta=V_RW+eta b, b a supporto compatto esterno
al picco, eta piccolo per non cambiare il massimo dominante. Confrontare
E_M, Lambda3, errore WKB e frequenza indipendente, includendo perturbazioni
mai usate per calibrare il rapporto. Un esito negativo è informativo;
non predisporre la deformazione per garantire una conclusione desiderata.

### 2. Il diagnostico locale non è invariante sotto riscalatura master

Per Atilde=e^g A, con g reale regolare e Q_M=-eps^2 A''/A,

    Qtilde_M-Q_M=-eps^2 [g''+2g' A'/A+(g')^2].

Anche A costante e g localmente lineare danno un cambiamento non nullo.
Si può scegliere g a supporto compatto, così la trasformazione è invertibile
e non modifica le condizioni asintotiche, se si trasforma anche l'operatore.
Il nuovo operatore può contenere una derivata prima: non sono due diverse
soluzioni dello stesso operatore di Schrödinger. È un controesempio alla
pretesa di indipendenza dalla variabile master arbitraria, non alla covarianza
di un funzionale geometrico con metrica, misura e densità trasformate insieme.
Inoltre l'ampiezza u^(-1/2) complessa non è automaticamente una densità di
probabilità reale positiva, e i profili QNM non sono normalmente integrabili.

### 3. Condizionamento: conservare entrambe le direzioni

Sul ramo asintotico Q_M=-(Im omega)^2/L^2: il fattore 2 è quello della
variazione proporzionale delta omega=delta omega_base; per perturbazioni
complesse arbitrarie in norma relativa il fattore è 2|omega|/|Im omega|.
Il manoscritto conserva ancora la ritrattazione eccessiva del secondo numero.
Il resoconto precedente sul condizionamento va recepito, non cancellato.

## Calcoli avviati dopo l'audit

`kerr_madelung_analytic.py`: mantenuta la correzione geometrica precedente;
rimossa un'integrazione del profilo storico usata soltanto per recuperare
il raggio del picco; sostituita la quadratura uniforme singolare di H/Delta
con la differenza analitica della coordinata tortoise non estremale.
L'estremo finale è ora realmente r=60 entro errore dell'integrazione.

Aggiunto bordo entrante di primo ordine mediante z=psi_*/psi:

    v z_r=-q-z^2,
    z0=-i(omega-m Omega_H),
    z1=-q_r(r_plus)/(v_r(r_plus)+2z0).

q_r è stimato con differenza centrale, passo 1e-5; non lo chiamiamo bordo
esatto. Lo sweep confronta ordini 0/1, offset 1e-4/1e-6, griglie 12001/24001,
frequenze eikonali storiche/WKB3 autoconsistenti, modi (28,19),(61,41) a mu=2/3.
La frequenza eikonale storica è soltanto un baseline: il suo riferimento
angolare finito e la derivata in omega richiedono ancora un audit dedicato.
Un piccolo residuo WKB3 verifica quel sistema approssimato, non una frequenza
QNM esatta alla precisione del residuo.

I risultati del nuovo sweep sono in `kerr_corrected_profile_results_2026-09-12.md`.
Questa è la rigenerazione del riferimento d'ampiezza, non ancora dell'intera
tabella degli errori del §9.6, che richiede anche una previsione coerente.

## Ordine dei prossimi passi e soglia per il paper

1. Chiudere convergenza del riferimento Kerr (bordi, griglia, tolleranza),
   poi usare frequenze indipendenti dalla condizione WKB3 per l'errore fisico.
2. Test fuori campione su deformazioni di barriera, separando diagnostico
   integrale e mediana pesata: non sono lo stesso funzionale.
3. Riscrivere il core come risultato limitato sostenuto dai test, evitando
   "esatto", "universale", "invariante" dove mancano prove.
4. Vaidya: confrontare prima il caso a dot M costante con Capuano; poi una
   storia M(v) liscia con evoluzione PDE causale e segnale a raggio fissato.
   Verificare che la correzione migliori il residuo/errore all'ordine previsto
   e resti invariata sotto G -> c(M)G quando combinata con l'ampiezza temporale.
5. Solo se emerge una previsione quantitativa aggiuntiva o un limite generale
   non già coperto dai precedenti, promuoverlo a rivendicazione di originalità.

Non è stata completata una nuova simulazione Vaidya in questo turno.
Il manoscritto integrato resta una bozza da revisionare, non pronto all'invio.
