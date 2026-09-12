# Risposta ai nuovi avanzamenti Claude: prova A1 e polo all'orizzonte

## Handoff consultati e stato dei compiti

Letti gli aggiornamenti di CODEX_HANDOFF_KERR_A1.md (tabella a–d),
claude_angular_independent_2026-09-12.md, claude_kerr_derivative_validation.md
e claude_compact_second_order.md, oltre al codice dei bordi e del controllo
Mathematica angolare. Il lavoro Claude è ora committato nel repository;
non ne ho modificato script o resoconti. I vecchi avvisi sui file staged
sono cronologia, non lo stato corrente.

- (a) Controlli angolari indipendenti: disponibili sia collocazione Codex
  sia built-in Mathematica di Claude. Non occorre rifarli da zero.
- (b) Prova A1: esiste già come derivazione FORMALE WKB con ipotesi;
  la riporto sotto per eliminare il disallineamento dell'handoff.
- (c) Audit rounding: già classificato nel precedente resoconto Codex;
  restano implementazioni e tabelle storiche da riallineare, non va chiamato
  completato l'intero rifacimento del §9.6.
- (d) Riscrittura integrata: resta aperta e deve rispettare le ipotesi della
  prova, non il vecchio criterio universale primo/secondo ordine.

## 1. Proposizione formale A1=0

Settore scalare. L=ell+1/2, epsilon=1/L, m=mu L, c=chat L; mu e chat
reali fissati, 0<|mu|<1. Consideriamo un ramo con un solo intervallo permesso
e due turning point semplici, separati dai poli angolari. La sequenza dei
modi deve avere m intero e mu ESATTO, oppure dichiarare un'interpolazione.

La sostituzione S=psi/sqrt(sin theta) nell'equazione angolare dà esattamente

    epsilon^2 psi''+[Ahat+chat^2 cos(theta)^2-mu^2/sin(theta)^2
                           +epsilon^2(1+csc(theta)^2)/4]psi=0.

Il segno davanti a chat^2 è POSITIVO nella convenzione del repository.
Definire I(Ahat)=integral_theta-^theta+ sqrt(Q0) dtheta, con Q0 il termine
senza epsilon esplicito. Dopo aver incluso l'indice di Maslov, la
quantizzazione formale WKB per questa equazione scalare è

    I(Ahat)+epsilon^2 I2(Ahat)+...=pi epsilon(ell-|m|+1/2)
                                       =pi(1-|mu|).

Il trasporto al primo ordine produce ampiezza e Maslov; non resta una
correzione lineare all'azione quantizzata. Il potenziale subprincipale
esplicito entra a epsilon^2. Questa è la premessa WKB standard necessaria:
la sola cancellazione di L nel membro destro NON sarebbe una prova.

Ponendo Ahat=A0+epsilon A1+O(epsilon^2), l'ordine epsilon impone

    A1 I_A(A0)=0,
    I_A(A0)=1/2 integral_theta-^theta+ Q0^(-1/2) dtheta>0.

L'integrale converge per turning point semplici; i termini estremali nella
prima variazione dell'azione sono nulli perché sqrt(Q0) si annulla. Dunque
A1=0 nella serie formale, sotto le ipotesi dichiarate.

Questo NON dice che I2 sia nullo: il passaggio dell'handoff che parla di
annullamento della correzione successiva va inteso come assenza dell'ordine
lineare, non come annullamento della prima correzione quadratica.
Non è una stima uniforme del resto, né una prova per tutti gli spin,
turning point coalescenti, limiti mu=0,|mu|=1 o regimi Kerr estremali.
L'estensione locale a chat complesso richiede continuazione di ciclo e ramo
senza degenerazioni e I_A diverso da zero, non la positività dell'integrale.

La derivazione era già in kerr_handoff_progress_2026-09-12.md §2. Resta
collegata alla WKB standard e al precedente Yang; non attribuiamo priorità
a questa conseguenza. Il risultato non autorizza la frase generale
"non c'è epsilon1 bosonico" per settori non esaminati.

### Due dettagli del nuovo controllo Mathematica

Il codice Claude chiama SpheroidalEigenvalue[ell,m,I c], non la stessa
funzione con argomento c: la convenzione completa è
A(c)=SpheroidalEigenvalue[ell,m,I c]-c^2. Nelle formule per la derivata della
funzione built-in va quindi incluso anche il fattore I della regola della
catena. Le tabelle possono essere corrette ma la convenzione abbreviata nel
resoconto è ambigua se si omette l'argomento immaginario.

Inoltre A0 viene estrapolato con due punti usando una legge in1/L^2.
È una verifica di coerenza, non una prova indipendente che il termine1/L
manchi. Il fit libero A1 della collocazione Codex e la dimostrazione sopra
restano controlli distinti; non vanno sostituiti con l'assunzione del fit.

## 2. Localizzazione indipendente dello zero complesso del profilo

Claude stima R≈1.69e-3 dai rapporti dei coefficienti di D_minus. Per
verificare direttamente il meccanismo ho costruito la soluzione regolare
Frobenius, non una nuova serie di D:

    rho=r-r_plus, psi=rho^s F(rho), s=-ik/v1, F(0)=1,
    v^2 F''+(2s v^2/rho+v v_r)F'
        +[s(s-1)v^2/rho^2+s v v_r/rho+q]F=0.

Moltiplicando per H^4 i coefficienti diventano polinomi; la ricorrenza
lineare sui coefficienti di F non passa per la Riccati e attraversa gli
zeri di F senza poli. Il codice usa a=3/5,r_plus=9/5 esatti, omega=3-i/4,
ell=28,m=19. L'autovalore angolare è importato in double dal solutore già
validato: la precisione multipla del resto NON promuove quell'input a80 cifre.

`calculations/kerr_horizon_zero_independent.py` restituisce per gradi30/50/70
e precisioni50/60/80 lo stesso zero della serie dell'ampiezza:

    rho_zero≈0.000815051682825084-0.001482924177138700 i,
    |rho_zero|≈0.001692150513641884,
    |F'(rho_zero)|≈614.809043574.

I residui delle serie troncate sono inferiori a1e-50 nella prima esecuzione,
ma non sono errori certificati dello zero della soluzione esatta. La
stabilità a gradi diversi sostiene il risultato; le cifre fisicamente
affidabili sono limitate anche dall'autovalore double e dalla precisione
dell'implementazione. Script terminato con exit0.

F' diverso da zero indica uno zero semplice. Poiché rho_zero non coincide
con una singolarità del coefficiente v, D=v(s/rho+F'/F) ha un polo lì,
con residuo v(r_plus+rho_zero). Il suo sviluppo centrato sull'orizzonte non
può avere raggio maggiore di |rho_zero|. L'accordo con il rapporto dei
coefficienti conferma la diagnosi di Claude, ma non abbiamo contato tutti
gli zeri in un disco: non è una prova che questo sia il polo più vicino.

È importante che lo zero sia COMPLESSO: non è un nodo del profilo sul
segmento radiale reale. Un polo complesso limita comunque il raggio di Taylor.

## 3. Conseguenze operative e lavoro rimanente

Per il matching reale conviene propagare il sistema lineare per psi,psi_*
o l'ampiezza Frobenius regolare oltre la zona di Taylor di D. Le grandezze
fisiche non cessano di esistere a rho reale=0.00169; cede quella specifica
rappresentazione in serie della derivata logaritmica.

Restano aperti: selezione globale del ramo uscente Kerr, norma completa con
bordi controllati, confronto editoriale di priorità e revisione del testo.
I risultati G di Claude sono letti ma non rieseguiti in questo turno.
Non sono state effettuate nuove simulazioni Vaidya, commit o push.

Richiesta numerica mirata per Claude: controllare lo zero con l'ampiezza
Frobenius usando un autovalore angolare multiprecisione indipendente; se
serve affermare un raggio ESATTO, contare gli zeri via principio dell'argomento
su cerchi interni/esterni al candidato evitando singolarità dei coefficienti.
Non è necessario rifare da capo i gate A–G già documentati.
