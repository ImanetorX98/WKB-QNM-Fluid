# Secondo ordine: curvatura della risposta e zeri mobili del diagnostico

## Stato e convenzioni

Estensione analitica di `analytic_core_response_2026-09-12.md` e
`kerr_nonlinear_spectral_proof_2026-09-12.md`. Nella regione interna la
perturbazione esterna agisce solo tramite omega(eta); coordinate, peso,
finestra, parametro eps e altri parametri del potenziale sono fissati.
Supponiamo un ramo semplice sufficientemente regolare e assenza di zeri
di psi nella regione dove si usa la derivata logaritmica.

Definiamo d=omega_eta(0), e=omega_eta_eta(0). Sono DERIVATE: l'espansione
è omega=omega0+eta d+eta^2 e/2+..., non con coefficiente quadratico e.
Le formule sono conseguenze della differenziazione dell'ODE e di integrali
con zeri mobili, non rivendicazioni di un nuovo formalismo QNM.

## 1. Curvatura del ramo spettrale

Per un residuo spettrale F(omega,eta) con F_omega non nullo:

    d=-F_eta/F_omega,
    e=-(F_eta_eta+2F_omega_eta d+F_omega_omega d^2)/F_omega.

Tutte le derivate sono valutate alla radice imperturbata e includono la
dipendenza dei bordi. Derivare un residuo con condizioni al bordo congelate
non equivale a derivare il problema QNM. La formula segue derivando due volte
F(omega(eta),eta)=0. È invariata sotto moltiplicazione di F per un fattore
analitico non nullo se TUTTE le derivate vengono trasformate coerentemente.

## 2. Seconda risposta del profilo interno

Posto z=psi'/psi, k=z_omega, j=z_omega_omega, dalla Riccati segue

    k'+2zk=-q_omega,
    j'+2zj=-q_omega_omega-2k^2,
    k(a)=D_minus,omega, j(a)=D_minus,omega_omega.

Quindi

    z_eta=kd,
    z_eta_eta=ke+jd^2.

La forma integrale per j è

    j(x)=psi(x)^(-2)[D_minus,omega_omega psi(a)^2
                  -integral_a^x (q_omega_omega+2k^2)psi^2 dt].

Nel benchmark compatto q=omega^2-V0: q_omega=2omega,q_omega_omega=2,
D_minus=-iomega, dunque k(a)=-i,j(a)=0. Per Kerr occorre invece anche
la seconda derivata dell'autovalore sferoidale: non porre q_omega_omega=2.

## 3. Seconda risposta di Q_M

Per p=Im z e Q=eps^2(Re q-p^2), definiamo

    p1=Im(kd), p2=Im(ke+jd^2),
    Q1=eps^2[Re(q_omega d)-2p p1],
    Q2=eps^2[Re(q_omega e+q_omega_omega d^2)-2p1^2-2p p2].

Q1 e Q2 sono derivate rispetto a eta, non i coefficienti senza fattoriali.
La derivazione usa la regola della catena, q privo di dipendenza esplicita
da eta nella finestra. Se la perturbazione entra nella regione interna,
compaiono q_eta,q_omega_eta,q_eta_eta e la fattorizzazione non vale.

### Rango: due al primo ordine non significa due a perturbazione finita

La risposta lineare appartiene allo spazio generato dalle due componenti
reali di d. Al secondo ordine entrano le due componenti di e e i tre
monomi (Re d)^2,(Re d)(Im d),(Im d)^2: il rango delle colonne di seconda
risposta è AL PIÙ CINQUE, su un medesimo background/finestra.
Non si afferma che sia esattamente cinque; degenerazioni possono ridurlo.
Il profilo completo continua a dipendere da due parametri reali: la crescita
del rango dello sviluppo riflette la curvatura di quella famiglia, non
l'apparizione di nuove variabili fisiche indipendenti.

Per differenze finite simmetriche, (Q(eta)-Q(-eta))/(2eta)=Q1+O(eta^2).
Una terza singolare non nulla può quindi essere troncamento cubico e/o
errore numerico: va studiata al variare di eta, non interpretata subito
come fallimento della proposizione lineare.

## 4. Seconda derivata dell'integrale: il termine degli zeri mobili

Sia U(eta)=integral_J w(x)|Q(x,eta)|dx con w e J fissati. Supponiamo che
gli zeri x_j di Q(x,0) interni a J siano finiti e semplici (Q_x non nullo),
separati dagli estremi, e che non si creino/coalescano zeri per piccoli eta.
La funzione implicita dà x_j'(0)=-Q1(x_j)/Q_x(x_j).

Spezzando l'integrale in tratti di segno costante, la prima derivata non
contiene termini agli zeri perché |Q|=0. Nella seconda derivata, invece,
l'integrando della prima derivata ha un salto. Il risultato è

    U1=integral_J w sign(Q) Q1 dx,
    U2=integral_J w sign(Q) Q2 dx
           +2 sum_j w(x_j) Q1(x_j)^2/|Q_x(x_j)|.

Il termine aggiuntivo è non negativo per w non negativo. Equivalentemente
segue da partial_eta sign(Q)=2delta(Q)Q1 e dalla regola per delta(Q(x)),
ma la suddivisione dell'integrale evita di fondare la prova solo su una
notazione distribuzionale. Per zeri multipli o che attraversano gli estremi
questa formula non è automaticamente valida: vanno riesaminate regolarità
e differenziabilità. Il peso non deve dipendere da eta in questa forma.

Controesempio esatto all'omissione del termine: Q=x-eta,w=1,J=[-1,1].
Per |eta|<1, U=1+eta^2 e U2=2, mentre Q2=0. L'intero contributo viene
dallo zero mobile. Questo controllo non coinvolge approssimazioni QNM.

## 5. Rapporto integrale del benchmark compatto

Nel benchmark eps=1 poniamo

    D=integral_J w (|omega|^2+|V0|+p^2) dx,
    D1=integral_J w [2Re(conj(omega)d)+2p p1] dx,
    D2=integral_J w [2|d|^2+2Re(conj(omega)e)+2p1^2+2p p2] dx.

Per E=U/D, D>0:

    E1=(U1-E D1)/D,
    E2=(U2-E D2-2E1 D1)/D.

La mediana pesata non è un integrale di |Q|: queste formule NON vanno
applicate automaticamente al quantile discreto. Anche l'eventuale dipendenza
del peso/larghezza della finestra dal potenziale deve essere derivata.

## 6. Verifiche e limiti effettivi di questo turno

`verify_second_order_response.wl`, eseguito con WolframKernel con licenza:
**4/4 controlli True**, exit0. Controlla espansione del quadrato della
Riccati, Q2, regola del quoziente e controesempio esatto dello zero mobile.
Le ipotesi di regolarità e la prova generale sugli zeri sono quelle sopra,
non un teorema certificato automaticamente dal CAS.

Ho letto anche la parte finale del resoconto Claude: riporta una risposta
del diagnostico ma debole correlazione con l'errore WKB sul suo campione.
Non ho rieseguito quegli script qui. Una correlazione debole non prova un
no-go universale e una correlazione negativa non è di per sé "sbagliata"
senza una relazione monotona postulata in anticipo. Il test limita la pretesa
di usare quei diagnostici come stimatori monotoni crescenti dell'errore.

Non è stato calcolato numericamente E2 sul benchmark: il prossimo controllo
è prescritto per Claude, includendo esplicitamente i termini degli zeri.
Nessuna nuova verifica globale Kerr o simulazione Vaidya in questo turno.
Il core metodologico è più preciso; originalità e pubblicabilità restano aperte.
