# D_minus e D_plus per Claude: definizione, coefficienti e ricorrenze

## Quale passo è bloccato?

Il passo E1–E3 del brief controlla A_c e q_omega a ell,m,a fissati: NON
richiede D_plus/minus, integrazione radiale o una nuova frequenza QNM esatta.
La norma spettrale globale Kerr e il controllo F richiedono invece i D ai
bordi finiti. La nota precedente ne dava la definizione ma non una ricorrenza
operativa: questa nota colma la lacuna, senza chiamare esatte serie troncate.

## 1. Convenzioni e definizione esatta

M=1, |a|<1, tempo exp(-i omega t), ell,m fissati. Qui D è la derivata
logaritmica della variabile master psi=sqrt(r^2+a^2) R, NON della radiale R.

    H=r^2+a^2, Delta=r^2-2r+a^2, v=Delta/H, partial_*=v partial_r,
    lambda=A_lm(a omega)+a^2 omega^2-2am omega,
    G=v[v_r r/H+v a^2/H^2],
    q=(omega-ma/H)^2-Delta lambda/H^2-G,
    D_minus(r,omega)=v partial_r log(psi_in),
    D_plus(r,omega)=v partial_r log(psi_out).

psi_in è il ramo entrante all'orizzonte, psi_out quello uscente all'infinito.
Sono funzioni esterne definite anche FUORI dalla frequenza QNM: servono
proprio a costruire il residuo e le sue derivate. Entrambe soddisfano

    v D_r+D^2+q=0.

Se si usa R invece di psi, occorre aggiungere v r/H alla sua derivata
logaritmica. Una diversa variabile master cambia i coefficienti al bordo.

## 2. Orizzonte: coefficienti espliciti

    r_plus=1+sqrt(1-a^2), r_minus=1-sqrt(1-a^2),
    gap=r_plus-r_minus, H_plus=r_plus^2+a^2,
    Omega_H=a/H_plus, k=omega-m Omega_H,
    rho=r-r_plus, v1=gap/H_plus.

Il ramo entrante ha psi_in proporzionale a rho^(-i k/v1) per una serie
regolare. La sua derivata logaritmica ammette, lontano da risonanze della
ricorrenza e da zeri del profilo,

    D_minus=-i k+h1 rho+O(rho^2),
    h1=-q1/(v1-2i k),
    q1=4ma r_plus k/H_plus^2-gap lambda/H_plus^2
           -gap^2 r_plus/H_plus^3.

Il termine finale di q1 viene da G_r all'orizzonte e non va omesso.
Per la derivata in omega, con a,ell,m,r fissati:

    lambda_omega=a A_c+2a^2 omega-2am,
    q1_omega=4ma r_plus/H_plus^2-gap lambda_omega/H_plus^2,
    h1_omega=-q1_omega/(v1-2ik)-2i q1/(v1-2ik)^2,
    D_minus,omega=-i+h1_omega rho+O(rho^2).

### Ricorrenza a ordine arbitrario

Sviluppare v=sum_{j>=1}v_j rho^j, q=sum_{j>=0}q_j rho^j,
D_minus=sum_{n>=0}h_n rho^n, h0=-ik. I coefficienti sono quelli di Taylor
(includono divisione per il fattoriale), non le derivate grezze.

    h_n=-[q_n+sum_{i=1}^{n-1}h_i h_{n-i}
              +sum_{j=2}^{n}v_j(n-j+1)h_{n-j+1}]/(n v1+2h0), n>=1.

Somme vuote uguali a zero. Se un denominatore si annulla o è molto piccolo,
non dividere alla cieca: ramo Frobenius/risonanze richiedono analisi dedicata.
Il limite estremale v1=0 non è coperto. L'espansione si usa vicino
all'orizzonte, poi si propaga fino al punto di matching se necessario.

## 3. Infinito: coefficienti espliciti

Nella coordinata tortoise r_*, non nel solo r:

    q=omega^2-B/r^2+(2lambda-2)/r^3+O(r^-4),
    B=A_lm(a omega)+a^2 omega^2=lambda+2am omega,
    v=1-2/r+O(r^-3).

Per omega non nulla, il ramo uscente dà la serie formale

    D_plus=i omega+t2/r^2+t3/r^3+O(r^-4),
    t1=0, t2=B/(2i omega),
    t3=[2t2-(2lambda-2)]/(2i omega).

Non compare il termine1/r per D relativo a r_* e a psi. Non usare queste
formule per partial_r log R: lì le potenze e i termini geometrici differiscono.

    B_omega=a A_c+2a^2 omega,
    t2_omega=B_omega/(2i omega)-B/(2i omega^2),
    t3_omega=[2t2_omega-2lambda_omega]/(2i omega)
                 -[2t2-(2lambda-2)]/(2i omega^2),
    D_plus,omega=i+t2_omega/r^2+t3_omega/r^3+O(r^-4).

### Ricorrenza

Scrivere v=sum_{j>=0}v_j r^-j, q=sum_{j>=0}q_j r^-j,
D_plus=sum_{n>=0}t_n r^-n, con v0=1,q0=omega^2,t0=i omega.

    t_n=[sum_{k=1}^{n-1} k t_k v_{n-k-1}
              -sum_{k=1}^{n-1} t_k t_{n-k}-q_n]/(2i omega), n>=1.

Questa è una SERIE ASINTOTICA, non una promessa di convergenza per n->infinito.
Le derivate t_n,omega si ottengono differenziando la ricorrenza e A insieme.

## 4. Come ottenere le mappe a bordi finiti

Le espressioni esatte sono i rapporti delle soluzioni esterne specificate
nella §1; non sono in generale funzioni elementari di r e omega. Le serie
sono inizializzazioni approssimate, da controllare separatamente da A_c.

Procedura operativa: costruire coefficienti a più ordini, inizializzare vicino
a r_plus e a un raggio esterno R, e integrare le due soluzioni verso il matching.
Propagare anche K=D_omega mediante

    v K_r+2D K=-q_omega.

Usare i dati di K derivati dalle STESSE serie, non congelare A. In presenza
di poli di D usare il sistema lineare per psi,psi_* e relative variazioni,
oppure cambiare punto di matching.

Attenzione specifica ai QNM con Im omega<0: il ramo uscente cresce lungo
l'asse reale e una contaminazione entrante è esponenzialmente subdominante.
La sola serie in1/r, o il solo limite D->i omega, non seleziona univocamente
il ramo in quel semipiano. La prescrizione fisica è la continuazione del
ramo uscente dalla regione Im omega>0 (lontano da tagli/singolarità), oppure
un metodo equivalente di Leaver/Jost o un contorno complesso giustificato.
La stabilità cambiando R e ordine non prova da sola l'assenza di contaminazione
oltre tutti gli ordini. Per precisione QNM globale serve confronto indipendente.

Non promettiamo quindi una norma Kerr esatta usando soltanto +/-i omega
o tre termini della serie a raggi finiti. Per E1–E3 questo problema non
interviene: Claude può completare subito i controlli angolari e di q_omega.
Per D_eta occorre specificare la perturbazione eta anche nelle regioni
esterne: non c'è un'espressione universale senza tale scelta.

## 5. Verifica eseguita e controlli prescritti

`verify_kerr_boundary_series.wl`: WolframKernel con licenza, **3/3 True**.
Controllati q fino a r^-3, residuo Riccati uscente fino a r^-3 e q1
all'orizzonte con Delta(r_plus)=0. Non sono test numerici delle mappe globali.

Per Claude: E1–E3 non sono bloccati; per il successivo problema ai bordi
implementare le ricorrenze, confrontare ordini4/6/8 e posizioni del bordo,
verificare residui locali coerenti con il primo termine omesso e controllare
K con differenze in omega sul MEDESIMO ramo. Escludere inizialmente a≈1,
omega≈0 e denominatori Frobenius piccoli. Non dichiarare chiusa la norma
globale senza il controllo della selezione uscente descritto sopra.
