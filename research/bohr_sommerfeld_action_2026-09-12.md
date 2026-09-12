# Quantizzazione angolare di Bohr–Sommerfeld: riferimento indipendente

## Perché studiarla

È la condizione che collega l'equazione angolare al numero di nodi, non un
accessorio della prova A1=0. Letto il passaggio di Yang et al., eq.2.14–2.18
nel testo locale1207.4253: il matching ai turning point produce L-|m| con
L=ell+1/2. È un precedente diretto; non rivendichiamo come nuova questa
quantizzazione. Qui ne calcoliamo l'azione senza un solutore di autovalori.

## Azione principale e regolarizzazione della quadratura

Nel settore scalare a mu=m/L e c_hat fissati,

    Q0=A0+c_hat^2 cos(theta)^2-mu^2/sin(theta)^2,
    I(A0,c_hat,mu)=integral_theta-^theta+ sqrt(Q0) dtheta=pi(1-|mu|).

Ipotesi: parametri reali,0<|mu|<1, due turning point semplici, un intervallo
permesso che attraversa l'equatore. Per c_hat>0, x=cos(theta),

    P(x)=(A0+c_hat^2 x^2)(1-x^2)-mu^2
         =c_hat^2(t_plus-x^2)(x^2-t_minus),
    t_plus/minus=(c_hat^2-A0 +/- sqrt((A0+c_hat^2)^2-4c_hat^2 mu^2))/(2c_hat^2).

In questo ramo 0<t_plus<1 e t_minus<0. Ponendo x=sqrt(t_plus)sin(u),

    I=2c_hat t_plus integral_0^(pi/2)
        cos(u)^2 sqrt(t_plus sin(u)^2-t_minus)/(1-t_plus sin(u)^2) du,
    I_A=integral_0^(pi/2) du/[c_hat sqrt(t_plus sin(u)^2-t_minus)],
    I_c=2 integral_0^(pi/2) t_plus sin(u)^2/sqrt(t_plus sin(u)^2-t_minus) du,
    partial_c A0=-I_c/I_A.

La sostituzione rimuove la singolarità integrabile degli estremi. A c_hat=0,
I=pi(sqrt(A0)-|mu|), da cui A0=1: controllo analitico esatto.

## Risultato numerico dell'azione, non del problema angolare finito

`calculations/angular_bohr_sommerfeld_action.py` usa solo mpmath e le formule
sopra, senza importare autovalori dal progetto. A mu=2/3,c_hat=3/5:

    A0 = 0.897478121218106117251403583,
    partial_c A0 = -0.34922245866756854380202.

Le cifre stampate coincidono a precisioni35/55. Il residuo dell'azione è
4.5e-36/8.2e-56. Una differenza centrale con passo1e-5 della soluzione
dell'azione verifica la derivata con scarto assoluto2.17e-12.
Questa stabilità non è un bound rigoroso a intervalli.

È un riferimento di testa indipendente dai fit. L'A0=0.8974781250072938
estrapolato da Claude con due L finiti differisce di circa3.79e-9: non è
un disaccordo della teoria, ma il termine residuo dell'estrapolazione va
considerato prima di attribuire a L^2(Ahat-A0) cifre significative arbitrarie.

## A1 e il prossimo coefficiente A2

La forma di Liouville esatta è epsilon^2 psi''+(Q0+epsilon^2 Q2)psi=0,
Q2=(1+csc(theta)^2)/4. Dopo Maslov, la quantizzazione FORMALE è

    I(Ahat)+epsilon^2 I2(Ahat)+...=pi(1-|mu|).

Perciò A1 I_A=0 e A1=0; al passo successivo A2=-I2(A0)/I_A(A0).
Non basta la condizione principale per conoscere I2, e non è corretto
porlo zero. La ricorsione del momento p=sqrt(Q0) fornisce l'integrando
formale del secondo periodo:

    p2=Q2/(2p)-p''/(4p^2)+3p'^2/(8p^3),

con derivate in theta a Ahat fissato. I singoli termini non sono integrabili
ordinariamente fino ai turning point. Serve un periodo su ciclo complesso
che circonda il taglio (con normalizzazione che restituisca I al primo
ordine), oppure una regolarizzazione equivalente derivata con il matching.
NON dare a un integratore reale questa formula tagliando arbitrariamente
gli estremi: il cutoff diventerebbe un parametro spurio del coefficiente.

Il calcolo controllato di I2 resta aperto in questa nota; non lo abbiamo
sostituito con un fit. Il limite c_hat=0 deve dare Ahat=1-epsilon^2/4,
quindi A2=-1/4, poiché A=ell(ell+1) esatto. È un test vincolante di segni,
normalizzazione del ciclo e termini ai poli angolari.

## Prescrizione numerica aggiornata

Usare il nuovo A0 da azione nel confronto con autovalori a mu esatto;
non estrapolare A0 assumendo già assente A1. Ripetere il fit libero A1,A2
e studiarne la stabilità. Questo può procedere mentre viene definito il
secondo periodo regolarizzato. Non promettere un teorema uniforme nei limiti
di turning point coalescenti o mu=0,|mu|=1: non sono trattati qui.

Nessuna nuova simulazione Kerr radiale o Vaidya in questo turno.
