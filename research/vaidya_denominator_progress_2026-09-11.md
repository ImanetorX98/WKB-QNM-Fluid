# Vaidya: denominatore verificato e primo candidato di trasporto

11 settembre 2026. Seguito di kerr_fixed_mu_and_vaidya_transport_2026-09-11.md.
Unità M=1; modo fondamentale scalare ell=2; normalizzazione del ramo
entrante G(2)=1. Nessuna previsione osservativa è dichiarata in questa nota.

## 1. Nuovo calcolo indipendente del denominatore

Calcoliamo

~~~text
P = Reg integral mu G G_r dr,
mu = exp(-2i omega r_*), G = exp(i omega r_*) R.
~~~

Il codice usa lo stesso profilo fattorizzato del numeratore, ma un integrando
e primitive propri: al bordo interno il prodotto h h' viene sviluppato in
Frobenius e integrato termine a termine; al bordo esterno

~~~text
P_integrand ~ C^2 exp(2i omega r_*) Wp,
Wp = u [2i omega u/f + u'].
~~~

La primitiva esterna è costruita con algebra di serie e sottratta usando gli
stessi estremi della quadratura. Nessuna differenza finita sull'integrando.

Separatamente si calcola la norma generalizzata su un intervallo finito,

~~~text
Norm = integral_a^b R^2 dr/f
       + [Dplus_omega R(b)^2 - Dminus_omega R(a)^2]/(2 omega).
~~~

D è la derivata logaritmica rispetto a r* del ramo selezionato. Dminus
proviene dalla serie entrante di Frobenius; Dplus dalla serie uscente.
La derivata in frequenza di Dplus è calcolata derivando analiticamente la
ricorrenza dei coefficienti. Usare differenze finite per questo termine dava
un errore crescente con il cutoff, amplificato da R(b)^2.
Per Dminus si mantiene una differenza centrata di passo 1e-6 |omega|, sufficiente
ai livelli di accuratezza misurati qui; non ne rivendichiamo precisione arbitraria.

## 2. Identità che lega i due calcoli

In coordinate r* l'integrando di P è R R' + i omega R^2. Se il potenziale
Schwarzschild è indipendente dalla frequenza, il ramo radiale soddisfa

~~~text
D' + D^2 + omega^2 - V = 0,
(D_omega)' + 2 D D_omega = -2 omega.
~~~

Quindi le primitive locali sono

~~~text
F_norm = -D_omega R^2/(2 omega),
F_P = (1-i D_omega) R^2/2.
~~~

Derivandole rispetto a r* si ottengono rispettivamente R^2 e
R R'+i omega R^2. Le due identità sono verificate con Mathematica.
Applicando i termini di bordo dei rami entrante e uscente in modo coerente,
le parti R^2/2 si cancellano e si ottiene

~~~text
P = i omega Norm.
~~~

Per il fondamentale la primitiva interna di P tende a zero all'orizzonte;
l'integrale P è integrabile lì. La norma invece richiede il suo termine
di bordo. È importante non trasferire l'integrabilità di un integrando
all'altro senza controllare il fattore 1/f e le derivate.

Questa derivazione è per il problema radiale congelato a potenziale
indipendente da omega. Non va applicata senza modifiche al potenziale Kerr
dipendente dalla frequenza e dall'autovalore angolare.

## 3. Risultati numerici

Il raccordo interno è a r=2.5, il cambio di variabile esterna a r=25;
le primitive asintotiche usano 20 termini, Frobenius 60. Il profilo esterno
è integrato fino a 95, con C valutato a r circa 85.

| Cutoff esterno | Re P | Im P | abs(P/(i omega Norm)-1) |
|---:|---:|---:|---:|
| 40 | -3.7327390302 | -3.3753995761 | 8.27e-11 |
| 50 | -3.7327390308 | -3.3753995767 | 5.05e-11 |
| 60 | -3.7327390314 | -3.3753995747 | 2.09e-10 |
| 70 | -3.7327389456 | -3.3753996403 | 1.06e-8 |

A cutoff 40:

~~~text
Norm = -8.1951601953 + 6.0784091969 i,
P    = -3.7327390302 - 3.3753995761 i.
~~~

Spostando width da 0.2 a 1.0, P a cutoff 40 varia di meno di 3e-11 in
modulo assoluto; l'identità con la norma resta entro circa 1e-10 relativo.
Le cifre della tabella servono per riprodurre il calcolo, non sono tutte
certificazioni dell'accuratezza assoluta del modo o del problema dinamico.
L'aumento dello scarto a 70 è coerente con la sottrazione di termini crescenti.

## 4. Primo rapporto di trasporto, ancora condizionale

Combinando N del programma precedente e P alla stessa convenzione e cutoff 40:

~~~text
K_candidate = N/(2P) = 1.1642938694 + 4.3489002268 i.
~~~

Non è la frequenza fisica corretta. Questo rapporto è il candidato per
alpha-i delta_omega = Mdot K solo se la proiezione della sorgente temporale
completa si riduce a N, senza contributi di bordo ulteriori rispetto allo
schema usato. L'identità del denominatore non dimostra quest'ultima ipotesi.

La trasformazione di normalizzazione G -> c(M)G manda K -> K-c_M/c.
La combinazione candidata Xi=K+partial_M ln G elimina questa dipendenza.
Per la famiglia Schwarzschild G(r;M)=h(r/M), valutata a M=1,
partial_M ln G=-r h'/h. Integrando h e h' direttamente fino ai punti indicati:

| r | partial_M ln G | Xi_candidate |
|---:|---|---|
| 3 | -1.2038328312 - 4.3770561878 i | -0.0395389619 - 0.0281559610 i |
| 5 | -1.5460408646 - 6.4460792573 i | -0.3817469952 - 2.0971790305 i |
| 10 | -2.4519025947 - 11.4206464507 i | -1.2876087253 - 7.0717462239 i |

Sono valori diagnostici, non una tabella di frequenze previste. La forte
dipendenza radiale richiama il ruolo della propagazione e della coordinata
temporale: la fase uscente include il tempo di viaggio. Non è di per sé
memoria non locale. La piccola differenza a r=3 è una cancellazione misurata,
non una legge esatta o un annullamento dimostrato.

## 5. Riproduzione e verifiche

~~~sh
python3.13 calculations/vaidya_transport_denominator.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_transport_boundary.wl
~~~

Da calculations:

~~~sh
python3.13 -m unittest test_vaidya_transport_denominator.py test_vaidya_numerator_factored.py
~~~

Tre nuovi test: accordo P=i omega Norm e cutoff, variazione del raccordo
interno, derivata analitica Dplus_omega contro differenze centrate.
Si aggiungono ai dodici test del numeratore. Le due identità delle primitive
sono confermate dal kernel Mathematica (2/2 True).
Esito della suite selezionata: **15 test superati**.

## 6. Prossimo passo necessario

Il denominatore congelato ha ora una verifica analitica e numerica incrociata.
Restano da derivare le condizioni al bordo della correzione temporale e da
verificare che il numeratore regolarizzato impiegato sia proprio il funzionale
Fredholm completo della sorgente. Dopo questo passaggio si potrà interpretare
K e Xi come coefficienti di un'espansione controllata, e confrontarli con una
soluzione nel tempo a massa lentamente variabile. Non abbiamo ancora costruito
quella soluzione o dimostrato originalità/pubblicabilità.
