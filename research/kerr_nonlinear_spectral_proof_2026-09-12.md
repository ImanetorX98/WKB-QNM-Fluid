# Estensione analitica al problema spettrale dipendente da frequenza

## Avanzamento e attribuzione

Proseguimento di `analytic_core_response_2026-09-12.md`. Letti anche i passaggi
§2.2–2.3, eq.2.19–2.21 di Leung et al., dal testo già disponibile nel progetto.
La loro eq.2.21 contiene già la risposta della derivata logaritmica della
funzione d'onda: il nostro k e la conseguente fattorizzazione non sono un
nuovo formalismo perturbativo. La formulazione del test di rango può essere
utile operativamente; la sua priorità non è stabilita.

Questa nota generalizza esplicitamente l'algebra a q dipendente da omega e
la specializza alla radiale scalare Kerr. Non certifica numericamente la
risposta dei QNM Kerr né una nuova previsione fisica. La derivazione è una
specializzazione del metodo standard del Wronskiano/perturbazione spettrale.

## 1. Formula con bordi e operatore non lineare nel parametro spettrale

In coordinata x fissata consideriamo

    psi'' + q(x,omega,eta) psi = 0,
    psi'(a)=D_minus(omega,eta) psi(a),
    psi'(b)=D_plus(omega,eta) psi(b).

I D sono le derivate logaritmiche delle soluzioni esterne esatte, non
necessariamente +/-i omega a bordi FINITI. I punti a,b sono fissati, e le
funzioni e il ramo QNM sono differenziabili; omega è una radice semplice.
Per il seguito richiediamo anche denominatore non nullo e rappresentazioni
D prive di poli nei punti scelti. Si può cambiare punto di matching se serve.

Con phi=partial_eta psi, d=partial_eta omega, il Wronskiano dà

    (psi phi'-psi'phi)'=-(q_eta+q_omega d)psi^2.

Al bordo destro W=(D_plus,eta+D_plus,omega d)psi(b)^2, e analogamente
al sinistro. Pertanto

    d = [-integral_a^b q_eta psi^2 dx
         -D_plus,eta psi(b)^2+D_minus,eta psi(a)^2] / N,

    N = integral_a^b q_omega psi^2 dx
        +D_plus,omega psi(b)^2-D_minus,omega psi(a)^2.

La formula precedente a potenziale compatto segue da q=omega^2-V0-eta B,
D_plus=i omega,D_minus=-i omega. Non si può sostituire q_omega con2omega
se il potenziale effettivo contiene omega. Non si può usare la stessa
formula senza termini aggiuntivi se si cambiano coordinata o estremi con eta.
Su domini non compatti bisogna mantenere i D esatti a bordi finiti, oppure
giustificare una regolarizzazione equivalente: gli integrali QNM divergenti
non si trattano come ordinari integrali L2.

## 2. Specializzazione alla radiale scalare Kerr

M=1, a spin fissato, ell,m fissati, x=r_* con D_*=(Delta/H)partial_r,
H=r^2+a^2, Delta=r^2-2r+a^2. Nella convenzione del repository:

    q=(omega-ma/H)^2-Delta/H^2 [A_lm(a omega)+a^2 omega^2-2am omega]-G(r),
    G=D_*^2 sqrt(H)/sqrt(H).

A dipende da omega sul ramo angolare scelto. La derivata TOTALE a r,a,ell,m
fissati è

    q_omega=2(omega-ma/H)
       -Delta/H^2 [a A_c(a omega)+2a^2 omega-2am].

G non contribuisce perché qui a e la coordinata tortoise sono fissati.
Al limite a=0 si recupera q_omega=2omega. Per perturbazioni della geometria
che cambiano a o M non basta aggiungere un q_eta ingenuo: bisogna specificare
identificazione delle coordinate e condizioni ai bordi.

All'orizzonte il comportamento asintotico è D_minus=-i(omega-m Omega_H),
quindi D_minus,omega=-i nel limite; all'infinito D_plus,omega=+i nel limite.
Ai raggi di matching FINITI restano le correzioni dei profili esterni.
Un bordo troncato al primo ordine non è un D esatto a precisione arbitraria.

### Derivata angolare senza differenze finite dell'autovalore

L'equazione è

    (sin(theta) S')' + sin(theta)[A+c^2 cos(theta)^2-m^2/sin(theta)^2] S=0.

Per una famiglia angolare regolare e differenziabile, deriviamo in c,
moltiplichiamo per S e sottraiamo l'equazione originale moltiplicata per S_c.
Integrando su [0,pi] il concomitante al bordo si annulla per le condizioni
regolari e m fissato. Se integral S^2 sin(theta) non è nullo:

    A_c=-2c [integral_0^pi S^2 cos(theta)^2 sin(theta) dtheta] /
               [integral_0^pi S^2 sin(theta) dtheta].

Per c complesso è un pairing BILINEARE senza coniugazione. È la formula
Hellmann–Feynman per questo operatore formalmente simmetrico; non è nuova.
In una base armonica ortonormale la matrice è complessa simmetrica, non
hermitiana: usare il trasposto, non automaticamente il trasposto coniugato.
Vicino a degenerazioni/auto-ortogonalità questa rappresentazione va riesaminata.

## 3. Fattorizzazione del diagnostico con code non compatte

L'argomento di unicità non richiede un potenziale compatto. È sufficiente che
in una regione contenente il bordo sinistro e la finestra J l'operatore e
la condizione logaritmica esatta D_minus dipendano da eta SOLO tramite omega.
Normalizzata la soluzione al punto a:

    z_eta(x)=z(x;omega_eta),
    k'+2zk=-q_omega,
    k(a)=D_minus,omega,
    k(x)=psi(x)^(-2)[D_minus,omega psi(a)^2
                     -integral_a^x q_omega psi^2 dt].

Questo sostituisce il precedente integrando2omega. Il rango reale al più
due per i diagnostici reali differenziabili resta valido sotto le stesse
ipotesi; il carattere complesso di A non aumenta i parametri indipendenti,
finché A è una funzione del solo omega sul ramo fissato.

Ciò copre un problema matematico di scattering scalare Kerr con deformazione
esterna della radiale e settore interno/angulari invariati. NON prova che
una tale deformazione corrisponda a una soluzione di Einstein separabile:
la realizzabilità gravitazionale è un problema aggiuntivo. Nessuna pretesa
di aver esteso qui il risultato a perturbazioni arbitrarie della metrica.

## 4. Perché Vaidya non è un corollario

Nel testo locale di Capuano et al.2407.06009, §II eq.6–10, il caso a tasso
di massa costante è reso conformalmente statico mediante x=r/(2M(w)) e
W=integral dw/(2M(w)), con f(x)=1-1/x-4|M'|x.
Cambiano coordinata temporale, radiale e bordi del problema trasformato.
La frequenza coniugata al tempo trasformato non si confronta direttamente
con una frequenza locale rispetto al tempo avanzato originario.

Perciò la formula Kerr sopra non dimostra un risultato Vaidya: servono
la PDE evolutiva oppure un confronto nel caso a tasso costante dopo avere
allineato trasformazione del campo, tempo e bordi. La precedente connessione
K che cambia sotto G->c(M)G non è di per sé uno shift spettrale osservabile.
Non sono stati eseguiti nuovi conti numerici Vaidya in questo turno.

## 5. Verifiche e prossima soglia

`calculations/verify_kerr_spectral_response.wl` controlla la derivata di q,
il limite Schwarzschild e l'algebra della formula con i bordi. L'integrazione
per parti angolare, la regolarità dei rami e la validità dei D esterni sono
argomenti analitici da verificare nelle ipotesi, non certificazioni CAS.

Eseguito con WolframKernel con licenza: **3/3 controlli True, exit0**.

La prossima soglia numerica mirata è confrontare A_c bilineare e differenze
finite complesse del medesimo ramo, poi q_omega totale e derivata diretta di q.
Questo può controllare anche il termine omesso nel vecchio calcolo eikonale
che differenziava omega tenendo il riferimento angolare fisso. Non è stata
ancora corretta né rimisurata quella routine in questo turno.

Il contributo dimostrativo è dunque più generale e coerente con Kerr, ma
rimane basato su metodi standard. Non rivendichiamo una nuova teoria QNM;
la pubblicabilità richiede ancora un risultato discriminante validato.
