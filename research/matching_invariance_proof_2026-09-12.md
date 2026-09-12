# Indipendenza dal matching e semplicità della risonanza

## Contesto e ipotesi

Completiamo `kerr_nonlinear_spectral_proof_2026-09-12.md` per
psi''+q(x,omega,eta)psi=0 in una coordinata fissata. D_minus,D_plus
sono le derivate logaritmiche delle soluzioni esterne esatte. Consideriamo
punti senza poli di D e una radice con soluzione non nulla ai bordi.
Le derivate in omega sono complesse sul ramo analitico scelto.
La normalizzazione globale della soluzione è fissata durante lo spostamento
dei punti di matching; se viene cambiata, gli integrali bilineari riscalano.

## 1. Il denominatore non dipende dal punto di matching

Indichiamo con

    N(a,b)=integral_a^b q_omega psi^2 dx
            +D_plus,omega(b)psi(b)^2-D_minus,omega(a)psi(a)^2.

La Riccati esterna D_x+D^2+q=0 implica

    (D_omega)_x+2D D_omega=-q_omega.

Al QNM la soluzione interna è proporzionale alla soluzione esterna, quindi
psi'=D psi nella regione in cui si sposta il punto. Ne segue

    partial_x(D_omega psi^2)=-q_omega psi^2.

Al bordo destro la derivata dell'integrale è +q_omega psi^2 e quella del
termine superficiale è -q_omega psi^2: si cancellano. Al bordo sinistro
i segni si invertono e si cancellano nuovamente. Dunque partial_a N=partial_b N=0.

Questo non significa che N sia una norma positiva: è una quantità bilineare
complessa. Inoltre l'indipendenza vale con condizioni esterne esatte;
usare un D asintotico troncato a raggio finito lascia un difetto misurabile.
Non è lecito spostare il bordo oltre una perturbazione ignorandola nel D.

## 2. Anche il numeratore è indipendente

Per

    B(a,b)=-integral_a^b q_eta psi^2 dx
             -D_plus,eta(b)psi(b)^2+D_minus,eta(a)psi(a)^2,

la Riccati differenziata a omega fissata dà
(D_eta)_x+2D D_eta=-q_eta. La medesima cancellazione prova
partial_a B=partial_b B=0. Quindi d omega/d eta=B/N è indipendente dal
matching. Se psi viene moltiplicata per C, sia B sia N si moltiplicano per C^2
e il rapporto resta invariato. Non serve una normalizzazione probabilistica.

Questa dimostrazione rende esplicito, per q dipendente da omega, il
meccanismo già presente nel formalismo di Leung et al., eq.2.14–2.16.
Non è una rivendicazione di originalità.

## 3. Relazione con la derivata del residuo spettrale

Fissiamo psi(a;omega)=1 e psi'(a;omega)=D_minus(omega), integrando a destra.
Sia F(omega)=psi'(b;omega)-D_plus(omega)psi(b;omega).
Al QNM F(omega0)=0. Per phi=partial_omega psi,

    W=psi phi'-psi'phi, W'=-q_omega psi^2,
    W(a)=D_minus,omega psi(a)^2,
    W(b)=psi(b)F_omega+D_plus,omega psi(b)^2.

Integrando segue esattamente

    N=-psi(b) F_omega.

Poiché psi(b) non è nulla quando il D scelto è finito, N non nullo equivale
a F_omega non nullo: una radice semplice di QUESTO residuo. Se F è olomorfo
in omega e regolare in eta, il teorema della funzione implicita giustifica
la continuazione locale e la derivata perturbativa. La derivata del residuo
dipende dalla normalizzazione adottata; non confrontarla fra codici senza
allineare tale scelta. Per un Wronskiano normalizzato diverso da F compare
il corrispondente fattore di normalizzazione.

Questa relazione chiude una precedente ipotesi lasciata separata: non occorre
postulare indipendentemente semplicità della radice e N non nullo in questa
rappresentazione regolare. Non autorizza a dividere per N ai punti eccezionali.

## 4. Conseguenza per l'error budget

Vicino a una radice semplice, a primo ordine delta omega≈-delta F/F_omega.
Un residuo piccolo da solo non certifica un errore spettrale piccolo: occorre
conoscere F_omega e l'errore numerico nella costruzione di F. La stima è locale,
non un bound rigoroso in presenza di altri zeri vicini o termini non lineari.
Una derivata grande/piccola dipende dalla scalatura di F: riportare insieme
la normalizzazione del residuo e quella del suo errore.

## 5. Collegamento con i nuovi file di Claude

Durante questo turno sono comparsi sei file Claude già staged, inclusa
`research/claude_compact_validation.md`. Ho letto il resoconto fino ai test
di rango: riporta un controllo di spostamento dei bordi e quadratura spezzata.
Questi esiti sono risultati riportati da Claude, non rieseguiti qui.
La frase "nessun raffinamento uniforme la salva" non va interpretata come
teorema di non convergenza: il difetto di regolarità riduce l'efficienza e
l'ordine, non impedisce di per sé la convergenza di quadrature appropriate.
Non ho modificato né incluso nel mio commit i sei file staged da Claude.

## Verifiche

`calculations/verify_matching_invariance.wl` controlla le cancellazioni
algebriche delle derivate superficiali e la relazione finale fra N e F_omega.
Il terzo controllo è solo algebra della relazione: la sua derivazione
tramite Wronskiano è quella scritta sopra, non una prova automatica CAS.
La revisione bibliografica di priorità e la validazione globale Kerr/Vaidya
restano aperte. Non è stata eseguita una nuova simulazione Vaidya.

Verifica pre-commit: **12 test Python superati** (angolare indipendente,
scaling radiale a mu esatto, geometria Kerr, barriera compatta); **17 controlli
CAS superati** nei cinque script Wolfram del lavoro, incluso quello nuovo.
Non è l'intera suite del repository e non include la riesecuzione dei file Claude.
