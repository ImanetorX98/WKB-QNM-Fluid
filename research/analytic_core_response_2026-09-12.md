# Risposta analitica del diagnostico: teoremi, ipotesi e limiti

## Stato

Questo documento estende `core_proof_scope_2026-09-12.md`. Distingue tre
oggetti: il getto locale del potenziale, la frequenza QNM determinata ai due
bordi e un funzionale della soluzione in una finestra. Non sono equivalenti.
La novità bibliografica delle conseguenze qui formulate resta da verificare.

Letture pertinenti: §2.1 eq.2.6–2.16 del testo locale di Leung et al. (1998),
e §I.A eq.1.2 del testo locale Delabaere–Dillinger–Pham (1997). Quest'ultimo
esplicita la parte pari del momento formale Riccati e la parte dispari che
determina l'ampiezza. Non dimostra equivalenza fra soluzione esatta e una
singola correzione WKB troncata. Non sono state lette integralmente in questo
turno le altre fonti Kerr/Vaidya; non dichiariamo chiusa la priorità.

## 1. Risposta spettrale: risultato già noto, qui specializzato

Usiamo psi''+(omega^2-V_eta)psi=0, V_eta=V0+eta b reale e liscio,
a supporto compatto fra a e c. Convenzione temporale exp(-i omega t).
I bordi esatti sono psi'(a)=-i omega psi(a), psi'(c)=i omega psi(c).
Supponiamo omega0 non nulla, una risonanza semplice e un ramo differenziabile
in eta; normalizziamo psi(a)=1. Indichiamo con d=partial_eta omega|0.
Il Wronskiano fra psi e partial_eta psi dà

    d = B/N,
    B = integral_a^c b psi^2 dx,
    N = 2 omega integral_a^c psi^2 dx + i[psi(a)^2+psi(c)^2].

N deve essere non nullo. Non si usa |psi|^2 e non si omettono i bordi.
Per risonanze multiple/punti eccezionali la derivazione non autorizza questa
formula. La derivazione e l'esistenza di bump con B non nullo sono nel
resoconto precedente: si tratta della perturbazione logaritmica di Leung,
non di una nuova legge spettrale.

Se b è nulla vicino al massimo dominante, il suo getto resta invariato a
tutti gli ordini, mentre può essere d diverso da zero. Questo riguarda la
classe C-infinito: per funzioni analitiche su un dominio connesso, identità
di tutti i coefficienti di Taylor ha conseguenze diverse. Non trasferire
silenziosamente il controesempio a una classe di potenziali analitici.

## 2. Proposizione di fattorizzazione del diagnostico interno

Ipotesi aggiuntiva: b=0 su tutto [a,d0] e la finestra J è contenuta in tale
intervallo. a è nel tratto libero a sinistra; il supporto perturbato è a
destra di d0. In [a,d0] la soluzione è il problema a dati iniziali

    u''+(omega^2-V0)u=0, u(a;omega)=1, u'(a;omega)=-i omega.

Per unicità dell'ODE,

    psi_eta(x)=u(x;omega_eta),  x in [a,d0].

Con altra normalizzazione si introduce solo un fattore costante in x;
z=psi'/psi e Q_M ne sono indipendenti. Dunque ogni diagnostico interno
costruito da z, omega, V0 e una finestra/peso fissati soddisfa esattamente

    E_eta=F(Re omega_eta, Im omega_eta).

Questo è un risultato condizionato, non un'affermazione universale sui QNM:
fallisce se la perturbazione cambia il potenziale nella regione da a a J,
se cambiano il bordo sinistro o la prescrizione di misura. F è in generale
una funzione reale di due variabili, NON una funzione olomorfa di omega.

**Significato:** una perturbazione esterna può cambiare E anche a Lambda3
fissa; tuttavia in questa famiglia la risposta di E passa interamente per
la frequenza esatta già inserita. Non costituisce un nuovo dato indipendente
da quella frequenza e dal problema interno. Non esclude utilità pratica di
un indicatore più economico costruito con una frequenza approssimata.

## 3. Risposta locale esplicita e riduzione a due componenti reali

In assenza di zeri di psi nella regione usata, z=u'/u è olomorfa in omega.
Derivando la Riccati z'+z^2+omega^2-V0=0 rispetto a omega:

    k'+2zk=-2omega, k(a)=-i,
    k(x)=psi(x)^(-2)[-i psi(a)^2-2omega integral_a^x psi(t)^2 dt].

Ne segue h=partial_eta z=d k in J. Per eps fissato, ampiezza reale
A=|psi| e p=Im z:

    Q_M=-eps^2 A''/A=eps^2[Re(omega^2-V0)-p^2],
    partial_eta Q_M=2eps^2[Re(omega d)-p Im(k d)].

Ponendo d=d_R+i d_I, la risposta è una combinazione di sole due funzioni
reali di x, con coefficienti d_R,d_I. Per un insieme di perturbazioni
esterne, la matrice delle risposte linearizzate dei diagnostici interni ha
rango reale al più due. È una conseguenza di unicità/differenziazione;
non è un limite sul numero di gradi di libertà del campo completo.
La matrice di cambiamenti FINITI contiene termini quadratici e non deve
avere rango esattamente due. Un test deve estrarre la derivata a eta=0.

## 4. Derivata del rapporto integrale

Per il benchmark eps=1, J=[-0.8,0.8], w=exp(-x^2):

    U=integral_J w |Q_M| dx,
    D=integral_J w [|omega|^2+|V0|+p^2] dx, E=U/D,
    dot U=integral_J w sign(Q_M) dot Q_M dx,
    dot D=integral_J w [2Re(conj(omega)d)+2p Im(kd)] dx,
    dot E=(dot U D-U dot D)/D^2.

Richiediamo D>0, psi priva di zeri su J e regolarità uniforme locale in eta.
Se Q_M ha zeri isolati (o un insieme di misura nulla), la formula per dot U
segue dalla derivazione quasi ovunque e convergenza dominata. La derivata
ha salti in x ai cambi di segno: quadrature lisce globali possono convergere
male. Con Q_M nullo su un intervallo la differenziabilità va riesaminata.

La mediana pesata NON segue questa formula. Può non essere differenziabile
in presenza di plateaux, quantili non unici o cambi dell'ordinamento discreto.
Non presentare la derivata dell'integrale come prova per la mediana.

## 5. Controlli eseguiti e problema numerico ancora aperto

`verify_core_response.wl`, WolframKernel con licenza: 4/4 identità algebriche
True (Wronskiano, k, Q_M, sua prima variazione). Non certifica le ipotesi
analitiche né la novità dei teoremi. Il primo controllo di Q_M usava confronto
strutturale di espressioni equivalenti; è stato corretto a Simplify della
differenza. Nessuna modifica della formula fisica è stata necessaria.

`compact_indicator_linear_response.py` usa la risonanza precedente e calcola
la derivata da quadrature del profilo non perturbato, confrontandola con
soluzioni a eta diverso. Non fitta la formula analitica sui risultati.

| punti di quadratura | d omega/d eta | dE/d eta, Simpson |
|---|---|---|
| 8001 | 0.01036181821255 - 0.17599916946698 i | 0.00189086052182 |
| 16001 | 0.01036181821249 - 0.17599916946700 i | 0.00189891597344 |
| 32001 | 0.01036181821249 - 0.17599916946700 i | 0.00188899954578 |

La norma spettrale converge, ma la derivata integrale non è ancora convergente
alla precisione desiderata con questa quadratura globale. Una prima versione
spostava anche lievemente l'estremo sinistro selezionando nodi >=-0.8; ora
gli estremi della finestra sono inclusi esattamente. I numeri sopra sono quelli
della versione corretta. Il limite della quadratura resta visibile.

Usando ESATTAMENTE la quadratura trapezoidale a2001 punti del diagnostico
originale, la derivata analitica del funzionale discreto è0.00187362426521.
Le differenze centrali del medesimo funzionale danno:

| passo eta | derivata discreta FD | disaccordo relativo |
|---|---|---|
| 1e-3 | 0.00188825938837 | 7.81e-3 |
| 3e-4 | 0.00187364631662 | 1.18e-5 |
| 1e-4 | 0.00187364561741 | 1.14e-5 |

L'accordo circa1e-5 è un controllo utile, NON la convergenza al funzionale
continuo. Il plateau residuo va diagnosticato, non attribuito a priori a
una sola fonte. Possibili contributi: quadratura di k, tolleranza della
risonanza e non regolarità introdotta dal valore assoluto. I cambi di segno
vanno localizzati prima di richiedere ulteriori cifre a Wolfram.

## Conclusione per il core

Ora abbiamo una catena analitica: getto invariato non implica spettro
invariato; nella famiglia esterna il diagnostico interno fattorizza attraverso
omega; la sua risposta si calcola e si può falsificare con test indipendenti.
Questo dà un nucleo metodologico coerente, non ancora una novità pubblicabile.
Non prova una nuova previsione fisica su Kerr/Vaidya e non promette un
predittore migliore di WKB3. Il protocollo per Claude separa questi livelli.

Comandi:

    python3.13 calculations/compact_indicator_linear_response.py
    /Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_core_response.wl
