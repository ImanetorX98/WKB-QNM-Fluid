# Kerr e Vaidya: audit del nucleo dell'articolo

11 settembre 2026. Controllo richiesto dopo il consolidamento del condizionamento.
Non sono stati modificati i solutori di produzione né le tesi del manoscritto.

## Verdetto

**Kerr:** l'indipendenza dal vecchio §5 non basta a certificare il §9.
È emersa una contaminazione di ordine 1/L proprio nella misura che sostiene
il 'termine geometrico di ordine uno'. Questo risultato centrale va sospeso.
Il solutore WKB3 autoconsistente e la struttura dell'equazione separata non
sono invalidati da questa diagnosi, ma non sono stati nuovamente certificati
qui. Il condizionamento di Q_M su Kerr resta da ricontrollare.

**Vaidya:** il numeratore fondamentale regolarizzato è riprodotto e i test
passano. Il nucleo tecnico sopravvive al problema del §5; l'interpretazione
come nuova dinamica o memoria non è ancora stabilita.

## Kerr: il rapporto mu non è fisso nel test originale

calculations/kerr_eikonal_order_test.py, eikonal_fit, impone

~~~python
L = ell + 0.5
m = int(round(mu * L))
~~~

Poi interpreta A/L^2 come sviluppo a mu fissato. In realtà
mu_eff=mu+delta_m/L, delta_m=m-mu L. Di conseguenza, anche se il vero
sviluppo a mu costante fosse A0(mu)+A2(mu)/L^2, il dato campionato conterrebbe

~~~text
A/L^2 = A0(mu) + delta_m partial_mu A0(mu)/L + O(L^-2).
~~~

Per mu=0.5 e gli ell pari del test, delta_m=-0.25 costante: la contaminazione
ha una pendenza pulita pari a 1. Non produce rumore o fit scadenti.
Il controllo sferico non la vede perché a c=0 l'autovalore non dipende da m.
Il test di una singola ridefinizione di L non la esclude, perché partial_mu A0
dipende dai parametri. Anche il test di c reale non la elimina.

### Tre controlli riproducibili

~~~sh
python3.13 calculations/audit_kerr_fixed_mu.py
~~~

Con chat=c/L=0.4 reale:

| Test | A1 stimato |
|---|---:|
| Vecchio mu nominale 0.5, ell=40,60,80,120,160,240 | -0.01928806498 |
| Interpolazione cubica in m al medesimo mu=0.5 | +2.6341e-7 |
| Termine previsto dal rounding: -0.25 partial_mu A0 | -0.01929739954 |
| mu=2/3 esattamente, ell=40,61,82,121,160,241 | -3.0502e-8 |

L'interpolazione è un diagnostico, non una soluzione fisica con m non intero.
Il controllo mu=2/3 non interpola: m è intero e mu esattamente costante per
tutti gli ell scelti. La previsione tramite derivata è un'approssimazione
numerica che riproduce il vecchio A1 a circa 5e-4 relativo: non un'identità
esatta a tutte le cifre.

Questi dati rendono non sostenibile citare il precedente -0.0193 come prova
di un A1 geometrico a mu fissato. Non dimostrano invece assenza di ogni ordine
1/L per ogni settore o traiettoria QNM complessa: la dipendenza di chat da L
può introdurre termini di quel tipo e va separata dall'effetto qui individuato.

### Anche la misura radiale è esposta

In kerr_radial_order_profile.py sia angular_leading sia exact_scaled_potential
arrotondano m. Il potenziale di testa usa mu nominale, quello finito usa m/L:
il termine delta_m partial_mu Q0/L entra nel confronto. Inoltre A0 è stimato
a ell_reference=400 e contiene un errore finito, non è il limite estrapolato.
Non abbiamo ancora ricalcolato tutta la tabella radiale: prima occorre
una sequenza coerente a mu fissato e un limite di testa convergente.

**Conseguenza editoriale:** sospendere la curva Kerr di Figura 1, il confronto
'Kerr bosonico ordine uno come Dirac' e l'enunciato generale che vi si appoggia.
Non cancellare i dati storici: indicare la contaminazione e rifare la misura.
Il testo A=A0+A1/L+... del §9.1 è inoltre notazionalmente incoerente con il
codice, che espande A/L^2: correggere la normalizzazione nella riscrittura.

## Vaidya: ciò che è stato ricontrollato

Eseguito:

~~~sh
python3.13 -m unittest calculations.test_vaidya_numerator_factored
~~~

**12 test superati**: riferimento numerico, taglio esterno, raccordo interno,
punto di raccordo, modi di spin diversi, Frobenius e controlli della
fattorizzazione. Per ell=2,s=0,n=0, con cutoff 40,50,60,70,80:

~~~text
N = 20.6665453911 - 40.3265380859 i
spread = 1.6023e-7
~~~

È un valore nel fissato schema di normalizzazione/regolarizzazione del codice,
non un osservabile già indipendente da tali scelte. Non abbiamo rivalidato
indipendentemente il denominatore o una derivata della frequenza.

La struttura locale I~x^(-4i omega)B(x), B analitica, spiega l'integrabilità
del fondamentale quando 4 Im omega > -1. È indipendente dall'argomento di
condizionamento di Schwarzschild. La riduzione del residuo adiabatico a
2 partial_r partial_M Z resta una derivazione della PDE, non dipende dal §5.

### Quello che ancora non segue

N diverso da zero non dimostra una componente complementare del profilo
Z1 diversa da zero, né memoria osservabile. La componente risonante deve
essere assorbita tramite la condizione di solvibilità, distinguendo trasporto
di ampiezza/fase e frequenza. La normalizzazione G -> c(M)G cambia la sorgente
partial_M G; il campo completo, con ampiezza trasformata inversamente, resta
invariato. Per questo il numeratore isolato non è ancora il risultato fisico.

Passi necessari: fissare le convenzioni, validare il rapporto/proiezione
completo con i bordi appropriati, costruire la risposta non adiabatica e
confrontarla con evoluzione nel tempo e limite a massa costante. Partire dal
fondamentale: aumentare gli overtoni prima di chiarire l'osservabile non
risolve il problema scientifico.

## Originalità: distinzione dalla validità numerica

Non è stata effettuata una nuova ricerca esaustiva di priorità. I riferimenti
primari pertinenti confermano che gli ambiti generali sono già studiati:
[Yang et al. 2012, geometria e QNM eikonali Kerr](https://arxiv.org/abs/1207.4253)
e [Capuano et al. 2024, perturbazioni Vaidya a tasso di massa costante](https://arxiv.org/abs/2407.06009).
Non attribuiamo loro il difetto del nostro codice né una specifica formula A1
senza verifica della convenzione e della traiettoria di scala.

Per un possibile nucleo pubblicabile, Vaidya oggi offre un calcolo tecnico
più solido della classificazione Kerr proposta, ma richiede ancora il ponte
verso una quantità fisica o un vantaggio metodologico confrontabile con i
precedenti. Una raccolta di correzioni di bug non basta come novità fisica.
