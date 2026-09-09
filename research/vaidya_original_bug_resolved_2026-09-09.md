# Risolto il preciso errore del riproduttore di Claude

9 settembre 2026. Questo risultato supera lo stato 'causa non accertata' dei
due audit precedenti. Il nuovo file repro_numerator_discrepancy.py permette
ora di riprodurre direttamente il numero storico, senza ipotizzarne i parametri.

## Causa verificata: una cella mancante

Il codice originale confrontava:

~~~python
analytic = antiderivative[j] - antiderivative[i]
numeric = simpson(integrand[i:j], x=r[i:j])
~~~

La slice Python i:j termina al campione j-1. La primitiva, invece, è valutata
al campione j. Quindi numeric integra da r[i] a r[j-1], analytic da r[i] a r[j].
Non è un errore dell'asintotica, né un errore della formula di Simpson: i
due calcoli hanno domini diversi. Trapezio e Simpson concordavano proprio
perché integravano entrambi lo stesso intervallo incompleto.

Due correzioni indipendenti confermano la diagnosi:

~~~python
# A: includere l'estremo j nella quadratura
numeric = simpson(integrand[i:j+1], x=r[i:j+1])
analytic = antiderivative[j] - antiderivative[i]

# B: mantenere la slice e correggere l'estremo della primitiva
numeric = simpson(integrand[i:j], x=r[i:j])
analytic = antiderivative[j-1] - antiderivative[i]
~~~

## Numeri ottenuti con il comando esatto

~~~sh
python3.13 calculations/repro_numerator_discrepancy.py
~~~

Configurazione originale: ell=2, s=0, n=0; 100001 punti in r*;
r iniziale 2.0001, r finale 90; normalizzazione C al campione di r circa 70.
Nessun cambio di frequenza, integratore, dati al bordo, normalizzazione,
derivate o troncamento è necessario per risolvere questo difetto.

| Finestra nominale | Modulo rapporto originale, Simpson | Errore complesso originale | Errore complesso corretto A | Errore complesso corretto B |
|---|---:|---:|---:|---:|
| 40–50 | 0.999833883 | 1.046253e-3 | 3.262983e-7 | 3.265347e-7 |
| 50–60 | 0.999835624 | 1.038754e-3 | 2.616845e-7 | 2.616671e-7 |
| 60–70 | 0.999836606 | 1.032944e-3 | 2.100179e-7 | 2.100175e-7 |
| 70–80 | 0.999837069 | 1.028569e-3 | 2.174845e-7 | 2.174852e-7 |

Errore complesso significa abs(rapporto-1), non abs(abs(rapporto)-1).
La larghezza della cella omessa è circa 0.001108–0.001125.
La quasi costanza dello scarto è coerente con una cella quasi costante nella
regione dove la crescita esponenziale domina.

## Secondo difetto del controllo: si stampava solo il modulo

La tabella della discrepanza usava abs(trap) e abs(simp), non il rapporto
complesso. Il controllo puntuale dell'integrando, invece, stampava davvero
entrambe le componenti. Ad esempio, sulla finestra 50–60:

~~~text
rapporto originale = 0.9998350983831241 - 0.001025580974801166 i
rapporto corretto  = 0.9999997537709034 + 0.000000088600242681 i
~~~

Quindi il famoso 1.65e-4 misurava solo lo scarto del **modulo** dall'unità;
nascondeva un errore di fase più grande. I controlli separati non erano in
contraddizione: il confronto integrale usava estremi diversi e la stampa
perdeva l'informazione di fase.

## Modifiche e verifica

Il riproduttore conserva la sezione storica per documentare il bug e aggiunge
una sezione finale con entrambe le correzioni e gli errori complessi. Il suo
docstring ora dichiara esplicitamente la diagnosi. Non cambiamo setup per
evitare di confondere questa riparazione con gli altri miglioramenti numerici.

Nuovi test in calculations/test_repro_numerator_discrepancy.py:

1. Riproduzione del numero originale e verifica di entrambe le correzioni
   sulle quattro finestre, soglia sull'errore complesso corretto 5e-7.
2. Controllo indipendente su un polinomio con primitiva esatta, senza QNM.

Da calculations:

~~~sh
python3.13 -m unittest test_repro_numerator_discrepancy.py test_vaidya_quadrature.py test_outgoing_asymptotics.py
~~~

**13 test superati.** Verifiche di questa diagnosi eseguite con Python/SciPy;
non richiedono nuove identità simboliche o ricerca bibliografica.

## Stato scientifico dopo la correzione

**Il preciso paradosso numerico del handoff è risolto.** Per attribuire la causa
di questo numero non occorrono ulteriori analisi né altro codice originale.
Restano i residui numerici di ordine 1e-7, studiati separatamente nei precedenti
audit, e soprattutto la verifica del numeratore regolarizzato completo con
entrambi i bordi. Risolvere questo bug non dimostra ancora memoria fisica,
uno shift osservabile, né originalità/pubblicabilità del progetto.

---

## Seguito Claude, stessa data: che cosa sblocca davvero la correzione

Verifica indipendente della diagnosi. La cella mancante spiega **tutto** il numero:

| grandezza, finestra 50–60 | valore |
|---|---|
| larghezza della cella \(r_j-r_{j-1}\) | 1.1153e-3 |
| \(\lvert\)cella\(\rvert\) | 6.6956e+4 |
| \(\lvert\)integrale\(\rvert\) | 6.4455e+7 |
| cella / totale | **1.0388e-3** |
| errore complesso originale | **1.0387e-3** |

Coincidono a quattro cifre: nessun residuo da spiegare. La cella pesa 1e-3 pur
essendo 1/11000 dell'intervallo perché l'integrale è **cancellativo**
(\(\int\lvert I\rvert = 2.45\times10^8\) contro \(\lvert\int I\rvert=6.4\times10^7\))
e l'integrando cresce esponenzialmente verso il bordo destro.

Nessun altro punto del codice usa lo stesso schema: `grep searchsorted` trova
solo `robust_indicator_test.py:141`, che è una mediana, non una quadratura. Il
difetto era confinato al diagnostico.

### Il numeratore regolarizzato con la quadratura corretta

\(N(L_+)=\int_4^{L_+} I\,dr - F(L_+)\), sottrazione della primitiva asintotica:

| \(L_+\) | \(\lvert F\rvert\) | \(\lvert N\rvert\) | \(\lvert N/F\rvert\) |
|---|---|---|---|
| 40 | 6.95e+5 | 44.82 | 6.4e-5 |
| 50 | 6.58e+6 | 44.51 | 6.8e-6 |
| 60 | 5.88e+7 | 67.48 | 1.1e-6 |
| 70 | 5.06e+8 | 132.2 | 2.6e-7 |
| 80 | 4.22e+9 | 1390 | 3.3e-7 |
| 85 | 1.21e+10 | 3982 | **3.3e-7** |

\(\lvert N/F\rvert\) **scende fino a un pavimento 3.3e-7 e lì si ferma**: è
esattamente l'accordo integrando-numerico contro serie (2–3e-7) misurato dal
riproduttore. Oltre \(L_+\simeq60\) si sta integrando quel pavimento, non il
modo.

**Non è tolleranza dell'ODE.** Con `rtol=1e-13, atol=1e-16` la tabella è
identica cifra per cifra. Il pavimento viene dalla serie asintotica.

**È in parte la costante.** Lo scarto residuo è quasi un fattore moltiplicativo
costante, perché \(C\) è fittata su \(u\) mentre l'integrando usa la serie \(W\),
troncata diversamente. Rifittando \(C^2\) ai minimi quadrati contro \(I\) su
\(45<r<75\) — correzione \(1-1.97\times10^{-7}+9.0\times10^{-8}i\) — il pavimento
scende a 1.1e-7 e il plateau si allunga:

| \(L_+\) | 35 | 40 | 50 | 60 | 70 | 80 |
|---|---|---|---|---|---|---|
| \(\operatorname{Re}N\) | 43.2 | 43.8 | 44.5 | 52.3 | 24.8 | 31.8 |
| \(\operatorname{Im}N\) | 8.0 | 8.9 | 11.2 | 17.3 | −31.9 | 476 |

### Stato onesto

Il numeratore **esiste ed è misurabile**: \(N\approx44\) sul plateau
\(35\le L_+\le50\). Ma converge al **10%**, non allo 0.2% del denominatore, e il
limite non è più un bug: è il condizionamento. Un errore relativo 1e-7 sulla
serie asintotica moltiplica una primitiva che cresce di quattro ordini fra
\(L_+=40\) e \(L_+=85\).

Per scendere sotto l'1% servirebbe aritmetica multiprecisione (`mpmath`) su ODE
*e* serie insieme — restringere una sola delle due non muove nulla, come mostra
il test a `rtol=1e-13`.

Vale la pena notare che il pavimento è la **stessa** patologia del §2.2 del
riassunto di sessione: un funzionale che amplifica l'errore su ciò che vorrebbe
misurare. Qui \(\sim10^{2}\) su \(Q_M\), lì \(\sim10^{4}\) sul numeratore.
