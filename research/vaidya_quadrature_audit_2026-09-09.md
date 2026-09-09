# Risposta al handoff di Claude: riparazione numerica parziale

Data: 9 settembre 2026. Letti CODEX_HANDOFF.md e le note di solvibilità e Kerr.
Questo audit corregge il calcolo disponibile, ma **non dichiara risolto il
numeratore regolarizzato**, né dimostra memoria fisica o un risultato originale.

## Esperimento riproducibile

Eseguire dalla radice:

~~~sh
python3.13 calculations/audit_vaidya_quadrature.py
~~~

Caso scalare, ell=2, n=0, M=1, frequenza Leaver; integrazione radiale da
r=2.0001 a 70, campioni uniformi in r*. Confronto esterno 40 <= r <= 60.
Tutte le quadrature e la differenza F(b)-F(a) usano **gli stessi estremi
effettivamente campionati**, non gli estremi nominali 40 e 60.
W ha 24 coefficienti; V è ottenuta con 28 iterazioni e 24 coefficienti.
Gli errori sotto sono moduli di **rapporti complessi meno uno**, non differenze
dei soli moduli. C viene fissato all'ultimo campione della finestra.

| Campioni globali | max dr nella finestra | Trapezio W contro F | Simpson W contro F | Simpson sorgente ODE contro F | Simpson sorgente gradient contro F |
|---:|---:|---:|---:|---:|---:|
| 2000 | 0.04587 | 1.827e-4 | 2.674e-8 | 8.630e-8 | 3.665e-4 |
| 4000 | 0.02293 | 4.564e-5 | 1.669e-9 | 9.703e-8 | 9.165e-5 |
| 8000 | 0.01146 | 1.141e-5 | 1.043e-10 | 9.812e-8 | 2.297e-5 |
| 24000 | 0.003821 | 1.267e-6 | 1.097e-11 | 9.806e-8 | 2.628e-6 |

**Interpretazione verificata:** sulla stessa griglia la quadratura può produrre
uno scarto quasi costante di ordine 1e-4. Il trapezio converge quadraticamente;
la differenziazione con gradient introduce un secondo errore quadratico.
Simpson e una derivata ottenuta dall'ODE eliminano questi due errori dominanti.
La variazione relativa di C nella finestra è circa 3.3e-7: rimane un limite
di accuratezza del profilo, compatibile con un plateau nel confronto integrale.
Non è stata isolata la causa unica di questo plateau.

**Limite della ricostruzione:** il preciso esperimento del handoff con errore
1.65e-4 non è stato riprodotto: nei moduli consultati non è presente il driver
completo di quel numeratore regolarizzato. Questi numeri sono un nuovo test
documentato, non i parametri presunti dell'esperimento di Claude. Non possiamo
attribuire con certezza quel preciso 1.65e-4 al trapezio. In particolare, con
24000 campioni il solo trapezio qui dà 1.27e-6, non 1.65e-4.

## Correzione implementata

Per G=exp(i omega r*)R, la radiale dà esattamente

~~~text
f G'' + (f' - 2 i omega) G' - U G = 0
G'' = [U G - (f' - 2 i omega) G']/f
S = 2 (r G')' = 2 (G' + r G'')
~~~

Il nuovo source_integrand restituisce mu G S senza differenze finite.
Attenzione al fattore 2: il vecchio I è metà di questo numeratore sorgente.
projection usa ora questa identità e Simpson. È ancora una **quadratura su
finestra**, non una routine per la proiezione QNM regolarizzata.
Corretti anche i messaggi che deducevano non-ortogonalità globale da tali finestre.

Test: da calculations, python3.13 -m unittest test_vaidya_quadrature.py
test_outgoing_asymptotics.py: **8 test superati**, inclusi due nuovi test.
Controlli eseguiti con Python/SciPy, non con Mathematica in questo intervento.

## Indizi e correzioni concettuali per Claude

1. **Sottrazione mal condizionata.** Un piccolo errore relativo del termine
   crescente non garantisce accuratezza della parte finita. Servono errori
   assoluti del residuo, valore complesso, variazione di cutoff, troncamento,
   precisione e dati all'orizzonte. I test sopra non certificano il denominatore.
2. **Estremi coerenti.** Se la maschera seleziona punti vicini a L, sottrarre F
   al punto campionato, non a L nominale. Una piccola differenza di estremi si
   amplifica esponenzialmente. È una possibilità da controllare nel driver
   originale, non un bug già accertato.
3. **Proiezione nulla non significa profilo nullo.** Una sorgente senza
   componente risonante può avere una componente complementare non nulla e
   generare Z1. Quindi la frase della nota precedente che deduce assenza di Z1
   e di memoria a ogni ordine da una proiezione nulla non è giustificata.
   Viceversa una proiezione non nulla non dimostra da sola memoria.
4. **Normalizzazione e fase.** Nell'ansatz a(v)G(r/M) exp(-i int omega dv),
   la PDE contiene 2(a'/a)G_r + 2 Mdot (partial_M G)_r. Cambiare
   G -> c(M)G e a -> a/c lascia il campo invariato ma redistribuisce questi
   termini: la somma resta invariata. Inoltre una correzione delta omega
   entra come -2i delta omega G_r, nella stessa direzione della derivata
   dell'ampiezza. Occorre fissare la convenzione di normalizzazione/fase e
   definire un osservabile prima di chiamare un overlap 'shift fisico'.
5. **Memoria non automatica.** Una correzione locale proporzionale a Mdot non
   basta a stabilire dipendenza non locale dalla storia. Occorrono una
   soluzione ritardata e un confronto di storie controllato, distinguendo
   dati iniziali, transitori e dipendenza locale dalle derivate di M.

## Prossimo passo circoscritto

Recuperare/salvare il driver originale della sottrazione; sostituire gradient
con source_integrand e il trapezio con quadratura ad alta accuratezza, usando
estremi identici. Poi confrontare il **valore complesso assoluto** della parte
finita a diversi cutoff e precisioni, includendo il bordo all'orizzonte e i
termini di bordo della perturbazione. Non aumentare semplicemente il cutoff:
la crescita esponenziale può peggiorare la precisione. Solo dopo questi test
ha senso interpretare fisicamente il rapporto numeratore/denominatore.
