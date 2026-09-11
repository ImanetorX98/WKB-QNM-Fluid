# Verifica diretta della correzione forzata e del bordo uscente

11 settembre 2026. Seguito di vaidya_denominator_progress_2026-09-11.md.
Caso: ell=2, s=0, n=0, M=1. Questa verifica riguarda il problema radiale
congelato forzato nella prescrizione asintotica indicata, non la PDE Vaidya
completa con dati iniziali fisici.

## 1. Equazione integrata indipendentemente

Per evitare di confermare il numeratore soltanto con un'altra quadratura,
integriamo direttamente, insieme al modo G,

~~~text
L G = 0,
L H = S - 2 K G_r,
S = 2(r G_r)_r,
L = f partial_r^2 + (f_r-2i omega) partial_r - U.
~~~

H è la correzione del profilo per unità di Mdot, in una convenzione con
coefficiente di trasporto K. La scelta H(2)=0 fissa l'addizione del modo
omogeneo G. La soluzione H è analitica all'orizzonte: i dati iniziali sono
ottenuti da una ricorrenza di Frobenius **inhomogenea**, non imponendo H=H'=0
a distanza finita dal bordo.

La ricorrenza usa 60 coefficienti e viene valutata a r=2.5. Il programma
integra G,G_r,H,H_r con DOP853. Non integra il numeratore come variabile di
stato e non usa l'integrale cumulativo per costruire H.

## 2. Identità di bordo

Introduciamo il Wronskiano bilineare pesato

~~~text
J = mu f (G H_r - G_r H), mu=exp(-2i omega r_*).
~~~

Le due ODE danno esattamente

~~~text
J_r = mu G (S-2K G_r).
~~~

L'identità è verificata simbolicamente con Mathematica. Il ramo H analitico
ha J -> 0 all'orizzonte nel caso fondamentale in esame.
Nella regione esterna la primitiva particolare uscente è

~~~text
J_out = F_N - 2K F_P,
~~~

con le primitive asintotiche già definite per N e P. Ne segue che la
costante residua J-J_out è formalmente N-2KP, nella medesima prescrizione
di parte finita. Il test confronta questa previsione con il Wronskiano
calcolato dalla soluzione H dell'ODE, non con la quadratura di N.

Per G~C exp(2i omega r*) e una componente omogenea entrante della correzione
H_in~B, si ha J_in -> -2i omega C B. La costante residua è dunque un
diagnostico del coefficiente entrante rispetto al particolare uscente scelto.
L'aggiunta di un multiplo di G a H non cambia J e non può cancellarla.
La separazione di componenti asintotiche presuppone la scelta del ramo e
della continuazione analitica: non è una dimostrazione globale sulla PDE
non stazionaria o sulle ambiguità oltre tutti gli ordini delle serie.

## 3. Risultato a K=N/(2P)

Usiamo il candidato precedente

~~~text
K0 = 1.1642938693836598 + 4.348900226776741 i.
~~~

Residuo assoluto |J-J_out| con rtol=1e-12:

| r | Residuo |
|---:|---:|
| 30 | 3.94e-7 |
| 35 | 6.10e-7 |
| 40 | 2.42e-6 |
| 45 | 8.91e-6 |

Il residuo cresce all'esterno per cancellazione fra quantità esponenzialmente
grandi. Non misuriamo un errore relativo rispetto a zero e non dichiariamo
convergenza uniforme fino all'infinito.

### Controllo di tolleranza

A r=40:

| rtol richiesta | Residuo assoluto |
|---:|---:|
| 1e-10 | 2.70e-4 |
| 1e-12 | 2.42e-6 |
| 2e-14 | 1.45e-7 |

Nell'ultima riga SciPy porta automaticamente rtol al minimo macchina
2.220446049250313e-14 e avvisa: non si tratta di precisione arbitraria.
A r=30 il residuo non migliora monotonamente al diminuire della tolleranza;
vi contribuiscono anche troncamento asintotico, normalizzazione e accuratezza
del valore K0. Non è corretto attribuire tutti gli errori al solo integratore.

## 4. Test negativo: spostare K deve introdurre un residuo

Ripetiamo l'integrazione di H cambiando K, senza modificare il valore previsto
del residuo per accomodare la soluzione. La previsione è

~~~text
J-J_out = -2 deltaK P.
~~~

| deltaK | Residuo complesso previsto |
|---|---|
| 0.01 | 0.0746547806 + 0.0675079915 i |
| 0.01 i | -0.0675079915 + 0.0746547806 i |

I test verificano entrambe le componenti con errore assoluto inferiore a
1e-5 ai punti r=30,35,40, con rtol=1e-12. La scala del segnale introdotto è
circa 0.1, nettamente superiore al residuo numerico. Questo controllo
distingue una vera condizione selettiva da un azzeramento numerico automatico.

## 5. Limite fisico importante: espansione non uniforme a grande r

Se scriviamo H=G y, l'equazione esterna dominante diventa, con k=2i omega,

~~~text
y'' + k y' ~ 2 k^2 r + O(1),
y ~ k r^2 + termini subdominanti.
~~~

Quindi Mdot H/G cresce come Mdot k r^2 nel regime formale esterno. Anche
con coefficiente entrante nullo, l'espansione G+Mdot H non è uniformemente
piccola quando r -> infinito a Mdot fissato. I limiti adiabatico e di grande
distanza non possono essere scambiati senza un raccordo o una risommazione.
Il programma riporta |H/G| ai punti di controllo per rendere visibile il problema.

Questa osservazione non dimostra una nuova memoria. Indica che il prossimo
confronto deve fissare osservatore, coordinata temporale e scala di variazione
di M, oppure costruire un'espansione raccordata nel settore lontano.

## 6. Stato della verifica dei bordi

**Verificato:** entro il problema radiale congelato forzato e la prescrizione
di particolare uscente impiegata, N/(2P) seleziona la soluzione senza costante
Wronskiana residua, al livello di errore misurato. La sottrazione del numeratore
ha quindi un controllo indipendente tramite l'equazione della correzione.

**Ancora aperto:** dimostrare che questa prescrizione coincida con la
condizione uscente del problema Vaidya a massa variabile, trattando la regione
lontana e i transitori; validare la forma d'onda e la combinazione Xi.
Non sono stati calcolati nuovi overtoni o un'evoluzione nel tempo.

## 7. Riproduzione

~~~sh
python3.13 calculations/vaidya_forced_boundary.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt -script calculations/verify_forced_boundary.wl
~~~

Da calculations:

~~~sh
python3.13 -m unittest test_vaidya_forced_boundary.py
~~~

**3 test superati**: cancellazione al K candidato, due perturbazioni complesse
di controllo, miglioramento con tolleranza. Un'identità Mathematica confermata.
I tre test contengono più punti di misura; non è stata rieseguita l'intera suite.
