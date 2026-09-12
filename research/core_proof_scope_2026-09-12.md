# Fin dove arriva una dimostrazione del core — 12 settembre 2026

## Risposta breve

Non posso dimostrare il core forte attuale perché contiene affermazioni false
o più forti dell'evidenza. Posso dimostrare un limite preciso dei diagnostici
locali e verificare un controesempio on-shell alla generalizzazione della
ridondanza. La sua base teorica è già nota: non chiamarlo nuovo teorema fisico.
Per ora non è dimostrato un nucleo originale sufficiente alla pubblicazione.

## Proposizione: il getto al picco non determina lo spettro QNM

Considerare psi''+[omega^2-V_eta(x)]psi=0, V_eta=V0+eta b, con V0 e b reali,
lisci, a supporto compatto in (a,bound). Il supporto della perturbazione è
disgiunto da un intorno del massimo dominante x_p. Supporre una risonanza
semplice omega0 non nulla, un ramo differenziabile e il denominatore sotto
non nullo. I bordi sono esattamente uscenti:

    psi'(a)=-i omega psi(a), psi'(bound)=+i omega psi(bound).

Per ogni k, V_eta^(k)(x_p)=V0^(k)(x_p). Quindi qualsiasi formula WKB locale
di ordine finito valutata al picco, inclusa Lambda3 a overtone fissato, resta
invariata (finché si segue lo stesso massimo). Invece

    d omega/d eta = integral b(x) psi(x)^2 dx /
      [2 omega integral psi(x)^2 dx + i(psi(a)^2+psi(bound)^2)].

Gli integrali sono da a a bound, senza coniugazione complessa.

### Dimostrazione

Posto phi=partial_eta psi e W=psi phi'-psi' phi, differenziando l'ODE:

    W'=(b-2 omega omega_eta)psi^2.

Differenziando anche le condizioni uscenti:

    W(bound)=i omega_eta psi(bound)^2,
    W(a)=-i omega_eta psi(a)^2.

Integrando si ottiene la formula. Non è lecito omettere i termini al bordo:
sono proprio la parte superficiale della norma generalizzata QNM.

Scegliere un punto esterno al picco con psi non nulla e una bump reale
non negativa sufficientemente stretta attorno a quel punto. Per continuità,
l'integrale è vicino a psi(x_c)^2 integral b, dunque non nullo per una
bump abbastanza stretta. Ne segue che esistono perturbazioni che lasciano
immutato il getto ma spostano la risonanza. La prova richiede le ipotesi di
regolarità del ramo/simplezza, non vale indiscriminatamente ai punti eccezionali.

**Attribuzione:** è una conseguenza della perturbazione logaritmica dei QNM,
non una formula originale. Letta direttamente la §2.1, eq.2.6–2.16 di
[Leung et al. 1998](https://arxiv.org/abs/physics/9712037), disponibile localmente
in `Schw-QNM-WKB-Fluid/papers/leung1998.txt` e nel PDF associato.

## Test nuovo nel repository, con frequenza determinata ai due bordi

Definire B(x)=exp(1-1/(1-x^2)) per |x|<1, zero altrimenti, e

    V_eta(x)=10 B(x)+eta B((x-2)/0.4).

È C-infinito ma non analitico agli estremi dei supporti. La perturbazione è
zero in un intero intorno del picco x=0, non solo fino all'ordine sei.
Per i piccoli eta esaminati il massimo dominante resta quello centrale.
Questa è una famiglia di barriere modello, NON una soluzione di Einstein
né una dimostrazione sulla sola famiglia Schwarzschild del manoscritto.

Solutore indipendente: shooting DOP853 da -1 a 2.4, psi(-1)=1,
psi'(-1)=-i omega; si azzera psi'(2.4)-i omega psi(2.4).
Le condizioni ai bordi sono esatte perché il potenziale è nullo fuori dal
dominio; non vi sono code troncate. Si segue una risonanza per continuazione,
senza pretendere di catalogare tutto lo spettro o dimostrarla fondamentale.

Diagnostico integrale di tipo §6, eps=1, w=exp(-x^2), finestra [-0.8,0.8]:

    E = integral w |Q_M| / integral w (|omega|^2+|V|+(Im z)^2),
    z=psi'/psi, Q_M=Re(omega^2-V)-(Im z)^2.

La formula di Q_M usa l'ODE, non differenze finite di |psi|. Questa E non
è la mediana pesata del §6(iv): non confondere i due stimatori.

| eta | omega | E |
|---|---|---|
| -0.01 | 3.4231972718943964 - 0.605642918496682 i | 0.0608365709632 |
| -0.001 | 3.4232851815842094 - 0.6072117476227428 i | 0.0608506237872 |
| 0 | 3.4232954881855626 - 0.6073875920723604 i | 0.0608524987777 |
| 0.001 | 3.423305905508026 - 0.6075637463353818 i | 0.0608544002675 |
| 0.01 | 3.423404782399676 - 0.6091632508703276 i | 0.0608745678288 |

Il residuo uscente massimo del campione è 5.95e-14: residuo del problema
discretizzato, non certificazione delle cifre dell'autovalore.
La formula perturbativa dà omega_eta=0.01036180412-0.17599918163i;
la differenza centrale con passo1e-4 dà 0.01036182635-0.17599918370i,
scarto relativo1.27e-7. La quadratura della norma usa20001 punti ed è
potenzialmente soggetta a cancellazione: non assumere precisione arbitraria.

Stringendo rtol da2e-11 a2e-12, a eta0/0.01 le differenze in omega sono
1.8e-12/1.45e-10, e quelle in E sono2.0e-15/2.1e-11. Il cambiamento di E
fra0 e0.01 è2.21e-5, largamente risolto rispetto a questi controlli.
Resta da studiare anche convergenza della quadratura dell'indicatore e un
secondo solutore spettrale: non è un calcolo a intervalli certificato.

### Che cosa dimostra e che cosa no

- La proposizione prova che il dato locale non determina generalmente la
  risonanza. Il test mostra anche una risposta dell'indicatore integrale.
- Non dimostra che l'indicatore predica bene l'errore WKB. Non confronta
  ancora le prestazioni con predittori più economici su un campione ampio.
- Non smentisce una correlazione approssimata ristretta a Schwarzschild,
  né prova un risultato per la mediana pesata.
- Smentisce l'inferenza universale "Madelung è un'identità, dunque il suo
  funzionale esatto contiene soltanto Lambda3". WKB troncata e soluzione
  esatta sono oggetti diversi.
- Un effetto remoto è compatibile con il carattere asintotico/locale WKB:
  questo non prova l'invalidità della WKB nel proprio dominio di applicabilità.

## Letture necessarie, senza download duplicati

Sono già presenti PDF e testi estratti per Leung1998, Yang1207.4253,
Capuano2407.06009 e Yoo2510.25062. Non serve scaricarli nuovamente per iniziare.
Per il confronto di priorità serve leggere e confrontare formule, non solo
gli abstract. In questo turno ho letto la derivazione pertinente di Leung;
NON dichiaro completata una lettura integrale dei quattro lavori.

Priorità: Leung §2–4 per norma, bordi e sensibilità; Yang §II e App.A per
il risultato Kerr; Capuano per il benchmark Vaidya a tasso costante; Yoo per
definire frequenza osservata e confronto temporale. DDP1997, anch'esso locale,
serve a delimitare ciò che è già WKB esatta. Verificare versione/data dei
PDF locali prima di un confronto editoriale finale.

## Comandi

    python3.13 calculations/compact_barrier_core_test.py

Da calculations:

    python3.13 -m unittest test_compact_barrier_core.py

Eseguito: **2 test superati**; lo script numerico completo termina con exit0.
Questo turno usa SciPy; non è stato eseguito un nuovo controllo Mathematica
del problema spettrale. Il secondo metodo è proposto nell'incarico a Claude.

## Decisione sul core

Il vecchio core non va salvato mediante una dimostrazione forzata. L'ipotesi
di lavoro sostenibile è: stabilire limiti e utilità residua dei diagnostici
d'ampiezza con validazione fuori campione. Il test qui chiude la lacuna
"solo dati iniziali" del precedente audit, ma non chiude originalità o
pubblicabilità. Per una nuova previsione fisica Vaidya serve ancora la PDE.
