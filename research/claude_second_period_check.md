# Controllo H — verifica indipendente del secondo periodo

Verifica di `research/angular_A2_closed_2026-09-12.md`. **La formula regge su
tutti i controlli richiesti**, compreso quello che il brief segnala come
delicato: il moto degli estremi.

## 1. La catena algebrica, verificata simbolicamente

`calculations/claude_verify_second_period_chain.py`. Con \(p=\sqrt{Q_0}\), la
parte derivativa di \(p_2\) si semplifica in

\[
\frac{5Q_0'^2-4Q_0Q_0''}{32\,Q_0^{5/2}} .
\]

I due passi della riduzione sono **esatti**, con residuo identicamente nullo:

| passo | affermazione | coefficiente della derivata totale | residuo |
|---|---|---|---|
| 1 | parte derivativa \(\equiv -Q_0'^2/(32Q_0^{5/2})\) | \(-\tfrac18\,\partial_\theta[Q_0'/Q_0^{3/2}]\) | **0** |
| 2 | \(-Q_0'^2/(32Q_0^{5/2})\equiv-Q_0''/(48Q_0^{3/2})\) | \(+\tfrac1{48}\,\partial_\theta[Q_0'/Q_0^{3/2}]\) | **0** |
| 3 | \(-Q_0''/(48Q_0^{3/2})=\tfrac1{24}\partial_A[Q_0''/\sqrt{Q_0}]\) | — (usa \(\partial_AQ_0=1\)) | **True** |

La riduzione
\(I_2=\tfrac12\int Q_2/\sqrt{Q_0}+\tfrac1{24}\partial_A\!\int Q_0''/\sqrt{Q_0}\)
è quindi confermata al livello dell'algebra, indipendentemente
dall'implementazione.

## 2. Il moto degli estremi è indispensabile, e si misura

Il brief avverte: «la derivata in A include il moto degli estremi». Ecco la
differenza, a \(\mu=2/3\), \(\hat c=0.6\), \(A=A_0\):

**Derivando l'integrale intero**, con i turning point ricalcolati a ogni \(A\):

| passo | \(10^{-6}\) | \(10^{-7}\) | \(10^{-8}\) |
|---|---|---|---|
| \(\partial_A\!\int Q_0''/\sqrt{Q_0}\) | −16.66214049 | −16.66214060 | **−16.66214054** |

Stabile a otto cifre.

**Derivando sotto il segno**, cioè \(-\tfrac12\int Q_0''/Q_0^{3/2}\) con un
taglio \(\epsilon\) agli estremi:

| cutoff | \(10^{-2}\) | \(10^{-3}\) | \(10^{-4}\) |
|---|---|---|---|
| valore | 847.7 | 8625.6 | 86405.9 |

**Diverge come \(1/\epsilon\)**: un fattore dieci per ogni decade. Il cutoff
diventerebbe un parametro spurio del coefficiente, esattamente come avvertito.
Non è una sfumatura: le due strade differiscono di ordini di grandezza.

## 3. Il limite sferico

\(\hat c=0\) dà \(I_2=0.3926990816987241548078\) e \(A_2=-0.25\) **esatto**, a
precisione 30 e 45. Riproduce \(A=\ell(\ell+1)=L^2-\tfrac14\), cioè il termine di
Langer del §8 del manoscritto, non soltanto il segno.

## 4. Contro gli autovalori: **nessun parametro adattato**

Fissando \(A_0\) e \(A_2\) entrambi dal lato analitico — zero regressioni — e
usando `SpheroidalEigenvalue` a precisione 70 su \(\mu=2/3\), \(\hat c=3/5\),
\(m\) intero esatto:

| \(L\) | 40.5 | 67.5 | 97.5 | 127.5 | 160.5 |
|---|---|---|---|---|---|
| \(L^4(\hat A-A_0-A_2/L^2)\) | −0.049070 | −0.048900 | −0.048850 | −0.048831 | **−0.048821** |

Converge a \(A_4\simeq-0.04882\), con la deriva residua compatibile con
\(A_6/L^2\). È il test più severo possibile: se \(A_2\) fosse sbagliato anche
solo nella settima cifra, \(L^4\) lo amplificherebbe.

## 5. Fit con \(A_1\) **libero**

Come richiesto, \(A_0\) fissato dall'azione e \(A_1\) lasciato libero insieme ad
\(A_2\) e ai termini superiori, variando grado e intervallo:

| intervallo | grado | \(A_1\) | \(A_2\) | residuo |
|---|---|---|---|---|
| \(\ell\) 40–160 | 3 | \(-1.5\times10^{-7}\) | −0.2310308 | 2.4e−10 |
| \(\ell\) 40–160 | 4 | \(1.2\times10^{-9}\) | −0.2310655 | 6.3e−13 |
| \(\ell\) 40–160 | 5 | \(-2.6\times10^{-10}\) | −0.2310651 | 4.6e−14 |
| \(\ell\) 40–97 | 4 | \(2.6\times10^{-9}\) | −0.2310658 | 3.1e−13 |
| \(\ell\) 97–160 | 4 | \(\mathbf{6.4\times10^{-11}}\) | −0.2310652 | 3.9e−16 |

\(A_1\) **scende di quattro ordini** man mano che grado e intervallo migliorano,
senza mai stabilizzarsi su un valore non nullo: si comporta come un coefficiente
che è zero e sta assorbendo l'errore di modello.

\(A_2\) converge al valore in forma chiusa
\(-0.2310651878076874644\): lo scarto migliore è \(3\times10^{-8}\), ed è il
limite del fit, non un disaccordo.

## Esito

Tutti i controlli del brief passano. La formula del secondo periodo è verificata
per tre strade indipendenti — algebra simbolica, derivata con estremi mobili,
confronto con autovalori multiprecisione — e il limite sferico è esatto.

Resta quello che Codex dichiara: le cifre indicano stabilità numerica, non una
certificazione a intervalli, e non è un teorema uniforme del resto.

## Comandi

```sh
python3.13 calculations/claude_verify_second_period_chain.py
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt \
  -script calculations/claude_verify_second_period.wl
```
