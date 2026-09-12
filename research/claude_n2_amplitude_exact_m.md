# N2 — l'audit del §9.6 su \(m\) intero esatto

**Esito: la conclusione qualitativa regge, i numeri no** — e il limite non è
quello che credevo. Non è la frequenza: è la **condizione al bordo**.

## 1. La cura non richiedeva codice nuovo

Con \(\mu=2/3\) e \(\ell\equiv1\pmod3\) si ha \(m=(2\ell+1)/3\) **intero esatto**,
e `int(round(...))` non arrotonda nulla. Anche il riferimento angolare di
default, \(\ell=400\), soddisfa la condizione — \(400\equiv1\pmod3\).

Il difetto del §9.3 era quindi evitabile **scegliendo i parametri**, non
riscrivendo i moduli. Lo script verifica l'esattezza con un `assert` invece di
assumerla: scarto \(0.0\times10^{0}\) su tutti gli \(\ell\) usati.

## 2. La tabella corretta

\(\mu=2/3\), \(m\) esatto, derivate analitiche. Errore relativo mediano su
\(20<r<50\):

| \(a\) | \(\ell\) | \(\omega\) eikonale | \(\omega\) autoconsistente |
|---|---|---|---|
| 0.0 | 40 | 2.3e−3 | 7.3e−4 |
| 0.0 | 100 | 2.7e−3 | 4.5e−3 |
| 0.3 | 40 | 3.4e−2 | 8.1e−4 |
| 0.3 | 100 | 1.3e−1 | 5.0e−3 |
| 0.6 | 40 | 7.4e−2 | 1.6e−3 |
| 0.6 | 100 | 2.8e−1 | 9.8e−3 |
| 0.9 | 40 | 6.9e−1 | 3.2e−2 |
| 0.9 | 100 | 9.7e−1 | 2.0e−1 |

**Il confronto fra colonne regge**: la frequenza autoconsistente migliora di uno
o due ordini a ogni spin e multipolo.

**Ma i valori assoluti crescono con \(\ell\)**, dove dovrebbero calare. Nella
vecchia tabella calavano — ed era l'arrotondamento a farli sembrare tali.

## 3. Di che cosa è fatto quel limite

Regola di sempre: variare qualcosa che non deve contare. Qui l'offset \(\delta\)
con cui la condizione entrante è imposta a \(r=r_++\delta\), a \(a=0\),
frequenza autoconsistente:

| \(\ell\) | \(\delta=10^{-4}\) | \(10^{-5}\) | \(10^{-6}\) | \(10^{-7}\) |
|---|---|---|---|---|
| 40 | 1.2e−2 | 3.0e−3 | 7.3e−4 | 1.8e−4 |
| 61 | 2.8e−2 | 6.9e−3 | 1.7e−3 | 4.1e−4 |
| 100 | 7.7e−2 | 1.8e−2 | 4.5e−3 | 1.1e−3 |

Un fattore \(\simeq4\) per ogni decade di \(\delta\), a ogni \(\ell\).

> Le cifre del §9.6 misurano la **condizione al bordo troncata**, non
> l'accuratezza della previsione di \(Q_M\). E questo spiega anche perché
> crescono con \(\ell\): a multipolo più alto la soluzione è più sensibile al
> bordo.

## 4. E adesso è riparabile

N1 ha mostrato che il bordo interno di Kerr ha una rappresentazione
**convergente** in `HeunC`, con residuo della Riccati \(10^{-13}\) fino a
\(\rho=1\). Usarla come condizione iniziale al posto della troncatura
rimuoverebbe esattamente il termine che domina questa tabella.

Non l'ho fatto in questo passaggio: il §9.6 è una verifica di contorno — serve a
dire che nulla di patologico accade all'ampiezza — e il suo contenuto è il
confronto fra colonne, che non cambia. Ma è il seguito naturale, e il manoscritto
ora dichiara sia il limite sia la via d'uscita.

## 5. Aggiornato nel manoscritto

* tabella sostituita con quella a \(m\) esatto;
* aggiunto il riquadro che dichiara che cosa misura la colonna di destra, con la
  scala in \(\delta\);
* menzionata la forma di Heun confluente come rimedio non ancora applicato;
* tabella di riproducibilità aggiornata con i tre script.

## Comando

```sh
python3.13 calculations/claude_kerr_amplitude_exact_m.py
```

Circa 100 s: risolve la frequenza autoconsistente e integra il profilo per
dodici combinazioni.
