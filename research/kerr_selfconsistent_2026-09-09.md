# Kerr: frequenza autoconsistente e chiusura del passo B

**Data:** 9 settembre 2026.
**Natura:** **misurato**, con validazione contro un risultato indipendente.
**Riproduzione:** `python3.13 calculations/kerr_wkb3_selfconsistent.py`;
test in `calculations/test_kerr_wkb3_selfconsistent.py` (6 test).

---

## 1. Il problema che chiude

Il confronto fra previsione analitica di \(Q_M\) (Appendice C del manoscritto) e
valore estratto dall'ampiezza integrata riusciva su Schwarzschild
(\(1.9\times10^{-4}\) a \(\ell=100\)) e falliva su Kerr, con errori fino a
\(2.3\times10^{-1}\) e non monotoni in \(\ell\).

La diagnosi in [kerr_madelung_sensitivity](kerr_madelung_sensitivity_2026-09-09.md)
aveva escluso risoluzione, quasi-zeri di \(\psi\) e troncamento dell'autovalore
al reale, e imputato la colpa alla **frequenza eikonale**, il cui errore per
\(a\neq0\) è \(\sim10^{-3}\) mentre \(Q_M\) amplifica di \(\sim10^2\).

## 2. Costruzione

A \(\omega\) fissata l'autovalore sferoidale è un **numero**, quindi

\[
q(r)=\left(\omega-\frac{ma}{H}\right)^2-\frac{\Delta\lambda}{H^2}-\frac{1}{h}\frac{d^2h}{dr_*^2}
\]

è esplicita in \(r\): le derivate rispetto a \(r_*\) si ottengono
simbolicamente iterando \(D_*=(\Delta/H)\,d/dr\), senza differenze finite di
sesto ordine. La condizione di barriera di Iyer–Will, scritta in termini di
\(q\) anziché di \(V\),

\[
\frac{q_0}{\sqrt{2q_0''}}=\Lambda_2-i\left(n+\tfrac12\right)(1+\Lambda_3),
\]

si risolve per \((\omega,r_0)\) **complessi simultaneamente**, essendo
\(\lambda=A_{\ell m}(a\omega)+a^2\omega^2-2am\omega\) dipendente da \(\omega\).

I \(\Lambda_j\) sono costruiti dai rapporti \(q_4/q_2\), \((q_3/q_2)^2\), … che
sono invarianti sotto \(q\leftrightarrow-V\): è la nota di convenzione
dell'Appendice A, qui usata attivamente.

## 3. Validazioni

| controllo | esito |
|---|---|
| \(a=0\) contro il WKB3 scalare di `core/` | \(3\times10^{-11}\) |
| residuo del sistema su Kerr | \(10^{-14}\)–\(10^{-12}\) |
| indipendenza da \(\mu\) per \(a=0\) | \(<10^{-9}\) |
| spostamento rispetto all'eikonale, decrescente in \(\ell\) | verificato |
| \(\operatorname{Re}\omega\) cresce con \(a\) (frame dragging) | verificato |

Lo spostamento rispetto alla frequenza eikonale va da \(1.5\times10^{-3}\) a
\(1.4\times10^{-2}\): **esattamente l'ordine che la diagnosi aveva previsto**.
È la terza conferma indipendente della stessa diagnosi.

## 4. Risultato

Errore relativo mediano fra previsione analitica e \(Q_M\) numerico, finestra
\(20<r<50\):

| \(a\) | \(\ell\) | con \(\omega\) eikonale | con \(\omega\) autoconsistente | guadagno |
|---|---|---|---|---|
| 0.0 | 30 | \(1.6\times10^{-3}\) | \(3.4\times10^{-4}\) | 5× |
| 0.3 | 30 | \(3.6\times10^{-2}\) | \(3.4\times10^{-4}\) | 107× |
| 0.3 | 70 | \(2.1\times10^{-2}\) | \(3.2\times10^{-4}\) | 65× |
| 0.6 | 30 | \(7.8\times10^{-2}\) | \(4.4\times10^{-4}\) | 179× |
| 0.6 | 70 | \(1.2\times10^{-3}\) | \(3.2\times10^{-5}\) | 37× |
| 0.6 | 100 | \(3.5\times10^{-2}\) | \(3.9\times10^{-4}\) | 90× |
| 0.9 | 30 | \(4.4\times10^{-3}\) | \(7.8\times10^{-5}\) | 57× |
| 0.9 | 70 | \(2.3\times10^{-1}\) | \(1.5\times10^{-3}\) | 161× |
| 0.9 | 100 | \(1.9\times10^{-1}\) | \(1.4\times10^{-3}\) | 134× |

Tutti i valori finiscono nella banda \(3\times10^{-5}\)–\(1.5\times10^{-3}\),
comparabile al controllo statico, e il comportamento erratico scompare.

**Conclusione.** La previsione dell'Appendice C, costruita e validata su
Pöschl–Teller, vale anche in presenza di rotazione. Non c'è nulla di patologico
nell'ampiezza su Kerr: la struttura di ordine del §9.2 resta l'unico contenuto.
Il passo B non è più parziale, e il manoscritto non ha più sezioni che
dichiarino «non stabilito».

## 5. Limite dichiarato

Il caso \(a=0.9\), \(\mu=0.9\) **non è incluso**: il bordo interno richiesto per
sopprimere la contaminazione dal ramo riflesso avvicina l'orizzonte al punto in
cui l'integratore esaurisce la precisione di macchina. Il regime quasi estremale
richiede un trattamento dedicato, che non affrontiamo e che è comunque oggetto di
letteratura propria (Yang *et al.* sui modi ramificati near-extremal;
Hatsuda–Shiga e Lo Chiatto *et al.* con exact WKB).

## 6. Effetto sul manoscritto

Nuovo **§9.6**, con la tabella prima/dopo e il limite dichiarato. Kerr è ora
trattato completamente, su potenziale (§9.2) e su ampiezza (§9.6). Appendice A
aggiornata con lo script.

Test totali: 39 (19 in `calculations`, 20 in `core`).
