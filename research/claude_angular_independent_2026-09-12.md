# Verifica indipendente dell'autovalore angolare: compito (a) chiuso

**Non è servito scrivere un solutore.** `SpheroidalEigenvalue` è built-in in
Mathematica: implementazione esterna, autori esterni, metodo diverso dalla
nostra matrice in base armonica sferica. Il compito (a) del
`CODEX_HANDOFF_KERR_A1.md` — «non fidatevi del nostro numero perché è nostro» —
si chiude senza lavoro analitico.

## 1. La convenzione, identificata e verificata

I due valori non coincidono, e lo scarto è **esattamente \(c^2\)**:

| \((\ell,m)\) | \(c\) | nostro | Mathematica | differenza | \(c^2\) |
|---|---|---|---|---|---|
| (28,19) | 0 | 812.0000000000 | 812.0000000000 | 0 | 0 |
| (28,19) | 5.7 | 802.9513373617 | 835.4413373617 | 32.49 | **32.49** |
| (28,19) | 11.4 | 775.4855918968 | 905.4455918968 | 129.96 | **129.96** |
| (28,19) | 8.55+2.85i | 793.9079−13.6952i | 858.8879+35.0398i | 64.98+48.735i | **64.98+48.735i** |

Le cifre decimali coincidono (…3617, …8968): è una convenzione, non un
disaccordo. Quindi \(A_{\rm nostro}=\lambda_{\rm MMA}-c^2\) e di conseguenza

\[
\partial_c A_{\rm nostro} = \partial_c\lambda_{\rm MMA} - 2c .
\]

## 2. \(A_c\): quattordici cifre

Formula bilineare \(A_c=-2c\,(v^{T}M_2v)/(v^{T}v)\) contro differenziazione di
`SpheroidalEigenvalue` a precisione 40:

| \((\ell,m)\) | \(c\) | Mathematica \(-2c\) | nostro |
|---|---|---|---|
| (28,19) | 5.7 | −3.18456480678527482 | −3.18456480678528 |
| (28,19) | 11.4 | −6.47871147148527310 | −6.47871147148529 |
| (61,41) | 12.3 | −6.87404773688933829 | −6.87404773688934 |
| (61,41) | 24.6 | −13.98312312105877890 | −13.98312312105877 |
| (28,19) | 8.55+2.85i | −4.7919915721 − 1.6437501787i | −4.791991572 − 1.643750179i |
| (61,41) | 18.45+6.15i | −10.3435931225 − 3.5474379816i | −10.343593122 − 3.547437982i |

Quattordici cifre sui reali, dieci sui complessi.

## 3. \(A_1=0\) rifatto **interamente** in Mathematica

Nessun nostro codice nella catena. \(\mu=2/3\), \(\ell\equiv1\pmod 3\) così che
\(m=(2\ell+1)/3\) sia intero **esatto**, \(\hat c=0.6\), precisione 50.

\(A_0 = 0.89747812500729378719\) (il nostro solutore dava 0.89747812).

\[
L^2\!\left(\frac{A}{L^2}-A_0\right):
\]

| \(L\) | 40.5 | 43.5 | 46.5 | 49.5 | 52.5 | 55.5 |
|---|---|---|---|---|---|---|
| | −0.2311013 | −0.2310983 | −0.2310960 | −0.2310945 | −0.2310934 | −0.2310927 |

**Costante a cinque cifre**, con la deriva residua compatibile con il termine
\(A_3/L^3\). Il nostro valore era \(A_2\approx-0.2311\): identico.

> \(A_1=0\) è confermato da un'implementazione con cui non abbiamo nulla in
> comune. Il compito (a) è chiuso.

## 4. Che cosa resta

Il compito **(b)** — dimostrare \(A_1=0\) dalla quantizzazione di
Bohr–Sommerfeld angolare, dove né il membro sinistro né il destro contengono
\(1/L\) — resta aperto ed è **analitico**. È l'unico pezzo che trasformerebbe un
fatto numerico ben stabilito in un enunciato.

## Comando

```sh
/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noprompt \
  -script calculations/claude_verify_angular_independent.wl
```
