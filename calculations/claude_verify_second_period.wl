prec = 70; chat = 6/10;
A0 = SetPrecision[8974781212181061172514/10^22, 50];
A2 = SetPrecision[-2310651878076874644327/10^22, 50];
ells = Range[40, 160, 3];
data = Table[Module[{L, m, c, lam},
    L = ell + 1/2; m = (2 ell + 1)/3; c = chat L;
    lam = SpheroidalEigenvalue[ell, m, I SetPrecision[c, prec]];
    {L, (lam - c^2)/L^2}], {ell, ells}];

Print["TEST SENZA PARAMETRI ADATTATI:  L^4 (Ahat - A0 - A2/L^2)  deve essere costante"];
Print["   L        L^4 (Ahat - A0 - A2/L^2)"];
Do[Module[{L, v},
   L = data[[i, 1]]; v = data[[i, 2]];
   Print["  ", SetPrecision[L, 5], "     ",
    SetPrecision[L^4 (v - A0 - A2/L^2), 12]]],
  {i, {1, 5, 10, 15, 20, 30, 41}}]

Print[""];
Print["FIT CON A1 LIBERO, A0 fissato dall'azione, al variare di grado e intervallo"];
freefit[range_, deg_] := Module[{sub, m2, rhs, sol, res},
   sub = data[[range]];
   m2 = Table[Table[d[[1]]^(-k), {k, 1, deg}], {d, sub}];
   rhs = Table[d[[2]] - A0, {d, sub}];
   sol = LeastSquares[N[m2, 40], N[rhs, 40]];
   res = Max[Abs[N[m2, 40].sol - N[rhs, 40]]];
   Print["  ell ", sub[[1, 1]] - 1/2, "-", sub[[-1, 1]] - 1/2, ", grado ", deg,
     ":  A1 = ", SetPrecision[sol[[1]], 6],
     "   A2 = ", SetPrecision[sol[[2]], 12],
     "   res = ", SetPrecision[res, 3]]];
freefit[All, 3];
freefit[All, 4];
freefit[All, 5];
freefit[1 ;; 20, 4];
freefit[20 ;; 41, 4];
Print[""];
Print["  A2 di riferimento (forma chiusa) = ", SetPrecision[A2, 12]];
