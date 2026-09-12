prec = 60; chat = 6/10;
A0 = SetPrecision[89747812121810611725140358/10^26, 45];
ells = Range[40, 100, 3];  (* tutti = 1 mod 3 *)
data = Table[Module[{L, m, c, lam},
    L = ell + 1/2; m = (2 ell + 1)/3; c = chat L;
    lam = SpheroidalEigenvalue[ell, m, I SetPrecision[c, prec]];
    {L, L^2 ((lam - c^2)/L^2 - A0)}], {ell, ells}];
Print["Codex: la quantizzazione ha solo potenze PARI di epsilon -> A1 = A3 = 0."];
Print["Test: adattare  A2 + b/L + c/L^2  e guardare se b e' zero.\n"];
fit[basis_, label_] := Module[{m2, rhs, sol, res},
   m2 = Table[Table[d[[1]]^(-k), {k, basis}], {d, data}];
   rhs = data[[All, 2]];
   sol = LeastSquares[N[m2, 30], N[rhs, 30]];
   res = Max[Abs[N[m2, 30].sol - N[rhs, 30]]];
   Print[label, "  coefficienti: ", SetPrecision[sol, 10]];
   Print["      residuo max = ", SetPrecision[res, 4]];
   sol];
fit[{0, 1, 2}, "A2 + b/L + c/L^2   (b libero, test di A3)"];
fit[{0, 2}, "A2 + c/L^2         (solo pari)                "];
fit[{0, 2, 4}, "A2 + c/L^2 + d/L^4 (solo pari, due termini)"];
fit[{0, 1, 2, 3}, "A2 + b/L + c/L^2 + d/L^3                  "];
