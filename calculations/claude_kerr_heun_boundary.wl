a = 3/5; mm = 19; om = 3 - I/4;
Amma = SpheroidalEigenvalue[28, 19, I SetPrecision[a om, 40]];
lam = (Amma - (a om)^2) + a^2 om^2 - 2 a mm om;
rp = 1 + Sqrt[1 - a^2]; rm = 1 - Sqrt[1 - a^2]; b = rm - rp;
kap = (rp^2 + a^2) om - a mm;
p = -I kap/(rp - rm); qe = SetPrecision[0.0625 - 6.375 I, 30]; s = -I om b;
qH = SetPrecision[746.823871198005 + 16.110026570174174 I, 30];
al = SetPrecision[0.8 + 9.6 I, 30]; ga = 2 p + 1; de = 2 qe + 1; ep = 2 s;

(* q corretto: v derivata SIMBOLICAMENTE, non con D0 gia' valutato *)
vSym[x_] := (x^2 - 2 x + a^2)/(x^2 + a^2);
qT[r0_] := Module[{S = r0^2 + a^2, D0 = r0^2 - 2 r0 + a^2, vv, vr, G},
   vv = D0/S; vr = D[vSym[x], x] /. x -> r0;
   G = vv (vr r0/S + vv a^2/S^2);
   (om - mm a/S)^2 - D0 lam/S^2 - G];
Sig[r0_] := r0^2 + a^2; v[r0_] := (r0^2 - 2 r0 + a^2)/Sig[r0];
Dm[r0_] := Module[{z0 = (r0 - rp)/b, h, hp},
   h = HeunC[qH, al, ga, de, ep, z0]; hp = HeunCPrime[qH, al, ga, de, ep, z0];
   v[r0] (r0/Sig[r0] + (1/b) (p/z0 + qe/(z0 - 1) + s + hp/h))];

Print["residuo della Riccati con q CORRETTO"];
Print["   rho        |v D' + D^2 + q|"];
Do[Module[{r0 = rp + rho, dd, der, hh = 10^-7},
   dd = Dm[r0]; der = (Dm[r0 + hh] - Dm[r0 - hh])/(2 hh);
   Print["  ", N[rho, 3], "     ", SetPrecision[Abs[v[r0] der + dd^2 + qT[r0]], 4]]],
  {rho, {10^-4, 10^-3, 2 10^-2, 1/10, 1/2, 1}}]
Print[""];
Print["D_minus a raggi dove la serie di Frobenius diverge:"];
Do[Print["  rho = ", N[rho, 3], "   D = ", SetPrecision[Dm[rp + rho], 12]],
  {rho, {2 10^-2, 1/10, 1/2, 1}}]
