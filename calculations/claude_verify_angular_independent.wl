(* A_1 = 0 verificato interamente con SpheroidalEigenvalue di Mathematica.
   mu = 2/3 ed ell = 1 mod 3  ->  m = (2 ell + 1)/3 intero ESATTO.
   Convenzione: A_nostro = lambda_MMA - c^2, verificata con lo scarto c^2. *)
prec = 50; chat = 6/10;
ells = {40, 43, 46, 49, 52, 55, 58, 61};
data = Table[
   Module[{L, m, c, lam, A},
    L = ell + 1/2; m = (2 ell + 1)/3; c = chat L;
    lam = SpheroidalEigenvalue[ell, m, I SetPrecision[c, prec]];
    A = lam - c^2;
    {N[L, 20], SetPrecision[A/L^2, 25]}], {ell, ells}];
(* A0 per estrapolazione di Richardson sui due L piu grandi *)
A0 = ((data[[-1, 1]]^2 data[[-1, 2]] - data[[-2, 1]]^2 data[[-2, 2]])/
     (data[[-1, 1]]^2 - data[[-2, 1]]^2));
Print["A0 = ", SetPrecision[A0, 18]];
Print["L^2 (A/L^2 - A0):  se costante, A_1 = 0"];
Do[Print["  L = ", SetPrecision[data[[i, 1]], 6], "   ",
   SetPrecision[data[[i, 1]]^2 (data[[i, 2]] - A0), 14]], {i, 1, 6}]
