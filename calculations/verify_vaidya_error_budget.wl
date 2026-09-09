(* Local symbolic checks for the numerical error-budget audit. *)
ClearAll[d, w, h1, h2, h, f, u, residual, a, b, z, e];
f = 1 - 2/(2+d);
u = 6/(2+d)^2 + 2/(2+d)^3;
h1 = (7/4)/(1/2-2 I w);
h2 = ((7/4+1/2) h1-15/8)/(2-4 I w);
h = 1+h1 d+h2 d^2;
residual = f D[h,{d,2}] + (D[f,d]-2 I w) D[h,d] - u h;
checks = {
 Simplify[SeriesCoefficient[residual,{d,0,0}]] === 0,
 Simplify[SeriesCoefficient[residual,{d,0,1}]] === 0,
 Simplify[Limit[((z/2) Coth[z/2]-1)/z^2,z->0]] === 1/12,
 Simplify[D[Exp[a e],e]/.e->0] === a
};
Print["Vaidya error budget symbolic checks: ", checks];
If[And@@checks, Exit[0], Exit[1]];
