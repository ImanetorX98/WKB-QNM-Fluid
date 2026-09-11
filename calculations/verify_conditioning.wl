ClearAll[a,b,ll,delta,q,kappa];
q=-b^2/ll^2;
assumptions=Element[{a,b,ll},Reals] && b<0 && ll>0;
radial=FullSimplify[(a D[q,a]+b D[q,b])/q,assumptions];
kappa=FullSimplify[Sqrt[a^2+b^2] Sqrt[D[q,a]^2+D[q,b]^2]/Abs[q],assumptions];
checks={radial===2, D[q,a]===0,
 FullSimplify[b D[q,b]/q,assumptions]===2,
 FullSimplify[kappa-2 Sqrt[a^2+b^2]/Abs[b],assumptions]===0};
Print["Conditioning checks: ",checks];
If[And@@checks,Exit[0],Exit[1]];
