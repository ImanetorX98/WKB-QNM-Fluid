ClearAll[x,ep,amp,g,z,h,u,t];
qOriginal=-ep^2 amp''[x]/amp[x];
qTransformed=-ep^2 D[Exp[g[x]]amp[x],{x,2}]/(Exp[g[x]]amp[x]);
c1=Simplify[qTransformed-qOriginal+ep^2(g''[x]+2 g'[x]amp'[x]/amp[x]+g'[x]^2)]===0;
(* Linearized Riccati at fixed frequency: z'+z^2+q=0. *)
c2=Expand[Coefficient[D[z[x]+t h[x],x]+(z[x]+t h[x])^2+t u[x],t]]===
 Expand[h'[x]+2 z[x]h[x]+u[x]];
(* Ingoing first-order logarithmic boundary coefficient. *)
ClearAll[v1,z0,z1,q1,d];
c3=Simplify[Coefficient[v1 d z1+q1 d+2 z0 z1 d,d]/.
 z1->-q1/(v1+2 z0)]===0;
Print["Diagnostic scope and boundary identities: ",{c1,c2,c3}];
If[And[c1,c2,c3],Exit[0],Exit[1]];
