(* Exact local normalization covariance; NOT a boundary/Fredholm proof. *)
ClearAll[r,m,c,g,alpha,mdot,dw,nn,pp];
source=-2 D[g[r,m],r,m];
transformed=-2 D[c[m] g[r,m],r,m];
sourceCheck=Simplify[transformed-(c[m] source-2 c'[m] D[g[r,m],r])]===0;
original=2 (alpha-I dw) D[g[r,m],r]-mdot source;
changed=2 (alpha-mdot c'[m]/c[m]-I dw) D[c[m] g[r,m],r]-mdot transformed;
covariance=Simplify[changed-c[m] original]===0;
newN=c[m]^2 nn-2 c[m] c'[m] pp;
newP=c[m]^2 pp;
connection=Simplify[newN/(2 newP)-(nn/(2 pp)-c'[m]/c[m])]===0;
observable=Simplify[newN/(2 newP)+D[c[m] g[r,m],m]/(c[m] g[r,m])
 -(nn/(2 pp)+D[g[r,m],m]/g[r,m])]===0;
checks={sourceCheck,covariance,connection,observable};
Print["Vaidya normalization checks (local): ",checks];
If[And@@checks,Exit[0],Exit[1]];
