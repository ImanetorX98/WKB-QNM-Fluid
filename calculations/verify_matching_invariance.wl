ClearAll[x,p,d,dw,de,qw,qe,n,fw,pb,ip];
(* d=exterior logarithmic derivative; dw,de are parameter derivatives. *)
surfaceW=dw[x]p[x]^2;
c1=Simplify[(D[surfaceW,x]+qw[x]p[x]^2)/.
 {p'[x]->d[x]p[x],dw'[x]->-2d[x]dw[x]-qw[x]}]===0;
surfaceE=de[x]p[x]^2;
c2=Simplify[(D[surfaceE,x]+qe[x]p[x]^2)/.
 {p'[x]->d[x]p[x],de'[x]->-2d[x]de[x]-qe[x]}]===0;
(* Integrated Wronskian at a root: pb F_omega + N = 0. *)
c3=Simplify[(pb fw+n)/.fw->-n/pb]===0;
Print["Matching invariance algebra: ",{c1,c2,c3}];
If[And[c1,c2,c3],Exit[0],Exit[1]];
