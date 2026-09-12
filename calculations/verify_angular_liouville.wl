ClearAll[t,x,aa,cc,mm,ll,eps,mu,chat,ahat,y,psi];
ss=psi[t]/Sqrt[Sin[t]];
op=D[ss,{t,2}]+Cot[t] D[ss,t]+(aa+cc^2 Cos[t]^2-mm^2/Sin[t]^2)ss;
expected=D[psi[t],{t,2}]+(aa+1/4+cc^2 Cos[t]^2-(mm^2-1/4)/Sin[t]^2)psi[t];
c1=FullSimplify[Sqrt[Sin[t]]op-expected,0<t<Pi]===0;
scaled=eps^2 (aa+1/4+cc^2 Cos[t]^2-(mm^2-1/4)/Sin[t]^2)/.
 {aa->ahat/eps^2,cc->chat/eps,mm->mu/eps};
c2=Simplify[scaled-(ahat+chat^2 Cos[t]^2-mu^2/Sin[t]^2+eps^2 (1+Csc[t]^2)/4)]===0;
c3=Simplify[(ll-mu ll) /ll-(1-mu)]===0;
(* Counterexample to classifying order from differential order alone. *)
c4=Simplify[(ll^2 y[x]^2+ll y'[x])/ll^2-(y[x]^2+y'[x]/ll)]===0;
Print["Angular Liouville identities: ",{c1,c2,c3,c4}];
If[And[c1,c2,c3,c4],Exit[0],Exit[1]];
