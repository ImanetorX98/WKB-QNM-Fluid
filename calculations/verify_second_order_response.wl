ClearAll[t,z,k,j,d,e,q,qw,qww,p,pr,ps,eps,u,v];
zz=z+t k d+t^2(k e+j d^2)/2;
c1=Simplify[D[zz^2,{t,2}]/.t->0]===2 k^2 d^2+2 z (k e+j d^2);
c1=Simplify[(D[zz^2,{t,2}]/.t->0)-(2 k^2 d^2+2 z (k e+j d^2))]===0;
qm=eps^2(q+t qw+t^2 qww/2-(p+t pr+t^2 ps/2)^2);
c2=Simplify[(D[qm,{t,2}]/.t->0)-eps^2(qww-2pr^2-2p ps)]===0;
c3=Simplify[D[u[t]/v[t],{t,2}]-(u''[t]-(u[t]/v[t])v''[t]-
 2 D[u[t]/v[t],t]v'[t])/v[t]]===0;
(* Moving simple zero: exact integral of |x-t| on [-1,1] is 1+t^2. *)
moving=Integrate[t-x,{x,-1,t}]+Integrate[x-t,{x,t,1}];
c4=Simplify[D[moving,{t,2}]-2]===0;
Print["Second-order response algebra: ",{c1,c2,c3,c4}];
If[And[c1,c2,c3,c4],Exit[0],Exit[1]];
