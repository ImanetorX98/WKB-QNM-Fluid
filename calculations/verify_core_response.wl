ClearAll[x,p,f,v,b,om,dw,j,k,q,zr,zi,t,hr,hi,dq,eps];
w=p[x] f'[x]-p'[x] f[x];
c1=Simplify[(D[w,x]/.{p''[x]->(v[x]-om^2)p[x],
 f''[x]->(v[x]-om^2)f[x]+(b[x]-2om dw)p[x]})-
 (b[x]-2om dw)p[x]^2]===0;
k=(-I-2om j[x])/p[x]^2;
c2=Simplify[(D[k,x]+2p'[x]/p[x] k+2om)/.j'[x]->p[x]^2]===0;
(* Real amplitude: A''/A=Re(-q-z^2)+(Re z)^2. *)
c3=Simplify[ComplexExpand[-eps^2(Re[-q-(zr+I zi)^2]+zr^2)]-eps^2(q-zi^2)]===0;
c4=Simplify[Coefficient[eps^2((q+t dq)-(zi+t hi)^2),t]-
 eps^2(dq-2zi hi)]===0;
Print["Core response identities: ",{c1,c2,c3,c4}];
If[And[c1,c2,c3,c4],Exit[0],Exit[1]];
