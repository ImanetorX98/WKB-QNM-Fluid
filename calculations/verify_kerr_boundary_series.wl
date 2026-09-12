ClearAll[r,a,m,om,lam,del,h,v,g,q,t,rp,rm,gap,k,b,d2,d3];
h=r^2+a^2; del=r^2-2r+a^2; v=del/h;
g=v(D[v,r]r/h+v a^2/h^2);
q=(om-m a/h)^2-del lam/h^2-g;
qinf=Normal[Series[q/.r->1/t,{t,0,3}]];
b=lam+2a m om;
c1=Simplify[qinf-(om^2-b t^2+(2lam-2)t^3)]===0;
d2=b/(2I om); d3=(2d2-(2lam-2))/(2I om);
z=I om+d2 t^2+d3 t^3;
c2=Simplify[Normal[Series[-t^2(v/.r->1/t)D[z,t]+z^2+qinf,{t,0,3}]]]===0;
(* Express horizon parameters through rp; M=1 gives a^2=2rp-rp^2. *)
q1=D[q,r]/.r->rp;
hp=rp^2+a^2; gap=2rp-2; k=om-m a/hp;
expect=4m a rp k/hp^2-gap lam/hp^2-gap^2 rp/hp^3;
c3=FullSimplify[q1-expect,Assumptions->a^2==2rp-rp^2 && rp>1]===0;
Print["Kerr boundary coefficients: ",{c1,c2,c3}];
If[And[c1,c2,c3],Exit[0],Exit[1]];
