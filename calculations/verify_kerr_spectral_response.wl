ClearAll[r,a,m,om,aa,bb,hh,dd,qq,geom,iw,ie,pm,pp,dmw,dpw,dme,dpe,dw];
hh=r^2+a^2; dd=r^2-2r+a^2;
qq=(om-m a/hh)^2-dd (aa[a om]+a^2 om^2-2a m om)/hh^2-geom[r];
expected=2(om-m a/hh)-dd/hh^2 (a aa'[a om]+2a^2 om-2a m);
c1=Simplify[D[qq,om]-expected]===0;
c2=Simplify[expected/.a->0]===2om;
(* Integrated Wronskian: endpoint difference + integral q_eta+q_omega dw=0. *)
balance=(dpe+dpw dw)pp-(dme+dmw dw)pm+ie+iw dw;
answer=(-ie-dpe pp+dme pm)/(iw+dpw pp-dmw pm);
c3=Simplify[balance/.dw->answer]===0;
(* Finite-domain factorization k'+2zk=-q_omega is separate from HF proof. *)
Print["Kerr and nonlinear spectral identities: ",{c1,c2,c3}];
If[And[c1,c2,c3],Exit[0],Exit[1]];
