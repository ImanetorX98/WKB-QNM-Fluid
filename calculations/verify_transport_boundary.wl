(* Derivative along r_*; R'=D R, d(D_omega)/dr_*=-2omega-2D D_omega. *)
ClearAll[rr,dd,dw,w];
derivative[expr_]:=D[expr,rr] dd rr + D[expr,dw] (-2 w-2 dd dw);
fp=(1-I dw) rr^2/2;
fn=-dw rr^2/(2 w);
checks={Simplify[derivative[fp]-rr^2 (dd+I w)]===0,
        Simplify[derivative[fn]-rr^2]===0};
Print["Boundary primitive identities: ",checks];
If[And@@checks,Exit[0],Exit[1]];
