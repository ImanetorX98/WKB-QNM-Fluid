ClearAll[mu,f,fp,k,g,gp,z,zp,u,forcing];
gpp=(u g-(fp-k)gp)/f;
zpp=(u z-(fp-k)zp+forcing)/f;
jprime=mu (fp-k)(g zp-gp z)+mu f (g zpp-gpp z);
check=Simplify[jprime-mu g forcing]===0;
Print["Forced Wronskian identity: ",check];
If[check,Exit[0],Exit[1]];
