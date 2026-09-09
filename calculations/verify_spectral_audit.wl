(* Independent symbolic checks for research/spectral_audit_2026-09-08.md. *)
ClearAll["Global`*"];
Print["Wolfram kernel: ", $Version];
checks = 0;
assertZero[label_String, expression_, assumptions_: True] := Module[{value},
  value = FullSimplify[expression, Assumptions -> assumptions];
  If[!TrueQ[value === 0],
    Print["[FAIL] ", label, ": ", InputForm[value]]; Exit[1]];
  checks++; Print["[ok] ", label]
];

(* Fundamental amplitude; its overall normalization cancels. *)
amp = Sqrt[Cosh[y]];
assertZero["fundamental amplitude curvature",
  D[amp, {y, 2}]/amp - (1 + Sech[y]^2)/4, Element[y, Reals]];
phaseMomentum = Sqrt[1 - eps^2/4] Tanh[y];
omegaHat0 = Sqrt[1 - eps^2/4] - I eps/2;
assertZero["fundamental normalized frequency modulus",
  Abs[omegaHat0]^2 - 1, 0 < eps < 2];
assertZero["fundamental denominator density",
  1 + Sech[y]^2 + phaseMomentum^2 - (2 - eps^2 Tanh[y]^2/4)];
assertZero["fundamental integrated indicator",
  (eps^2 (2 - moment)/4)/(2 - eps^2 moment/4)
    - eps^2 (2 - moment)/(8 - eps^2 moment)];

(* N is reserved in Wolfram Language; use overtoneN = n + 1/2. *)
bb = overtoneN^2 + 1/4;
omegaHatExact = Sqrt[1 - eps^2/4] - I overtoneN eps;
omegaHat1 = Sqrt[1 - 2 I overtoneN eps];
omegaHat3 = Sqrt[1 - bb eps^2 - 2 I overtoneN eps (1 - eps^2/8)];
assertZero["WKB1 leading relative complex error",
  Normal@Series[omegaHat1/omegaHatExact - 1, {eps, 0, 2}] - bb eps^2/2];
assertZero["WKB3 leading relative complex error",
  Normal@Series[omegaHat3/omegaHatExact - 1, {eps, 0, 5}]
    + I overtoneN eps^5/128];
(* Independent reduction of the even-barrier Iyer-Will coefficients. *)
v2 = -2 barrierL^2; v4 = 16 barrierL^2; v6 = -272 barrierL^2;
lambda2 = (v4/v2) bb/(8 (2 barrierL));
lambda3 = ((v4/v2)^2 (67 + 68 overtoneN^2)/2304
  - (v6/v2) (5 + 4 overtoneN^2)/288)/(4 barrierL^2);
assertZero["WKB3 second-order coefficient", 2 barrierL lambda2 + bb];
assertZero["WKB3 third-order coefficient", lambda3 + 1/(8 barrierL^2)];

nuAtExact = I (omegaHatExact^2 - 1)/(2 eps) - 1/2;
nuMismatch = overtoneN (Sqrt[1 - eps^2/4] - 1) - I bb eps/2;
assertZero["quadratic index mismatch",
  nuAtExact - (overtoneN - 1/2) - nuMismatch];
assertZero["quadratic index imaginary part",
  Im[nuAtExact] + bb eps/2, 0 < eps < 2 && overtoneN > 0];
assertZero["parabolic WKB1 quantization",
  I (omegaHat1^2 - 1)/(2 eps) - 1/2 - (overtoneN - 1/2)];

(* Reduction to Gauss' equation at arbitrary trial frequency, not at a QNM. *)
aa = 1/2 + I spectralLambda - I omega;
ab = 1/2 - I spectralLambda - I omega;
ac = 1 - I omega;
hz = 2 z (1 - z);
lp = (-I omega/2) (1/z - 1/(1 - z));
coefficientF1 = (2 hz^2 lp + hz D[hz, z])/(4 z (1 - z));
coefficientF0 = (hz^2 (D[lp, z] + lp^2) + hz D[hz, z] lp
  + omega^2 - 4 barrierL^2 z (1 - z))/(4 z (1 - z));
assertZero["Gauss first-derivative coefficient", coefficientF1 - (ac - (aa + ab + 1) z)];
assertZero["Gauss zeroth-derivative coefficient", coefficientF0 + aa ab,
  spectralLambda^2 == barrierL^2 - 1/4];
assertZero["Gauss right incoming exponent", ac - aa - ab - I omega];
assertZero["positive-frequency QNM from gamma zero",
  (aa /. omega -> spectralLambda - I (nn + 1/2)) + nn];
cin = Gamma[1 - I omega] Gamma[-I omega]/(Gamma[aa] Gamma[ab]);
logSlope = -I (PolyGamma[0, 1 - I omega] + PolyGamma[0, -I omega]
  - PolyGamma[0, aa] - PolyGamma[0, ab]);
assertZero["Jost incoming logarithmic derivative", D[cin, omega]/cin - logSlope];

(* Local polar residual and its exact connection to a boundary Wronskian. *)
trialPsi = amplitude[x] Exp[I phase[x]/eps];
polarResidual = (qr[x] - phase'[x]^2 + eps^2 amplitude''[x]/amplitude[x]
  + I (qi[x] + eps (phase''[x] + 2 phase'[x] amplitude'[x]/amplitude[x])));
assertZero["full polar residual includes real and transport defects",
  (eps^2 D[trialPsi, {x, 2}] + (qr[x] + I qi[x]) trialPsi)/trialPsi
    - polarResidual];
wronskian = testChi[x] trialF'[x] - testChi'[x] trialF[x];
lagrangeIdentity = (D[wronskian, x]
  - testChi[x] (eps^2 trialF''[x] + q[x] trialF[x])/eps^2);
assertZero["boundary Wronskian sourced by complex residual",
  lagrangeIdentity /. testChi''[x] -> -q[x] testChi[x]/eps^2];
Print["All ", checks, " spectral audit checks passed."];
