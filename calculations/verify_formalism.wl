(* Independent Wolfram Language checks for the WKB-QNM-Madelung programme. *)

ClearAll["Global`*"];
Print["Wolfram kernel: ", $Version];

assertZero[label_String, expression_] := Module[{reduced},
  reduced = FullSimplify[expression];
  If[TrueQ[reduced === 0],
    Print["[ok] ", label],
    Print["[FAIL] ", label, ": ", InputForm[reduced]];
    Exit[1]
  ]
];

(* 1. Exact one-dimensional phase-amplitude closure. *)
psi = u[x]^(-1/2) Exp[I theta[x]/eps];
ode = Expand[
    eps^2 D[psi, {x, 2}]/psi + q[x]
  ] /. {
    Derivative[1][theta][x] -> u[x],
    Derivative[2][theta][x] -> u'[x]
  };
closure = q[x] - u[x]^2 +
  eps^2 D[u[x]^(-1/2), {x, 2}]/u[x]^(-1/2);
assertZero["exact 1D Madelung closure", ode - closure];

(* 2. Even WKB hierarchy through epsilon^4. *)
u0 = Sqrt[q[x]];
u2 = (5 q'[x]^2 - 4 q[x] q''[x])/(32 q[x]^(5/2));
u4 = (64 q[x]^3 q''''[x] - 448 q[x]^2 q'[x] q'''[x]
      - 304 q[x]^2 q''[x]^2 + 1768 q[x] q'[x]^2 q''[x]
      - 1105 q'[x]^4)/(2048 q[x]^(11/2));
uSeries = u0 + eps^2 u2 + eps^4 u4;
seriesResidual = Normal@Series[
    q[x] - uSeries^2 +
      eps^2 D[uSeries^(-1/2), {x, 2}]/uSeries^(-1/2),
    {eps, 0, 5}
  ];
assertZero["even WKB recursion through epsilon^4", seriesResidual];

(* 3. Exact Langer split of the Regge-Wheeler potential. *)
fSchw = 1 - 2/x;
vReggeWheeler = fSchw ((1/eps^2 - 1/4)/x^2 + 2 (1 - spin^2)/x^3);
vEikonal = fSchw/x^2 +
  eps^2 fSchw (2 (1 - spin^2)/x^3 - 1/(4 x^2));
assertZero["Regge-Wheeler Langer split", eps^2 vReggeWheeler - vEikonal];

(* 4. Exact scalar Kerr radial reduction and its slow-rotation limit. *)
ClearAll[rr, aa, mm, massK, om, sepK, psiK];
deltaK = rr^2 - 2 massK rr + aa^2;
h2K = rr^2 + aa^2;
hK = Sqrt[h2K];
dstarK[expression_] := deltaK/h2K D[expression, rr];
kK = h2K om - aa mm;
qK = (om - aa mm/h2K)^2 - deltaK sepK/h2K^2 -
  dstarK[dstarK[hK]]/hK;
radialK = D[deltaK D[psiK[rr]/hK, rr], rr] +
  (kK^2/deltaK - sepK) psiK[rr]/hK;
schrodingerK = dstarK[dstarK[psiK[rr]]] + qK psiK[rr];
assertZero[
  "exact scalar Kerr radial reduction",
  schrodingerK - deltaK radialK/hK^3
];

fK0 = 1 - 2 massK/rr;
qSchwarzschild = om^2 - fK0 (ellTerm/rr^2 + 2 massK/rr^3);
assertZero[
  "scalar Schwarzschild limit of Kerr",
  (qK /. {aa -> 0, sepK -> ellTerm}) - qSchwarzschild
];

qSlowKerr = Normal@Series[qK /. sepK -> ellTerm - 2 aa mm om,
  {aa, 0, 1}];
assertZero[
  "linear slow-Kerr frame dragging",
  qSlowKerr - (qSchwarzschild - 4 aa mm massK om/rr^3)
];

(* 5. Vaidya reduced Klein-Gordon equation and Madelung split. *)
ClearAll[amp, phase, f, pot, v, r];
psiV = amp[v, r] Exp[I phase[v, r]/eps];
reducedVaidya = eps^2 (
    2 D[psiV, v, r] + D[f[v, r] D[psiV, r], r] - pot[v, r] psiV
  )/psiV;
hjV = 2 D[phase[v, r], v] D[phase[v, r], r] +
  f[v, r] D[phase[v, r], r]^2;
ampV = (2 D[amp[v, r], v, r]
  + D[f[v, r] D[amp[v, r], r], r]
  - pot[v, r] amp[v, r])/amp[v, r];
continuityV = (
    D[amp[v, r]^2 D[phase[v, r], r], v]
    + D[amp[v, r]^2 (D[phase[v, r], v]
        + f[v, r] D[phase[v, r], r]), r]
  )/amp[v, r]^2;
assertZero[
  "Vaidya real/imaginary Madelung split",
  reducedVaidya - (-hjV + eps^2 ampV + I eps continuityV)
];

(* 6. Kodama-energy drift for ingoing Vaidya, K = partial_v. *)
mass = massFunction[v];
fV = 1 - 2 mass/r;
orbitMetric = {{-fV, 1}, {1, 0}};
waveVector = {wv, wr};
kodamaDrift = -1/2 waveVector . D[orbitMetric, v] . waveVector;
assertZero[
  "Kodama-energy drift",
  kodamaDrift + massFunction'[v] wv^2/r
];

Print["All Wolfram Language checks passed."];
