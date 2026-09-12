import sympy as sp
th, A = sp.symbols("theta A")
Q = sp.Function("Q0")(th)
p = sp.sqrt(Q)
Q2 = sp.Symbol("Q2")

p2 = Q2 / (2 * p) - sp.diff(p, th, 2) / (4 * p**2) + 3 * sp.diff(p, th) ** 2 / (8 * p**3)
derivative_part = sp.simplify(p2 - Q2 / (2 * p))
print("parte derivativa di p2:")
print("   ", sp.simplify(derivative_part))

total = sp.diff(sp.diff(Q, th) / Q ** sp.Rational(3, 2), th)   # d/dtheta [Q0'/Q0^{3/2}]

first = -sp.diff(Q, th) ** 2 / (32 * Q ** sp.Rational(5, 2))
print("\npasso 1:  parte derivativa  ==  -Q0'^2/(32 Q0^{5/2})  + derivata totale?")
print("   residuo con coefficiente -1/8 :",
      sp.simplify(derivative_part - first - sp.Rational(-1, 8) * total))

second = -sp.diff(Q, th, 2) / (48 * Q ** sp.Rational(3, 2))
print("\npasso 2:  -Q0'^2/(32 Q0^{5/2})  ==  -Q0''/(48 Q0^{3/2}) + derivata totale?")
print("   residuo con coefficiente +1/48:",
      sp.simplify(first - second - sp.Rational(1, 48) * total))

print("\npasso 3:  (1/24) d/dA [Q0''/sqrt(Q0)]  con dQ0/dA = 1")
rest = sp.Symbol("rest")
Qpp = sp.Symbol("Qpp")
expr = sp.Rational(1, 24) * sp.diff(Qpp / sp.sqrt(A + rest), A)
print("   =", sp.simplify(expr), "   cioe' -Q0''/(48 Q0^{3/2})  ->",
      sp.simplify(expr + Qpp / (48 * (A + rest) ** sp.Rational(3, 2))) == 0)
