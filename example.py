from p_solver import *

import sympy as sp

x = sp.Symbol("x")
e = sp.Symbol("epsilon")
#plhs = sp.sin(x)*e+x**5
#plhs = taylor(sp.sin(x*e), x, 3)+x**(105)
#plhs = 12*sp.log(2)**4 * sp.factor(2*taylor(sp.exp(x)/sp.ln(x+2), x, 3)+x**2+e*x)
plhs = taylor(sp.tan(x), x, 4) + x*e
rhs = -2

s = PPolySolver(x, e, plhs, 3, rhs)
pprint(sp.Eq(s.pLHS, rhs))
s.print_inf_series(expanded=True)
try:
    pprint(s.solve())
except KeyboardInterrupt:
    pprint(s.coeffs_inf)
pprint(s.coeffs_inf)