from p_solver import *

import sympy as sp

x = sp.Symbol("x")
e = sp.Symbol("epsilon")
plhs = sp.sin(x)*e+x**5
s = PPolySolver(x, e, taylor(plhs, x, 5), 3)
pprint(s.pLHS)
s.print_inf_series(expanded=True)
pprint(s.solve())
pprint(s.coeffs_inf)