from p_solver import *

import sympy as sp

x = sp.Symbol("x")
e = sp.Symbol("epsilon")
plhs = x**5 + taylor(sp.exp(x), 0, 3) + x - e
s = PPolySolver(x, e, plhs, 3)
print(s)
s.print_inf_series(expanded=True)
pprint(s.solve())
pprint(s.coeffs_inf)