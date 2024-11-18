from p_solver import *

import sympy as sp

x = sp.Symbol("x")
e = sp.Symbol("epsilon")
plhs = sp.sin(e*x)
s = PTranscendentalSolver(x, e, plhs, 3, 4)
pprint(s.pLHS)
#s.print_inf_series(expanded=True)
#pprint(s.solve())
#pprint(s.coeffs_inf)