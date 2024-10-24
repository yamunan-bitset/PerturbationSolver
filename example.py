from p_solver import *

import sympy as sp

x = sp.Symbol("x")
e = sp.Symbol("epsilon")
lhs = x**5 + x
plhs = e*x**5 + x**7
s = PSolver(x, e, lhs, plhs, 3)
print(s)
s.print_perturbed()
s.print_inf_series(expanded=True)
pprint(s.solve())