import sympy as sp
from sympy.abc import n

def pprint(x):
    print()
    sp.pprint(x)

class PSolver:
    def __init__(self, x, e, LHS, perturbed_LHS, power_of_epsilon):
        self.LHS = LHS # RHS is always = 1.
        self.x = x
        self.e = e
        self.pLHS = perturbed_LHS
        self.a = sp.IndexedBase(sp.Symbol("a"))
        self.poe = power_of_epsilon
        self.inf_series = sp.Sum(self.a[n]*self.e**n, (n, 0, self.poe))
        self.coeffs_inf = {}

    def solve(self):
        self.inf_sol = self.inf_series
        eqnL = self.pLHS.subs(self.x, self.inf_series.doit().expand()).expand()
        eq = []
        eq_sL = []
        for j in range(self.poe + 1):
            if j == 0:
                eq.append(sp.Eq(eqnL.coeff(self.e, j), 1))
                if eq == [False]:
                    print("Unable to find perturbation. No constant term in perturbed LHS!")
                    exit(-1)
            else:
                eq.append(sp.Eq(eq_sL[j - 1], 0))
            self.coeffs_inf["a" + str(j)] = sp.solveset(eq[j], self.a[j], sp.Reals)
            tmpL = eqnL.coeff(self.e, j + 1)
            for k in range(self.poe):
                try:
                    tmpL = tmpL.subs(self.a[k], self.coeffs_inf["a" + str(k)].args[0])
                except KeyError:
                    pass
            eq_sL.append(tmpL)

            self.inf_sol = self.inf_sol.doit().subs(self.a[j], self.coeffs_inf["a" + str(j)].args[0])

        self.inf_sol = self.inf_sol.subs(self.e, 1)
        self.sol_poly = sp.Eq(self.x, self.inf_sol.n(20))
        return self.sol_poly

    def __repr__(self):
        pprint(sp.Eq(self.LHS, 1))
        return ""

    def __str__(self):
        pprint(sp.Eq(self.LHS, 1))
        return ""
    
    def print_perturbed(self):
        pprint(sp.Eq(self.pLHS, 1))
    
    def print_inf_series(self, expanded=False):
        if not expanded: 
            pprint(sp.Eq(self.x, self.inf_series))
        else: 
            pprint(sp.Eq(self.x, self.inf_series.doit().expand()))
