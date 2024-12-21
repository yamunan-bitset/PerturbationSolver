import sympy as sp
from inspect import getsource
from sympy.abc import n

def pprint(x):
    print()
    sp.pprint(x)

def taylor(pLHS, x, pot):
    # compute taylor expansion symbolically 
    i = 0
    taylor_exp = sp.Integer(0)
    while i <= pot:
        taylor_exp = taylor_exp + (sp.diff(pLHS, x, i).subs(x, 0))/(sp.factorial(i))*(x-0)**i
        i += 1
    return taylor_exp

class PPolySolver:
    '''Perturbation Polynomial Solver'''
    def __init__(self, x, e, perturbed_LHS, power_of_epsilon, rhs=1):
        self.x = x
        self.e = e
        self.pLHS = perturbed_LHS
        self.rhs = rhs
        self.a = sp.IndexedBase(sp.Symbol("a"))
        self.poe = power_of_epsilon
        self.inf_series = sp.Sum(self.a[n]*self.e**n, (n, 0, self.poe))
        self.coeffs_inf = {}

    def solve(self, verbose=False):
        '''Group terms of same powers and solve for a[n]'''
        self.inf_sol = self.inf_series
        '''Compute perturbation equation by plugging in infinite series for x'''
        eqnL = self.pLHS.subs(self.x, self.inf_series.doit().expand()).expand()
        eq = []
        eq_sL = []
        for j in range(self.poe + 1):
            if j == 0:
                eq.append(sp.Eq(eqnL.coeff(self.e, j), self.rhs))
                if eq == [False]:
                    print("FIRST: Unable to find perturbation. Try changing perturbed LHS")
                    exit(-1)
            else:
                eq.append(sp.Eq(eq_sL[j - 1], 0))
            self.coeffs_inf["a" + str(j)] = sp.solveset(eq[j], self.a[j], sp.Reals)
            '''Group like powers'''
            tmpL = eqnL.coeff(self.e, j + 1)
            for k in range(self.poe):
                try:
                    tmpL = tmpL.subs(self.a[k], self.coeffs_inf["a" + str(k)].args[0])
                except KeyError:
                    pass
                except IndexError:
                    print("SECOND: Unable to find perturbation. Try changing perturbed LHS")
                    exit(-1)
            eq_sL.append(tmpL)

            self.inf_sol = self.inf_sol.doit().subs(self.a[j], self.coeffs_inf["a" + str(j)].args[0])
            pprint(self.coeffs_inf)

        self.inf_sol = self.inf_sol.subs(self.e, 1)
        self.sol_poly = sp.Eq(self.x, self.inf_sol.n(20))
        return self.sol_poly

    def __repr__(self):
        pprint(sp.Eq(self.pLHS, 1))
        return ""

    def __str__(self):
        pprint(sp.Eq(self.pLHS, 1))
        return ""
    
    def print_inf_series(self, expanded=False):
        if not expanded: 
            pprint(sp.Eq(self.x, self.inf_series))
        else: 
            pprint(sp.Eq(self.x, self.inf_series.doit().expand()))

class PTranscendentalSolver(PPolySolver):
    '''Perturbation Transcendental Solver'''
    def __init__(self, x, e, perturbed_LHS, power_of_epsilon, power_of_taylor):
        PPolySolver.__init__(self, x, e, perturbed_LHS, power_of_epsilon)
        self.power_of_taylor = power_of_taylor
        self.pLHS = self.convert_taylor()
    
    def convert_taylor(self):
        return taylor(self.pLHS, self.x, self.power_of_taylor)
