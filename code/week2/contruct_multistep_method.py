
"""
Given the coefficients a0, a1, ..., am
and b0, b1, ..., bm, use the order condition
to construct a s-step methods of order p that converges.
Furtheremore check for convergence
"""
from sympy import *

def am(s):
    """
    return am array
    """
    return ['a_{}'.format(i) for i in range(s+1)]
    
    
def bm(s):
    """
    return bm array
    """
    return ['b_{}'.format(i) for i in range(s+1)]
    
    
def solve_system(am,bm,*args):
    """
    construct the k systems, 
    k = 1, ..., p from the order 
    condition where p is the order of the multi
    step method
    """
    # make the a's and b's symbolic
    p = 3
    y = " ".join(am(*args))
    x = " ".join(bm(*args))
    am = symbols(y)
    bm = symbols(x)
    sys = []
    for k in range(1, p+1):
        Am = [am[i]*(i**k) for i in range(len(am)) ]
        Bm = [bm[i]*(i**(k-1)) for i in range(len(bm)) ]
        sys.append(sum(Am)-k*sum(Bm))
    
    #variables = am+bm
    equations = [Eq(equation) for equation in sys]
    solutions = solve(equations)
    return solutions
    
    


s = 3
print solve_system(am,bm,s)