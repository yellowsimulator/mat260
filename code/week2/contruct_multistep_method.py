
"""
Given the coefficients a0, a1, ..., am
and b0, b1, ..., bm, use the order condition
to construct a s-step methods.
Furtheremore check for convergence
"""

from sympy import *
def am(s):
    return ['a_{}'.format(i) for i in range(s+1)]
    
sum_am = am(3)
print sum_am