"""
sympy solution of problem 3.8 from book
"""

from sympy import *
from sympy import Rational as r
c1,c2,c3,t,cl = symbols('c1,c2,c3,t,cl')
c = [r(1./4), r(1./2), r(3./4)]
q = (t-c1)*(t-c2)*(t-c3)
q1 = q/(t-c1)
q2 = q/(t-c2)
q3 = q/(t-c3)

#for each j stage compute all a_(j,i)
def a(j,i,g,x,c):
    d = g.evalf(subs={x:x})
    f = g/( d)
    return integrate(f, (t, 0, c[j]) )

def b(j,g,x,c):
    d = g.evalf(subs={x:x})
    f = g/( d)
    return integrate(f, (t, 0, 1) )
    
    
A = {}
B = {}
Q = [q1,q2,q3]
for j in range(3):
    for i in range(3):
        A['A_({},{})'.format(j,i)]= a(j,i,Q[i],c[i],c)
        B['A_({},{})'.format(j,i)]= b(j,Q[j],c[j],c)
        
        
print A
print B
        

