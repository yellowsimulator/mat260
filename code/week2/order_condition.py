"""
implement the order condition
and return the order of the method
"""
from math import factorial


def first_condition(p,am,bm):
    """
    check first conditon
    """
    if (float(sum(am))==0):
        return True
    else:
        return False

def second_condition(p,am,bm):
    """
    check_first condition
    """
    values = []
    for k in range(1, p+1):
        Am = [am[i]*(i**k) for i in range(len(am)) ]
        Bm = [bm[i]*(i**(k-1)) for i in range(len(bm)) ]
        condition1 = float(sum(Am)-k*sum(Bm))
        values.append(condition1)
    
    generator = [x == 0. for x in values]
    return any(generator)
    
    
def third_condition(p,am,bm):
    """
    check second condition
    """
    values = []
    A1m = [am[i]*(i**(p+1)) for i in range(len(am)) ]
    B1m = [bm[i]*((i**p)) for i in range(len(bm)) ]
    condition2 = float(sum(A1m)-(p+1)*sum(B1m))
    if condition2 == 0:
        return False
    else:
        return True
    
        
def return_order(am,bm):
    """
    return order if it satisfies 
    the order condition
    """
    order = range(1, 11)
    for p in order:
        C1 = first_condition(p,am,bm)
        C2 = second_condition(p,am,bm)
        C3 = third_condition(p,am,bm)
        if C1 and C2 and C3:
            return p
            
    
    
am = [-1,0,0,1]
bm = [3/8., 9/8.,9/8.,3/8.]

print "order of convergence: ", return_order(am,bm)
