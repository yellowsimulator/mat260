"""
implement the Runge-Kutta method.
input vector c, b and A, with yn and f
"""
import numpy as np
def runge_kutta(b):
    """
    usage
    """
    mu = len(b)
    #compute ksi
    """
    for n in range(0, N):
        for m in range(0, M):
            for i in range(1,mu):
                ksi[m][i] = y[m][n] + h*sum(a[r][i]*f(t[n]+c[i]*h,ksi[m][i]))
    
    for n in range(0, N):
        for m in range(0,M):
            for j in range(1,mu+1)
                y[m][n+1] = y[m][n] + h*sum(b[j]*f(t[n]+c[j]*h,ksi[j]))[m][n]
    return y

    """
    print b
    
b = np.array([1./6,1/3.,1/3.,1/6.])
runge_kutta(b)