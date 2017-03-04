"""
solver for oblig1 in mat260
"""
import numpy as np
def RKMethod(b,y,f,M,N,h,t):
    """
    Lower oder RK method
    """
    for n in range(0, N):
        for m in range(0, M):
            k1 = f(t[n],y)
            y1 = y + (1./3)*h*k1
            k2 = f(t[n]+c1*h,y1)
            y2 = y + (1./6)*h*k1 + (1./6)*h*k2
            k3 = f(t[n]+c2*h,y2)
            y3 = y + (1./8)*h*k1 + (0.)*h*k2 + (3./8)*h*k3
            k4 = f(t[n]+c3*h,y3)
            y4 = y + (1./2)*h*k1 + 0.*h*k2 + (-3./2)*h*k3 + 2*h*k4
            k5 = f(t[n]+c4*h,y4)
            y[m][n+1] = y[m][n]+h*(b[0]*k1[m][n] +\
            b[1]*k2[m][n]+b[2]*k3[m][n]+b[3]*k4[m][n]+b[4]*k5[m][n])
    return y
     

        
            
def RKstepSizeAlgorithm(t0,T,c,tol,h,y,f,M,h,t):
    """
    variable step size algorithm
    """ 
    N = int(T/h);
    # higher order method coefficients
    b_higher = [1/6.,0,0,2./3,1/6.]
    #lower order coefficients
    b_lower = [1/2.,0,-3./2,2,0]
    K = max(np.abs(RK_higer-RK_lower))
    t0 = 0
    while(t<T):
        N = int(T/h)
        u = RKMethod(b_higher,y,f,M,N,h,t)
        v = RKMethod(b_lower,y,f,M,N,h,t) 
        K = max(np.abs(RK_higer-RK_lower))
        
        while(K > (tol*h)/(T-t0) ):
             h = c*h*( tol/((T-t0)*K) )**(1./(p+1))
             u = RKMethod(b_higher,y,f,M,N,h,t)
             v = RKMethod(b_lower,y,f,M,N,h,t) 
             K = max(np.abs(RK_higer-RK_lower))
             
        h = h*( tol/((T-t0)*K) )**(1./(p+1))
        t0 += h
    return u,v
        
             
             
        
        
        
        
        
        
        
        
        
        
        