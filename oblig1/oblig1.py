"""
solver for oblig1 in mat260
"""
import numpy as np
import matplotlib.pylab as plt
from parameters import *

def RKMethod(b,y,f,M,i,h,t):
    """
    Embeded Runge-Kutta method.
    """
    c1, c2, c3, c4, c5 = 0,1./3,1./3,1/2.,1.
    for n in range(i, i+1):
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
            #print "{}".format(i)
    print " "
    print "{}".format(i)
    print y
    return y
     
     
def plot_solutions(method,b_higher,b_lower,y,f,t,N,M,h,T,c,tol,se,ie):
    #e_h,i_h,r_h,N_h,e_v,i_v,N_v = RKstepSizeAlgorithm(method,b_higher,b_lower,T,c,tol,y,f,M,h,t)
    s,i = RKstepSizeAlgorithm(method,b_higher,b_lower,y,f,t,N,M,h,T,c,tol,se,ie)
    plt.plot(t,i,t,ie)
    #plt.legend(['{} exact '.format(sol),'{} numerical'.format('0')], loc='best')
    #plt.title('exact solution vs numerical solution for h = {}, T= {}'.format(h,T))
    plt.ylim([-0.2,1.2])
    #plt.savefig("latex_files/figures/%s.png"%fig, transparent = True)
    plt.show()
        
            
def RKstepSizeAlgorithm(method,b_higher,b_lower,y,f,t,N,M,h,T,c,tol,se,ie):
    """
    variable step size algorithm for Runge-Kutta method
    questions: whuch types of norm and how
    to compute the norm?
    is N_v and N_h constants?
    what are the initial conditions?
    """ 
    t0 = 0
    p = 3;
    N = int(T/h)
    print " N ", N
    for i in range(0,N):
        Y = method(b_higher,y,f,M,i,h,t)
        X = method(b_lower,y,f,M,i,h,t) 
        
        #K = max(error)
        #print " error in for ", error
        #while(K > (tol*h)/(T-t0) ): 
            #h = c*h*( tol/((T-t0)*K) )**(1./(p+1))
            #Y = method(b_higher,y,f,M,i,h,t)
            #X = method(b_lower,y,f,M,i,h,t)
            #error = np.sqrt(e1**2 + e2**2 + e3**2 + e4**2 + e5**2 + e6**2 + e7**2)
            #K = max(error)  
        #if K!=0:  
            #h = h*( tol/((T-t0)*K) )**(1./(p+1))
    return X
        
if __name__ == "__main__":
    h = 0.1
    #y,f,t,N,M,h,T = MalariaModel()
    y, f, se, ie, t, N, M, h,T = SIModel(h)
    b_higher = [1/6.,0,0,2./3,1/6.]
    b_lower = [1/2.,0,-3./2,2,0]
    c,tol = 0.5, 0.01
    plot_solutions(RKMethod,b_higher,b_lower,y,f,t,N,M,h,T,c,tol,se,ie)
    #plot_solutions(RKMethod,b_higher,b_lower,y,f,t,N,M,h,T,c,tol) 
              
             
        
        
        
        
        
        
        
        
        
        
        