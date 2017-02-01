
import numpy as np
import matplotlib.pylab as plt
import sys
from scipy.integrate import odeint
from math import log
    
def EulersMethod(y,f,M,N,h,t):
    """
    General eulers method
    M is the number of variables defining the system
    N is the number of grid points
    """
    for n in range(0, N):
        for m in range(0,M):
            y[m][n+1] = y[m][n] + h*f(t[n],y)[m][n]
    return y
            
        
def call_solver(h, M, T):
    """
    define your sytem here
    """
    N=int(T/h)
    t = np.linspace(0, T, N+1)
    x1 = np.zeros(N+1);x2 = np.zeros(N+1)
    x3 = np.zeros(N+1);y1 = np.zeros(N+1)
    y2 = np.zeros(N+1)
    x1[0] = 1; x2[0] = 0; x3[0] = 0; y1[0] = 1; y2[0] = 1
    y = np.array([x1,x2,x3,y1,y2])
    
    def f(t,y):
        f1=y[1];f2=y[2];f3=t*t*y[0]*y[2] - y[0]*y[4]
        f4=y[4];f5=t*y[3]*y[4]+4*y[1]
        return np.array([f1,f2,f3,f4,f5])
        
    return EulersMethod(y,f,M,N,h,t),t
 
  
def test_convergence_rate():
    """
    test convergence rate for Eulers method
    for a simple problems:
    y' = -2y
    y(0) = 1
    y is the exact solution while ya is the approximate
    solution
    """
    E , E_half, p = [], [], []
    M=1; T=1
    H = [0.1*2**(-i) for i in range(10)]
    for i,h in enumerate(H):
        N = int(T/h)
        t = np.linspace(0,T,N+1)
        x1 = np.zeros(N+1)
        x1[0] = 1
        y = np.array([x1])
        
        def f(t,y):
            f1 = -2*y[0]
            return np.array([f1])
        
        def ye(t):
            """
            exact solution
            """
            return np.exp(-2*t)
    
        ya = EulersMethod(y,f,M,N,h,t)
        z1 = abs(ya[0]-ye(t))
        E.append(max(z1))
        
    for i in range(len(E)):
        p.append( log(E[i-1]/E[i])/log(H[i-1]/H[i])   )
        
    #plt.plot(range(len(E)), p)
    #plt.title("convergence rate p with error")
    #plt.xlabel("error")
    #plt.ylabel("convergence rate p")
    #plt.show()
    return [round(r,2) for r in p]
    
    
    

            
def plot_solution():
    """
    plot solution for different time step h
    """
    M = 5; T = 1
    H = [0.1,0.01,0.001]
    for h in H:
        Y,t = call_solver(h,M,T)
        u,v = Y[0], Y[3]
        
        N = int(1/h)
        plt.subplot(211)
        plt.plot(t,u)   
        plt.title("u as a function of t for different values of h", color='blue')
        plt.subplot(212)
        plt.title("v as a function of t for different values of h", color='blue')
        plt.plot(t,v)
    plt.legend(['h={}'.format(H[0]),'h={}'.format(H[1]),'h={}'.format(H[2])] ,loc='best')    
    plt.show()


    
if __name__ == "__main__":
    #plot_solution()
    p = test_convergence_rate()
    print " convergence rate ", p
    
    
    
    