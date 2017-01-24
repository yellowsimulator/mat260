
import numpy as np
import matplotlib.pylab as plt
import sys
from scipy.integrate import odeint

def eulers_method(h):
    """
    solve the system of ode:
    
    u''' = t*t*u*u'' - uv' 
    v''  = tv*v' + 4*u'
    u(0) = 1; u'(0) = 0; u''(0) = 0
    v(0) = v'(0) = 1
    
    by Euler's method
    Y[n+1] = Y[n] + h*F[n]
    """
    T = 1 # simulation time t in [0, T=1]
    N = int(T/h) 
    t = np.linspace(0, T, N+1)
    x1 = np.zeros(N+1);x2 = np.zeros(N+1)
    x3 = np.zeros(N+1);y1 = np.zeros(N+1)
    y2 = np.zeros(N+1)
    
    #initial conditions
    x1[0] = 1; x2[0] = 0; x3[0] = 0; y1[0] = 1; y2[0] = 1
    for n in range(0, N):
        x1[n+1] = x1[n] + h*x2[n]
        x2[n+1] = x2[n] + h*x3[n]
        x3[n+1] = x3[n] + h*(t[n]*t[n]*x1[n]*x3[n] - x1[n]*y2[n])
        y1[n+1] = y1[n] + h*y2[n]
        y2[n+1] = y2[n] + h*(t[n]*y1[n]*y2[n] + 4*x2[n])
    return x1, y1, t
    

def main():
    """
    call solver and plot solutions
    """
    
    H = [ 0.25, 0.1, 0.001]
    for h in H:
        u, v, t = eulers_method(h)
        N = int(1/h)
        plt.subplot(211)
        plt.plot(t,u)   
        plt.title("u as a function of t for different values of h", color='blue')
        plt.subplot(212)
        plt.title("v as a function of t for different values of h", color='blue')
        plt.plot(t,v)
    plt.legend(['h={}'.format(H[0]),'h={}'.format(H[1]),'h={}'.format(H[2])] ,loc='best')    
    plt.show()
    
# a class version
class EulersMethod(object):
    """
    class version: in progress
    """
    
    def __init__(self,y0,y,f,h,T,t,dim):
        self.y0 = np.array(self.dim)
        self.y = np.array(self.dim)
        self.f = f
        self.h = h
        self.T = T
        self.t = t
        
    def solver(self):
        """
        general solver for Euler's method.
        """ 
        N = int(self.T/self.h)
        for n in range(0,n):
            self.y[n+1] = self.y[n] + self.h*self.f(y[n],t[n])
        return y


    
if __name__ == "__main__":
    #h = float(sys.argv[1])
    main()