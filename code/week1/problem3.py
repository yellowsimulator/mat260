import numpy as np
import matplotlib.pylab as plt
import sys
from math import sin, cos, sqrt
from mpl_toolkits.mplot3d import Axes3D
def eulers_method(h,sz):
    """
    solve the system of ode:
    
    xz' = sz*(R*cos(w*t) - xz)*1/alpha
    yz' = sz*(R*sin(w*t) - yz)*1/alpha
    alpha = sqrt((R*cos(w*t)-xz)**2 + (R*sin(w*t)-yz)**2)
    xz(0) = a
    yz(0) = -a
    0< a < sqrt(0.5)
    by Euler's method
    Y[n+1] = Y[n] + h*F[n]
    """
    T = 8*np.pi # simulation time t in [0, T=1]
    N = int(T/h) 
    t = np.linspace(0, T, N+1)
    xz = np.zeros(N+1)
    yz = np.zeros(N+1)
    
    
    #initial conditions
    a = sqrt(0.5)/2
    w = 1; R = 1
    xz[0] = a; yz[0] = -a
    
    for n in range(0, N):
        alpha = sqrt( (R*cos(w*t[n])-xz[n])**2 + (R*sin(w*t[n])-yz[n])**2   )
        xz[n+1] = xz[n] + ((sz*h)/alpha)*(R*cos(w*t[n])-xz[n])
        yz[n+1] = yz[n] + ((sz*h)/alpha)*(R*sin(w*t[n])-yz[n])
        
    return xz, yz, t
    


def main(sz):
    """
    call solver and plot solutions
    """
    
    H = [0.001]
    zz = 0
    a = sqrt(0.5)/2.
    for h in H:
        xz, yz, t = eulers_method(h,sz)
        N = int(8*np.pi/h)
        plt.subplot(211)
        plt.plot(t,xz)   
        plt.title("xz,  sz={} a = {}".format(sz, a), color='blue')
        plt.subplot(212)
        plt.title("yz", color='blue')
        plt.plot(t,yz)
    #plt.legend(['h={}'.format(H[0])] ,loc='best')    
    
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(xz, yz, zz, label='Trajectory')
    plt.title(" Zombi path in the plane z = 0")
    plt.xlabel("xz")
    plt.xlabel("yz")
    plt.xlabel("z")
    #plt.plot(xz,yz)
    plt.show()
    
    
    
if __name__ == "__main__":
    sz = float(sys.argv[1])
    main(sz)