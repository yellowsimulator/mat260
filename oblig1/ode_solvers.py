"""
implemenation of common ode methods.
Usage of the solvers:
N=int(T/h)
t = np.linspace(t0, T, N+1)
yi = np.zeros(N+1), i = 0,1,...,N

yi[t0] = y0, i = 0,1,...,N
y = np.array([y1,y2,...,yn])
def f(t,y):
    f = np.array([f1,f2,...,fn])
    #where fi = fi(y[0],y[1],...,y[N]), i = 0,1,...,N
    return f

for example call EulersMethod like this:
M = 2 # system of two equations
T = 10
h = 0.1
N=int(T/h)
u = EulersMethod(y,f,M,N,h,t)
u1, u2 = u[0], u[1]
file-line-error

figure for latex

\begin{figure}
\centering
\begin{python}
#
from pyx import *

g = graph.graphxy(width=8)
g.plot(graph.data.function("y(x)=sin(x)/x", min=-15, max=15))
g.writePDFfile("function")
print r'\includegraphics{function}'
\end{python}
\caption{$y(x)=\frac{\sin(x)}{x}$}
\end{figure}


"""
from sympy import *
import numpy as np
from numpy.linalg import inv

def Jacobian(f):
    """
    compute the jakobian matrix
    """
    pass

def EulersMethod(y,f,M,N,h,t):
    """
    Explicite Euler's method for the system:
    y' = f(t,y); y(0) = y0, given by
    y[n+1] = y[n] + h*f(t[n], y[n]) where
    h is the step size.
    M is the number of equation in the system
    N is the number of grid points
    """
    for n in range(0, N):
        for m in range(0,M):
            y[m][n+1] = y[m][n] + h*f(t[n],y)[m][n]
    return y


def ERK2Method(y,f,M,N,h,t):
    """
    Explicite order 2 Runge-Kutta method for the system:
    y' = f(t,y); y(0) = y0, given by
    k1 = h*f(t[n], y[n])
    k2 = h*f(t[n], y[n] + 0.5*k1)
    y[n+1] = y[n] + h*k2.
    M is the number of equation in the system
    N is the number of grid points
    """ 
    for n in range(0, N):
        for m in range(0, M):
            k1 = f(t[n],y)
            y2 = y + k1*0.5*h
            k2 = f(t[n],y2)[m][n]
            y[m][n+1] = y[m][n] + h*k2
    return y

 
def ERK4Method(y,f,M,N,h,t):
    """
    Explicite order 4 Runge-Kutta method for the system
    y' = f(t,y); y(0) = y0, given by
    k1 = h*f(t[n], y[n])
    k2 = h*f(t[n], y[n] + 0.5*k1)
    k3 = h*f(t[n], y[n] + 0.5*k2)
    k4 = h*f(t[n], y[n] + k3)
    y[n+1] =  (1/6)*(k1 + 2*k2 + 2*k3 + k4).
    M is the number of equation in the system
    N is the number of grid points
    """ 
    for n in range(0, N):
        for m in range(0, M):
            k1 = f(t[n],y)
            y1 = y + 0.5*h*k1
            k2 = f(t[n] + h/2.,y1)
            y2 = y + 0.5*h*k2
            k3 = f(t[n] + h/2.,y2)
            y3 = y + h*k3
            k4 = f(t[n] + h/2.,y3)[m][n]
            y[m][n+1] = y[m][n] + (1./6)*h*(k1[m][n] + 2*k2[m][n] + 2*k3[m][n] + k4)
    return y



def YRKMethod(y,f,M,N,h,t):
    """
    an other RK method
    """
    c1, c2, c3, c4, c5 = 0,1./3,1./3,1/2.,1.
    b = [1/6.,0,0,2./3,1/6.]
    #b = [1/2.,0,-3./2,2,0]
    for n in range(0, N):
        for m in range(0, M):
            k1 = f(y)
            y1 = y + (1./3)*h*k1
            k2 = f(y1)
            y2 = y + (1./6)*h*k1 + (1./6)*h*k2
            k3 = f(y2)
            y3 = y + (1./8)*h*k1 + (0.)*h*k2 + (3./8)*h*k3
            k4 = f(y3)
            y4 = y + (1./2)*h*k1 + 0.*h*k2 + (-3./2)*h*k3 + 2*h*k4
            k5 = f(y4)
            y[m][n+1] = y[m][n]+h*(b[0]*k1[m][n] +\
            b[1]*k2[m][n]+b[2]*k3[m][n]+b[3]*k4[m][n]+b[4]*k5[m][n])
    return y
    
    

    
    
    
    
    
    
    
    
    
