"""
solver for oblig1 in mat260
"""
import numpy as np
import matplotlib.pylab as plt
from parameters import *
from math import sqrt
def RKMethod(b,y,f,M,n,h,t):
    """
    Embeded Runge-Kutta method.
    """
    c1, c2, c3, c4, c5 = 0,1./3,1./3,1/2.,1.
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
     
                 
def RKstepSizeAlgorithm(f,c,tol,y,h,M):
    """
    variable step size algorithm for Runge-Kutta method
    questions: whuch types of norm and how
    to compute the norm?
    is N_v and N_h constants?
    what are the initial conditions?
    """ 
    
    b_higher = [1/6.,0,0,2./3,1/6.]
    b_lower = [1/2.,0,-3./2,2,0]
    for n in range(0,N):
        #h = t[n+1]-t[n]
        Y = RKMethod(b_higher,y,f,M,n,h,t)
        X = RKMethod(b_lower,y,f,M,n,h,t)
       
        #update
        y, x = Y, X
    error = Y-X
    print " error ", error
    return Y
    
    
    
def plot_solutions(solution,t):
    u1,u2,u3,u4,u5,u6,u7 = solution
    plt.plot(t,u1,t,u2)
    #plt.legend(['{} exact '.format(sol),'{} numerical'.format('0')], loc='best')
    #plt.title('exact solution vs numerical solution for h = {}, T= {}'.format(h,T))
    plt.ylim([-0.2,1.2])
    #plt.savefig("latex_files/figures/%s.png"%fig, transparent = True)
    plt.show()
        
        
if __name__ == "__main__":
    T = 100; 
    
    h = 0.01
    N = int(T/h)
    t = np.linspace(0,T,N+1)
    
    x1 = np.zeros(N+1);x2 = np.zeros(N+1)
    x1[0] = 0.9; x2[0] = 0.1
    y1 = np.array([x1,x2])
    
    
    
    e_h=np.zeros(N+1);i_h=np.zeros(N+1);r_h=np.zeros(N+1);N_h = np.zeros(N+1);e_v=np.zeros(N+1)
    i_v=np.zeros(N+1);N_v = np.zeros(N+1)
    N_h[0] = 400.+10.+30.+0 
    N_v[0] = 1000.+100.+50.
    e_h[0] = 10./(400+10+30+0)
    i_h[0] = 10./(400+10+30+0)
    r_h[0] = 0
    e_v[0] = 100./(1000+100+50)
    i_v[0] = 10./(1000+100+50)
    y = np.array([e_h,i_h,r_h,N_h,e_v,i_v,N_v])
    
    print y
    #f = MalariaModel()
    def f1(t,y1):
        beta = 0.5
        f1 = -beta*y1[0]*y1[1]
        f2 = beta*y1[0]*y1[1]
        return np.array([f1,f2])
        
        
    
    def f2(t,y):
        " function f"
        #print " t = ",t
        sigma_v = 0.6; sigma_h = 18.; beta_hv = 2.*1e-2; nu_h = 8.333*1e-2; phi_h = 7.666*1e-5
        l_h = 3.285*1e-2; delta_h = 3.454*1e-4 ; gamma_h = 3.704*1e-3; rho_h = 1.460*1e-2
        mu_1h = 0.1429;  mu_2h = 1e-7; beta_bar_vh =  8.333*1e-2 ; beta_vh = 2.*1e-2
        nu_v = 0.1; phi_v = 0.4; mu_1v =0.1429 ;mu_2v = 1e-7
        e_h = y[0];i_h = y[1];r_h = y[2];N_h = y[3];e_v = y[4];i_v = y[5];N_v = y[6]
        #f2,f3,f4,f5,f6,f7
        N = len(e_h)-1; f1 = np.zeros(N+1);f2 = np.zeros(N+1);f3 = np.zeros(N+1);f4 = np.zeros(N+1)
        f5 = np.zeros(N+1);f6 = np.zeros(N+1);f7 = np.zeros(N+1)
        
        f1[t] = ((sigma_v*sigma_h*N_v[t]*beta_hv*i_v[t])/(sigma_v*N_v[t]+sigma_h*N_h[t]))*\
             (1-e_h[t]-i_h[t]-r_h[t])-(nu_h+phi_h+(l_h/N_h[t]) )*e_h[t] + delta_h*i_h[t]*e_h[t]
        f2[t] = nu_h*e_h[t]-(gamma_h+sigma_h+phi_h+(l_h/N_h[t]))*i_h[t]+sigma_h*i_h[t]*i_h[t]
        f3[t] = gamma_h*i_h[t]-(rho_h+phi_h+(l_h/N_h[t]))*r_h[t]+sigma_h*i_h[t]*r_h[t]
        f4[t] = l_h+phi_h*N_h[t]-(mu_1h+mu_2h*N_h[t])*N_h[t]-sigma_h*i_h[t]*N_h[t]
        f5[t] =( (sigma_h*sigma_h*N_v[t])/(sigma_v*N_v[t]+sigma_h*N_h[t]) )*(beta_vh*i_h[t]+beta_bar_vh*r_h[t])-\
               (nu_v+phi_v)*e_v[t]
        f6[t] = nu_v*e_v[t]-phi_v*i_v[t]
        f7[t] = phi_v*N_v[t]-(mu_1v+mu_2v*N_v[t])*N_v[t]
        ft = np.array([f1,f2,f3,f4,f5,f6,f7])
        #print "f7 :",f7
        return ft
    

        
    
    M = 7 
    b_higher = [1/6.,0,0,2./3,1/6.]
    b_lower = [1/2.,0,-3./2,2,0]
    c,tol = 0.5, 0.01
    #sol = RKstepSizeAlgorithm(f1,c,tol,y1,h,M)
    
    sol = RKstepSizeAlgorithm(f2,c,tol,y,h,M)
    plot_solutions(sol,t)
    #plot_solutions(RKMethod,b_higher,b_lower,y,f,t,N,M,h,T,c,tol) 
              
             
        
        
        
        
        
        
        
        
        
        
        