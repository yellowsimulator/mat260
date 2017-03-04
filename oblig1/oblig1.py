"""
solver for oblig1 in mat260
"""
import numpy as np
import matplotlib.pylab as plt
from parameters import *

def RKMethod(b,y,f,M,h,t):
    """
    Embeded Runge-Kutta method.
    """
    y_new = np.zeros(M)
    c1, c2, c3, c4, c5 = 0,1./3,1./3,1/2.,1.
    for m in range(0, M):
        k1 = f(t,y)
        y1 = y + (1./3)*h*k1
        k2 = f(t + c1*h,y1)
        y2 = y + (1./6)*h*k1 + (1./6)*h*k2
        k3 = f(t + c2*h,y2)
        y3 = y + (1./8)*h*k1 + (0.)*h*k2 + (3./8)*h*k3
        k4 = f(t + c3*h,y3)
        y4 = y + (1./2)*h*k1 + 0.*h*k2 + (-3./2)*h*k3 + 2*h*k4
        k5 = f(t + c4*h,y4)
        y_new[m] = y[m]+h*(b[0]*k1[m] +\
        b[1]*k2[m]+b[2]*k3[m]+b[3]*k4[m]+b[4]*k5[m])
    #print y_new
    return y_new
     
     
def plot_solutions(time, u):
    #beta = 0.5
    e_h,i_h,r_h,N_h,e_v,i_v,N_v = u[0],u[1],u[2],u[3],u[4],u[5],u[6]
    time = np.array(time)
    #ie = 1./(1+9*np.exp(-beta*time))
    #se = 1-ie
    #s,i = solution[0],solution[1]
    #x = (-1/beta)*np.log( ((1./i[2]) -1)/9  )
    #su = np.array(s)+np.array(i)
    #print len(e_h)
    plt.plot(time,N_h)
    #plt.xlim([0,1])
    plt.show()
        
            
def RKstepSizeAlgorithm(method,M,y,f,h,c,tol):
    """
    variable step size algorithm for Runge-Kutta method
    questions: whuch types of norm and how
    to compute the norm?
    is N_v and N_h constants?
    what are the initial conditions?
    """ 
    time = []
    solutionY = [[]*i for i in range(1,M+1)]
    solutionX = [[]*i for i in range(1,M+1)]
    
    #initialise solutions
    S_h0 = 400.; E_h0 = 10.; I_h0 = 30.; R_h0 = 0.; S_v0 = 1000.; E_v0 = 100.; I_v0 = 50.
    N_h0 = S_h0 + E_h0 + I_h0 + R_h0; N_v0 = S_v0 + E_v0 + I_v0  
    e_h0 = E_h0/N_h0; i_h0 = I_h0/N_h0; r_h0 = R_h0/N_h0; e_v0 = E_v0/N_v0; i_v0 = I_v0/N_v0
    
    solutionY[0].append(e_h0); solutionY[1].append(i_h0); solutionY[2].append(r_h0)
    solutionY[3].append(N_h0); solutionY[4].append(e_v0); solutionY[5].append(i_v0)
    solutionY[6].append(N_v0)
    solutionX[0].append(e_h0); solutionX[1].append(i_h0); solutionX[2].append(r_h0)
    solutionX[3].append(N_h0); solutionX[4].append(e_v0); solutionX[5].append(i_v0)
    solutionX[6].append(N_v0)
    
    p = 4; T = 50; n = 0; t0 = 0; t = 0; n1 = 0; h_new = h
    while (t < T):
        t = (n+1)*h_new
        y = method(b_higher,y,f,M,h_new,t)
        x = method(b_lower,y,f,M,h_new,t)
        error = np.array(y)-np.array(x)
        K = np.sqrt(error[0]**2 + error[1]**2 + error[2]**2 +\
        error[3]**2 +error[4]**2 + error[5]**2+ error[6]**2 )
        #print " error K : ", K
        """
        while(K > (tol*h)/(T-t0) ): 
            K_min = (tol*h)/(T-t0)
            print " in inner loop K: ", K
            print " in inner look K_min: ", K_min
            h_new = c*h*( tol/((T-t0)*K) )**(1./(p+1))
            t1 = (n1+1)*h_new
            time.append(t1)
            y = method(b_higher,y,f,M,h_new,t1)
            x = method(b_lower,y,f,M,h_new,t1)
            print " e_v ", x[0]
            error = np.array(y)-np.array(x)
            K = np.sqrt(error[0]**2 + error[1]**2 + error[2]**2 +\
            error[3]**2 +error[4]**2 + error[5]**2+ error[6]**2 )
            solutionX[0].append(x[0]);solutionX[1].append(x[1]);solutionX[2].append(x[2])
            solutionX[3].append(x[3]);solutionX[4].append(x[4]);solutionX[5].append(x[5])
            solutionX[6].append(x[6])
            
            n1+=1
            print " time t in inner loop :\n" , t1
            if t1 >=T:
                break
        """
        #if n1 < n:
        solutionX[0].append(x[0]);solutionX[1].append(x[1]);solutionX[2].append(x[2])
        solutionX[3].append(x[3]);solutionX[4].append(x[4]);solutionX[5].append(x[5])
        solutionX[6].append(x[6])
        time.append(t)
            
        h_new = h*( tol/((T-t0)*K) )**(1./(p+1))
        #n = n1-1
        #increment
        n+=1;     
    time.append(T)
    return time,solutionX
        
if __name__ == "__main__":
    h = 0.01
    y,f = MalariaModel()
    #y, f = SIModel()
    M = len(y)
    method = RKMethod
    b_higher = [1/6.,0,0,2./3,1/6.]
    b_lower = [1/2.,0,-3./2,2,0]
    c = 0.1
    tol = 0.025
    time, u = RKstepSizeAlgorithm(method,M,y,f,h,c,tol)
    plot_solutions(time, u)
             
        
        
        
        
        
        
        
        
        
        
        