"""
define parameters for a specific problem here
example: set vectors y, f for y' = f, with all parameters
"""
import numpy as np

    
def MalariaModel():
    """
    set parameters for the malaria model
    """

    def f(t,y):
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
        print "f7 :",f7
        return ft
    
    #f(0,y)
    return f
    
    
def SIModel(y):
    beta = 0.5
    f1 = -beta*y[0]*y[1]
    f2 = beta*y[0]*y[1]
    return np.array([f1,f2])   
    
    
if __name__ == "__main__":
    print " Ok"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    