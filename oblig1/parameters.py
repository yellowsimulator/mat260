"""
define parameters for a specific problem here
example: set vectors y, f for y' = f, with all parameters
"""
import numpy as np
def SIModel():
    """
    setting parameters for SI
    spreading of deseases model:
    (s,i)' = (f1,f2) where 
    f1 = -beta*s(t)*i(t)
    f2 = beta*s(t)*i(t)
    """
    #vector solution
    beta = 0.5
    #s = np.array([]); i = np.array([])
    s = 0.9; i = 0.1 # x1 = s, x2 = i
    y = np.array([s,i])
    #setting function f = (f1,f2)
    def f(t,y):
        s,i = y[0],y[1]
        beta = 0.5
        f1 = -beta*s*i # f for s
        f2 = beta*s*i # f for i
        return np.array([f1,f2])

    # setting analytical solutions
    #ie = 1./(1+9*np.exp(-beta*t))
    #se = 1-ie
    return y, f
    

def SIRModel(h,sigma,beta):
    """
    setting parameters for SI
    spreading of deseases model:
    """
    gamma = float(beta/sigma)
    M=3; T = 50
    N = int(T/h)
    t = np.linspace(0,T, N+1)
    s = np.zeros(N+1);i = np.zeros(N+1); r = np.zeros(N+1)
    s[0] = 0.8; i[0] = 0.2 ; r[0] = 1-s[0]-i[0]
    y = np.array([s,i,r])
    #setting function f = (f1,f2)
    def f(t,y):
        s = y[0]; i = y[1]; r = y[2] 
        f1 = -beta*s*i 
        f2 = beta*s*i - gamma*i
        f3 = gamma*i 
        return np.array([f1,f2,f3])
    return y, f, t, N, M, h,T
    
    

def MalariaModel():
    
    
    #initial conditions
    sigma_v = 0.6; sigma_h = 18.; beta_hv = 2.*1e-2; nu_h = 8.333*1e-2; phi_h = 7.666*1e-5
    l_h = 3.285*1e-2; delta_h = 3.454*1e-4 ; gamma_h = 3.704*1e-3; rho_h = 1.460*1e-2
    mu_1h = 4.212*1e-5;  mu_2h = 1e-7; beta_bar_vh =  8.333*1e-2 ; beta_vh = 0.8333
    nu_v = 0.1; phi_v = 0.4; mu_1v =0.1429 ;mu_2v = 2.279*1e-4
    S_h0 = 400.; E_h0 = 10.; I_h0 = 30.; R_h0 = 0.; S_v0 = 1000.; E_v0 = 100.; I_v0 = 50.
    N_h0 = S_h0 + E_h0 + I_h0 + R_h0; N_v0 = S_v0 + E_v0 + I_v0 
    N_h = N_h0; N_v = N_v0 
    e_h = E_h0/N_h0; i_h = I_h0/N_h0; r_h = R_h0/N_h0; e_v = E_v0/N_v0; i_v = I_v0/N_v0
    y = np.array([e_h,i_h,r_h,N_h,e_v,i_v,N_v])
    def f(t,y):
        f1 = ((sigma_v*sigma_h*N_v*beta_hv*i_v)/(sigma_v*N_v+sigma_h*N_h))*\
             (1-e_h-i_h-r_h)-(nu_h+phi_h+(l_h/N_h) )*e_h + delta_h*i_h*e_h
        f2 = nu_h*e_h-(gamma_h+delta_h+phi_h+(l_h/N_h))*i_h+delta_h*i_h*i_h
        f3 = gamma_h*i_h-(rho_h+phi_h+(l_h/N_h))*r_h+delta_h*i_h*r_h
        f4 = l_h+phi_h*N_h-(mu_1h+mu_2h*N_h)*N_h-delta_h*i_h*N_h
        f5 =( (sigma_v*sigma_h*N_v)/(sigma_v*N_v+sigma_h*N_h) )*(beta_vh*i_h+beta_bar_vh*r_h)*(1-e_v-i_v)-\
               (nu_v+phi_v)*e_v
        f6 = nu_v*e_v-phi_v*i_v
        f7 = phi_v*N_v-(mu_1v+mu_2v*N_v)*N_v
        return np.array([f1,f2,f3,f4,f5,f6,f7])
        
    return y,f
    
    
    
    
    
if __name__ == "__main__":
    print " Ok"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    