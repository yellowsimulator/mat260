"""
define parameters for a specific problem here
example: set vectors y, f for y' = f, with all parameters
"""
import numpy as np
def SIModel(h):
    """
    setting parameters for SI
    spreading of deseases model:
    (s,i)' = (f1,f2) where 
    f1 = -beta*s(t)*i(t)
    f2 = beta*s(t)*i(t)
    """
    #vector solution
    beta = 0.5
    M=2; T = 20
    N = int(T/h)
    t = np.linspace(0,T, N+1)
    s = np.zeros(N+1);i = np.zeros(N+1)
    s[0] = 0.9; i[0] = 0.1 # x1 = s, x2 = i
    y = np.array([s,i])
    #setting function f = (f1,f2)
    def f(t,y):
        s,i = y[0],y[1]
        beta = 0.5
        f1 = -beta*s*i # f for s
        f2 = beta*s*i # f for i
        return np.array([f1,f2])

    # setting analytical solutions
    ie = 1./(1+9*np.exp(-beta*t))
    se = 1-ie
    return y, f, se, ie, t, N, M, h,T
    
    
    
def MalariaModel():
    """
    set parameters for the malaria model
    """
    M = 7
    T = 300; h = 0.8
    N = int(T/h)
    t = np.linspace(0,T, N+1)
    e_h,i_h,r_h,N_h,e_v,i_v,N_v = np.zeros(N+1),np.zeros(N+1),np.zeros(N+1),\
    np.zeros(N+1),np.zeros(N+1),np.zeros(N+1),np.zeros(N+1)
    #initial conditions and model parameters
    IC1 = {'S_h':400,'E_h':10,'I_h':30,'R_h':0, \
                          'S_v':1000,'E_v':100,'I_v':50
                         }
    IC2 = {'S_h':700,'E_h':10,'I_h':30,'R_h':0, \
                          'S_v':1000,'E_v':100,'I_v':50
                         }
    p = {'lamda_h':3.285*1e-2,'phi_h':7.666*1e-5,'phi_v':0.4, \
             'beta_vh':0.8333,'beta_hv':2.*1e-2,'beta_bar_vh':8.333*1e-2,\
             'sigma_h':18.,'sigma_v':0.6,'nu_h':8.333*1e-2, 'nu_v':0.1, \
             'gamma_h':3.704*1e-3,'delta_h':3.454*1e-4,'rho_h':1.460*1e-2,\
             'mu_1h':4.212*1e-5,'mu_1v':0.1429,'mu_2h':1e-7,'mu_2v':2.279*1e-4
             }
    
    N_h[0] = IC1['S_h'] + IC1['E_h'] + IC1['I_h'] + IC1['R_h']
    N_v[0] = IC1['S_v'] + IC1['E_v'] + IC1['I_v']
    
    e_h[0] = (IC1['E_h']/N_h[0])
    i_h[0] = (IC1['I_h']/N_h[0])
    r_h[0] = (IC1['R_h']/N_h[0])
    e_v[0] = (IC1['E_v']/N_v[0])
    i_v[0] = (IC1['I_v']/N_v[0])
    y = np.array([e_h,i_h,r_h,N_h,e_v,i_v,N_v])

    def f(t,y):
        " function f"
        e_h = y[0]
        i_h = y[1]
        r_h = y[2]
        N_h = y[3]
        e_v = y[4]
        i_v = y[5]
        N_v = y[6]
        
        f1 = ( (p['sigma_v']*p['sigma_h']*N_v*p['beta_hv']*i_v)/(p['sigma_v']*N_v+p['sigma_h']*N_h) )* \
        (1-e_h-i_h-r_h)-(p['nu_h']+p['phi_h']+(p['lamda_h']/N_h )*e_h )+p['delta_h']*i_h*e_h
        
        f2 = p['nu_h']*e_h-(p['gamma_h']+p['sigma_h']+p['phi_h']+(p['lamda_h']/N_h))*\
           i_h+p['sigma_h']*i_h*i_h
           
        f3 = p['gamma_h']*i_h-(p['rho_h']+p['phi_h']+(p['lamda_h']/N_h))*r_h+\
             p['sigma_h']*i_h*r_h
             
        f4 = p['lamda_h']+p['phi_h']*N_v-(p['mu_1h']+p['mu_2h']*N_h)-p['sigma_h']*i_h*N_h
        
        f5 = ((p['sigma_v']*p['sigma_h']*N_h)/(p['sigma_v']*N_v+p['sigma_h']*N_h))*\
             (p['beta_vh']*i_h+p['beta_bar_vh']*r_h)*(1-e_v-i_v)-(p['nu_v']+p['phi_v'])*e_v
             
        f6 = p['nu_v']*e_v-p['phi_v']*i_v
        f7 = p['phi_v']*N_v-(p['mu_1v']+p['mu_2v']*N_v)*N_v
        f = np.array([f1,f2,f3,f4,f5,f6,f7])
        print f
        return f
    
    #f(0,y)
    return y,f,t,N,M,h,T
    
    
    
    
    
if __name__ == "__main__":
    print " Ok"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    