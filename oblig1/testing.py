from ode_solvers import *
import matplotlib.pylab as plt

def test():
    T = 300.
    h = 0.01
    N=int(T/h)
    t = np.linspace(0, T, N+1)
    N_v = np.zeros(N+1)
    N_v[0] = 1150;
   
    y = np.array([N_v])
    M = len(y)
    
    def f(y):
        phi_v = 0.4
        mu_1v = 0.1429
        mu_2v = 2.279*1e-4
        f1 = phi_v*y[0]-(mu_1v+mu_2v*y[0])*y[0]
        print " f ",f1
        return np.array([f1])
       
    return y, f, t, h ,M,N
    


def model():
    T = 300.
    h = 0.01
    N=int(T/h)
    t = np.linspace(0, T, N+1)
    e_h = np.zeros(N+1); i_h = np.zeros(N+1); r_h = np.zeros(N+1)
    N_h = np.zeros(N+1); e_v = np.zeros(N+1);N_v = np.zeros(N+1)
    i_v = np.zeros(N+1);
    
    
    S_h0 = 400.; E_h0 = 10.; I_h0 = 30.; R_h0 = 0.; S_v0 = 1000.; E_v0 = 100.; I_v0 = 50.
    N_v0 = S_v0 + E_v0 + I_v0 
    N_h0 = S_h0 + E_h0 + I_h0 + R_h0
    N_h[0] = S_h0 + E_h0 + I_h0 + R_h0; N_v[0] = S_v0 + E_v0 + I_v0 
    e_h[0] = E_h0/N_h0; i_h[0] = I_h0/N_h0; r_h[0] = R_h0/N_h0
    e_v[0] = E_v0/N_v0; i_v[0] = I_v0/N_v0
    
    print " "
    
    N_v[0] = 1150;
    
    y = np.array([e_h,i_h,r_h,N_h,e_v,i_v,N_v])
    M = len(y)
    
    def f(t,y):
        sigma_v = 0.6; sigma_h = 18.; beta_hv = 2.*1e-2; nu_h = 8.333*1e-2; phi_h = 7.666*1e-5
        l_h = 3.285*1e-2; delta_h = 3.454*1e-4 ; gamma_h = 3.704*1e-3; rho_h = 1.460*1e-2
        mu_1h = 4.212*1e-5;  mu_2h = 1e-7; beta_bar_vh =  8.333*1e-2 ; beta_vh = 0.8333
        nu_v = 0.1; phi_v = 0.4; mu_1v =0.1429 ;mu_2v = 2.279*1e-4
        with np.errstate(divide='ignore', invalid='ignore'):
            a = sigma_v*sigma_h*y[3]*beta_hv*y[5]
            b = sigma_v*y[6]+sigma_h*y[3]
            c = np.true_divide(a,b)
            c[c == np.inf] = 0
            c = np.nan_to_num(c)
            
            d = np.true_divide(l_h,y[3])
            c[c == np.inf] = 0
            d = np.nan_to_num(d)
            f1 = c*(1-e_h[0]-y[1]-y[2])-( nu_h+phi_h+(d) )*y[0] + delta_h*y[1]*y[0]
             
            f2 = nu_h*y[0]-(gamma_h+delta_h+phi_h+(d))*y[1]+delta_h*y[1]*y[1]
            f3 = gamma_h*y[1]-(rho_h+phi_h+(d))*y[2]+delta_h*y[1]*y[2]
        
            f4 = l_h+phi_h*y[3]-(mu_1h+mu_2h*y[3])*y[3]-delta_h*y[1]*y[3]
            
            e = (sigma_v*sigma_h*y[6])
            f = (sigma_v*y[6]+sigma_h*y[3])
            g = np.true_divide(e,f)
            g[g == np.inf] = 0
            g = np.nan_to_num(g)
            f5 =(g)*(beta_vh*y[1]+beta_bar_vh*y[2])*(1-y[4]-y[5])-\
               (nu_v+phi_v)*y[4]
            f6 = nu_v*y[4]-phi_v*y[5]
            f7 = phi_v*y[0]-(mu_1v+mu_2v*y[0])*y[0]
        return np.array([f1,f2,f3,f4,f5,f6,f7])
    return y, f, t, h ,M,N
        

def plot(t,u):
    plt.plot(t,u[0])
    plt.show()
    

y, f, t,h,M,N = test()
u = YRKMethod(y,f,M,N,h,t)

plot(t,u)