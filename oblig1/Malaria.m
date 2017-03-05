function [y0 ] = Malaria( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes her

% initial condition in order [Sh, Eh, Ih, Rh, Sv, Ev, Iv]
initial_cond_1 = [400,10,30,0,1000,100,50]; 
%initial_cond_2 = [700, 10, 30,0,1000,100,50];

%variables in order [lambda_h,phi_h,phi_v,beta_vh,beta_hv,beta_bar_vh,
%             sigma_h,sigma_v,nu_h, nu_v,gamma_h,delta_h,rho_h,
%             mu_1h,mu_1v,mu_2h,mu_2v]
lamda_h= 3.285*1e-2;phi_h= 7.666*1e-5;phi_v=0.4; 
beta_vh=0.8333;beta_hv=2.*1e-2;beta_bar_vh=8.333*1e-2;
sigma_h=18;sigma_v=0.6;nu_h=8.333*1e-2; nu_v=0.1;
gamma_h=3.704*1e-3; delta_h=3.454*1e-4;rho_h=1.460*1e-2;
mu_1h=4.212*1e-5;mu_1v=0.1429;mu_2h=1e-7;mu_2v=2.279*1e-4;

N_h = initial_cond_1(1) + initial_cond_1(2) +initial_cond_1(3) + initial_cond_1(4);
N_v = initial_cond_1(5) + initial_cond_1(6) + initial_cond_1(7);
 
e_h = initial_cond_1(2)/N_h(1);
i_h = initial_cond_1(3)/N_h(1);
r_h = initial_cond_1(4)/N_h(1);
e_v = initial_cond_1(6)/N_h(1);
i_v = initial_cond_1(7)/N_h(1);

y0 = [e_h;i_h;r_h;N_h;e_v;i_v;N_v];

        
end

