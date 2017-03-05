function [ x ] = f(y)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
lamda_h= 3.285*1e-2;phi_h= 7.666*1e-5;phi_v=0.4; 
beta_vh=0.8333;beta_hv=2.*1e-2;beta_bar_vh=8.333*1e-2;
sigma_h=18;sigma_v=0.6;nu_h=8.333*1e-2; nu_v=0.1;
gamma_h=3.704*1e-3; delta_h=3.454*1e-4;rho_h=1.460*1e-2;
mu_1h=4.212*1e-5;mu_1v=0.1429;mu_2h=1e-7;mu_2v=2.279*1e-4;




f1 = (sigma_v*sigma_h*y(7)*beta_hv*y(6))/(sigma_v*y(7)+sigma_h*y(4))...
    * (1-y(1)-y(2)-y(3))-(nu_h+phi_h+lamda_h/y(4))*y(1) +delta_h*y(2)*y(1);
        
f2 = nu_h*y(1)-(gamma_h+sigma_h+phi_h+(lamda_h/y(4)))*y(2)+sigma_h*(y(2))^2;
           
f3 = gamma_h*y(2)-(rho_h+phi_h+(lamda_h/y(4)))*y(3)+sigma_h*y(2)*y(3);
             
f4 = lamda_h+phi_h*y(7)-(mu_1h+mu_2h*y(4))*y(4)-sigma_h*y(2)*y(4);
        
f5 = ((sigma_v*sigma_h*y(7))/(sigma_v*y(7)+sigma_h*y(4)))*(beta_vh*y(2)+beta_bar_vh*y(3))*(1-y(5)-y(6))-(nu_v+phi_v)*y(5);
             
f6 = nu_v*y(5)-phi_v*y(6);
f7 = phi_v*y(7)-(mu_1v+mu_2v*y(7))*y(7);


x = [f1;f2;f3;f4;f5;f6;f7];


end
