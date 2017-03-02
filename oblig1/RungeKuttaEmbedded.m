function [ y ] = RungeKuttaEmbedded( h, tol, y0, T ,f) 
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


t = 0:h:T;


c=0.1;
p=4;
y = y0;

b_higher=[1/6,0,0,2/3,1/6];
b_lower=[1/2,0,-3/2,2,0];

n=1;
while t(n)<T
    if n~=1
        t = t(n):h:T;
        
    end
    
    x_higher = RungeKutta( b_higher, y(:,n), f, h, t, n);
    x_lower = RungeKutta( b_lower, y(:,n), f, h, t, n);
    
    K = norm(x_higher-x_lower);
    
    while K>(tol*h)/(T-t(0))
        h = c(tol/(T-t(0)*K))^(1/(p+1))*h;
        x_higher = RungeKutta( b_higher, y, f, h, t, n);
        x_lower = RungeKutta( b_lower, y, f, h, t, n);
    
        K = norm(x_higher-x_lower);
    end 
    h = h*(tol/(T-t(0)*K))^(1/(p+1));
    n = n+1;
    y = [y,x_lower];
  
end
  
end

