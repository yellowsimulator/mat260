function [ y ] = RungeKuttaEmbedded( h, tol, y0, T ,f) 
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

t1=0;
c=0.1;
p=4;
x_lower = y0;
x_higher = y0;
y=y0;



b_higher=[1/6,0,0,2/3,1/6];
b_lower=[1/2,0,-3/2,2,0];
t0 = h;
n=1;
while t0<T
    x_higher = RungeKutta( b_higher, x_higher, f, h);
    x_lower = RungeKutta( b_lower,x_lower, f, h);
    K = norm(x_higher-x_lower);
%     fprintf('K = %d\n',K)
%     
%     while K>(tol*h)/(T-t1)
%         h = c*(tol/((T-t1)*K))^(1/(p+1))*h;
%         x_higher = RungeKutta( b_higher, x_higher, f, h);
%         x_lower = RungeKutta( b_lower, x_lower, f, h);
%         K = norm(x_higher-x_lower);
%         fprintf('K in side = %d\n',K);
%     end 
%     h = h*( tol/(T-t1)*K )^(1/(p+1));
    
    t0 = t0 + h;
    %disp(t0);
    %fprintf('h1 = %d\n',h);
    
    y = [y , x_lower];
    
end
  
end

