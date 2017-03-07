function [tt, y ] = RungeKuttaEmbedded( h, tol, y0, T ,f) 
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


c=0.3;
p=4;
x_lower = y0;
x_higher = y0;
y=y0;

b_higher=[1/6,0,0,2/3,1/6];
b_lower=[1/2,0,-3/2,2,0];
t0 = 0;
tt = [t0];
while t0<T
    x_higher = RungeKutta( b_higher, x_higher, f, h);
    x_lower = RungeKutta( b_lower,x_lower, f, h);
    d = norm(x_higher-x_lower); %global error
    K = c*h^(p+1);
    error_bound = (d*h)/(T-t0);
    fprintf('global error = %d\n',d)
   
    if K > error_bound
    while K>(tol*h)/(T-t0)
        K_new = (d*h)/((T-t0)*K);
        h = c*(tol/((T-t0)*K_new))^(1/(p+1))*h;
        x_higher = RungeKutta( b_higher, x_higher, f, h);
        x_lower = RungeKutta( b_lower, x_lower, f, h);
        d = norm(x_higher-x_lower);
        t0 = t0 + h;
        tt = [tt,t0];
        y = [y , x_lower];
        K = K_new;
        %fprintf('K in side = %d\n',K);
    end %end while
    h = h*( tol/(T-t0)*K )^(1/(p+1));
    t0 = t0 + h;
    tt = [tt,t0];
    y = [y , x_lower];
    end %end if
   
    t0 = t0 + h;
    tt = [tt,t0];
    y = [y , x_lower];
    
  
end % en
end %while

