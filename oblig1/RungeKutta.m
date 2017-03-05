function [ x ] = RungeKutta( b, y, f, h)
% Embedded Runge Kutta

m=length(y);
x=zeros(m,1);



k1 = f(y);
k2 = f(y+h*((1/3)*k1));
k3 = f(y+h*((1/6)*k1 + (1/6)*k2));
k4 = f(y+h*((1/8)*k1 + (3/8)*k3));
k5 = f(y+h*((1/2)*k1 + (-3/2)*k3 + 2*k4));

x = y + h*(b(1)*k1 + b(2)*k2 + ...
           b(3)*k3 + b(4)*k4 + b(5)*k5);



end

