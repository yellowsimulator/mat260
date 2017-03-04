function [ x ] = RungeKutta( b, y, f, h, t, i )
% Embedded Runge Kutta

m=length(y);
x=zeros(m,1);

for j = 1:m
    k1 = f(t(i),y(j));
    k2 = f(t(i)+(1/3)*h, y(j)+h((1/3)*k1));
    k3 = f(t(i)+(1/3)*h, y(j)+h((1/6)*k1 + (1/6)*k2));
    k4 = f(t(i)+(1/2)*h, y(j)+h((1/8)*k1 + (3/8)*k3));
    k5 = f(t(i)+h, y(j)+h((1/2)*k1 + (-3/2)*k3 + 2*k4));

    x(j) = y(j) + h(b(1)*k1(j,i) + b(2)*k2(j,i) + ...
               b(3)*k3(j,i) + b(4)*k4(j,i) + b(5)*k5(j,i));
end


end

