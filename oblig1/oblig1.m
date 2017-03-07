%plot by choosing a number 
% between 1 and 7 corresponding 
% to the solution you want to plot
i = 6;

y0 = Malaria();
h = 0.1;  % initialize time step
T = 30;
c = 0.1; % safety factor
tol = 0.1;
[tt,y] = RungeKuttaEmbedded(h, tol, y0, T ,@f);
plot(tt,y(i,:))
if i ==2
axis([0,1,0,0.09])
end