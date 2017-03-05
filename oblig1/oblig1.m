y0 = Malaria();

h = 0.1;
T = 100;
c = 0.1;
tol = 0.1;
y = RungeKuttaEmbedded(h, tol, y0, T ,@f);

t = 0:h:T;
plot(t,y(7,:))