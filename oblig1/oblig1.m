y0 = Malaria();

h = 0.1;
T = 300;
c = 0.1;
tol = 0.1;
[tt,y] = RungeKuttaEmbedded(h, tol, y0, T ,@f);

%length(tt)
%length(y(7,:))
plot(tt,y(7,:))