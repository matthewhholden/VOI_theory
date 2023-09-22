n=50;
b1 = linspace(0,1,n);
b2 = linspace(0,1,n);

[B1 B2] = meshgrid(b1, b2);
Vstar = B1.*B2./(B1+B2);
Vstar(1)=0;

figure(1)
surf(B1, B2, Vstar)
zlabel('EVPI', FontSize=16)
xlabel('\beta_1', FontSize=16)
ylabel('\beta_2', FontSize=16)
title('p = p^*', FontSize=16)


%%%%%%%%%%%%%%%%%
%%%% Fig 3
%%%%%%%%%%%%%%%%%

h = figure(3)
b2 = .1
p = linspace(0,1,n);
b1 = linspace(0,1,n);
[P B1] = meshgrid(p,b1);
pstar = b2./(B1+b2)
V = zeros(n,n)
V(P<pstar) = P(P<pstar).*B1(P<pstar)
V(P>=pstar) = (1-P(P>=pstar))*b2

subplot(1,3,1)
surf(P, B1, V)
zlim([0 0.5])
zlabel('EVPI', FontSize=16)
xlabel('p', FontSize=16)
ylabel('\beta_1', FontSize=16)
title('\beta_2 = 0.1', FontSize=16)

b2 = .5
p = linspace(0,1,n);
b1 = linspace(0,1,n);
[P B1] = meshgrid(p,b1);
pstar = b2./(B1+b2)
V = zeros(n,n)
V(P<pstar) = P(P<pstar).*B1(P<pstar)
V(P>=pstar) = (1-P(P>=pstar))*b2

subplot(1,3,2)
surf(P, B1, V)
zlim([0 0.5])
zlabel('EVPI', FontSize=16)
xlabel('p', FontSize=16)
ylabel('\beta_1', FontSize=16)
title('\beta_2 = 0.5', FontSize=16)


b2 = 1
p = linspace(0,1,n);
b1 = linspace(0,1,n);
[P B1] = meshgrid(p,b1);
pstar = b2./(B1+b2)
V = zeros(n,n)
V(P<pstar) = P(P<pstar).*B1(P<pstar)
V(P>=pstar) = (1-P(P>=pstar))*b2

subplot(1,3,3)
surf(P, B1, V)
zlim([0 0.5])
zlabel('EVPI', FontSize=16)
xlabel('p', FontSize=16)
ylabel('\beta_1', FontSize=16)
title('\beta_2 = 1', FontSize=16)


% export_fig(h)
% set(h,'PaperSize',[20 10])