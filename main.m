clear 
close all
kb = 1; % For simplicity
p = [0 ; 0];
lambda = 2*pi/kb; %meter
ps = [lambda/2; 10*lambda];


x = 0:lambda/20:lambda;
y = 0:lambda/20:lambda;

figure(1)
hold on
xlim([-3*lambda 3*lambda])
grid on
set(gca, 'YDir','reverse')
scatter(ps(1),ps(2))
line([-lambda 2*lambda], [1.5*lambda 1.5*lambda])
annotation('textarrow',[0.4429,0.5179],[0.292857142857143,0.195238095238095],'String','\bf(\rho_{s})')
line([0 lambda],[lambda lambda])
yregion([0 lambda] , [0 lambda])
hold off

[X, Y] = meshgrid(x,y);
N = length(X)*length(X);

u_inc = (-1j/4)*besselh(0,2,kb*(sqrt((X - ps(1)).^2 + (Y - ps(2)).^2)));


figure(2)
hold on
title(['Real part of u_{inc} with k_b = ' num2str(kb), ' \rho_s = (' num2str(ps(1)) ',' num2str(ps(2)) ')'])
xlabel("x-axis")
imagesc(real(u_inc))
axis equal tight
set(gca, 'YDir','reverse')
clim([-0.1 0.1])
colorbar
hold off


figure(3)
hold on
title(['Imaginary part of u_{inc} with k_b = ' num2str(kb), ' \rho_s = (' num2str(ps(1)) ',' num2str(ps(2)) ')'])
imagesc(imag(u_inc))
axis equal tight
set(gca, 'YDir','reverse')
clim([-0.1 0.1])
colorbar
hold off

% k_rho = omega/c();
% chi = (k_rho/kb)^2 - 1; 