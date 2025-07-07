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
annotation('textarrow',[0.4429,0.5179],[0.292857142857143,0.195238095238095],'String','\bf(\rho_{s})')
yregion( 0 , lambda)
line([0 lambda],[lambda lambda])
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


%% Define the object (Diamond)
omega = 3;
c = 3;
radius = 3;
n = 21;
center = ceil(n/2);
[X_obj, Y_obj] = meshgrid(1:n, 1:n);
obj = (abs(X_obj - center) + abs(Y_obj - center)) <= radius;

% Convert to double
obj = double(obj);

k_rho = omega./(c*obj);
k_rho(isinf(k_rho)) = 0;
chi = (k_rho/kb).^2 - 1; 

figure(4)
hold on
title(['Contrast'])
imagesc(chi)
axis equal tight
set(gca, 'YDir','reverse')
colorbar
hold off

%% Define the Receiver domain
L = 3*lambda;
Rec_x = [-lambda;   2*lambda];
Rec_y = [1.5*lambda;  1.5*lambda];

figure(1)
hold on 
line(Rec_x, Rec_y)
hold off
%% Intgral formulation
upperBound = [2; 2];
lowerBound = [0; 0];
n_int = 100;
h = (upperBound(1) - lowerBound(1))/n;
x = 0;
for i = 1:n-1
    pos_int_x = lowerBound(1) + i*h;
    pos_int_y = lowerBound(2) + i*h;
    temp = besselh(0,2,kb*sqrt((X - pos_int_x).^2 + (Y - pos_int_y).^2))*chi;
    x = x + temp; 
end

int = ((-1*(kb*sqrt((X - ps(1)).^2+ (Y - ps(2)).^2)*kb^2))/(16))*(h*x);

figure(5)
imagesc(real(int))