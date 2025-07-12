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
% omega = 3;
% c = 3;
% radius = 5;
% n = 21;
% center = ceil(n/2);
% [X_obj, Y_obj] = meshgrid(1:n, 1:n);
% obj = (abs(X_obj - center) + abs(Y_obj - center)) <= radius;

n = 21;
line_length = 7;
obj = zeros(n, n);
center = ceil(n/2);
half_length = floor(line_length/2);
%obj(center, center-half_length:center+half_length) = 1;

obj = zeros(n, n);
obj(center-half_length:center+half_length, center) = 1;


% Convert to double
chi = double(obj);

% k_rho = omega./(c*obj);
% k_rho(isinf(k_rho)) = 0;
% chi = (k_rho/kb).^2 - 1; 

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
axis equal tight
hold off
%% Intgral formulation
% upperBound = [2; 2];
% lowerBound = [0; 0];
% n_int = 100;
% h = (upperBound(1) - lowerBound(1))/n;
% x = 0;
% % for i = 1:n-1
% %     pos_int_x = lowerBound(1) + i*h;
% %     pos_int_y = lowerBound(2) + i*h;
% %     temp = besselh(0,2,kb*sqrt((X - pos_int_x).^2 + (Y - pos_int_y).^2))*...
% %            besselh(0,2,(kb*sqrt((X - ps(1)).^2+ (Y - ps(2)).^2)*kb^2));
% %     x = x + temp; 
% % end


m = 40;
antenna_cordintes = [];
%u_sc = zeros([m 1]);
antenna_x = linspace(Rec_x(1) -2 , Rec_x(2)+2, m/2);
antenna_y = ones([1 m/2])*Rec_y(1);
%antenna_y2 = -1*ones([1 m/2])*Rec_y(1);
%antenna_y2 = linspace(-Rec_y(1) - 5, -Rec_y(2), m/2);
antenna_cordintes = [antenna_x ; antenna_y];
%antenna_cordintes = [antenna_cordintes, [antenna_x ; antenna_y2]];
% Calculate distances between x and y
dist_source2image = zeros(n);
dist_image2antenne = cell(m, 1);
for m = 1:length(antenna_cordintes)
    temp = zeros(n);
    for i = 1:length(X)
        for j = 1:length(Y)
            dist_source2image(i,j)  = sqrt((X(1,i) - ps(1))^2 + (Y(j,1) - ps(2))^2) ; 
            temp(i,j) = sqrt((X(1,i) - antenna_cordintes(1,m))^2 + (Y(j,1) - antenna_cordintes(2,m))^2);
        end
    end
    dist_image2antenne{m} = temp;
end
A = zeros([m N]);
for j = 1:m
    sum = 0;
    %for i = 1:n
        sum = sum + besselh(0,2,kb*dist_image2antenne{j})*besselh(0,2,kb*dist_source2image);
    %end    
%u_sc(j) = (-kb^2)/(16).*sum;
A(j,:) = vec((-kb^2)/(16).*sum)';

end    
%A = ((-1)/(16))*(h.*u_sc);

figure(5)
imagesc(real(A)) % Check documentation of imagesc for X and Y axis


figure(6)
clf(6)
hold on
ylabel('Y-axis')
xlabel('X-axis')
%axis equal tight
set(gca, 'YDir','reverse')
scatter(antenna_cordintes(1, 1:m/2),antenna_cordintes(2,1:m/2))
scatter(antenna_cordintes(1, m/2:end), antenna_cordintes(2, m/2:end))
hold off
%% Reconstruction
chi_vec = vec(chi);
%u_sc = vec(int);

u_sc = A*chi_vec;

image = pinv(A)*u_sc;
image = reshape(image, [n n]);
figure(7)
hold on
title(['Reconstruction using ', num2str(m), ' antennas'])
colorbar;
%clim([0 0.25])
imagesc(real(image))
axis equal tight
set(gca, 'YDir','reverse')
hold off

[~, S, ~] = svd(A);

figure(8)
clf(8)
hold on
title("Singular values of the A matrix")
plot(diag(S))
set(gca, 'YScale','log')
hold off