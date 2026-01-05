% A function to define the crypt domain using a parameterization
% in theta and z.  This file can also craete a partiall split or fully
% split crypt by varying the parameter c. 
clear; close all; clc;

% Mesh parameters:
N_theta = 100;
N_z     = 100;

% Splitting function witdh parametrs:
b = .25;
A = @(z,c) 1 - 1./(1+exp(-b*(z-c)));
g = @(theta,z,c) 1 - A(z,c) + 2*A(z,c).*sin(theta).^2;

% Radius function with parameters:
r_b = 41/2/pi;
r_t = 10/pi;
a   = .3;
L   = 78;
R = @(theta,z,c) g(theta,z,c).*(r_b*(1-exp(-a*z)) + r_t*exp(a*(z-L)));

% Mesh:
theta = linspace(eps,2*pi,N_theta);
z     = linspace(0,L,N_z);
[T,Z] = meshgrid(theta,z);

% Define splitting parameter, c: 
% c = 60;

for c = -2:.5:90

% Plot:
figure(1)
subplot(1,2,1)
surf(Z,T,R(T,Z,c))
set(gca,'fontsize',20)
set(gcf,'WindowState','maximized')
title('Plot of $r(\theta,z)$','fontsize',30,'interpreter','latex')
xlabel('$z$','fontsize',25,'interpreter','latex')
ylabel('$\theta$','fontsize',25,'interpreter','latex')
zlabel('$r(\theta,z)$','fontsize',25,'interpreter','latex')
view([1 -1 1])
grid on; grid minor
shading interp; camlight; lighting phong
subplot(1,2,2)
mesh(R(T,Z,c).*cos(T),R(T,Z,c).*sin(T),Z)
set(gca,'fontsize',20)
set(gcf,'WindowState','maximized')
title('Plot of Crypt Domain','fontsize',30,'interpreter','latex')
xlabel('$x$','fontsize',25,'interpreter','latex')
ylabel('$y$','fontsize',25,'interpreter','latex')
zlabel('$z$','fontsize',25,'interpreter','latex')
xlim([-30 30])
ylim([-30 30])
zlim([0 85])
view([-1 1 1])
grid on; grid minor
shading interp; camlight; lighting phong

if c == -2
    pause
else
    pause(.1)
end
end
