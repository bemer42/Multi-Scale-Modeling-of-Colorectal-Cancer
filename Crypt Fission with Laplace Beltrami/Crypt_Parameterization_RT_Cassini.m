% A function to define the crypt domain using a parameterization
% in r and theta with the Cassini splitting maechaniz:
close all, clear all, clc;

% Mesh parameters:
N_r     = 200;
N_theta = 200;

% Cassini with parametrs:
X = @(r,theta) r.*cos(theta); 
Y = @(r,theta) r.*sin(theta); 
rho_1 = @(r,theta,d) (X(r,theta)-d).^2 + Y(r,theta).^2; 
rho_2 = @(r,theta,d) (X(r,theta)+d).^2 + Y(r,theta).^2; 
Z = @(r,theta,d) rho_1(r,theta,d).*rho_2(r,theta,d)./(rho_1(r,theta,d) + rho_2(r,theta,d)); 
Z = @(r,theta,d) Z(r,theta,d).*(Z(r,theta,d)<=82) + 82.*(Z(r,theta,d)>82);

% Mesh:
r     = linspace(0,200,N_r);
theta = linspace(eps,2*pi,N_theta);
[R,T] = meshgrid(r,theta);

% Define splitting parameter, d: 
% d = 10;

for d = 0:.1:16

% Plot:
figure(1)
mesh(X(R,T),Y(R,T),Z(R,T,d))
set(gca,'fontsize',20)
set(gcf,'WindowState','maximized')
title('Plot of Crypt Domain','fontsize',30,'interpreter','latex')
xlabel('$x$','fontsize',25,'interpreter','latex')
ylabel('$y$','fontsize',25,'interpreter','latex')
zlabel('$z$','fontsize',25,'interpreter','latex')
xlim([-30 30])
ylim([-30 30])
zlim([0 100])
view([-1 1 -1])
grid on; grid minor
shading interp; camlight; lighting phong

if d == 0
    pause
else
    pause(.1)
end
end
