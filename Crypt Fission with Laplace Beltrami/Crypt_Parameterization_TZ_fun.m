function R = Crypt_Parameterization_TZ_fun(T,Z,c)
% This function is used to update the crypt domain as the cell density
% increases.  

% Splitting function with parametrs:
b = .25;
A = @(z,c) .95*(1 - 1./(1+exp(-b*(z-c))));
g = @(theta,z,c) 1 - A(z,c) + 2*A(z,c).*sin(theta).^2;

% Radius function with parameters:
r_b = 41/2/pi;
r_t = 10/pi;
a   = .3;
L   = 78;
R = @(theta,z,c) g(theta,z,c).*(r_b*(1-exp(-a*z)) + r_t*exp(a*(z-L)));

% Output the Radius: 
R = R(T,Z,c); 

