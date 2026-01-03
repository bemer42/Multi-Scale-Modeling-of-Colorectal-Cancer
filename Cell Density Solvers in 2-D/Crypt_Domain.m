function [x,y,z,radii] = Crypt_Domain(s,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
%   s - this is a vector of crypt level values.  The typical input would be
%   a vector from 1 to 82, spaced by 1.  
%   theta - this is a vector of values from 0 to 2pi with any number of
%   points.  
% OUTPUT: 
%   x - this is a matrix of size N_theta x N_s of all the x coordinates for
%   the entire crypt domain. 
%   y - this is a matrix of size N_theta x N_s of all the y coordinates for
%   the entire crypt domain. 
%   z - this is a matrix of size N_theta x N_s of all the z coordinates for
%   the entire crypt domain. z is the matrix form of s.
%   radii - this is a vector the same size as s which gives the radius of
%   the crypt at the crypt level s. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a meshgrid for the crypt domain: 
[S,Theta] = ndgrid(s,theta);

% Crypt Shape Parameters:
rb = 41/2/pi; 
rt = 10/pi;
L  = 82;

% Crypt Geometry:
R = @(s) (1/2/pi).*(s==0) + ...
         sqrt(rb^2-4/pi^2.*(s-pi/2*rb).^2).*(s>0).*(s<pi/2*rb) + ...
         rb.*(s>=pi/2*rb).*(s<L-pi/2*rt) + ...
         (rb + rt - sqrt(rt^2-4/pi^2*(s-(L-pi/2*rt)).^2)).*(s>=L-pi/2*rt).*(s<=L);

% Define outputs:
radii = R(S);
x = R(S).*cos(Theta);
y = R(S).*sin(Theta);
z = S;

end
     
     