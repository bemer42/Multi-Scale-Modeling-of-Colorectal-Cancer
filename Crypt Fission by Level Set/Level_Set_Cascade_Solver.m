function [t,CRYPT_ARC,CRYPT] = Level_Set_Cascade_Solver(N_xi,N_zeta,t,radii,Lag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   N_xi - this is the number of points in the xi-domain.
%   N_zeta - this is the number of points inthe zeta-domain. 
%   t - is the time vector for which the level set formulation will be
%   solved on. 
%   radii - this is the vector of crypt radius values at each crypt level.
%   This should be the R(s) output of Crypt_Domain.m
%   Lag - this is a 2 x 1 vector of parameter values the first input of
%   which is the arc_lag and the second is the speed.  Both parameters
%   together determine the speed at which the crypt splits.  arc_lag is a
%   number that tells the current level to start splitting based on the
%   arclength of the crypt level below it.  I've been using .92 for this
%   parameter.  speed is the velocity of the propagating surface.  I've
%   been using 1. 
% OUTPUT: 
%   t - the same time vector that is inputted.
%   CRYPT_ARC - this is a matrix of size N_t by N_s (here N_s is the number
%   of crypt levels or the length of radii) that gives the arclength of the
%   current level at time t.
%   CRYPT - this is a N_t x N_s structure with the following fields: 
%       XL - this gives the x-coordinates of the left half of the current
%       level at time t
%       XR - this gives the x-coordinates of the right half of the current
%       level at time t
%       Y  - this gives the y-coordinates of the current level at time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Number of Layers:
N_s = length(radii);

% Discretization of time variable:
N_t = length(t);
dt  = t(2) - t(1);

% Preallocation:
CRYPT_ARC = nan(N_t,N_s);
CRYPT_ARC(1,:) = 2*pi*radii;
PHI = nan(N_xi,N_zeta,N_s);

% Solve for Phi, ArcLength, and Contour of each layer at each time:
for k = 1:N_t-1
    for j = 1:N_s

    % Initial Condition:
    if k <= 1
        xi     = linspace(-2*radii(j),2*radii(j),N_xi)';
        zeta   = linspace(-2*radii(j),2*radii(j),N_zeta)';
        [XI,ZETA] = meshgrid(xi,zeta);
        Phi_0 = sparse(max(0,radii(j)-sqrt(XI.^2+ZETA.^2)));    
    else
        Phi_0 = PHI(:,:,j);
    end

    if j == 1
         vel = Lag(2);
    else
         vel = Lag(2).*(CRYPT_ARC(k,j-1)>Lag(1)*CRYPT_ARC(1,j-1));
    end    
    
    [t_layer,LAYER,COORDS] = Level_Set_Single_Layer_Solver(Phi_0,N_xi,N_zeta,[t(k) t(k)+dt],radii(j),vel);

    CRYPT(k,j).XL = COORDS(1).XLCoords;    
    CRYPT(k,j).XR = COORDS(1).XRCoords;   
    CRYPT(k,j).Y = COORDS(1).YCoords;
    CRYPT_ARC(k+1,j) = LAYER.Arclength(end);
    PHI(:,:,j) = LAYER.Phi(:,:,end);
     
    end 
k

end

end



