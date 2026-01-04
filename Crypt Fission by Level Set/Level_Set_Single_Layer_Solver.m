function [t,LAYER,COORDS] = Level_Set_Single_Layer_Solver(Phi_0,N_xi,N_zeta,t,radius,speed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   N_xi - this is the number of points in the xi-domain.
%   N_zeta - this is the number of points inthe zeta-domain. 
%   t - is the time vector for which the level set formulation will be
%   solved on. 
%   radius - the starting radius of the crypt level.
%   speed - the velocity of the propagating surface (this is the value c in 
%   the paper.  
% OUTPUT: 
%   t - this is the same time vector as the input. 
%   LAYER - a 1 x 1 structure with the following fields: 
%       Arclength - a vector of arclengths for at time t.
%       Phi - an Nxi x Nzeta x Nt array that stores the surface at time t.
%   COORDS - an Nt x 1 structure with the following fields:
%       XL - this gives the x-coordinates of the left half of the circles
%       at time t
%       XR - this gives the x-coordinates of the right half of the circles
%       at time t
%       Y  - this gives the y-coordinates of the circles at time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discretization of time variable:
N_t = length(t);
dt  = t(2) - t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radius of Crypt Layer with buffer:
eps     = 0.25;
rad     = (1+eps)*radius;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xi- and zeta-Space discretization:
xi     = linspace(-2*rad,2*rad,N_xi)';
dxi    = xi(2)-xi(1);
zeta   = linspace(-2*rad,2*rad,N_zeta)';
dzeta  = zeta(2)-zeta(1);

% Mesh Generation:
[XI,ZETA] = meshgrid(xi,zeta);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Matrices in xi- and zeta-space:
[Dxi_f,Dxi_b]     = diffper(N_xi);
[Dzeta_f,Dzeta_b] = diffper(N_zeta);
Dxi_f   = Dxi_f/dxi;        Dxi_b   = Dxi_b/dxi;
Dzeta_f = Dzeta_f/dzeta;    Dzeta_b = Dzeta_b/dzeta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation of Phi and Arc:
Phi = nan(N_xi,N_zeta,N_t);
Phi(:,:,1) = Phi_0;
Arc = nan(N_t,1);
Arc(1,:) = 2*pi*radius;

% Solve for Phi, ArcLength, and Contour for each time value:
for k = 1:N_t-1
    
        % Forcing Function:
        Fr      = @(XI,ZETA) max(0,(rad-sqrt((XI-rad).^2+ZETA.^2)));
        Fl      = @(XI,ZETA) max(0,(rad-sqrt((XI+rad).^2+ZETA.^2)));
        Phi_end = @(XI,ZETA) Fr(XI,ZETA) + Fl(XI,ZETA);
        F       = @(Phi,XI,ZETA) speed*(Phi-Phi_end(XI,ZETA));
            
        % Solving for Phi:
        Phi(:,:,k+1) = Phi(:,:,k) -...
        dt*(max(F(Phi(:,:,k),XI,ZETA),0).*sqrt((max(Dxi_b*Phi(:,:,k),0).^2+...
                                                 min(Dxi_f*Phi(:,:,k),0).^2)+...
                                                (max(Phi(:,:,k)*Dzeta_b',0).^2+...
                                                 min(Phi(:,:,k)*Dzeta_f',0).^2))+...
            min(F(Phi(:,:,k),XI,ZETA),0).*sqrt((max(Dxi_f*Phi(:,:,k),0).^2+...
                                                 min(Dxi_b*Phi(:,:,k),0).^2)+...
                                                (max(Phi(:,:,k)*Dzeta_f',0).^2+...
                                                 min(Phi(:,:,k)*Dzeta_b',0).^2)));

        % Finding the Arc-Length of zero level set: 
        C = contourc(xi,zeta,Phi(:,:,k),[eps*radius eps*radius]);

        x_mesh = [C(1,2:floor(end/2)) C(1,ceil(end/2+2):end)];
        y_mesh = [C(2,2:floor(end/2)) C(2,ceil(end/2+2):end)];

        xy_l_ind = find(x_mesh<=0);
        xy_r_ind = find(x_mesh>=0);

        x_l = x_mesh(xy_l_ind)';
        y_l = y_mesh(xy_l_ind)';

        x0_ind = find(x_l == 0);

        if isempty(x0_ind)==0
            endpoint = max(y_l(x0_ind));
        else
            endpoint = 0;
        end

        x_r = x_mesh(xy_r_ind)';
        y_r = y_mesh(xy_r_ind)';

        xy_ul_ind = find(y_l>=0); x_ul = x_l(xy_ul_ind); y_ul = y_l(xy_ul_ind);
        xy_ll_ind = find(y_l<=0); x_ll = x_l(xy_ll_ind); y_ll = y_l(xy_ll_ind);
        xy_ur_ind = find(y_r>=0); x_ur = x_r(xy_ur_ind); y_ur = y_r(xy_ur_ind);
        xy_lr_ind = find(y_r<=0); x_lr = x_r(xy_lr_ind); y_lr = y_r(xy_lr_ind);

        [x_ul,mix1] = sort(x_ul); y_ul = y_ul(mix1); 
        [x_ll,mix2] = sort(x_ll); y_ll = y_ll(mix2);
        [x_ur,mix3] = sort(x_ur); y_ur = y_ur(mix3);    
        [x_lr,mix4] = sort(x_lr); y_lr = y_lr(mix4); 

        Arc(k+1) = Crypt_Arclength(x_ul,y_ul) +...
                     Crypt_Arclength(x_ll,y_ll) +...
                     Crypt_Arclength(x_ur,y_ur) +...
                     Crypt_Arclength(x_lr,y_lr);
    
    % Store the coordinates of Phi for each time value in a structure:
    COORDS(k).XLCoords  = [x_ul; x_ll];       
    COORDS(k).XRCoords  = [x_ur; x_lr];       
    COORDS(k).XCoords  = [x_ul; x_ll; x_ur; x_lr];       
    COORDS(k).YCoords  = [y_ul; y_ll; y_ur; y_lr];
            
end

    % Store each curve and arclength in structure:
    LAYER.Arclength    = Arc(:);
    LAYER.Phi          = Phi(:,:,:);

end



