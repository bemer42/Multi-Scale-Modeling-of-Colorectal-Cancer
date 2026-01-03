% Solving spring model for cell movement on 2-D surface domain:
close all;
clear all;
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the time variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_t  = 200;
t_0   = 0;
t_end = 2;
t     = linspace(t_0,t_end,N_t);
dt    = t(2) - t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the crypt level variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_s   = 82;
s_0   = 0;
s_end = 1;
s     = linspace(s_0,s_end,N_s)';
ds    = s(2) - s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the circumference variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_c   = 82;
c_0   = 0;
c_end = 1;
c     = linspace(c_0,c_end,N_c)';
dc    = c(2) - c(1);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L       = 82;
gamma_s = .1667;
sp      = 27/L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Differentiation Matrices:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s,D_s,D_ss] = q_diffmat2(N_s,s);
[c,D_c,D_cc] = q_diffmat2(N_c,c);
D_c(1,:)     = [0 1/2 zeros(1,N_c-3) -1/2]/dc;
D_c(end,:)   = [1/2 zeros(1,N_c-3) -1/2 0]/dc;
D_c          = sparse(D_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Heaviside:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = @(S) (S>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Boundary and Grid:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Boundary = q_bndy(N_s,N_c);
bndy = Boundary.istop + Boundary.isbottom;

[S,C] = ndgrid(s,c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crypt Geometry:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,z,radii] = Crypt_Domain(L*s',2*pi*c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the BC functions and the RHS function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vec   = @(Q) Q(~bndy);
unvec = @(q) q_unpack(q,bndy,N_s,N_c);

RHS    = @(t,q) vec((1/gamma_s)*D_s*((1./unvec(q).^2).*(D_s*unvec(q))) + ...
                    (1/gamma_s)*(1./radii.^2).*((1./unvec(q).^2).*(unvec(q)*D_c'))*D_c' + ...
                     H(sp-S).*unvec(q));

% RHS    = @(t,q) vec((1/gamma_s)*D_s*((1./unvec(q).^2).*(D_s*unvec(q))) + ...
%                     (1/gamma_s)*(1./radii.^2).*((1./unvec(q).^2).*(unvec(q)*D_c'))*D_c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supply an initial condition and solve the IBVP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_0 = 1+.5*H(sp-S).*H(.65-C).*H(C-.35);
 
tic
[t,q] = ode15s(RHS,t,vec(Q_0));
toc
q = q';

%% Movie Plot of Surface: 

for i = 1:N_t
    
    figure(1) 
    surf(S,C,unvec(q(:,i)))
    title('Cell Density','fontsize',30,'interpreter','latex')
    xlabel('$s$','fontsize',25,'interpreter','latex')
    ylabel('$\theta$','fontsize',25,'interpreter','latex')
    zlabel('$q$','fontsize',25,'interpreter','latex')
    set(gca,'fontsize',20) 
    light 
    lightangle(5,30)
    lighting gouraud
    shading interp
    material shiny
    colormap default
    box on
    grid on
    grid minor
    xlim([s_0 s_end])
    ylim([c_0 c_end])
    zlim([.9 1.1*max(max(Q_0))])
%     view([0 0 1])
    if i == 1 
    pause
    end

end

%% Plot of the cell density on the crypt domain: 
for k = 1:N_t
    
    Q = unvec(q(:,k));

    markerSize = 100;
    figure(1)
    scatter3(x(:),y(:),z(:),markerSize,Q(:),'filled')
    xlim([-L/1.5 L/1.5])
    ylim([-L/1.5 L/1.5])
    zlim([-1 L+1])
    grid on
    grid minor
    box on
    colorbar
%     caxis([1 1.5])
    if k == 1 
        pause
    end
    hold off

end



