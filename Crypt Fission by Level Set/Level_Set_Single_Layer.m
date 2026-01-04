%%  This program solves the level set equations in (27)-(29) of the paper.
% The lower sections will create Figure 4 of the paper.  
clear all;
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radius of Crypt Layer:
radius = 1;
eps    = 0.25;
rad    = (1+eps)*radius;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time discretization:
N_t   = 2000;
t_0   = 0;
t_end = 10;
t     = linspace(t_0,t_end,N_t);
dt    = t(2)-t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xi-Space discretization:
xi_min = -2*rad;
xi_max = 2*rad;
N_xi   = 250;
xi     = linspace(xi_min,xi_max,N_xi)';

% zeta-Space discretization:
zeta_min = -2*rad;
zeta_max = 2*rad;
N_zeta   = 250;
zeta     = linspace(xi_min,xi_max,N_xi)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid Generation:
[XI,ZETA] = meshgrid(xi,zeta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Condition:
Phi_0 = sparse(max(0,rad-sqrt(XI.^2+ZETA.^2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the PDE:
[t,LAYER,COORDS] = Level_Set_Single_Layer_Solver(Phi_0,N_xi,N_zeta,t,radius,2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie Plot of Phi-Function:
scrsz = get(0,'ScreenSize');
hf = figure('Position',scrsz);
for n = 1:1:N_t
    
    Phi = LAYER.Phi(:,:,n);
    
    figure(1)
    surf(XI,ZETA,Phi)
    title('Surface Propagation','fontsize',34)
    xlabel('\xi','fontsize',30)
    ylabel('\zeta','fontsize',30)
    zlabel('\phi(\xi,\zeta,t)','fontsize',30)
    set(gca,'fontsize',28)
    xlim([min(min(XI)) max(max(XI))])
    ylim([min(min(XI)) max(max(XI))])
    zlim([0 max(max(max(Phi)))])
    light 
    lightangle(5,30)
    lighting gouraud
    shading interp
    box on
    grid on
    grid minor
        
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie Plot of Single Ring:
scrsz = get(0,'ScreenSize');
hf = figure('Position',scrsz);
dn = 6;
for n = 1:dn:N_t
    
    Phi = LAYER.Phi(:,:,n);
    Arc = LAYER.Arclength(n);

    % Plot of Surface and Single layer with Arclength:
    figure(1)
    subplot(1,2,2)
    plot(COORDS(n).XCoords,COORDS(n).YCoords,'ko','linewidth',4)
    title('Level Set','fontsize',34)
    xlabel('\xi','fontsize',30)
    ylabel('\zeta','fontsize',30)
    set(gca,'fontsize',28)
    axis([min(xi) max(xi) min(zeta) max(zeta)])
    text('string',['Arc Length = ' num2str(Arc)],'fontsize',30,...
      'units','norm','pos',[.15 .23])
    axis square
    grid on
    grid minor
    subplot(1,2,1)
    surf(XI,ZETA,Phi)
    title('Level Set Function, \Phi','fontsize',34)
    xlabel('\xi','fontsize',30)
    ylabel('\zeta','fontsize',30)
    zlabel('\Phi','fontsize',30)
    set(gca,'fontsize',28)
    set(gca,'Box','on')
    xlim([min(xi) max(xi)])
    ylim([min(zeta) max(zeta)])
    zlim([min(min(min(Phi))) max(max(max(Phi)))])
    light 
    lightangle(5,30)
    lighting gouraud
    shading interp
    colormap summer
    box on
    grid on
    grid minor
    hold off
    if n==1
        pause
    end
end


