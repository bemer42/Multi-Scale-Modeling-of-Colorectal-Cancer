%% This file essentially calls the level set cascade solver.
% The lower section creates Figure 5 of the paper.
clear all;
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split Lag and plotting parameters:
arc_lag = .92;
speed   = 1;
Lag     = [arc_lag speed];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Crypt domain and radii:
N_theta = 100;
N_s     = 82;
L       = 82;
s       = linspace(0,L,N_s);
theta   = linspace(0,2*pi,N_theta);
[x,y,z,radii] = Crypt_Domain(s,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial Points in xi-zeta space at each level:
N_xi   = 100;
N_zeta = 100;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time discretization:
N_t   = 1e3;
t_0   = 0;
t_end = 17;
t     = linspace(t_0,t_end,N_t);
dt    = t(2)-t(1);
 
tic
[t,CRYPT_ARC,CRYPT] = Level_Set_Cascade_Solver(N_xi,N_zeta,t,radii,Lag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the layers into a single plot:
ind_max = 0;
for k = 1:N_t-1
    Crypt_X = []; Crypt_Y = []; Crypt_Z = [];
  
    for j = 1:N_s
        if j <= ind_max
            CRYPT(k,j).XL = CRYPT(k,j).XL-rad_plot+radii(j); 
            CRYPT(k,j).XR = CRYPT(k,j).XR+rad_plot-radii(j);
        end
        Crypt_X = [Crypt_X; [CRYPT(k,j).XL; CRYPT(k,j).XR]];
        Crypt_Y = [Crypt_Y; CRYPT(k,j).Y];
        Crypt_Z = [Crypt_Z; s(j).*ones(size(CRYPT(k,j ).Y))];
    end
    
    Crypt_Coords(k).X = Crypt_X;
    Crypt_Coords(k).Y = Crypt_Y;
    Crypt_Coords(k).Z = Crypt_Z;
    
    ind_plot = find(abs(CRYPT_ARC(k,:)./(2*pi*radii))>1.01);
    ind_max  = max(ind_plot);
    
    if isempty(ind_plot) == 1
        rad_plot = 0;
        ind_max  = 0;
    else 
        rad_plot = 1.075*radii(ind_max);
    end
    
end
toc
   
%%
% Movie Plot of Stack of Rings:
scrsz = get(0,'ScreenSize');
hf = figure('Position',scrsz);
dk = 2;
Tc      = 15.1954;
for k = 1:dk:N_t-1
    
    % Plot of Surface and Single layer with Arclength:
    figure(1)
    plot3(Crypt_Coords(k).X,Crypt_Coords(k).Y,Crypt_Coords(k).Z,'k.','linewidth',3)
    hold on
    title(['Level Set Stack, $t = $' num2str(t(k)*Tc/log(2))],'fontsize',34,'interpreter','latex') 
    xlabel('$\xi$','fontsize',30,'interpreter','latex')
    ylabel('$\zeta$','fontsize',30,'interpreter','latex')
    zlabel('crypt level ($\hat{s}$)','fontsize',30,'interpreter','latex')
    set(gca,'fontsize',28)
    axis([-6*max(radii) 6*max(radii) -6*max(radii) 6*max(radii) -1 L+1])
    view([-.4 1 .2])
    box on
    grid on
    grid minor
    hold off
    if k==1
        pause
    end
%     elseif k == 300
%         pause
%     elseif k == 1130
%         pause
%     end
end



