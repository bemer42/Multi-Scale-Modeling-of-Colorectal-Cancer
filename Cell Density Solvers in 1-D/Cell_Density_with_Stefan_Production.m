%% This program will run the Cell Density PDE solver with the production 
% and shedding terms.  The very last section will create Figure 3 of the 
% paper. 
close all;
clear all;
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the time variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_t  = 1000;
t_0   = 0;
t_end = 6.4;
t     = linspace(t_0,t_end,N_t);
dt    = t(2) - t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the space variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_s   = 41+1;
s_0   = 0;
s_end = 1;
s     = linspace(s_0,s_end,N_s)';
ds    = s(2) - s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lstar = 82;
alpha = 1840.3;
Tc    = 15.1954;
ell   = 1;
mu    = 5.428;
sp    = 27;
ss    = 77;

gamma = @(t,L) L.^2.*log(2)/ell^2/alpha/Tc/Lstar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Differentiation Matrix for the Second Derivative:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s,D_s,D_ss] = diffcheb(N_s-1,[s_0 s_end]);
% [s,D_s,D_ss] = diffmat2(N_s,[s_0 s_end]); s = s(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Heaviside:
%%%%%%%%%%%%%%%%%7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = @(s) (s>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the BC functions and the RHS function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_R = @(v) 1;
q_L = @(v) -D_s(1,2:end)*[v; 1]./D_s(1,1);

Chop   = @(q) q(2:end-1);
Extend = @(v) [q_L(v); v; 1];

RHS_Q  = @(t,v,L) Chop(-s.*(1./gamma(t,L)).*(D_s(end,:)*Extend(v)).*(D_s*Extend(v)) + ...
                       (1./gamma(t,L)).*(D_s*((1./Extend(v).^2).*(D_s*Extend(v)))) + ...
                       (H(sp-s.*L) - mu*H(s.*L-ss)).*Extend(v));
                       

RHS_L  = @(t,v,L) -1./gamma(t,L).*(D_s(end,:)*Extend(v)).*L;

RHS    = @(t,u) [RHS_Q(t,u(1:end-1),u(end)); 
                 RHS_L(t,u(1:end-1),u(end))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supply an initial condition and solve the IBVP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_0 = 1;
q_0 = ones(size(s));

[t,U] = ode15s(RHS,t,[Chop(q_0); L_0]);

V = U(:,1:end-1)';
L = U(:,end)';
Q = zeros(N_s,N_t);

for j = 1:N_t   
    Q(:,j) = Extend(V(:,j));  
end

% Surface Plot of solution:
[S,T] = meshgrid(s,t);

figure(1)
surf(S,T,Q')
hold on
plot3(s,zeros(size(s)),q_0,'k-','linewidth',4)
title('Solution to Cell Density Equation','fontsize',30,'interpreter','latex')
xlabel('crypt level ($s$)','fontsize',25,'interpreter','latex')
ylabel('time ($t$)','fontsize',25,'interpreter','latex')
zlabel('cell density ($q$)','fontsize',25,'interpreter','latex')
set(gca,'fontsize',20)
xlim([s_0 s_end])
ylim([t_0 t_end])
zlim([min(min(Q)) max(max(Q))])
light 
lightangle(5,30)
lighting gouraud
shading interp
box on
grid on
grid minor

for i = 1:50:N_t
    
    figure(1)
    plot3(s,t(i)*ones(size(s)),Q(:,i),'k-','linewidth',2)
    
end

% Plot of the length of the crypt:
figure(2)
plot(t,L,'k-','linewidth',5)
set(gca,'fontsize',20)
title('Length of the Crypt','fontsize',30,'interpreter','latex')
xlabel('time ($t$)','fontsize',25,'interpreter','latex')
ylabel('crypt length ($L$)','fontsize',25,'interpreter','latex')
xlim([t_0 t_end])
ylim([min(L) max(L)])
box on
grid on
grid minor

%% Movie Plot:
for k = 1:N_t
    
   v_s = -alpha.*Lstar./(Q(:,k).^3).*(D_s*Q(:,k));
    
   figure(3)
   plot(linspace(s_0,L(k),N_s),Q(:,k),'k','linewidth',5)
   hold on
   plot(linspace(s_0,L(k),N_s),v_s,'k-.','linewidth',5)
   plot([L(k) L(k)],[0 100],'k:','linewidth',5)
   set(gca,'fontsize',20)
   title('Solution to Cell Density Equation','fontsize',30,'interpreter','latex')
   xlabel('crypt level ($s$)','fontsize',25,'interpreter','latex')
   ylabel('cell density ($q$)','fontsize',25,'interpreter','latex')
   xlim([s_0 1.1*max(L)])
   ylim([0 1.5*max(max(Q))])
   grid on
   grid minor
   hold off
   if k == 1
       pause
   end
   
   k
   
end

%% Creating Figure 3 of Crypt Fission Paper:
% Four k values:
k1 = 1;
k2 = 499;
k3 = 790;
k4 = 1000; 

k = [k1 k2 k3 k4];

for i = 1:4
      
   v_s = -alpha.*Lstar./(Q(:,k(i)).^3).*(D_s*Q(:,k(i)));
   time   = round(100*Tc/log(2)*t(k(i)))/100;
   length = round(100*L(k(i)))/100;
    
   figure(4)
   subplot(4,1,i)
   plot(linspace(s_0,L(k(i)),N_s),Q(:,k(i)),'k','linewidth',5)
   hold on
   plot(linspace(s_0,L(k(i)),N_s),v_s,'k-.','linewidth',5)
   plot([length length],[0 100],'k:','linewidth',5)
   set(gca,'fontsize',20)
   if i == 1
   title('Cell Density and Velocity with Cell Production and Shedding','fontsize',35,'interpreter','latex')
   legend('$\hat{q}(\hat{s},\hat{\tau})$','$\hat{v_s}(\hat{s},\hat{\tau})$','$\hat{L}(\hat{\tau})$','interpreter','latex')
   ylabel('cells and cell lengths/hour ($\hat{q}, \hat{v_s}$)','fontsize',25,'interpreter','latex')
   end
   if i == 4
   xlabel('crypt level ($\hat{s}$)','fontsize',25,'interpreter','latex')
   end
%    if i ==1 
%    ylabel('cells and cell lengths/hour ($\hat{q}, \hat{v_s}$)','fontsize',25,'interpreter','latex')
%    end
   text(.6,.8,['$\hat{\tau} = $ ' num2str(time) ' hrs'],'units','normalized','interpreter','latex','fontsize',20)
   text(.6,.55,['$\hat{L} = $ ' num2str(length)],'units','normalized','interpreter','latex','fontsize',20)
%    if i ==1
%    legend('$\hat{q}(\hat{s},\hat{\tau})$','$\hat{v_s}(\hat{s},\hat{\tau})$','$\hat{L}(\hat{\tau})$','interpreter','latex')
%    end
   xlim([s_0 1.1*max(L)])
   ylim([0 1.5*max(max(Q))])
   grid on
   grid minor
   hold off

end



