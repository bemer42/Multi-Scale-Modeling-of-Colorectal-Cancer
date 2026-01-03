% Solving density model with no production:
close all;
clear all;
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the time variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_t  = 1e3;
t_0   = 0;
t_end = 2;
t     = linspace(t_0,t_end,N_t);
dt    = t(2) - t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize the space variable: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_s   = 100;
s_0   = 0;
s_end = 1;
s     = linspace(s_0,s_end,N_s)';
ds    = s(2) - s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = .1667;
sp    = 27/82;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Differentiation Matrix for the Second Derivative:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s,D_s,D_ss] = diffcheb(N_s-1,[s_0 s_end]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Heaviside:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = @(s) (s>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the BC functions and the RHS function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_R = @(v) 1;
q_L = @(v) -D_s(1,2:end)*[v; q_R(v)]./D_s(1,1);

Chop   = @(q) q(2:end-1);
Extend = @(v) [q_L(v); v; q_R(v)];


RHS    = @(t,v) Chop((1/gamma)*D_s*((1./Extend(v).^2).*(D_s*Extend(v))) + ...
                H(sp-s).*Extend(v));
% RHS    = @(t,v) Chop((1/gamma)*D_s*((1./Extend(v).^2).*(D_s*Extend(v))));          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supply an initial condition and solve the IBVP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_0 = 1 + H(sp-s);

[t,V] = ode15s(RHS,t,Chop(q_0));

V = V';
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
zlim([0 max(max(Q))])
light 
lightangle(5,30)
lighting gouraud
shading interp
box on
grid on
grid minor

for i = 1:25:N_t
    
    figure(1)
    plot3(s,t(i)*ones(size(s)),Q(:,i),'k-','linewidth',2)
    
end


%% Movie Plot:
for k = 1:N_t
    
   figure(2)
   plot(s,Q(:,k),'k-','linewidth',5)
   title('Solution to Cell Density Equation','fontsize',30,'interpreter','latex')
   xlabel('crypt level ($s$)','fontsize',25,'interpreter','latex')
   ylabel('cell density ($q$)','fontsize',25,'interpreter','latex')
   set(gca,'fontsize',20)
   xlim([s_0 s_end])
   ylim([.9*min(min(Q)) 1.1*max(max(Q))])
   grid on
   grid minor
   if k == 1
       pause
   end
   
end


