clear
close all
clc

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       V_PLC K_PLC ...
       Vpm Kpm ...
       alpha0 alpha1 ...
       delta gamma ct ...
       T_end;
   
% Parameter values

Kf = 10; Kc = 0.2; Kp = 0.2; Kb = 0.4;
tau_max = 100; K_tau = 0.1; Kh = 0.08;
Vs = 0.9; Kbar = 1.957e-5; Ks = 0.2;
tau_p = 1; V_PLC = 0.51; K_PLC = 0.1;
Vpm = 0.11; Kpm = 0.3;
alpha0 = 0.0027; alpha1 = 0.015;
delta = 2.5; gamma = 5.5; ct = 2;

% ODE solver

T_end = 500;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode45(@hep_ct,0:0.001:T_end,[0.270238,0,0.627277],opts);

% Plots

c = Y(:,1);
h = Y(:,2);
p = Y(:,3);

figure(1)
plot(T,Y(:,1),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
title('Simulation of hormone induced oscillations in a closed cell')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

% figure(2)
% plot(Y(:,3),Y(:,1),'r','LineWidth',1)
% xlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% title('Phase portrait of the closed system')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(3)
plot(Y(:,2),Y(:,1),'r','LineWidth',1)
xlabel('h','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
title('Phase portrait of the closed system')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

tilefigs