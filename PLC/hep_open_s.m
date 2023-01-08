clear all
close all
clc

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       V_PLC K_PLC ...
       Vpm Kpm ...
       alpha0 alpha1 ...
       delta gamma Vp_start...
       T_end;
   
% Parameter values

Kf = 10; Kc = 0.2; Kp = 0.2; Kb = 0.4;
tau_max = 100; K_tau = 0.1; Kh = 0.08;
Vs = 0.9; Kbar = 1.957e-5; Ks = 0.2;
tau_p = 1; V_PLC = 0.1; K_PLC = 0.1;
Vpm = 0.11; Kpm = 0.3;
alpha0 = 0.0027; alpha1 = 0.015;
delta = 2.5; gamma = 5.5; Vp_start = 300;

% ODE solver

T_end = 300;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode15s(@hep_open,0:0.001:T_end,[0.073622,16.642,0.58233,0],opts);

% Plots

c = Y(:,1);
ce = Y(:,2);
h = Y(:,3);
p = Y(:,4);

figure(1)
plot(T,Y(:,1),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i conc. \muM','fontsize',20,'fontweight','b','fontname','arial')
title('Simulation of hormone-induced Ca^{2+} oscillations')
axis([0 300 0 1])
set(gca,'FontSize',30,'fontweight','b','fontname','arial')

% figure(2)
% plot(T,Y(:,2),'LineWidth',1)

% 
% figure(3)
% plot(T,Y(:,3),'LineWidth',1)

figure(4)
plot(T,Y(:,4),'r','LineWidth',1)
ylabel('IP_3 conc. \muM','fontsize',20,'fontweight','b','fontname','arial')
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
title('Simulation of hormone-induced IP_3 oscillations')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')

% figure(5)
% plot(Y(11001:70001,1),Y(11001:70001,4),'k','LineWidth',2)
% ylabel('IP_3 conc. \muM','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i conc. \muM','fontsize',20,'fontweight','b','fontname','arial')
% % ylabel('Ca^{2+}_e \muM','fontsize',20,'fontweight','b','fontname','arial')
% title('Phase space of Ca^{2+}_i & IP_3 trajectory')
% set(gca,'FontSize',30,'fontweight','b','fontname','arial')
% 
% plot3(Y(:,1),Y(:,3),Y(:,4))
% xlabel('c')
% ylabel('h')
% zlabel('p')

m_a = c.^4./(Kc^4 + c.^4);
m_b = m_a;
h_a = Kh^4./(Kh^4 + c.^4);
h_inf = h_a;
b = p.^2./(Kp^2 + p.^2);
a = 1- b;
Beta = b.*m_b.*h;
Alpha = a.*(1-(m_a.*h_a));

Po = Beta./(Beta + Kb.*(Beta + Alpha));

J_ipr = Kf.*Po.*(ce - c);

Jserca = Vs.*((c.^2 - Kbar.*ce.^2)./(Ks^2 + c.^2));

nu_PLC = V_PLC.*c.^2./(K_PLC^2 + c.^2);

Jpm = Vpm.*c.^2./(Kpm^2 + c.^2);
Jin = alpha0 + alpha1*V_PLC;

figure(6)
plot(T,J_ipr,'LineWidth',1)
title('J_{IPR}')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')

figure(7)
plot(T,Jserca,'LineWidth',1)
title('J_{SERCA}')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')

figure(8)
plot(T,Jpm,'LineWidth',1)
title('J_{PM}')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')

tilefigs