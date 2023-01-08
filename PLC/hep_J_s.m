clear all
close all
clc

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       A K_PLC ...
       Vpm Kpm ...
       alpha0 alpha1 ...
       delta gamma pm_block in_block...
       T_end;
   
% Parameter values

Kf = 10; Kc = 0.2; Kp = 0.2; Kb = 0.4;
tau_max = 20; K_tau = 0.1; Kh = 0.08;
Vs = 0.9; Kbar = 1.957e-5; Ks = 0.2;
tau_p = 1; A = 0.5; K_PLC = 0.1;
Vpm = 0.11; Kpm = 0.3;
alpha0 = 0.0027; alpha1 = 0.02;
gamma = 5.5; pm_block = 200; in_block = 200;

% ODE solver

T_end = 200;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode15s(@hep_J,0:0.001:T_end,[0.0119,0.9965,7.12,0.007],opts);

% Plots

c = Y(:,1);
h = Y(:,2);
ce = Y(:,3);
p = Y(:,4);

ct = c + (1/gamma).*ce;

PLC = A.*c.^2./(K_PLC^2 + c.^2);

phi_c = c.^4./(c.^4 + Kc^4);
phi_p = p.^2./(p.^2 + Kp^2);
phi_p_down = Kp^2./(p.^2 + Kp^2);

Jpm = Vpm.*c.^2./(Kpm^2 + c.^2);
Jin = (alpha0 + alpha1*A);

h_inf = Kh^4./(Kh^4 + c.^4);
tau = tau_max*K_tau^4./(K_tau^4 + c.^4);

Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + Kb.*(beta + alpha));

Jipr = Kf.*Po.*(ce-c);

figure(1)
plot(T,Y(:,1),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
h=vline2(in_block,'r')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(2)
plot(T,Y(:,2),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(3)
plot(T,Y(:,3),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_e \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(4)
plot(T,Y(:,4),'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(5)
plot3(Y(:,1),Y(:,2),Y(:,4),'LineWidth',1)
ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
zlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(6)
plot3(Y(:,1),Y(:,3),Y(:,4),'LineWidth',1)
ylabel('Ca^{2+}_e \muM','fontsize',20,'fontweight','b','fontname','arial')
xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
zlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')

figure(7)
plot(T,Jipr,'lineWidth',1)
title('J_{IPR}')

figure(8)
plot(T,Jserca,'lineWidth',1)
title('J_{SERCA}')
hold on
% figure(9)
plot(T,Jpm,'lineWidth',1)
title('J_{PM}')
legend('serca','pm')

figure(7)
plot(T,ct)
title('C_{tot}')

% figure(10)
% plot(T,phi_c,'LineWidth',1)
% title('\phi_c')
% 
% figure(11)
% plot(T,phi_p,'LineWidth',1)
% title('\phi_p')
% 
% figure(12)
% plot(T,phi_p_down,'LineWidth',1)
% title('\phi_p_{down}')
% 
% figure(13)
% plot(T,h_inf,'LineWidth',1)
% title('h_{\infty}')
% 
% figure(14)
% plot(T,beta,'LineWidth',1)
% title('\beta')
% 
% figure(15)
% plot(T,alpha,'LineWidth',1)
% title('\alpha')

tilefigs