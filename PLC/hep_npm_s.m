clear all
close all
clc

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       A K_PLC ...
       gamma ct...
       T_end;
   
% Parameter values

Kf = 10; Kc = 0.2; Kp = 0.2; Kb = 0.4;
tau_max = 100; K_tau = 0.09; Kh = 0.08;
Vs = 0.9; Kbar = 1.957e-5; Ks = 0.2;
tau_p = 1; A = 0.26; K_PLC = 0.4;
gamma = 5.5; ct = 3.7;

% ODE solver

T_end = 500;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode15s(@hep_npm,0:0.001:T_end,[0.0713,0.614,0],opts);

% Plots

c = Y(:,1);
h = Y(:,2);
p = Y(:,3);


% PLC = A.*c.^2./(K_PLC^2 + c.^2);
% 
% phi_c = c.^4./(c.^4 + Kc^4);
% phi_p = p.^2./(p.^2 + Kp^2);
% phi_p_down = Kp^2./(p.^2 + Kp^2);
% 
% h_inf = Kh^4./(Kh^4 + c.^4);
% tau = tau_max*K_tau^4./(K_tau^4 + c.^4);
% 
% Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);
% 
% beta = phi_p.*phi_c.*h;
% alpha = phi_p_down.*(1-phi_c.*h_inf);
% 
% Po = beta./(beta + Kb.*(beta + alpha));
% 
% Jipr = Kf.*Po.*(ce-c);
% 
% h_inf = Kh^4./(Kh^4 + c.^4);
% tau = tau_max*K_tau^4./(K_tau^4 + c.^4);
% 
% Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);
% 
% beta = phi_p.*phi_c.*h;
% alpha = phi_p_down.*(1-phi_c.*h_inf);
% 
% Po = beta./(beta + Kb.*(beta + alpha));
% 
% Jipr = Kf.*Po.*(ce-c);

figure(1)
plot(T,Y(:,1),'k','LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
title('C_t = 3\muM, A = 0.4\muMs^{-1}')
h=vline2(in_block,'r');
set(gca,'FontSize',20,'fontweight','b','fontname','arial')
box off