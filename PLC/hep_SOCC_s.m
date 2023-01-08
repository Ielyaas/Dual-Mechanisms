clear all
close all
clc

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       A K_PLC ...
       Vpm Kpm ...
       tau_i alpha0 alpha1 Kce...
       delta gamma pm_block in_block...
       kappa_p p_s step1 step2 step3...
       T_end;
   
% Parameter values

Kf = 4; Kc = 0.2; Kp = 0.2; Kb = 0.4;
tau_max = 80; K_tau = 0.09; Kh = 0.08;
Vs = 0.9; Kbar = 1.5e-5; Ks = 0.2;
tau_p = 1; A = 0.2; K_PLC = 0.1;
Vpm = 0.07; Kpm = 0.3;
tau_i = 1; alpha0 = 0.004; alpha1 = 0.01; Kce = 14;
delta = 2.5; gamma = 5.5; pm_block = 500; in_block = 500;
kappa_p = 0; p_s = 0;

% Kf = 4; Kc = 0.2; Kp = 0.2; Kb = 0.4;
% tau_max = 80; K_tau = 0.09; Kh = 0.08;
% Vs = 0.9; Kbar = 1.5e-5; Ks = 0.2;
% tau_p = 1; A = 0.55; K_PLC = 0.1;
% Vpm = 0.07; Kpm = 0.3;
% tau_i = 1; alpha0 = 0.004; alpha1 = 0.01; Kce = 14;
% delta = 2.5; gamma = 5.5; pm_block = 500; in_block = 500;
% kappa_p = 0; p_s = 1;

% ODE solver

T_end = 400;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode15s(@hep_SOCC,0:0.001:T_end,[0.0862898 4.137186 0.6433274 p_s],opts);
% [T,Y] = ode15s(@hep_SOCC,0:0.01:200,Y(end,:),opts);

% Plots

c = Y(:,1);
h = Y(:,2);
ce = Y(:,3);
p = Y(:,4);

ct = c + (1/gamma).*ce;

% PLC = A.*c.^2./(K_PLC^2 + c.^2);

phi_c = c.^4./(c.^4 + Kc^4);
phi_p = p.^2./(p.^2 + Kp^2);
phi_p_down = Kp^2./(p.^2 + Kp^2);

Jpm = Vpm.*c.^2./(Kpm^2 + c.^2);
Jin = tau_i.*(alpha0 + alpha1.*(Kce^4./(Kce^4+ce.^4)));

h_inf = Kh^4./(Kh^4 + c.^4);
tau = tau_max*K_tau^4./(K_tau^4 + c.^4);

Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + Kb.*(beta + alpha));

Jipr = Kf.*Po.*(ce-c);

load('slice51260318.mat')

fig = figure;
figure(1)
left_color = [0 0 1];
right_color = [0 0.8 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(T,Y(:,1),'b','LineWidth',2)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('[Ca^{2+}]_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% axis([0 250 0 0.8])
hold on
yyaxis right
plot(T,Y(:,4),'color',[0 0.8 0],'LineWidth',2)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('[IP_3] \muM','fontsize',20,'fontweight','b','fontname','arial')
% h=vline2(50,'r--');
% h=vline2(100,'r--');
% h=vline2(200,'g--');
set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% axis([0 250 0 1.5])
ax=gca;
% ax.YTick = 0.3;
set(ax,'Linewidth',3)
ax.FontSize=40;
box off

figure(2)
plot(T,Y(:,4),'color',[0 0.8 0],'LineWidth',2)
% figure(2)
% plot(Y(:,1),Y(:,2))
% 
% figure(3)
% plot(T,Y(:,2))
% hold on
% h=vline2(126,'g--');
% 
% figure(4)
% plot(T,Y(:,4))

% figure(5)
% plot(T,Y(:,3))
% hold on
% h=vline2(126,'g--');
% figure(2)
% plot(T,Y(:,2),'LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% 
% figure(3)
% plot(T,Y(:,3),'LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('Ca^{2+}_e \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% 
% figure(4)
% plot(T,Y(:,4),'r','LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',30,'fontweight','b','fontname','arial')
% 
% figure(5)
% plot3(Y(:,1),Y(:,2),Y(:,4),'LineWidth',1)
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% zlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% 
% figure(6)
% plot3(Y(:,1),Y(:,3),Y(:,4),'LineWidth',1)
% ylabel('Ca^{2+}_e \muM','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% zlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% 
% figure(7)
% plot(T,Jipr,'lineWidth',1)
% title('J_{IPR}')
% 
% figure(8)
% plot(T,Jserca,'lineWidth',1)
% title('J_{SERCA}')
% hold on
% % figure(9)
% plot(T,Jpm,'lineWidth',1)
% title('J_{PM}')
% legend('serca','pm')
% 
% figure(9)
% plot(T,ct)
% title('C_{tot}')

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
% 
% figure(16)
% plot(T,PLC,'LineWidth',1)
% hold on
% plot(T,Y(:,1),'r','LineWidth',1)
% title('PLC')

tilefigs