clear
close all
clc

%% Parameter set

param.Kf = 40; param.Kc = 0.2; param.Kp = 0.3; param.Kb = 0.4;
param.tau_max = 200; param.K_tau = 0.09; param.Kh = 0.1;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.V_plc = 0.43; param.Kplc = 0.11;
param.Vpm = 0.07; param.Kpm = 0.3;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.delta = 2.5; param.gamma = 5.5; param.tau_c = 2;
param.k_i = 2; param.tau_p = 0;

%% ODE solver

bm_solve = @(x,t)hep_base_model(x,t,param);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T1,Y1] = ode15s(bm_solve,[0 20],[0.0805296 3.861011 0.70395 0],opts);
[T2,Y2] = ode15s(bm_solve,[20 150],[0.0805296 3.861011 0.70395 0.05],opts);
[T3,Y3] = ode15s(bm_solve,[150 250],[0.0805296 3.861011 0.70395 0.1],opts);
% [T4,Y4] = ode15s(bm_solve,[250 500],[0.0805296 3.861011 0.70395 0],opts);
% [T,Y] = ode113(bm_solve,[0 250],Y(end,:),opts);

%% Results/Plots
T = [T1;T2;T3];
Y = [Y1;Y2;Y3];

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

L = param.V_plc*c.^2./(c.^2 + param.Kplc^2);
tau = param.tau_max.*param.K_tau^4./(param.K_tau^4 + c.^4);

fig = figure;
figure(1)
left_color = [0 0 1];
right_color = [0 0.8 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(T,c,'LineWidth',2)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
hold on
hold on
yyaxis right
% figure(3)
plot(T,p,'color',[0 0.8 0],'LineWidth',2)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('[IP_3] \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',20,'fontweight','b','fontname','arial')
% axis([0 400 0 0.7])
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=40;
box off
% % hold on
% % h=vline2(50,'r--');
% % hold on
% % h=vline2(150,'r--');
% % hold on
% % h=vline2(250,'g--');
% % figure(2)
% % plot(T,L,'LineWidth',2)
% hold on
% hold on
% yyaxis right
% % figure(3)
% plot(T,p,'color',[0 0.8 0],'LineWidth',2)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('[IP_3] \muM','fontsize',20,'fontweight','b','fontname','arial')
% figure(4)
% plot(T,h,'LineWidth',2)

% figure(4)
% plot(T,c,'LineWidth',2)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% set(gca,'FontSize',20,'fontweight','b','fontname','arial')
% ax=gca;
% set(ax,'Linewidth',3)
% ax.FontSize=45;
% box off
% tilefigs