clear
close all
clc

%%   
% Parameter values

param_ct.Kf = 4; param_ct.Kc = 0.2; param_ct.Kp = 0.2; param_ct.Kb = 0.4;
param_ct.tau_max = 80; param_ct.K_tau = 0.09; param_ct.Kh = 0.08;
param_ct.Vs = 0.9; param_ct.Kbar = 1.5e-5; param_ct.Ks = 0.2;
param_ct.tau_p = 1; param_ct.R_act = 0.3; param_ct.K_PLC = 0.2;
param_ct.Vpm = 0.07; param_ct.Kpm = 0.3;
param_ct.alpha0 = 0.004; param_ct.alpha1 = 0.01; param_ct.Kce = 14;
param_ct.delta = 2.5; param_ct.gamma = 5.5; param_ct.tau_cdum = 2;

%%
% ODE solver

hep_1 = @(x,t)hep_SOCC_ct_tr(x,t,param_ct);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T,Y] = ode15s(hep_1,[0 300],[0.0862898,0.424888,4.137186,0],opts);
% [T,Y] = ode15s(hep_1,0:0.001:200,Y(end,:),opts);

%%
% Plots

c = Y(:,1);
h = Y(:,2);
ct = Y(:,3);
p = Y(:,4);

% Bifurcation diagram
% load('sliceA06.mat')

figure(1)
plot(T,c,'Color',[0.29 0.53 0.87],'LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
hold on
% plot(T,ct,'LineWidth',1)
% legend('Ca^{2+}_i','Ca^{2+}_t','Location','Best')
% title('R_{act} = 0.45')
% h=vline2(100,'r');
set(gca,'FontSize',20,'fontweight','b','fontname','arial')
% axis([0 500 0 1])
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=40;
box off
% 

%figure(17)
% plot(ct,c,'k','LineWidth',2)
% xlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% title('V_{PLC}=0.2231373')
% set(gca,'FontSize',20,'fontweight','b','fontname','arial')
% axis([0 4 0 1])
% ax=gca;
% set(ax,'Linewidth',3)
% ax.FontSize=30;
% box off

% figure(20)
% plot(T,c,'b','LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% yyaxis left
% ylabel('[Ca^{2+}]_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% % axis([0 250 0 0.8])
% hold on
% plot(T,h,'r','LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% yyaxis right
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% h=vline2(100,'r--');
% % h=vline2(100,'r--');
% % h=vline2(200,'g--');
% set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% axis([0 250 0 1.5])
% ax=gca;
% set(ax,'Linewidth',3)
% ax.FontSize=30;
% box off
% 
% figure(2)
% plot(T,h,'LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% % 
% figure(3)
% plot(T,ct,'LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% % 
% figure(4)
% plot(T,p,'r','LineWidth',1)
% xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',30,'fontweight','b','fontname','arial')
% % 
% figure(5)
% plot3(c,h,p,'LineWidth',1)
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% zlabel('IP_3 \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')
% % 
% figure(6)
% plot3(c,h,ct,'LineWidth',1)
% zlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')

tilefigs