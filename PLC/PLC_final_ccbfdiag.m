clear
close all
clc

%% Parameter values

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.1; param.Kplc = 0.11;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3; param.k_i = 2; 
param.test = 1; param.tau_c = 2; param.in = 1;
param.tscale = 0.10;

%% ODE solver

rm_solver = @(x,t)PLC_final(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-30);
[T,Y] = ode15s(rm_solver,[0 4700],[0.0805296 4.861011 0.70395 0],opts);
[T1,Y1] = ode113(rm_solver,[0 2000],Y(end,:),opts);

%% LOAD FILES

load('PLC1.dat')
load('PLC6.dat')

%% Results/Plots

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

c1 = Y1(:,1);
ct1 = Y1(:,2);
h1 = Y1(:,3);
p1 = Y1(:,4);

figure(1)
% plot(ct,c,':','Color',[0.635 0.078 0.184],'LineWidth',3)
% hold on
% plot(ct1,c1,'Color',[0.635 0.078 0.184],'LineWidth',4)
% hold on
plot(PLC1(1:243,1),PLC1(1:243,2),'k','LineWidth',4)
hold on
plot(PLC1(244:606,1),PLC1(244:606,2),'k--','LineWidth',4)
hold on
plot(PLC1(607:1244,1),PLC1(607:1244,2),'k','LineWidth',4)
hold on
plot(PLC1(1245:1349,1),PLC1(1245:1349,2),'--','Color',[0 0.447 0.741],'LineWidth',4)
hold on
plot(PLC1(1245:1349,1),PLC1(1245:1349,3),'--','Color',[0 0.447 0.741],'LineWidth',4)
hold on
plot(PLC1(1350:end,1),PLC1(1350:end,2),'Color',[0.466 0.674 0.188],'LineWidth',4)
hold on
plot(PLC1(1350:end,1),PLC1(1350:end,3),'Color',[0.466 0.674 0.188],'LineWidth',4)
hold on
% plot(PLC6(1:52,1),PLC6(1:52,2),'k','LineWidth',4)
% hold on
% plot(PLC6(53:138,1),PLC6(53:138,2),'k--','LineWidth',4)
% hold on
% plot(PLC6(139:449,1),PLC6(139:449,2),'k','LineWidth',4)
% hold on
% plot(PLC6(450:end,1),PLC6(450:end,2),'--','Color',[0 0.447 0.741],'LineWidth',4)
% hold on
% plot(PLC6(450:end,1),PLC6(450:end,3),'--','Color',[0 0.447 0.741],'LineWidth',4)
% hold on
% plot(0.1768,0.1165,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(0.5082,0.04501,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(0.1672,0.09684,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% text(0.52,0.065,'SN','FontSize',60)
% hold on
% text(0.06,0.09684,'SN','FontSize',60)
% hold on
% text(0.08,0.145,'HB','FontSize',60)
% hold on
% plot(3.623,0.2604,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(2.405,0.06924,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(1.341,0.1363,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(4.635,0.4668,'o','MarkerSize',20,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
hold on
text(2.45,0.08,'SN','FontSize',60)
hold on
text(0.97,0.1363,'SN','FontSize',60)
hold on
text(3.75,0.24,'HB','FontSize',60)
hold on
text(4.7,0.4668,'SNPO','FontSize',60)
% hold on
% text(0.03,0.755,'V_{PLC} = 0.6 (\muM/s)','FontSize',80)
xlabel('c_t \muM')
ylabel('c \muM')
% set(gca,'FontSize',20)
axis([0 6 0 0.8])
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=70;
box off
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
saveas(gcf,'1Dman_base','epsc')

%% FUNCTION

function out = PLC_final(~,in,param)

% Reduced PKC model, PLC, DAG and PKC assumed to be at quasi steady state.

c = in(1);
ct = in(2);
h = in(3);
p = in(4);


% Functions

    
ce = param.gamma.*(ct - c);

L = param.Vplc.*c.^2./(c.^2 + param.Kplc^2);
 
phi_c = c.^4./(c.^4 + param.Kc^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = param.tau_max*param.Ktau^4./(c.^4 + param.Ktau^4);

beta = phi_c.*phi_p.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jin = param.in.*(param.alpha0 + param.alpha1*param.Kce^4./(ce.^4 + param.Kce^4));
Jpm = param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = param.tscale.*(Jipr - Jserca + param.delta.*(Jin - Jpm))./param.tau_c;
out(2) = param.tscale.*param.delta.*(Jin - Jpm);
out(3) = param.tscale.*(h_inf - h)./tau;
out(4) = param.tscale.*param.test.*(L - param.k_i.*p);
out = out';

end