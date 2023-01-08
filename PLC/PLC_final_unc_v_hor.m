clear
close all
clc

%% Parameter values

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.5; param.Kplc = 0.11;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3; param.k_i = 2; 
param.test = 1; param.tau_c = 2; param.in = 1;
param.tau_p = 0.0027; param.ps = 0;
param.tscale = 0.40;

%% ODE solver

rm_solver = @(x,t)PLC_final(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-30);
[T1,Y1] = ode15s(rm_solver,[0 200],[0.0805296 3.861011 0.70395 0 0],opts);
[T2,Y2] = ode15s(rm_solver,[200 450],[0.0805296 3.861011 0.70395 0 0.025],opts);
[T3,Y3] = ode15s(rm_solver,[450 950],[0.0805296 3.861011 0.70395 0 0.04],opts);
[T4,Y4] = ode15s(rm_solver,[950 1500],[0.0805296 3.861011 0.70395 0 0],opts);
% [T,Y] = ode113(rm_solver,[0 2000],Y(end,:),opts);

%% Results/Plots

T = [T1;T2;T3;T4];
Y = [Y1;Y2;Y3;Y4];

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p1 = Y(:,4);
p2 = Y(:,5);

p = p1 + p2;

fig = figure;
figure(1)
l_color = [0 0.447 0.741];
r_color = [0.466 0.674 0.188];
set(fig,'defaultAxesColorOrder',[l_color; r_color])
yyaxis left
plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
hold on
text(210,1.15,'UV','FontSize',60)
text(460,1.15,'UV','FontSize',60)
text(960,1.15,'VP','FontSize',60)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
% ylim([0 0.7])
hold on
yyaxis right
% xlim([0 4700])
% xticks([0 1000 2000 3000 4000])
% ylim([0.056153 0.056154])
% yticks([])
plot(T,p,'color',[0.466 0.674 0.188],'LineWidth',3)
xlabel('time (s)')
ylabel('[IP_3] \muM')
hold on
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=70;
box off
% axis([xmin, xmax, ymin, ymax])
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
saveas(gcf,'figure6','epsc')

% figure(1)
% plot(T,c,'b','LineWidth',4)
% % text(10,0.75,'PKC activated','FontSize',70)
% % text(260,0.75,'PKC inhibited','FontSize',70)
% % text(402,0.77,'PKC = 0.1','FontSize',80)
% % ylim([0 0.85])
% % xlim([0 500])
% % h = vline2(2500,'r--');
% % text(40,0.95,'V_{PLC} = 0.02 (\muM/s)','FontSize',60)
% % text(2040,0.95,'V_{PLC} = 0.16 (\muM/s)','FontSize',60)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% % axis([0 4000 0 0.8])
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% hold off
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'Figure9c','epsc')

%% FUNCTION

function out = PLC_final(t,in,param)

% Reduced PKC model, PLC, DAG and PKC assumed to be at quasi steady state.

c = in(1);
ct = in(2);
h = in(3);
p1 = in(4);
p2 = in(5);


% % Functions
% 
if t < 950
    param.test = 0;
else
    param.test = 1;
end
% 
if t < 950
    param.Vplc = 0;
else
    param.Vplc = 0.6;
end

% if t < 2000
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end

% if t < 2500
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end

% if t < 2500
%     param.in = 1;
% else
%     param.in = 0.1;
% end

% if t < 2000
%     param.out = 1;
% else
%     param.out = 0;
% end
    
ce = param.gamma.*(ct - c);
p = p1 + p2;

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
out(4) = param.tscale.*param.test.*(L - param.k_i.*p1);
out(5) = param.tscale.*param.tau_p.*(param.ps - p2);
out = out';

end