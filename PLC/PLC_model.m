clear
close all
clc

%% Parameter values

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.6; param.Kplc = 0.11;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3; param.k_i = 2; 
param.test = 1; param.tau_c = 2; param.in = 1;
param.tscale = 1; param.out = 1;

%% ODE solver

rm_solver = @(x,t)PLC_final(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-30);
[~,Y] = ode15s(rm_solver,[0 5050],[0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(rm_solver,[0 500],Y(end,:),opts);

%% Results/Plots

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

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
Jpm = param.out.*(param.Vpm.*c.^2./(c.^2 + param.Kpm^2));

% fig = figure;
% figure(1)
% l_color = [0 0.447 0.741];
% r_color = [0.466 0.674 0.188];
% set(fig,'defaultAxesColorOrder',[l_color; r_color])
% yyaxis left
% plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% % text(30,0.75,'V_{PLC} = 0.07 (\muM/s)','FontSize',60)
% % text(2050,0.75,'V_{PLC} = 0.14 (\muM/s)','FontSize',60)
% % ylim([0 0.8])
% hold on
% yyaxis right
% % xlim([0 4700])
% % xticks([0 1000 2000 3000 4000])
% % ylim([0.04 0.06])
% % yticks(0.05)
% plot(T,p,'color',[0.466 0.674 0.188],'LineWidth',3)
% xlabel('time (s)')
% ylabel('[IP_3] \muM')
% % text(30,0.33,'V_{PLC} = 0.07 (\muM/s)','FontSize',60)
% % text(2050,0.33,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% hold on
% % xline(1000,'r--','LineWidth',3)
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=70;
% box off
% % axis([xmin, xmax, ymin, ymax])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'figure15d','epsc')

figure(1)
plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
% hold on
% plot(T,Jserca,'r','LineWidth',2)
% text(10,0.75,'PKC activated','FontSize',70)
% text(260,0.75,'PKC inhibited','FontSize',70)
% text(402,0.77,'PKC = 0.1','FontSize',80)
ylim([0 0.9])
% xlim([0 500])
% xline(750,'r--','LineWidth',3)
text(10,0.85,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% % text(2050,1.05,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
% axis([0 2000 0 0.9])
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=70;
box off
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Figure12b','epsc')

% figure(3)
% subplot(3,1,1);
% plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
% ylabel('c \muM')
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=50;
% box off
% hold on
% subplot(3,1,2);
% plot(T,ct,'Color',[0.635 0.078 0.184],'LineWidth',3)
% ylabel('c_t \muM')
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=50;
% box off
% hold on
% subplot(3,1,3);
% plot(T,p,'Color',[0.466 0.674 0.188],'LineWidth',3)
% ylabel('p \muM')
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=50;
% box off
% % text(10,0.75,'PKC activated','FontSize',70)
% % text(260,0.75,'PKC inhibited','FontSize',70)
% % text(402,0.77,'PKC = 0.1','FontSize',80)
% % ylim([0 0.9])
% % xlim([0 500])
% % xline(750,'r--','LineWidth',3)
% % text(10,0.85,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% % % text(2050,1.05,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% xlabel('time (s)','FontSize',70)
% % axis([0 2000 0 0.9])
% ax=gca;
% set(ax,'Linewidth',5)
% % set(gca,'XScale','log','YScale', 'log')
% % ylim([0 200])
% % ax.FontSize=50;
% box off
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'QT2','epsc')

% figure(4)
% plot(h,c,'Color',[0 0.447 0.741],'LineWidth',3)
% % hold on
% % plot(T,Jserca,'r','LineWidth',2)
% % text(10,0.75,'PKC activated','FontSize',70)
% % text(260,0.75,'PKC inhibited','FontSize',70)
% % text(402,0.77,'PKC = 0.1','FontSize',80)
% % ylim([0 0.9])
% % xlim([0 500])
% % xline(750,'r--','LineWidth',3)
% % text(10,0.85,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% % % text(2050,1.05,'V_{PLC} = 0.6 (\muM/s)','FontSize',60)
% xlabel('h')
% ylabel('[Ca^{2+}_i] \muM')
% % axis([0 2000 0 0.9])
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=70;
% box off
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'Figure12b','epsc')

%% FUNCTION

function out = PLC_final(~,in,param)

% Reduced PKC model, PLC, DAG and PKC assumed to be at quasi steady state.

c = in(1);
ct = in(2);
h = in(3);
p = in(4);


% % Functions
% 
% if t < 2500
%     param.test = 0.0027;
% else
%     param.test = 0.25;
% end
% % 
% if t < 2000
%     param.Vplc = 0.07;
% else
%     param.Vplc = 0.6;
% end

% if t < 750
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end
% 
% if t < 2500
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end

% if t < 800
%     param.in = 1;
% else
%     param.in = 0;
% end
% 
% if t < 1000
%     param.out = 1;
% else
%     param.out = 0.3;
% end
    
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
Jpm = param.out.*(param.Vpm.*c.^2./(c.^2 + param.Kpm^2));

% DE's of the model

out(1) = param.tscale.*(Jipr - Jserca + param.delta.*(Jin - Jpm))./param.tau_c;
out(2) = param.tscale.*param.delta.*(Jin - Jpm);
out(3) = param.tscale.*(h_inf - h)./tau;
out(4) = param.tscale.*param.test.*(L - param.k_i.*p);
out = out';

end