clear
close all
clc

%% PARAMETER FILE

param.Kf = 10; param.Kc = 0.2; param.Kp = 0.2; param.Kb = 0.4;
param.tau_max = 100; param.Ktau = 0.09; param.Kh = 0.08;
param.Vs = 0.9; param.Kbar = 1.957e-5; param.Ks = 0.2;
param.tau_p = 1; param.Vplc = 0.51; param.Kplc = 0.1;
param.Vpm = 0.11; param.Kpm = 0.3;
param.alpha0 = 0.0027; param.alpha1 = 0.015; param.Kce = 14;
param.delta = 2.5; param.gamma = 5.5;
param.cbuf = 2;

%% ODE sol

HEP_sol = @(x,t)HEP_test_fm(x,t,param);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[~,Y] = ode15s(HEP_sol,[0 200],[0.078 17.6328 0.5252 0],opts);
[T,Y] = ode15s(HEP_sol,[0 200],Y(end,:),opts);

c = Y(:,1);
ce = Y(:,2);
h = Y(:,3);
p = Y(:,4);

fig = figure;
figure(1)
left_color = [0 0 1];
right_color = [0 1 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(T,c,'b','LineWidth',2)
yyaxis right
plot(T,p,'g','LineWidth',2)
hold on
box off
ax = gca;
set(ax,'LineWidth',2)
ax.FontSize = 20;
xlabel('time (s)')
yyaxis left
ylabel('[Ca^{2+}] (\muM)')
yyaxis right
ylabel('[IP_3] (\muM)')

%% FUNCTION FILE

function out = HEP_test_fm(~,in,param)

c = in(1);
ce = in(2);
h = in(3);
p = in(4);

phi_c = c.^4./(c.^4 + param.Kc^4);
phi_p = p^2/(p^2 + param.Kp^2);
phi_p_down = param.Kp^2/(p^2 + param.Kp^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = param.tau_max*param.Ktau^4./(c.^4 + param.Ktau^4);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

L = param.Vplc*c.^2./(c.^2 + param.Kplc^2);
cdum = param.cbuf./(c + param.cbuf);

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);

Jin = param.alpha0 + param.alpha1*param.Kce^4./(ce.^4 + param.Kce^4);
Jpm = param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

out(1) = (Jipr - Jserca + param.delta.*(Jin - Jpm))./cdum;
out(2) = param.gamma.*(Jserca - Jipr);
out(3) = (h_inf - h)./tau;
out(4) = param.tau_p.*(L - p);
out = out';

end
