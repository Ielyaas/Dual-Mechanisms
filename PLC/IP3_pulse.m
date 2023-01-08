clear 
close all
clc

%% PARAMETER FILE

param.Vplc = 0; param.p_s = 1;
param.Kf = 40; param.Kc = 0.2; param.Kp = 0.3; param.Kb = 0.4;
param.tau_max = 200; param.K_tau = 0.09; param.Kh = 0.1;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.Kplc = 0.11;
param.Vpm = 0.07; param.Kpm = 0.3;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.delta = 2.5; param.gamma = 5.5; param.tau_c = 2;
param.k_i = 2;
param.tau_p = 0.0027;
param.kL = 1;
param.tscale = 10;

%% ODE SOLVER

solver = @(x,t)PLC_model(x,t,param);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[T1,Y1] = ode15s(solver,[0 200],[0.0805296 3.861011 0.70395 0],opts);
[T2,Y2] = ode15s(solver,[200 1500],[0.0805296 3.861011 0.70395 0.03],opts);
[T3,Y3] = ode15s(solver,[1500 2500],[0.0805296 3.861011 0.70395 0.045],opts);
[T4,Y4] = ode15s(solver,[2500 5000],[0.0805296 3.861011 0.70395 0.045],opts);


% [T,Y] = ode113(solver,[0 500],Y(end,:),opts);

%% OUTPUT/PLOTS

T = [T1;T2;T3;T4];
Y = [Y1;Y2;Y3;Y4];

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

fig = figure;
figure(1)
l_color = [0 0 1];
r_color = [0 0.6 0];
set(fig,'defaultAxesColorOrder',[l_color; r_color])
yyaxis left
plot(T,c,'LineWidth',4)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
% title('Halving V_{PLC}')
% xlim([0 270])
% set(gca,'yticklabel',[])
hold on
yyaxis right
plot(T,p,'color',[0 0.6 0],'LineWidth',4)
xlabel('time (s)')
ylabel('[IP_3] \muM')
% set(gca,'yticklabel',[])
hold on
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',6)
ax.FontSize=50;
box off
% axis([xmin, xmax, ymin, ymax])
hold off
set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Testfig12','epsc')

%% FUNCTION FILE

function out = PLC_model(t,in,param)

c = in(1);
ct = in(2);
h = in(3);
p = in(4);

if t < 250
    param.Vplc = 0;
else
    param.Vplc = 0.6;
end

if t < 250
    param.tau_p = 0.0027;
else
    param.tau_p = 1;
end

ce = param.gamma.*(ct - c);

phi_c = c.^4./(c.^4 + param.Kc^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);

h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = param.tau_max.*param.K_tau^4./(c.^4 + param.K_tau^4);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

L = param.Vplc.*c.^2./(c.^2 + param.Kplc^2);

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);

Jin = param.alpha0 + param.alpha1*param.Kce^4./(ce.^4 + param.Kce^4);
Jpm = param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's

out(1) = param.tscale.*(Jipr - Jserca + param.delta.*(Jin - Jpm))/param.tau_c;
out(2) = param.tscale.*param.delta.*(Jin - Jpm);
out(3) = param.tscale.*(h_inf - h)/tau;
out(4) = param.tscale.*param.tau_p.*(param.kL.*L - param.k_i.*p);
out = out';


end