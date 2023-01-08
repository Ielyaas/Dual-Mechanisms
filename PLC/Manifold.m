clear
close all
clc

%% PARAMETER VALUES

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.6; param.Kplc = 0.11;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3; param.k_i = 2; 
param.test = 1; param.tau_c = 2; param.in = 1;
param.tscale = 1; param.out = 1; param.epsilon = 0.1;

%% 2D MANIFOLD

ce = @(c,ct) param.gamma.*(ct - c);
p = @(c) (param.Vplc.*c.^2./(param.Kplc^2 + c.^2))./param.k_i; % p=PLC

phi_c = @(c) c.^4./(c.^4 + param.Kc.^4);
phi_p = @(c) p(c).^2./(p(c).^2 + param.Kp^2);
phi_p_down = @(c) param.Kp^2./(p(c).^2 + param.Kp^2);
beta = @(c,h) phi_p(c).*phi_c(c).*h;
alpha = @(c,h) phi_p_down(c).*(1-phi_c(c)*h);
Po = @(c,h) beta(c,h)./(beta(c,h) + param.Kb.*(beta(c,h) + alpha(c,h)));

Jipr = @(c,h,ct) param.Kf.*Po(c,h).*(ce(c,ct)-c);
Jserca = @(c,ct) param.Vs.*(c.^2 - param.Kbar*ce(c,ct).^2)./(c.^2 + param.Ks^2);

Jin = @(c,ct) param.alpha0 + param.alpha1.*(param.Kce^4./...
      (ce(c,ct).^4 + param.Kce^4));
Jpm = @(c) param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

f = @(ct,h,c) Jipr(c,h,ct) - Jserca(c,ct);
interval = [0.0001 2.5 0.0001 0.5 0.0001 1];

fs = fimplicit3(f,interval,'k','lines','none');
fs.FaceAlpha = 0.2;
fs.LineWidth = 0.01;
fs.MeshDensity = 100;
zlabel('c \muM')
xlabel('c_t \muM')
ylabel('h')
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=70;
box off
hold on
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Figure12b','epsc')

%% 1D MANIFOLD

load('1DMan.dat')

plot3(X1DMan(:,4),X1DMan(:,8),X1DMan(:,7),'Color',[0.466 0.674 0.188],'LineWidth',4)


%% ODE solver

rm_solver = @(x,t)PLC_final(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-30);
[~,Y] = ode15s(rm_solver,[0 5000],[0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(rm_solver,[0 2000],Y(end,:),opts);

%% Results/Plots

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

plot3(ct,h,c,'Color',[0.635 0.078 0.184],'LineWidth',3)

% 'Color',[0.96 0.73 0.91]
xlim([0 2.5])
ylim([0 0.5])
% saveas(gcf,'2Dmanifold','epsc')

%% FUNCTION

function out = PLC_final(~,in,param)

% Reduced PKC model, PLC, DAG and PKC assumed to be at quasi steady state.

c = in(1);
ct = in(2);
h = in(3);
p = in(4);

    
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

out(1) = param.tscale.*(Jipr - Jserca + param.epsilon*param.delta.*(Jin - Jpm))./param.tau_c;
out(2) = param.tscale.*param.epsilon.*param.delta.*(Jin - Jpm);
out(3) = param.tscale.*param.epsilon.*(h_inf - h)./tau;
out(4) = param.tscale.*param.test.*(L - param.k_i.*p);
out = out';

end