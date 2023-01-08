clear
close all
clc

%% Parameter values

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.62; param.Kplc = 0.11;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3; param.k_i = 2; 
param.tau_c = 2;
param.tscale = 0.10;

%% ODE solver

rm_solver = @(x,t)PLC_final(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-30);
[~,Y] = ode15s(rm_solver,[0 1000],[0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(rm_solver,[0 200],Y(end,:),opts);

%% Results/Plots

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

figure(1)
plot(T,c,'b','LineWidth',2)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
% axis([0 4000 0 0.8])
ax=gca;
set(ax,'Linewidth',4)
ax.FontSize=20;
box off
hold off
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Figure9c','epsc')

[pks,locs] = findpeaks(c);
locs = locs(pks>0.5);
hold on
plot(T(locs),Y(locs,1),"*")

per_Y = Y(locs(end-1):locs(end),:);
per_T = T(locs(end-1):locs(end));

figure(2)
plot(per_T,per_Y)

for i=1:4
    per(:,i) = interp1(per_T,per_Y(:,i),...
    linspace(T(locs(end-1)),T(locs(end)),10000),'spline');
end

T_per = linspace(T(locs(end-1)),T(locs(end)),10000);

figure(3)
plot(T_per,per)

Z = [T_per',per];
fileID = fopen('L62.dat','w');
fprintf(fileID,'%8d %8d %8d %8d %8d\n',Z');
fclose(fileID);

tilefigs

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
Jin = param.alpha0 + param.alpha1*param.Kce^4./(ce.^4 + param.Kce^4);
Jpm = param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = (Jipr - Jserca + param.delta.*(Jin - Jpm))./param.tau_c;
out(2) = param.delta.*(Jin - Jpm);
out(3) = (h_inf - h)./tau;
out(4) = L - param.k_i.*p;
out = out';

end