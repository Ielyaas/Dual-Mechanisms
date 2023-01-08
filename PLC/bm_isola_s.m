clear
close all
clc


%% Parameter set

param.Kf = 4; param.Kc = 0.2; param.Kp = 0.2; param.Kb = 0.4;
param.tau_max = 80; param.K_tau = 0.09; param.Kh = 0.08;
param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;
param.V_plc = 0.525; param.Kplc = 0.1;
param.Vpm = 0.07; param.Kpm = 0.3;
param.alpha0 = 0.004; param.alpha1 = 0.01; param.Kce = 14;
param.delta = 2.5; param.gamma = 5.5; param.tau_c = 2;

%% ODE solver

bm_solve = @(x,t)hep_base_model(x,t,param);
opts = odeset('RelTol',1e-13,'AbsTol',1e-30);
[~,Y] = ode15s(bm_solve,[0 500],[0.0862898 4.137186 0.424888 0],opts);
[T,Y] = ode113(bm_solve,[0 300],Y(end,:),opts);

%% Results/Plots

c = Y(:,1);
ct = Y(:,2);
h = Y(:,3);
p = Y(:,4);

figure(1)
plot(T,c,'k','LineWidth',2)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
set(gca,'FontSize',20,'fontweight','b','fontname','arial')
% axis([0 500 0 1])
ax=gca;
set(ax,'Linewidth',3)
ax.FontSize=20;
box off

[pks,locs] = findpeaks(c);
locs = locs(pks>0.4);
hold on
plot(T(locs),Y(locs,2),"*")

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
fileID = fopen('bm_per525.dat','w');
fprintf(fileID,'%8d %8d %8d %8d %8d\n',Z');
fclose(fileID);

tilefigs