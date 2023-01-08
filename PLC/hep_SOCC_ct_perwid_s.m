clear
close all
clc

%% Parameter values

param_ct.Kf = 10; param_ct.Kc = 0.2; param_ct.Kp = 0.2; param_ct.Kb = 0.4;
param_ct.tau_max = 100; param_ct.K_tau = 0.09; param_ct.Kh = 0.08;
param_ct.Vs = 0.9; param_ct.Kbar = 1.957e-5; param_ct.Ks = 0.2;
param_ct.tau_p = 1; param_ct.R_act = 0.4; param_ct.K_PLC = 0.1;
param_ct.Vpm = 0.11; param_ct.Kpm = 0.3;
param_ct.alpha0 = 0.0027; param_ct.alpha1 = 0.015; param_ct.Kce = 14;
param_ct.delta = 2.5; param_ct.gamma = 5.5; param_ct.tau_cdum = 2;

%% ODE solver

hep_1 = @(x,t)hep_SOCC_ct_perwid(x,t,param_ct);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[~,Y] = ode15s(hep_1,[0 300],[0.078,0.525247,3.28396,0],opts);
[T,Y] = ode15s(hep_1,0:0.01:200,Y(end,:),opts);

c = Y(:,1);
h = Y(:,2);
ct = Y(:,3);
p = Y(:,4);

[pks,locs] = findpeaks(c);
CP = c(locs(end-10):locs(end));
Time = T(locs(end-10):locs(end));

figure()
[pks,locs,w,p] = findpeaks(CP,Time,'Annotate','extent','WidthReference','halfheight','MinPeakProminence',0);
findpeaks(CP,Time,'Annotate','extent','WidthReference','halfheight','MinPeakProminence',0)

per = mean(diff(locs));
ind = (c>0.1);
width = find(ind(2:end)-ind(1:end-1));