clear
close all
clc

%%   
% Parameter values

param_ct.Kf = 10; param_ct.Kc = 0.2; param_ct.Kp = 0.2; param_ct.Kb = 0.4;
param_ct.tau_max = 100; param_ct.K_tau = 0.09; param_ct.Kh = 0.08;
param_ct.Vs = 0.9; param_ct.Kbar = 1.957e-5; param_ct.Ks = 0.2;
param_ct.tau_p = 1; param_ct.R_act = 0.4; param_ct.K_PLC = 0.1;
param_ct.Vpm = 0.11; param_ct.Kpm = 0.3;
param_ct.alpha0 = 0.0027; param_ct.alpha1 = 0.015; param_ct.Kce = 14;
param_ct.delta = 2.5; param_ct.gamma = 5.5; param_ct.tau_cdum = 2;

%%
% ODE solver

R_act = 0.05:0.005:0.22;
% 
for j = 1: length(R_act)
    hep_1 = @(x,t)hep_SOCC_ct_new(x,t,param_ct,R_act(j));
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [~,Y] = ode15s(hep_1,[0 300],[0.078,0.525247,3.28396,0],opts);
    [T,Y] = ode15s(hep_1,0:0.01:200,Y(end,:),opts);
    
    c = Y(:,1);
    h = Y(:,2);
    ct = Y(:,3);
    p = Y(:,4);
    
    [pks,locs] = findpeaks(c);
    CP = c(locs(end-2):locs(end));
    Time = T(locs(end-2):locs(end));
    [~,~,w]=findpeaks(CP,Time);
    Width(j,1) = w;
    Period(j,1) = T(locs(end))-T(locs(end-1));
end

% figure(1)
% plot(R_act,Period,'LineWidth',3)
% xlabel('R_{act} (\muM)')
% ylabel('Period (sec)')
% ax=gca;
% ax.LineWidth=1.5;
% ax.FontSize=20;
% box off
% 
% 
% figure(2)
% plot(R_act,Width,'LineWidth',3)
% xlabel('R_{act} (\muM)')
% ylabel('Width (sec)')
% ax=gca;
% ax.LineWidth=1.5;
% ax.FontSize=20;
% box off


R = [0.36 .37 .38 .39 .4 .41 .42 .43 .44 .45 .46 .47 .48 .49 .5 .51];
wid_VP = [13.81 12.01 11.73 11.53 11.4 11.27 11.2 11.12 11.08 11.05...
         11.05 11.1 11.19 11.41 11.76 16.71];
per_VP = [32.79 30.73 30.14 29.67 29.25 28.87 28.51 28.17 27.88 27.6...
         27.35 27.13 26.97 26.91 27.4 31.83];

ave_PE_wid = mean(Width);
ave_VP_wid = mean(wid_VP);
nam = categorical({'PE','VP'});
comp = [ave_PE_wid ave_VP_wid];
bar(nam,comp)
% xlabel('R_{act} (\muM)')
ylabel('mean width (sec)')
ax=gca;
ax.LineWidth=1.5;
ax.FontSize=20;
box off

% figure(3)
% plot(R,wid_VP,'LineWidth',3)
% xlabel('R_{act} (\muM)')
% ylabel('Width (sec)')
% ax=gca;
% ax.LineWidth=1.5;
% ax.FontSize=20;
% box off
% 
% figure(4)
% plot(R,per_VP,'LineWidth',3)
% xlabel('R_{act} (\muM)')
% ylabel('Period (sec)')
% ax=gca;
% ax.LineWidth=1.5;
% ax.FontSize=20;
% box off
% 
% 
% [pks,locs] = findpeaks(c);
% 
% 
% CP = c(locs(end-4):locs(end));
% 
% 
% % 
% Time = T(locs(end-4):locs(end));
% 
% 
% % 
% [~,~,w]=findpeaks(CP,Time);
% % 
% figure()
% findpeaks(CP,Time,'Annotate','extent','WidthReference','halfheight','MinPeakProminence',0)
% % hold on
% findpeaks(CP,Time,'Annotate','extent')
% [pks,locs,w,p] = findpeaks(c,'Annotate','extent');
% %figure(17)
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

tilefigs
