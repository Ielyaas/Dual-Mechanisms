clear 
close all
clc

%% Function file


%% Parameter file for open cell model with ct as a varaiable

param_ct.Kf = 10; 
param_ct.Kc = 0.2; 
param_ct.Kp = 0.2; 
param_ct.Kb = 0.4;
param_ct.tau_max = 100; 
param_ct.K_tau = 0.09; 
param_ct.Kh = 0.08;
param_ct.Vs = 0.9; 
param_ct.Kbar = 1.957e-5; 
param_ct.Ks = 0.2;
param_ct.tau_p = 1; 
param_ct.R_act = 0.25; 
param_ct.K_PLC = 0.06;
param_ct.Vpm = 0.11; 
param_ct.Kpm = 0.3;
param_ct.alpha0 = 0.0027; 
param_ct.alpha1 = 0.015; 
param_ct.Kce = 14;
param_ct.delta = 2.5; 
param_ct.gamma = 5.5;
param_ct.epsilon = 1;

%% ODE script

hep_ct = @(x,t)hep_SOCC_ct1(x,t,param_ct);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
% for ct_i = 2.8:0.1:4.2
%     ct_i
[~,Y] = ode15s(hep_ct,[0 1000],[0.078,0.525273,3.283961,0],opts);
% [m,i] = max(Y(:,3));
[T,Y] = ode15s(hep_ct,0:0.1:200,Y(end,:),opts);

%% Plots

c = Y(:,1);
h = Y(:,2);
ct = Y(:,3);
p = Y(:,4);

load('slice51260318.mat')

figure(1)
hold on
plot(T,c,'k','LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% h=vline2(in_block,'r');
set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% axis([0 500 0 1])
box off
% hold on
% plot(T,Y(:,3),'k','LineWidth',1)

% figure(2)
% hold on
% plot(Y(:,3),Y(:,1),'k','LineWidth',1)
% xlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% % axis([0 4 0 1])
% box off
% hold on
% plot(slice51260318((1:126),1),slice51260318((1:126),2),'r','LineWidth',1)
% hold on
% plot(slice51260318((127:286),1),slice51260318((127:286),2),'k--','LineWidth',1)
% hold on
% plot(slice51260318((287:987),1),slice51260318((287:987),2),'r','LineWidth',1)
% axis([0 4 0 1])
% 
% figure(3)
% plot3(Y(:,1),Y(:,2),Y(:,3),'LineWidth',1)
% zlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
% xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
% set(gca,'FontSize',18,'fontweight','b','fontname','arial')

% figure(4)
% plot(T,Y(:,2),'k','LineWidth',1)

tilefigs

%% Manifold plot

% Create the mesh

load('mesh2.mat')  
    

%[c,ct]=meshgrid(linspace(0.01,1,50),linspace(0.01,1,50));
[c,ct]=meshgrid(cv,ctv);

Cond = c < ct;

c = c(Cond);
ct = ct(Cond);
triangulation = delaunay(c,ct);

% Functions

ce = param_ct.gamma.*(ct - c);

PLC = param_ct.R_act.*c.^2./(param_ct.K_PLC^2 + c.^2);
p = PLC;

phi_c = c.^4./(c.^4 + param_ct.Kc^4);
phi_p = p.^2./(p.^2 + param_ct.Kp^2);
phi_p_down = param_ct.Kp^2./(p.^2 + param_ct.Kp^2);

Jpm = param_ct.Vpm.*c.^2./(param_ct.Kpm^2 + c.^2);
Jin = param_ct.alpha0 + param_ct.alpha1.*(param_ct.Kce^4./(param_ct.Kce^4+ce.^4));

h_inf = param_ct.Kh^4./(param_ct.Kh^4 + c.^4);
tau = param_ct.tau_max.*param_ct.K_tau^4./(param_ct.K_tau^4 + c.^4);

Jserca = param_ct.Vs.*(c.^2 - param_ct.Kbar.*ce.^2)./(c.^2 + param_ct.Ks^2);

alpha = phi_p_down.*(1-phi_c.*h_inf);

f = Jserca - param_ct.epsilon.*param_ct.delta.*(Jin - Jpm);
h = (param_ct.Kb.*f.*alpha)./(phi_p.*phi_c.*(param_ct.Kf.*(ce-c) - 1.4*f));
figure(4)
%trisurf(triangulation,c,ct,h,'lines','none','facecolor',[1 0.1 0.1],'facealpha',0.3)
trisurf(triangulation,c,ct,h);
xlabel('c');
ylabel('ct');
zlabel('h');
% hold on
% plot3(Y(:,3),Y(:,2),Y(:,1),'LineWidth',1)

%% OPEN cell

% param_ct.epsilon = 1;
hep_ct = @(x,t)hep_SOCC_ct1(x,t,param_ct);
[T,Y] = ode15s(hep_ct,[0 200],[0.078,0.525273,3.5,0],opts);
[m,i] = max(Y(:,3));
[T,Y] = ode15s(hep_ct,0:0.001:1000,Y(end,:),opts);
figure(1)
hold on
plot(T,Y(:,1),'k','LineWidth',1)
xlabel('time (s)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
% h=vline2(in_block,'r');
set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% axis([0 500 0 1])
box off

figure(2)
hold on
plot(Y(:,3),Y(:,1),'k','LineWidth',1)
xlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',40,'fontweight','b','fontname','arial')
% axis([0 4 0 1])
box off

figure(6)
hold on
plot3(Y(:,1),Y(:,2),Y(:,3),'LineWidth',1)
zlabel('Ca^{2+}_t \muM','fontsize',20,'fontweight','b','fontname','arial')
xlabel('Ca^{2+}_i \muM','fontsize',20,'fontweight','b','fontname','arial')
ylabel('h','fontsize',20,'fontweight','b','fontname','arial')
set(gca,'FontSize',18,'fontweight','b','fontname','arial')
