close all
clear all
clc

load('psop.mat')
load('psopnp.mat')
load('hepcc.mat')

% plot(psop(:,4),psop(:,2),'k','LineWidth',1)
% hold on
% plot(psopnp(:,4),psopnp(:,2),'b','LineWidth',1)
% % hold on
plot(hepcc((778:1019),1),hepcc((778:1019),2),'r','LineWidth',1)
hold on
plot(hepcc((1020:1332),1),hepcc((1020:1332),2),'r--','LineWidth',1)
hold on
plot(hepcc((1333:1633),1),hepcc((1333:1633),2),'r','LineWidth',1)
hold on
plot(hepcc((1634:1879),1),hepcc((1634:1879),2),'r','LineWidth',1)
hold on
plot(hepcc((1880:2191),1),hepcc((1880:2191),2),'r--','LineWidth',1)
hold on
plot(hepcc((2192:2492),1),hepcc((2192:2492),2),'r','LineWidth',1)
hold on
plot(hepcc((2992:3131),1),hepcc((2992:3131),2),'g','LineWidth',2)
hold on
plot(hepcc((2992:3131),1),hepcc((2992:3131),3),'g','LineWidth',2)