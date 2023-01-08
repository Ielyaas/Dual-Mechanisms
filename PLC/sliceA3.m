close all
clear all
clc

load('sliceA3.mat')

plot(sliceA3((1:139),1),sliceA3((1:139),2),'r','LineWidth',1)
hold on
plot(sliceA3((140:320),1),sliceA3((140:320),2),'k--','LineWidth',1)
hold on
plot(sliceA3((321:885),1),sliceA3((321:885),2),'r','LineWidth',1)
hold on
plot(sliceA3((886:end),1),sliceA3((886:end),2),'b--','LineWidth',2)
hold on
plot(sliceA3((886:end),1),sliceA3((886:end),3),'b--','LineWidth',2)
axis([0 6 0 0.6])
box off
xlabel('C_t \muM')
ylabel('Ca^{2+}_i \muM')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')