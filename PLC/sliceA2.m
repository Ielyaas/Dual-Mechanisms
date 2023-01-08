close all
clear
clc

load('sliceA2.mat')

plot(sliceA2((1:160),1),sliceA2((1:160),2),'r','LineWidth',1)
hold on
plot(sliceA2((161:362),1),sliceA2((161:362),2),'k--','LineWidth',1)
hold on
plot(sliceA2((363:880),1),sliceA2((363:880),2),'r','LineWidth',1)
hold on
plot(sliceA2((881:1028),1),sliceA2((881:1028),2),'b--','LineWidth',2)
hold on
plot(sliceA2((881:1028),1),sliceA2((881:1028),3),'b--','LineWidth',2)
hold on
plot(sliceA2((1029:1176),1),sliceA2((1029:1176),2),'g','LineWidth',2)
hold on
plot(sliceA2((1029:1176),1),sliceA2((1029:1176),3),'g','LineWidth',2)
hold on
plot(sliceA2((1177:end),1),sliceA2((1177:end),2),'b--','LineWidth',2)
hold on
plot(sliceA2((1177:end),1),sliceA2((1177:end),3),'b--','LineWidth',2)
axis([0 6 0 0.8])
box off
xlabel('C_t \muM')
ylabel('Ca^{2+}_i \muM')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')