close all
clear
clc

load('sliceA06.mat')

plot(sliceA06((1:256),1),sliceA06((1:256),2),'r','LineWidth',1)
hold on
plot(sliceA06((257:522),1),sliceA06((257:522),2),'k--','LineWidth',1)
hold on
plot(sliceA06((523:671),1),sliceA06((523:671),2),'r','LineWidth',1)
hold on
plot(sliceA06((672:677),1),sliceA06((672:677),2),'b--','LineWidth',2)
hold on
plot(sliceA06((672:677),1),sliceA06((672:677),3),'b--','LineWidth',2)
hold on
plot(sliceA06((678:905),1),sliceA06((678:905),2),'g','LineWidth',2)
hold on
plot(sliceA06((678:905),1),sliceA06((678:905),3),'g','LineWidth',2)
hold on
plot(sliceA06((906:end),1),sliceA06((906:end),2),'b--','LineWidth',2)
hold on
plot(sliceA06((906:end),1),sliceA06((906:end),3),'b--','LineWidth',2)
axis([0 6 0 0.6])
box off
xlabel('C_t \muM')
ylabel('Ca^{2+}_i \muM')
set(gca,'FontSize',30,'fontweight','b','fontname','arial')