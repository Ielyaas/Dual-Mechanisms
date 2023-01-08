clear
close all
clc

%% LOAD FILE

load('VplcKplc_final.dat')

%% PLOT

Kplc = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.31 0.32 0.33 0.34];
Vplc = [0.285 0.2 0.17 0.19 0.21 0.245 0.28 0.29 0.297 0.307 0.316];

plot(VplcKplc_final(1:41,1),VplcKplc_final(1:41,2),'k','LineWidth',4)
hold on
plot(VplcKplc_final(47:90,1),VplcKplc_final(47:90,2),'k','LineWidth',4)
hold on
plot(VplcKplc_final(162:456,1),VplcKplc_final(162:456,2),'color','[0 0.447 0.741]','LineWidth',4)
hold on
plot(VplcKplc_final(921:1091,1),VplcKplc_final(921:1091,2),'color','[0 0.447 0.741]','LineWidth',4)
hold on
plot(VplcKplc_final(1092:1160,1),VplcKplc_final(1092:1160,2),'color','[0 0.447 0.741]','LineWidth',4)
hold on
plot(VplcKplc_final(1161:1224,1),VplcKplc_final(1161:1224,2),'color','[0 0.447 0.741]','LineWidth',4)
hold on
plot(Vplc,Kplc,'--sr','LineWidth',4)
text(0.1,0,'a','FontSize',70,'color','[0.4660 0.6740 0.1880]');
hold on
text(0.1,0.1,'b','FontSize',70,'color','[0.4660 0.6740 0.1880]');
hold on
text(0.6,0.11,'c','FontSize',70,'color','[[0.4660 0.6740 0.1880]');
hold on
text(0.6,0.3,'d','FontSize',70,'color','[0.4660 0.6740 0.1880]');
hold on
text(0.13,0.16,'i','FontSize',70);
hold on
text(0.8,0.4,'ii','FontSize',70);
hold on
text(0.3,0.4,'HB','FontSize',70,'color','[0 0.447 0.741]')
hold on
text(0.58,0.17,'HB','FontSize',70,'color','[0 0.447 0.741]')
hold on
text(0.78,0.1,'SNPO','FontSize',70,'color','[0 0 0]')
% plot(0.4,0.3,'o','MarkerSize',40,'MarkerEdgeColor','r','MarkerFaceColor','r')
axis([0 1 0 0.7])
box off
% grid minor
xlabel('V_{PLC} (\muM/s)')
ylabel('K_{PLC} (\muM)')
set(gca,'FontSize',70)
set(gca,'GridLineStyle','--')
ax=gca;
set(ax,'LineWidth',6)
ax.FontSize=70;
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'figure14','epsc')