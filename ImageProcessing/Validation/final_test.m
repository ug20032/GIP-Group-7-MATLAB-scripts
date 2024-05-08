clc
clear
close all

formatting

predicted = [0.119 0.144];


iter3a = [0.102 0.125 0.099];
iter3b = [0.149 0.090];

experimental = [mean(iter3a) mean(iter3b)];

err_pos = experimental - [max(iter3a) max(iter3b)];
err_neg = experimental - [min(iter3a) min(iter3b)];


dx = 0.1;
x = [0+dx 1+dx;
     0-dx 1-dx];

figure


hold on

bar(x(2,:),experimental,'FaceColor',colour(1,:),'EdgeColor','none','FaceAlpha',0.75,'BarWidth',0.25)

errorbar(x(2,1),experimental(1),err_neg(1),err_pos(1),'k','Linewidth',1.5)
errorbar(x(2,2),experimental(2),err_neg(2),err_pos(2),'k','Linewidth',1.5)

bar(x(1,:),predicted,'FaceColor',colour(3,:),'EdgeColor','none','FaceAlpha',0.75,'BarWidth',0.25)

legend('Experimental','','','Predicted','location','northeastoutside')
xticks([0 1])
xticklabels({'3a','3b'})
xlabel('Iteration')
ylabel('Mass flow rate of water (g/min)')