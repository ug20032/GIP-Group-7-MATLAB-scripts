clc
clear
close all


formatting

data =  [3.7 1.1 2.8;
         3.7 0.8 1.9;
         3.9 1.1 2.8;
         3.8 0.8 1.9;
         3.9 0.8 2.8;
         2.9 0.3 1.6];

dx = 0.25;
x = [0:5;
     (0+dx):(5+dx);
     (0-dx):(5-dx);];

fig = figure;
fig.Position = [100 100 800 450];
hold on
for i = 1:3
    bar(x(i,:),data(:,i),'FaceColor',colour(i,:),'EdgeColor','none','FaceAlpha',0.75,'BarWidth',0.3)
end
xticks(0:5)
xlim([-1 6])
xticklabels({'1a','1b','2a','2b','3a','3b'})
ylabel('Major semi-axis $\bar{r}_a$ (mm)')
legend('Experimental','Initial','Updated','location','northeastoutside','Interpreter','latex')
xlabel('Iteration')

%%
print('validation_size', '-dpng', '-r300');

