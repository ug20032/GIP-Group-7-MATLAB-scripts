clc
clear
close all

%% Velocity validation
iter = '1b';
formatting
fig = figure;
fig.Position = [50 200 600 450];
hold on
load(['iter' iter '_vel.mat'])

scatter(t*1e3,z_c,'MarkerEdgeColor',colour(1,:),'Marker','.')
load(['iter' iter '_vel_sim1.mat'])
plot(tlist*1e3,z*1e3,'Color',colour(2,:),'Linewidth',1.5)
load(['iter' iter '_vel_sim2.mat'])
plot(tlist*1e3,z*1e3,'Color',colour(3,:),'Linewidth',1.5)

xlabel('Time (ms)')
ylabel('Displacement $z_c$ (mm)')
ylim([0 100])
xlim([0 400])
legend('Experimental','Initial model','Updated model','location','SE','Interpreter','latex')

%% Size Validation

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
fig.Position = [700 200 800 450];
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