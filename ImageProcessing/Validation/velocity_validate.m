clc
clear
close all

iter = '1b';
formatting
fig = figure;
fig.Position = [100 100 600 450];
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

print(['validation_vel_' num2str(iter)], '-dpng', '-r300');
