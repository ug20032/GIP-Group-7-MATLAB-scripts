clc
clear
close all

formatting


%% Compare iteration 3a and 3b actual
load("iteration3a_veldata.mat")
figure
hold on
scatter(t*1e3,y_pos,'MarkerEdgeColor',colour(1,:),'Marker','.')
plot(t*1e3,yfit,'Color',colour(2,:))
load("iteration3b_veldata.mat")
scatter(t*1e3,y_pos,'MarkerEdgeColor',colour(3,:),'Marker','.')
plot(t*1e3,yfit,'Color',colour(4,:))
legend('3a experiment','3a fit','3b experiment','3b fit','Location','SouthEast')

xlabel('Time (ms)')
ylabel('Displacement of highest bubble (mm)')

%% Compare iteration 3b actual, old and new

% Experiment
load("iteration3a_veldata.mat")
figure
hold on
scatter(t*1e3,y_pos,'MarkerEdgeColor',colour(1,:),'Marker','.')
disp(['Actual mean velocity (m/s)' num2str(y_pos(end) ./ t(end)/1000)])

% Old
load("Iteration3a_sim_old.mat")
y_sim_old = 0;
for i = 1:length(tlist)-1
    y_sim_old(i+1) = y_sim_old(i) + dx_dt(i) * (tlist(i+1)-tlist(i));
end
plot(tlist*1e3,y_sim_old*1e3,'Color',colour(2,:))
disp(['Simulated mean velocity [OLD] (m/s)' num2str(max(y_sim_old) ./ max(tlist))])

% New
load("Iteration3a_sim.mat")
y_sim_new = 0;
for i = 1:length(tlist)-1
    y_sim_new(i+1) = y_sim_new(i) + dx_dt(i) * (tlist(i+1)-tlist(i));
end
plot(tlist*1e3,y_sim_new*1e3,'Color',colour(3,:))
disp(['Simulated mean velocity (m/s)' num2str(max(y_sim_new) ./ max(tlist))])

% Graph formatting
legend('Experimental','Initial','Updated','Location','SE','interpreter','Latex')
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('Time (ms)')
ylabel('Displacement of Bubble $z_c$ (mm)')
xlim([0 round(max(t)*1e1)*1e2])

%% 

sizedata = [3.06 0.79 2.96;
            2.21 0.33 1.69];

x = [0 1];
bwidth = 0.25;
alpha = 0.95;

figure
hold on
bar(x,sizedata(:,1),'BarWidth',bwidth,'FaceColor',colour(1,:),'EdgeColor','none','FaceAlpha',alpha)
bar(x-0.2,sizedata(:,2),'BarWidth',bwidth,'FaceColor',colour(2,:),'EdgeColor','none','FaceAlpha',alpha)
bar(x+0.2,sizedata(:,3),'BarWidth',bwidth,'FaceColor',colour(3,:),'EdgeColor','none','FaceAlpha',alpha)

legend('Experimental','Initial','Updated','Location','NE','interpreter','Latex')
xlabel('Time (ms)')
ylabel('Major semi-axis of bubble $r_a$ (mm)')
ylim([0 4])
xticks([0 1])
set(gca, 'TickLabelInterpreter', 'latex');
xticklabels(gca, {'3.a','3.b'})
xlabel('Iterations')



