clc
clear
close all

%% User selection
sim_version = 2; % either initial (1) or updated (2)
iter = "3b"; % either iteration 1a, 1b, 2a, 2b, 3a or 3b
T_water_c = 85; % Temperature of Water (deg C) - either 25 for bubble dynamics validation or 85 for mass/heat transfer validation

%% Design Parameters 

if iter == "1a"
    n_holes_list = [8 3]; % total number of holes vs actual number of holes working
    D_mm = 1;
elseif iter == "1b"
    n_holes_list = [31 7];
    D_mm = 0.7;
elseif iter == "2a"
    n_holes_list = [8 3];
    D_mm = 1;
elseif iter == "2b"
    n_holes_list = [24 7];
    D_mm = 0.7;
elseif iter == "3a"
    n_holes_list = [8 3];
    D_mm = 0.7;
elseif iter == "3b"
    n_holes_list = [24 10];
    D_mm = 0.3;
end
n_holes = n_holes_list(sim_version);

Q_L = 0.5 / n_holes ; % Flow Rate of N2 (L/min)
H_mm = 84; % Height of Water in Flask (mm)
theta = 90; % Angle of injection - vertical = 90, radial = 0

% Design parameters and fluid properties
Des = conv2SI(D_mm,H_mm,Q_L,T_water_c,theta);
fluidproperties
        
%% Initial Conditions
% Initial conditions
[R_bubble] = bubble_formation(Des.Q,Des.R_pipe,Prop,sim_version);
d_break = 0;
D_eq = 2*R_bubble;
[P,IC] = initial_conditions(Des,Prop,d_break,D_eq,sim_version);  

%% Solve system of ODEs (use Euler ODE variable time step)
dt_0 = 1e-4; % baseline timestep
[tlist,z,dz_dt,x,dx_dt,D_eq,AH,T_air,F,r_a] = sysODE(0,dt_0,Des,Prop,IC);

%% Performance Parameters
[~, AH_sat_max] = water_uptake(0,0,0,T_water_c,P); % Calculate maximum humidity
eta = AH(end) / AH_sat_max;
m_dot_water = AH(end) * Q_L;

v_av = max(z) / max(tlist);

disp(['Mass Flow Rate of Water [per hole]: ' num2str(round(m_dot_water,3)) ' g/min'])
disp(['Mass Flow Rate of Water [total]: ' num2str(round(m_dot_water*n_holes,3)) ' g/min'])
disp(newline)
disp(['Average major semi-axis: ' num2str(round(mean(r_a)*1e3,1)) ' mm'])
disp(['Average velocity: ' num2str(round(v_av*100,1)) ' cm/s'])

%% Plots
formatting

% Dynamics graph ----------------------------------------------------------
figure
hold on
plot(tlist,dz_dt*1e2,'k')
xlabel('Time (s)')
ylabel('Vertical velocity [$\dot{z}_c$] (cm/s)')
ylim([0 100])
yyaxis right
plot(tlist(1:end-1),D_eq*1e3,'Color',colour(1,:))
ylabel('Diameter of Bubble [$D_{eq}$] (mm)')
set(gca,'YColor',colour(1,:))

% Mass and Heat transfer graph --------------------------------------------
figure
plot(tlist,T_air-273,'Color','k')
xlabel('Time (s)')
ylabel('Bubble Temperature [$T_B$] ($^\circ$C)')
ylim([0 100])

yyaxis right
plot(tlist,AH,'Color',colour(1,:)) % 1g/L = 1kg/m3
hold on
yline(AH_sat_max,'Color',colour(1,:),'Linestyle','--')
xlabel('Time (s)')
ylabel('Humidity [$H_{s,B}$] (g$_{H_2O}$/L$_{N_2}$)') % nominal litres -> i.e. at room temperature
legend('','Bubble Humidity','Saturation Humidity','location','NW')
set(gca,'YColor',colour(1,:))

% Forces graph ------------------------------------------------------------
figure
subplot(1,2,1)
hold on
plot(tlist(1:end-1),-F.dragx,'Color',colour(1,:))
plot(tlist(1:end-1),F.buoyancyx,'Color',colour(2,:))
plot(tlist(1:end-1),-F.bassetx,'Color',colour(3,:))
set(gca,'xscale','log')
grid on
legend('Drag','Buoyancy','Basset History')
xlabel('Time (s)')
ylabel('Force (N)')
title('X-direction')

subplot(1,2,2)
hold on
plot(tlist(1:end-1),-F.dragy,'Color',colour(1,:))
plot(tlist(1:end-1),F.buoyancyy,'Color',colour(2,:))
plot(tlist(1:end-1),-F.bassety,'Color',colour(3,:))
set(gca,'xscale','log')
grid on
legend('Drag','Buoyancy','Basset History')
xlabel('Time (s)')
ylabel('Force (N)')
title('Y-direction')

% Trajectory Graph --------------------------------------------------------
if theta ~= 90
    figure
    plot(x,z,'Color',colour(1,:))
    axis equal
    title('Trajectory')
end




