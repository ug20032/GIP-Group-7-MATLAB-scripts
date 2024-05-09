clc
clear
close all

%% User selection
sim_version = 1; % either initial (1) or updated (2)
n_holes = 8; % either 8, 24 or 31

%% Design Parameters 
formatting

Q_L = 0.5/n_holes; % Flow Rate of N2 (L/min)
T_water_c = 85; % Temperature of Water (deg C)
theta = 90; % angle of injection

n_list = 100;
D_mm_list = logspace(-1,1,n_list); % Diameter of Pipe (mm) 
H_mm_list = [84]; % Height of Water in Flask (mm)

%% Iterate through different diameters and water levels
for jj = 1:length(H_mm_list)
    for ii = 1:length(D_mm_list)
    
        disp([num2str(jj) ',' num2str(ii) ' / ' num2str(length(H_mm_list)) ',' num2str(length(D_mm_list))])
        D_mm = D_mm_list(ii); % Diameter of Pipe (mm) 
        H_mm = H_mm_list(jj); % Height of Water in Flask (mm)

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
        m_dot_water(ii,jj) = AH(end) * Q_L;
        D_av(ii,jj) = mean(D_eq);
        
    end
end
disp('----------------- SIMULATION COMPLETE -----------------')

%% Plots
formatting
figure
hold on; grid on;
for i = 1:jj
    plot(D_mm_list,m_dot_water(:,i),'Color','k')
end
xlabel('Diameter of tube $D_t$ (mm)')
ylabel('Mass flow rate of water $\dot{m}_w$ (g/min)')
lgd = legend(split(num2str(H_mm_list)));
lgd.Title.String = 'Height of water (mm)';
lgd.AutoUpdate = 'off';
lgd.Visible = 'off';
yline(AH_sat_max*Q_L,'k--')
xticks([0.2 0.5 1 2 5]); xlim([0.1 8]);
set(gca,'xscale','log'); 

yticks(0:0.005:0.025)

yyaxis right
plot(D_mm_list,D_av(:,1)*1e3,'Color',colour(1,:))
ylabel('Average bubble diameter $\bar{D}_{eq}$ (mm)')
set(gca,'YColor',colour(1,:))
ylim([0 10])
