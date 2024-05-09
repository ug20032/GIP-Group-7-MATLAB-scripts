clc
clear
close all

%% User selection
sim_version = 2; % either initial (1) or updated (2)

%% Determine the optimum diameter tube
nn_list = [1 2 3 5 7 10 20 30 50 70 100 200 300];
Q_L_list = 0.5 ./ nn_list; % Flow Rate of N2 (L/min)
T_water_list = [80 85 90]; % Temperature of Water (deg C)

count = 0;
for bb = 1:length(T_water_list)
    for aa = 1:length(nn_list)
        count = count + 1;
        Q_L = Q_L_list(aa);
        T_water_c = T_water_list(bb);
        [m_dot_water,D_mm_list] = main_OPTIMISE_func(Q_L,T_water_c,sim_version);
        if mod(count,1) == 0
            disp([' ---------- ' num2str(bb) '/' num2str(length(T_water_list)) ' and ' num2str(aa) '/' num2str(length(nn_list)) ' ---------- '])
        end
        m_dot_max(aa,bb) = max(m_dot_water);
        d_opt_list = D_mm_list(m_dot_max(aa,bb) == m_dot_water);
        d_opt(aa) = max(d_opt_list); % more than one optimum - but choose largest hole for ease of manufacture
            
    end
end

%% Plots
formatting
figure
plot(nn_list,d_opt,'k')
xlabel('Number of tubes $n_t$')
ylabel('Optimum diameter of tubes $D_t$ (mm)')
grid on
set(gca,'xscale','log')
xticks([1 2 5 10 20 50 100 200])
xlim([1 200])
ylim([0 2.5])
yticks([0:0.5:2.5])

yyaxis right
hold on
for ii = 1:length(T_water_list)
    plot(nn_list,m_dot_max(:,ii).*nn_list','Color',colour(ii,:),'LineStyle','-')
end
ylabel('Total mass flow rate $\dot{m}_{w,tot}$ (g/min)')

set(gca,'ycolor',colour(1,:))
lgd = legend([''; string(split(num2str(T_water_list)))],'location','NE');
lgd.Title.String = 'Temperature $T_l$ ($^\circ$C)';
lgd.Title.Interpreter = 'latex';
ylim([0.05 0.3+1e-5])
yticks([0:0.05:1])

%% Functions
function [m_dot_water,D_mm_list] = main_OPTIMISE_func(Q_L,T_water_c,sim_version)

    % sim_version - either initial (1) or updated (2)
    % Q_L - Flow Rate of N2 (L/min)
    % T_water_c - Temperature of Water (deg C)
    theta = 90; % angle of injection
    
    n_list = 200;
    D_mm_list = logspace(-1,1,n_list); % Diameter of Pipe (mm) 
    H_mm_list = [84]; % Height of Water in Flask (mm)
    
    % Iterate through different diameters and water levels
    for jj = 1:length(H_mm_list)
        for ii = 1:length(D_mm_list)
        
            D_mm = D_mm_list(ii); % Diameter of Pipe (mm) 
            H_mm = H_mm_list(jj); % Height of Water in Flask (mm)
    
            % Design parameters and fluid properties
            Des = conv2SI(D_mm,H_mm,Q_L,T_water_c,theta);
            fluidproperties
            
            % Initial Conditions
            % Initial conditions
            [R_bubble] = bubble_formation(Des.Q,Des.R_pipe,Prop,sim_version);
            d_break = 0;
            D_eq = 2*R_bubble;
            [P,IC] = initial_conditions(Des,Prop,d_break,D_eq,sim_version);       
            
            % Solve system of ODEs (use Euler ODE variable time step)
            dt_0 = 1e-3; % baseline timestep
            [~,~,~,~,~,D_eq,H_sb,~,~,~] = sysODE(0,dt_0,Des,Prop,IC);
            
            % Performance Parameters
            [~, H_sat_max] = water_uptake(0,0,0,T_water_c,P); % Calculate maximum humidity
            eta = H_sb(end) / H_sat_max;
            m_dot_water(ii,jj) = H_sb(end) * Q_L;
            D_av(ii,jj) = mean(D_eq);
            
        end
    end

end