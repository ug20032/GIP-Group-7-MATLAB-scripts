%% Design Parameters 
sim_version = 2;

theta = 90; % angle of injection
D_mm = 1; % Diameter of Pipe (mm) 
H_mm = 84; % Height of Water in Flask (mm)

n_holes = 8;
Q_list = [0:0.1:0.5]/n_holes;
T_list = [22:90];

%% Iterate through different diameters and water levels
disp(' --------- Calculating Bubble Humidifier WVFR --------- ')
for jj = 1:length(Q_list)
    disp([num2str(jj) '/' num2str(length(Q_list))])
    for ii = 1:length(T_list)
    
        Q_L = Q_list(jj); % Flow Rate of N2 (L/min)
        T_water_c = T_list(ii); % Temperature of Water (deg C)

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
        dt_0 = 5e-4; % baseline timestep
        [tlist,z,dz_dt,x,dx_dt,D_eq,AH,T_air,F,r_a] = sysODE(0,dt_0,Des,Prop,IC);
        
        %% Performance Parameters
        [~, AH_sat_max] = water_uptake(0,0,0,T_water_c,P); % Calculate maximum humidity
        eta = AH(end) / AH_sat_max;
        m_dot_water(ii,jj) = AH(end) * Q_L;
        D_av(ii,jj) = mean(D_eq);
        
    end
end
m_dot_water_tot = m_dot_water * n_holes;
clearvars -except Q_list T_list m_dot_water_tot n_holes t_hold chem_sim_version
