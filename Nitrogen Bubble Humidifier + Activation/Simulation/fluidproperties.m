%% Fluid Properties 
Prop.rho_g = 1.2; Prop.rho_l = 1000; 
Prop.delta_rho = Prop.rho_l - Prop.rho_g;
Prop.mu_l = 1e-3; Prop.mu_g = 1.5e-5;
Prop.g = 9.81;
Prop.gamma = 7.3e-2; % N/m
Prop.Cp_air = 1005; Prop.Cp_vapour = 4200; % J/kg/K
Prop.R_air = 287; Prop.R_vapour = 461.5;
Prop.P_atm = 1.01e5;
Prop.k_l = 0.6; Prop.k_g = 0.025;

Prop.T_air_c = 20; % Inlet Temperature of Air (deg C)