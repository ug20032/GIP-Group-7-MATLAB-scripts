% Water Vapour Uptake function
function [m_dot,AH_sat] = water_uptake(S_B,V_B,m_water,T,P)

    % https://www.engineeringtoolbox.com/maximum-moisture-content-air-d_1403.html

    % Mass Transfer Coefficient (Liss, 1972)
    K_c_cm = 1000; % cm/hr
    K_c_m = K_c_cm / 1e2 / 3600; % m/s
    
    % Property Table
    properties=[15 1.7057;
                20 2.3392 
                25 3.1698 
                30 4.2469 
                35 5.6291
                40 7.3851 % T (deg C), P_sat (kPa)
                45 9.5953 
                50 12.352 
                55 15.763 
                60 19.947 
                65 25.043 
                70 31.202 
                75 38.597 
                80 47.416
                85 57.868 
                90 70.183];
    P_sat = interp1(properties(:,1),properties(:,2),T); % Sat. Pressure (kPa)
    P_atm = P / 1e3; % Pa -> kPa
    
    R_air = 287; R_water = 461; % J/mol/K
    SH = (R_air/R_water) * P_sat / (P_atm - P_sat); % Specific Humidity: mass of water / mass of air
    rho_air = ((P_atm - P_sat)*1e3) / (R_air*(273 + T)); % density of air (kg/m3)

    % Absolute humidity of water vapour
    AH_sat = SH*rho_air; % (kg of water / m3 air)

    % Absolute humidity of Bubble
    AH_B = m_water / V_B; % (kg of water / m3 air)

    m_dot = K_c_m*S_B*(AH_sat - AH_B);
end
