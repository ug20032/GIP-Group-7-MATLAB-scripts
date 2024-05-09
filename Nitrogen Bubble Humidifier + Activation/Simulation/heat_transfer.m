function [T_B] = heat_transfer(Des,Prop,IC,dt,dx_dt,D_eq,T_B,m_v,S_B)

    % Updates the temperature of the bubble (T_B) from the velocity (dx_dt),
    % equivalent spherical diameter (D_eq), mass of water vapour in bubble (m_v)
    % and total surface area (S_B)


    Re_l = Prop.rho_l.*dx_dt.*D_eq./Prop.mu_l;
    Pr_l = Prop.Cp_vapour.*Prop.mu_l./Prop.k_l;
    Pe_l = Re_l.*Pr_l;
    
    % Sideman, 1966 - Direct Contact Between Immiscible Liquids
    Nu1 = 1.13*Pe_l^(1/2); % eq. 18
    h_c = Prop.k_l/D_eq * Nu1; % Heat transfer across interface

    Re_g = Prop.rho_g.*dx_dt.*D_eq./Prop.mu_g;
    Pr_g = Prop.Cp_air.*Prop.mu_g./Prop.k_g;
    Pe_g = Re_g.*Pr_g;

    Nu2 = 0.0375*Pe_g; % eq. 28
    h_d = Prop.k_g/D_eq * Nu2; % Heat transfer inside bubble

    h = (1/h_c + 1/h_d)^(-1);

    Q_dot = h*S_B*(Des.T_water - T_B); % heat transfer rate (W)
    delta_T = (Q_dot*dt) / (IC.m_air*Prop.Cp_air + m_v*Prop.Cp_vapour); % change in temperature (K)
    T_B = T_B + delta_T;

end