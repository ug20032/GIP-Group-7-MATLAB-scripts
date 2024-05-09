function [P,IC] = initial_conditions(Des,Prop,d_break,D_eq,sim_version)

        if sim_version == 1
            V_exit = Des.Q / (pi*Des.R_pipe^2);
        elseif sim_version == 2
            V_exit = sqrt(2*Prop.g*D_eq);
        end

        IC.z = d_break;
        IC.dz_dt = V_exit*sin(Des.theta_rad) + 1e-3; % dx_dt_0 cannot be 0
        IC.x = 0;
        IC.dx_dt = V_exit*cos(Des.theta_rad) + 1e-3;
       
        IC.d2z_dt2 = 0;
        IC.d2x_dt2 = 0;
        
        % Initial Gas Conditions (Ideal Gas Law)
        IC.T_air = Des.T_air_c + 273;
        P = Prop.P_atm + Prop.rho_l*Prop.g*Des.H_m;
        IC.V_B_air = 4/3 * pi * (D_eq/2)^3;
        IC.m_air = (P*IC.V_B_air) / (Prop.R_air*IC.T_air); % mass of air in bubble
        IC.m_vapour = 0; % mass of water vapour in bubble

end