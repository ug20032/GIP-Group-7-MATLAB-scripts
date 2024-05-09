function [R_B] = bubble_formation(Q,r_t,Prop,sim_version)
                
        % Calculates the initial equivalent spherical radius (R_B) of a bubble
        % based on the volumetric flow rate of nitrogen (Q) and radius of
        % tube (r_t)

        % Chesters, 1977 - LOW gas flow bubble formation
        L_cap = sqrt(Prop.gamma/(Prop.delta_rho*Prop.g));
        V_low = 2*pi*r_t*L_cap.^2 * (1 + (3*r_t/(2*L_cap))^(2/3));
        R_b_low = (3/(4*pi)*V_low)^(1/3);
        
        % Davidson, 1977 - MEDIUM gas flow bubble formation
        g = 9.81;
        V_int = 1.14 * (Q^2/g)^(3/5);
        R_b_int = (3/(4*pi)*V_int)^(1/3);
        
        % Eggers, 2008 - HIGH gas flow bubble formation
        V_high = 40.84*r_t^3;
        R_b_high = (3/(4*pi)*V_high)^(1/3);
       
        R_B = max(R_b_low,R_b_int); % Radius of (spherical) bubble at given flow rate    

        if sim_version == 1 % jetting was simulated in the original model
            R_B = min(R_B,R_b_high);
        end

end