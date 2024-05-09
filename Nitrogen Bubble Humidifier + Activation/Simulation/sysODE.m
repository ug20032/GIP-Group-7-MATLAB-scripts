function [t,z,dz_dt,x,dx_dt,D_eq,AH,T_air,F,r_a] = sysODE(t,dt_0,Des,Prop,IC)

    % solves the displacement, heat and mass transfer system of ODEs

    % Set initial conditions
    z = IC.z; dz_dt = IC.dz_dt; d2z_dt2 = IC.d2z_dt2;
    x = IC.x; dx_dt = IC.dx_dt; d2x_dt2 = IC.d2x_dt2;
    T_air = IC.T_air; m_vapour = IC.m_vapour;

    i = 1;
    while z(i) < Des.H_m

        % Time
        dt = dt_0 / dz_dt(i);
        t(i+1) = t(i) + dt;
    
        % Pressure
        P = Prop.P_atm + Prop.rho_l*Prop.g*(Des.H_m - z(i));
        
        % Bubble Volume (air & water composition)
        V_B_air = idealgas(IC.m_air,Prop.R_air,T_air(i),P);
        V_B_vapour = idealgas(m_vapour(i),Prop.R_vapour,Des.T_water,P);
        V_B = V_B_air + V_B_vapour;
        D_eq(i) = (6*V_B / pi)^(1/3); % Equivalent Diameter
   
        % Non-dimensional Numbers
        [Re_x,Eo] = nondimensional_numbers(Prop,dz_dt(i),D_eq(i));
        [Re_y,~] = nondimensional_numbers(Prop,dx_dt(i),D_eq(i));

        % Calculate the cross-sectional and total surface area
        [E,A_B,S_B,r_a(i)] = bubble_area(V_B,Eo);
    
        % Temperature of Bubble
        T_air(i+1) = heat_transfer(Des,Prop,IC,dt,dz_dt(i),D_eq(i), T_air(i),m_vapour(i),S_B);
    
        % Water Uptake
        [m_dot,~] = water_uptake(S_B,V_B_air,m_vapour(i),T_air(i)-273,P);
        m_vapour(i+1) = m_vapour(i) + m_dot*dt;
        
        % Coefficient of Drag
        Cd_x(i) = TOMIYAMA_DragCoeff(Re_x,Eo,1);
        Cd_y(i) = TOMIYAMA_DragCoeff(Re_y,Eo,1);
    
        % Eccentricity
        [delta_A(i),delta_H(i)] = eccentricity_correction(E);
    
        % Force Balance
        [F.dragx(i), F.buoyancyx(i), F.bassetx(i)] = force_balance(i,t,dt,dz_dt(i),d2z_dt2,D_eq(i),A_B,V_B,Cd_x(i),Prop,1);
        F_rx = F.buoyancyx(i) - F.dragx(i) - delta_H(i)*F.bassetx(i); % resultant force

        [F.dragy(i), F.buoyancyy(i), F.bassety(i)] = force_balance(i,t,dt,dx_dt(i),d2x_dt2,D_eq(i),A_B,V_B,Cd_y(i),Prop,0);
        F_ry = - F.dragy(i) - delta_H(i)*F.bassety(i); % resultant force

        % Newton's 2nd Law
        m = (IC.m_air + m_vapour(i)) + delta_A(i)*Prop.rho_l/2*V_B; % mass of bubble plus added mass effect
        d2z_dt2(i) = F_rx / m; d2x_dt2(i) = F_ry / m;
    
        % Numerical Integration
        [z(i+1),dz_dt(i+1)] = EulerDoubleIntegration(z(i),dz_dt(i),d2z_dt2(i),dt);
        [x(i+1),dx_dt(i+1)] = EulerDoubleIntegration(x(i),dx_dt(i),d2x_dt2(i),dt);
    
        % Increment Count
        i = i + 1;
    
    end

    % Outputs
    AH = m_vapour/V_B_air;

end

%% Auxiliary functions
function [x_2,dx_dt_2] = EulerDoubleIntegration(x,dx_dt,d2x_dt2,dt)
    dx_dt_2 = dx_dt + d2x_dt2*dt;
    x_2 = x + dx_dt*dt;
end

function V = idealgas(m,R,T,P)
    V = m*R*T/P;
end
