function [F_drag, F_buoyancy, F_basset] = force_balance(i,t,dt,dx_dt,d2x_dt2,D_eq,A_B,V_B,Cd,Prop,isbuoyant)

    % Calculates the drag, buoyancy and Basset history forces. Note: added
    % mass is included in the mass term of Newton II. x is the displacement
    % of the bubble. 
    % Isbuoyant = 1 if calculating forces in vertical direction, or = 0 if in the radial direction. 
    
    % Buoyancy Force
    F_buoyancy = Prop.delta_rho * Prop.g * V_B * isbuoyant;

    % Drag Force
    F_drag = 1/2*Prop.rho_l*dx_dt^2*A_B*Cd;
    
    % Basset History Force
    F_basset = 3/2*D_eq^2*sqrt(pi*Prop.rho_l*Prop.mu_l)*sum(1./sqrt(t(i) - t(1:i-1)) .* d2x_dt2(1:i-1) .* dt);

end