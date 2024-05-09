
function Des = conv2SI(D_mm,H_mm,Q_L,T_water_c,theta)

    % Converts all design variables into SI units

    Des.R_pipe = D_mm/2/1e3; % Diameter -> Radius
    Des.Q = Q_L * 1000/(60*1e6); % L/min -> m3/s
    % Note: T_air is an initial condition
    Des.T_water = T_water_c + 273;
    Des.H_m = H_mm/1e3; % mm -> m
    Des.T_air_c = 20;

    Des.theta_rad = deg2rad(theta);

end