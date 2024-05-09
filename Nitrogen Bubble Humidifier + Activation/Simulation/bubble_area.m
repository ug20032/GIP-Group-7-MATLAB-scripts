function [E,A_B,S_B,a] = bubble_area(V_B,Eo)

    % Calculates the cross-sectional (A_B), total surface area (S_B),
    % aspect ration (E) and major semi-axis (a) of an oblate spheroidal
    % bubble based on the volume (V_B) and Eotvos number.

    E = (1 + 0.163*Eo^0.757)^-1; % aspect ratio
    e = 1 - E^2; % eccentricity
    
    R_eq = (3*V_B / (4*pi))^(1/3);
    a = R_eq / E^(1/3); % major semi-axis
    c = E * a; % minor semi-axis
    
    A_B = pi*R_eq^2; % cross-sectional area of sphere (for drag coefficient)
    S_B = 2*pi*a^2 + pi*c^2/e*log((1+e)/(1-e)); % cross-sectional area (spheroid)
    
end
        