function [delta_A,delta_H] = eccentricity_correction(E)

    % Calculates the added mass (delta_A) and Basset history (delta_H)
    % eccentricity correction from the aspect ratio (E)

    delta_A = 2*(E*acos(E)-sqrt(1-E^2)) / (E^2*sqrt(1-E^2) - E*acos(E)); % added mass
    delta_H = ((1.2*(4+E))/6)^2; % history term

end