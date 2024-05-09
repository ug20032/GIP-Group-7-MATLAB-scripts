function [Re,Eo] = nondimensional_numbers(Prop,u,L)
    
    %Calculates Reynold's and Eotvos numbers
    Re = Prop.rho_l*u*L/ Prop.mu_l;
    Eo = Prop.delta_rho*Prop.g*L^2/Prop.gamma;

end