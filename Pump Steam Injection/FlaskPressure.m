clc
clear
close all

%% User selection
Param.klist.valve = 10; % valve (either open = 10 or closed = inf)


%% Parameters - design/calculation variables

% ODE time variables
Param.dt = 1e-3; % timestep
Param.tmax = 20; % maximum time
tlist = 0:Param.dt:Param.tmax;

% Mass flow Rate
Param.m_in_av = 1.5 /1000/60; % mass flow rate (kg/s)

% Note: if t_interval = t_drip then constant flow rate
Param.t_interval = 3; % time interval between drips (s)
Param.t_drip = 0.5; % time for each drip to completely boil in flask (s)

Param.m_in = periodicdrips(Param.m_in_av,Param.t_interval,Param.t_drip,tlist);

% Minor Loss coefficients
Param.klist.entrance = 0.5; % flask to 1/4" tube
Param.klist.exit = 1; % 1/4" tube to tee-connector

%% Properties - NOT design variables
% Water vapour
Prop.rho_v = 0.59; % density of water vapour (kg/m3) https://www.engineeringtoolbox.com/saturated-steam-properties-d_457.html
Prop.mu_v = 0.013e-3; % dynamic viscosity of water vapour (Pa s) https://www.engineeringtoolbox.com/steam-viscosity-d_770.html
Prop.R = 461; % specific gas constant of water vapour (J/kg/K) https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html

% Pipe
Prop.D = 6.35e-3; % inner diameter of pipe (m)
Prop.A = pi*Prop.D^2/4; % cross-sectional area of pipe (m2)
Prop.L = 0.25; % length of pipe (m)

% Flask
Prop.T = 100 + 273; % temperature (K)
Prop.V = 1e-3; % volume (m3) - 1L conical flask

%% Calculations
P = zeros(1,length(tlist)); v = P; m = P; % pre-allocate arrays
m(1) = 1e-99; v(1) = 1e-99; % ensure non-zero ICs otherwise get inf errors

% Solve ODE
for i = 1:length(tlist)

    P(i) = idealgas(m(i),Prop);
    k(i) = pipelosses(v(i),Prop,Param);
    v(i+1) = bernoulliflow(P(i),Prop,k(i));
    m(i+1) = massflow(m,v,Prop,Param,i);

end
m = m(1:end-1); v = v(1:end-1); % truncate output to same length as tlist

%% Graphs
fig = figure;
fig.Position = [50 200 700 500];
plot(tlist,Param.m_in*60*1000,'k')
xlabel('Time (s)')
ylabel('Mass flow rate (g/min)')
ylim([0 ceil(max(Param.m_in*60*1000)*2)])
yline(Param.m_in_av*60*1000,'k--')
lgd = legend('Transient','Average');
lgd.AutoUpdate = 'off';

yyaxis right
plot(tlist,P/1e3,'r')
ylim([0 ceil(max(P/1e3)*1.1)])
ylabel('Gauge Pressure (kPa)')
set(gca,'YColor','r')

fig = figure;
fig.Position = [820 200 700 500];
plot(tlist,v)
xlabel('Time (s)')
ylabel('Velocity through pipe (m/s)')

%% Functions
function m_in = periodicdrips(m_in_av,t_interval,t_drip,tlist)

    %%% Generate the inlet mass flow rate curve based on a periodic square wave. %%%

    % create square wave (amplitude 0-1) with correct frequency
    s = 0.5*(square(tlist*(2*pi)/t_interval,t_drip/t_interval*100) + 1); 

    % alter amplitude so average mass flow rate is correct
    s_av = mean(s);
    m_in = s * m_in_av / s_av;
end

function P = idealgas(m,Prop)

    %%% Use the ideal gas law to calculate the pressure in the flask at any given time. %%%
    P = m*Prop.R*Prop.T/Prop.V;
end

function k = pipelosses(v,Prop,Param)

    %%% Calculate the major and minor (pressure) losses in the pipe from the
    %%% flask to nitrogen pipe. %%%

    % https://learn-eu-central-1-prod-fleet01-xythos.content.blackboardcdn.com/60e83182c0bd4/8706771?X-Blackboard-S3-Bucket=learn-eu-central-1-prod-fleet01-xythos&X-Blackboard-Expiration=1712156400000&X-Blackboard-Signature=zsvYwS336DtzWpXV2pW%2FFQS0nOc%2BDvmM4Rkeufhnf3g%3D&X-Blackboard-Client-Id=113292&X-Blackboard-S3-Region=eu-central-1&response-cache-control=private%2C%20max-age%3D21600&response-content-disposition=inline%3B%20filename%2A%3DUTF-8%27%27Turbulent_flows_notes.pdf&response-content-type=application%2Fpdf&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEoaDGV1LWNlbnRyYWwtMSJGMEQCIHSM6cINveky6pO827dSPNVlPUGXLZhd53dNZEkVdpWDAiB1qkD7QfTWS%2Bl1RW0jM58VSFfTvKX5bfWX04Z3%2B9Bj4Sq9BQhzEAMaDDYzNTU2NzkyNDE4MyIMsVt5AcI%2BmuH%2FPTcEKpoFKvMs%2BoE7W6rZmUy6mlTWtZbvHwqKwe926Md9clbsCbka3noN07chWUbW%2BfEZfmVU%2FmOLHMnXfaiGvF7fBkoJMD9nCD9g78tx%2FfAP5e10bfkGrSuggj%2BK8Zgc%2FXrR6NMNUIsXLDhpqtK7vSglE1mVeIhnnhlRz76HE9Ks7GytVHUoUgfKQqN1izRd4A1ej014%2BfdDrKG52sOBIJGHATJ4H8Hgynr%2Fx0eC2%2BsImgd76Q%2B9r3qKfGjbZzq0FFPUGk8nF4Hb0ggqG5uwAh%2B0IQVaugVmGmkpzwejyXJ5YVFHa%2FRgM5tScFQj3WSaO6HErG9ctf9r9Ty6WsokSJy2JLrDJqtLWbNbReT4MAAjGI4PrVOxf0uQuj9TqMTFu%2BVAuL8sd7neoEN18r3ob8vwCIqZEBP%2BgGukIyapc38y%2B68w1EYiMRvFVtLUniaSbNCIY3kISq4ihZYabXxFztlPcn30yP1kmKa%2Bx8ouYLBciicqI3pVC7%2BGsifNyJARIC0juq9klzcLN%2BNG4wfMMukDTRfEysKJOJx8V%2BarV2mKXhiKjQ%2BhtPwICgNfRN3Q7oeoet3hZd%2Bs%2B5LQqtDYcKVaInVcP9vhp3J0haQojhkJ94K5mTDQ9r7bp3Gohq5HTinulwDxg2MvLFShV7A1Uy0VE8JwNGyCpVj0YuaiJjyttQv4fEtLKcUefE19SdF3BJ0v4EpDQIosfA9iHxco95%2BBbe6kfEFfOCTaCw%2Fcb8goVnqsEzW1KxqWosGY7UaMErh8Z6fsFqu9a8NzgAkQ7%2BJJ4R3AmGUpc9n0psQBXeMD4VHgFlHqZNURCVZU5DtZzajyTopPKJINWbv1zgbsT%2BLTNVH7uhOUWUBP1bkenDExQV20Vt2b2PLo0oMzbR5HMLXatLAGOrIBqMu1ok2Cv%2BlUvEGTUP11fIDsR8NoSGFeD7VqAbROAEayfoRkO%2Fsb0qyojPA8zpU%2BsDVFihZ%2Fd3hNg%2Bfhsrx9Zr2A2CwL8aXEKZgbD300DWAu609%2BOJX9%2FiyInDhvKv6TmridfjD44xZe1q7R9psWnSGolAGW%2BtfmIXxXcQfV4FdbDcfShFoiDmZJtOH16UsbrJ0bv9rlM4M141H%2FcfGYz%2F5XwQAiKdIjiEVq6FG0%2FiWWFA%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240403T090000Z&X-Amz-SignedHeaders=host&X-Amz-Expires=21600&X-Amz-Credential=ASIAZH6WM4PLUSUGC4MY%2F20240403%2Feu-central-1%2Fs3%2Faws4_request&X-Amz-Signature=5f0189a83eb8aef486b48dc3f5898da31e678c69663291617ffa962dda374526

    % Major Losses
    Re = Prop.rho_v*abs(v)*Prop.D/Prop.mu_v;
    if Re < 1000 % laminar flow
        f = 64 / Re; 
    else % turbulent flow
        f = 0.316 * Re^(-1/4); % Blasius equation
    end
    k_major = f*Prop.L/Prop.D; % Darcy-Weisbach Equation: 1/2*rho*v^2*f*L/D

    % Minor Losses
    k_minor = Param.klist.entrance + Param.klist.valve + Param.klist.exit;  
    k = k_minor + k_major;
end

function v = bernoulliflow(P,Prop,k)

    %%% Calculate the velocity of water vapour in the pipe based on
    %%% Bernoulli's principle, assuming the velocity in the flask is zero
    %%% and the hydrostatic pressure is negligible. %%% 

    % P1 + 1/2*rho*v1^2 + rho*g*z1 = P2 + 1/2*rho*v2^2 + rho*g*z2 + delta_p
    % v1 = z1 = z2 = 0
    % P1 = 1/2*rho*v2^2 + delta_p = 1/2*rho*v2^2*(1 + k)

    v = sqrt((2*P)/(Prop.rho_v*(1+k)));

end

function m2 = massflow(m,v,Prop,Param,i)

    %%% Calculate the mass of water vapour in the flask based on conservation of mass. %%% 
    v2 = max(0,v(i+1));
    
    % dm/dt = (dm/dt)_in - (dm/dt)_out
    mdot_out = Prop.rho_v*Prop.A*v2; % Continuity 
    dm_dt = Param.m_in(i) - mdot_out;

    m2 = m(i) + dm_dt*Param.dt; % Euler forward integration
    m2 = max(m2,1e-99); % ensure mass of water > 0 
end

