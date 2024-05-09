clc
clear
close all

%% User selection
t_hold = 60; % either 17 or 60
chem_sim_version = 2; % either without bypass (1) or with bypass factor (2)

%% Bubble Humidifier Simulation
main_TOTAL_chemical
formatting

t_flask = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 21+1e-3 21+t_hold];
Q_flask = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.5 0.5]/n_holes;
if t_hold == 17
    tmax = 60;
    T_flask = [22.8 23.1 24.8 26.6 30.2 34.0 38.1 43.5 49.0 54.5 59.1 63.8 68.1 71.9 75.0 77.7 80.1 82.2 84.0 85.0 85.0 85.0 85.0 85];
    actual_burnoff = 0.64;
elseif t_hold == 60
    tmax = 120;
    T_flask = [22.1 23.0 25.1 26.7 28.2 32.2 36.9 41.8 46.6 51.3 55.9 60.3 64.6 68.7 72.5 75.7 78.4 80.9 82.7 84.1 84.8 85.0 85.0 85];
    actual_burnoff = 0.78;
end

% Carbon boat dimensions
w = 2e-2; % m (width of carbon boat)
L = 10e-2; % m (length of carbon boat)
D = 7.5e-2; % m
A = pi*D^2/4; % m2
V = A*L;
if chem_sim_version == 1
    bypass_ratio = 1;
    char_burnoff = 0;
elseif chem_sim_version == 2
    %bypass_ratio = w / (pi*D);
    bypass_ratio = 0.06;
    char_burnoff = 0.5985;
end

% get the water vapour and nitogen flow rate
[f_H2O,f_H2O_gmin,f_N2] = flowrates(t_hold,t_flask,bypass_ratio,Q_list,T_list,m_dot_water_tot,Q_flask,T_flask);

%% Furnace Environment 
% initial concentration (mol)
mC_i = 2; % g
nC_initial = mC_i / 12; % moles of C before carbonisation (g -> mol)
nC_i = nC_initial * (1 - char_burnoff); % moles of C after carbonisation
nCO2_i = 0; nCO_i = 0; nH2O_i = 1e-99; nH2_i = 0;

t = [linspace(0,1,1e4) linspace(1,10,2e3) linspace(10,140*60,2e3)]; % increase number of datapoints at start
t = unique(t);

T = 900 + 273; % Temperature of furnace (K)
P_tot = 1; % Temperature of furnace (atm)

%% Solve Ensemble of ODEs
ens_var = 5;
ensemble = create_ensemble(ens_var);
disp('---------- Simulating chemical reactions ----------')
for i = 1:length(ensemble)
    [nC(i,:),nH2O(i,:),nCO2(i,:),nCO(i,:),nH2(i,:),r1(i,:),r2(i,:),r3(i,:)] = solveRateEqs(nC_i,nH2O_i,nCO2_i,nCO_i,nH2_i,f_H2O,f_N2,t,V,P_tot,T,t_hold,ensemble(i,:));

    disp([num2str(i) '/' num2str(length(ensemble))])
end
disp('---------- Chemical Reaction Simulation Complete ----------')
%% Plots
% WVFR plot ---------------------------------------------------------------
figure
plot(t_flask,T_flask,'k','LineWidth',1)
xlabel('Time (mins)'); ylabel('Temp ($^\circ$C) or N$_2$ Flow Rate (cL/min)');

xlim([0 tmax]); ylim([10 90]); yticks([10:20:90])
hold on
plot(t_flask,Q_flask*n_holes*100,'Color',colour(1,:),'LineWidth',1)
lgd = legend('Temperature','Flow Rate','location','best');
lgd.AutoUpdate = 'off';
yyaxis right; set(gca,'YColor',colour(2,:));
plot(f_H2O_gmin(:,1),f_H2O_gmin(:,2),'Color',colour(2,:),'LineWidth',1)
ylabel('Water Vapour Flow Rate (g/min)');
ylim([-1e-3 0.2]);


% Concentration of reactants/products plot --------------------------------
figure
hold on
p_C = plot(t/60,nC(1,:),'Color',colour(1,:));
p_H2O = plot(t/60,nH2O(1,:),'Color',colour(2,:));
p_CO2 = plot(t/60,nCO2(1,:),'Color',colour(3,:));
p_CO = plot(t/60,nCO(1,:),'Color',colour(4,:));
p_H2 = plot(t/60,nH2(1,:),'Color',colour(5,:));
p_burn = yline(mC_i/12 * (1 - actual_burnoff),'k--');
xline(max(t_flask),'k:','LineWidth',2)
xlabel('Time (mins)'); ylabel('Concentration (mol)');

if chem_sim_version == 1
    xlim([0 tmax]); ylim([0 0.2]); yticks(0:0.04:0.2)
elseif chem_sim_version == 2
    xlim([0 tmax]); ylim([0 0.08]); yticks(0:0.02:0.2)
end

if chem_sim_version == 2 && t_hold == 17
    lgd = legend('C','H$_2$O','CO$_2$','CO','H$_2$','Location','SE','Interpreter','latex');
else
    lgd = legend('C','H$_2$O','CO$_2$','CO','H$_2$','Location','NE','Interpreter','latex');
end

lgd.AutoUpdate = 'off';
lgd.Visible = 'on';
% Ensemble plotting
nskip = 1e2;
fill([t(1:nskip:end)/60, fliplr(t(1:nskip:end)/60)], [min(nC(:,1:nskip:end)), fliplr(max(nC(:,1:nskip:end)))],colour(1,:),'FaceAlpha',0.1,'EdgeColor','none')
fill([t(1:nskip:end)/60, fliplr(t(1:nskip:end)/60)], [min(nH2O(:,1:nskip:end)), fliplr(max(nH2O(:,1:nskip:end)))],colour(2,:),'FaceAlpha',0.05,'EdgeColor','none')
fill([t(1:nskip:end)/60, fliplr(t(1:nskip:end)/60)], [min(nCO2(:,1:nskip:end)), fliplr(max(nCO2(:,1:nskip:end)))],colour(3,:),'FaceAlpha',0.05,'EdgeColor','none')
fill([t(1:nskip:end)/60, fliplr(t(1:nskip:end)/60)], [min(nCO(:,1:nskip:end)), fliplr(max(nCO(:,1:nskip:end)))],colour(4,:),'FaceAlpha',0.05,'EdgeColor','none')
fill([t(1:nskip:end)/60, fliplr(t(1:nskip:end)/60)], [min(nH2(:,1:nskip:end)), fliplr(max(nH2(:,1:nskip:end)))],colour(5,:),'FaceAlpha',0.05,'EdgeColor','none')
set(gca, 'XGrid', 'off');
set(gca, 'YGrid', 'on');

% Performance parameters --------------------------------------------------
burnoff = (1-nC(1,end)/nC_initial)*100;
burnoff_min = (1-max(nC(:,end))/nC_initial)*100;
burnoff_max = (1-min(nC(:,end))/nC_initial)*100;
disp(['Burn-off: ' num2str(burnoff) '%'])

%% Functions

function r = Boudouard(nCO2,nCO,T,P_tot)
    % (Jalan, 1977)
    % nX: moles/second
    % T: K
    % P_tot: Pa

    % r: mol/cm3 atm s

    R = 8.31;
    n_tot = nCO2 + nCO;
    P_tot = P_tot / 1e5;

    P_CO2 = nCO2/n_tot * P_tot;
    k_0 = exp(5.471 - 48020/(2.303*R*T));

    if n_tot == 0 % rate is zero if no moles
        r = 0;
    else 
        r = k_0*P_CO2;
    end
    r = max(r,0);
end

function r = carbon_steam(pH2O,pH2,T,W_r)

    % pX:atm
    % T: K
    % P_tot: Pa

    % Data from Johnstone, 1952
    % Custom curve fit
    k1_coeff = [141.09711 -5768.4370 1.139389 0.662035 0.0047877];
    k2_coeff = [1.3134177e-07 16976.162043 0.405285 0.155373 0.122469];
    k3_coeff = [3.0e-06 11365.14345 0.626765 0.0700145 -0.467086];

    T = T - 273; % K -> C

    k1 = k1_coeff(1).*exp(k1_coeff(2)./T)*(W_r+k1_coeff(3))^k1_coeff(4) + k1_coeff(5);
    k2 = k2_coeff(1).*exp(k2_coeff(2)./T)*(W_r+k2_coeff(3))^k2_coeff(4) + k2_coeff(5);
    k3 = k3_coeff(1).*exp(k3_coeff(2)./T)*(W_r+k3_coeff(3))^k3_coeff(4) + k3_coeff(5);

    r = k1*pH2O / (1 + k2*pH2 + k3*pH2O); % Langmuir type

    r = max(r,0);
end

function r = WGS(pCO,pH2O,T)

    % De Jong, 1980
    % pX: partial pressure (atm)
    % r_main: mol/gcat/s
    % r: mol/s

    gcat = 0.02e-5; % fudge factor -> would like to change

    R = 8.31;
    r_main = 25.9e3*exp(-1.6e4/(R*T)) * (pCO*pH2O)/(1 + 127*pCO*pH2O + 26*pCO);

    r = r_main*gcat; % mol/s
    r = max(r,0);
end


function [nC,nH2O,nCO2,nCO,nH2,r1,r2,r3] = solveRateEqs(nC,nH2O,nCO2,nCO,nH2,f_H2O,f_N2,t,V,P_tot,T,t_hold,Ensemble)

    for i = 1:length(t)-1
    
        dt = t(i+1)-t(i);
    
        % molar ratio (of gases)
        n_tot = nH2O(i) + nCO2(i) + nCO(i) + nH2(i);
        RH2O = nH2O(i) / n_tot;
        RCO2 = nCO2(i) / n_tot;
        RCO = nCO(i) / n_tot;
        RH2 = nH2(i) / n_tot;
    
        % partial pressure
        pH2O = RH2O*P_tot;
        pCO2 = RCO2*P_tot;
        pCO = RCO*P_tot;
        pH2 = RH2*P_tot;
    
        % Eq. 1: C + H2O -> CO + H2
        W_r = (nC(1)-nC(i))/nC(1); % fraction of C reacted
        r1(i) = -carbon_steam(pH2O,pH2,T,W_r) * nC(i) * Ensemble(1);
    
        % Eq. 2: CO + H2O -> CO2 + H2
        r2(i) = -WGS(pCO,pH2O,T) * Ensemble(2) ;
    
        % Eq. 3: C + CO2 -> 2CO
        r3(i) = -Boudouard(nCO2(i),nCO(i),T,P_tot) / (V*100^3) * Ensemble(3);

        % Flow Rate
        f1 = interp1(f_H2O(:,1),f_H2O(:,2),t(i));
        f2 = f_N2;

        if t(i) < (21+t_hold)*60
            nC(i+1) = nC(i) + r1(i)*dt + r3(i)*dt; % number of moles varies based on reaction (r) & flow (f)
            nH2O(i+1) = nH2O(i) + r1(i)*dt + r2(i)*dt + f1*dt - RH2O*f1*dt;
            nCO2(i+1) = nCO2(i) + r3(i)*dt - r2(i)*dt - RCO2*f1*dt;
            nCO(i+1) = nCO(i) + r2(i)*dt - r1(i)*dt - 2*r3(i)*dt - RCO*f1*dt;
            nH2(i+1) = nH2(i) - r1(i)*dt - r2(i)*dt - RH2*f1*dt;
        else
            nC(i+1) = nC(i) + r1(i)*dt + r3(i)*dt; % number of moles varies based on reaction (r) & flow (f)
            nH2O(i+1) = nH2O(i) + r1(i)*dt + r2(i)*dt - RH2O*f2*dt;
            nCO2(i+1) = nCO2(i) + r3(i)*dt - r2(i)*dt - RCO2*f2*dt;
            nCO(i+1) = nCO(i) + r2(i)*dt - r1(i)*dt - 2*r3(i)*dt - RCO*f2*dt;
            nH2(i+1) = nH2(i) - r1(i)*dt - r2(i)*dt - RH2*f2*dt;
        end

        nC(i+1) = max(0,nC(i+1));
        nH2O(i+1) = max(0,nH2O(i+1));
        nCO2(i+1) = max(0,nCO2(i+1));
        nCO(i+1) = max(0,nCO(i+1));
        nH2(i+1) = max(0,nH2(i+1));

    end
    r1 = [0 r1];
    r2 = [0 r2];
    r3 = [0 r3];

end

function ensemble = create_ensemble(ens_var)

    ensemble = [];
    for i = 1:1e3
        
        ens_list = [1/ens_var 1 ens_var];
        rand_indx = [ceil(3*rand()) ceil(3*rand()) ceil(3*rand())];
        ensemble = [ensemble;
                    ens_list(rand_indx(1)) ens_list(rand_indx(2)) ens_list(rand_indx(3))];
    end
    ensemble = unique(ensemble,'rows');
    ensemble = [1 1 1;
                ensemble];
end

function [f_H2O,f_H2O_gmin,f_N2] = flowrates(t_hold,t_flask,bypass_ratio,Q_list,T_list,m_dot_water_tot,Q_flask,T_flask)

    f_H2O_gmin = t_flask';
    f_H2O_gmin(:,2) = interp2(Q_list,T_list,m_dot_water_tot,Q_flask,T_flask);
    
    f_H2O_gmin = [f_H2O_gmin; 
            21+t_hold+1e-3 1e-6
            150 1e-6]; % append zero flow rate afterwards
    
    
    f_H2O(:,1) = f_H2O_gmin(:,1) * 60; % mins -> s
    f_H2O(:,2) = f_H2O_gmin(:,2) * 1/60/18; % flow rate - X g/min -> mol/s
    f_H2O(:,2) = f_H2O(:,2) * bypass_ratio;
    
    f_N2 = 0.5 / 60 / 14; % mol/s

end