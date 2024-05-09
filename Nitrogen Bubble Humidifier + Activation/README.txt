
---------------------------------------------------------------------------
SIMULATION

Contains the scripts for simulating the dynamics, heat transfer and mass transfer of the bubbles in the NBH.

main_SINGLE - simulation for a single tube diameter. Used in Table 6 and 7.

main_TOTAL - simulation of the mass flow rate of water vapour in bubbles for varying diameter tubes. 
(i.e. run main_SINGLE for different tube diameters). Used in Fig. 9a and 13a.

main_OPTIMISE - determines the optimum diameter for a varying number of tubes 
(i.e. run main_TOTAL for different flow rates per tube). Used in Fig. 9b and 13b.

main_RateOfReaction - simulates the chemical rate of reactions, based on the 
estimated water flow rate from main_TOTAL_chemical. Used in Fig. 17 & 18.

main_TOTAL_chemical - simulation of the mass flow rate of water vapour for varying liquid water temperatures and 
nitrogen flow rates. (i.e. run main_SINGLE for different flow rates and temperatures).


User Inputs

sim_version = 2; % either initial (1) or updated (2) bubble dynamics simulation
iter = "3b"; % either iteration 1a, 1b, 2a, 2b, 3a or 3b
T_water_c = 85; % Temperature of the liquid Water (deg C) -  bubble dynamics validation (25) or WVFR (85)
n_holes = 8; % either 8, 24 or 31 holes in the diffuser

t_hold = 60; % either 17 or 60 hold-time at 85 deg C during activation
chem_sim_version = 2; % either without bypass (1) or with bypass factor (2)

---------------------------------------------------------------------------
VALIDATION / IMAGE_PROCESSING

videoprocess_velocity - determines the displacement (and therefore velocity) of the centre
of the highest bubble rising into quiescent fluid. Used in Table 6 and Fig. 12b.

videoprocess_velocity - determines the maximum distance between any two points in a bubble. 
Used in Table 6 and Fig. 12a & 12 c.

---------------------------------------------------------------------------
VALIDATION / EXPERIMENTS_VS_PREDICTIONS

Validation - creates Fig. 12b & 12c

WVFR_test - creates WVFR figure in presentation
