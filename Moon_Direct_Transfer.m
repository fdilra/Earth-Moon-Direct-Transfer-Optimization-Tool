close all
clear all

tic

%%% GLOBAL VARIABLES AND KERNEL INITIALIZATION %%%
global GM  DEP_oe timestep rad_arr dist_moon_min target_switch
global tsol1 xsol1 DV1 DV2 ind_dist_min inclination_clos_app inc_arr 
global rad_tol inc_tol

% Load kernels
cspice_furnsh( 'metakr.tm' );


%--------------------------------------------------------------------------
%%% INITIAL GUESS & USER-DEFINED PARAMETERS %%%
% Launch date in UTC
timstr = '2024 21 AUG 12:00:00';  
% Spacecraft mass
mass_sc = 1000; %[kg]

% Departure orbital elements (with respect to J2000/ICRF. Units: [km, rad])
earth_radius = 6371; % [km]
DEP_oe(1) = earth_radius + 250; % radius of parking orbit/departure orbit pericenter
DEP_oe(2) = 0.965; % eccentricity of departure orbit
DEP_oe(3) = 0.3; % inclination of parking orbit/departure orbit
DEP_oe(4) = 1.5; % longitude of ascending node of parking orbit/departure orbit
DEP_oe(5) = 1.6; % argument of pericenter of parking orbit/departure orbit
% Mean anomaly of departure orbit - DO NOT CHANGE
DEP_oe(6) = 0;          

% Lunar orbit radius and inclination
moon_radius = 1737;          % [km]
rad_arr = moon_radius + 100; % [km]
inc_arr = pi/2;              % [rad]
% Tolerances on distance and inclination. Smaller tolerances can make
% convergence more difficult.
rad_tol = 0.5;   %[km]
inc_tol = 0.001; %[rad]

% Integration step in seconds
timestep = 60; 


%--------------------------------------------------------------------------
%%% PARAMETERS %%%
% Launch date in seconds past J2000
t0 = cspice_str2et(timstr);  

%%% Gravitational parameters %%%
GM    = zeros(1,4);
% Spacecraft parameter
GM(1) = mass_sc * 6.67139e-20; 
% Celestial bodies parameters
GM(2) = cspice_bodvrd('SUN', 'GM', 1);
GM(3) = cspice_bodvrd('MOON', 'GM', 1);
GM(4) = cspice_bodvrd('EARTH', 'GM', 1);


%--------------------------------------------------------------------------
%%% TARGETING AND OPTIMIZATION %%%
% Initial guess array 
S0(1) = DEP_oe(2); % eccentricity
S0(2) = DEP_oe(3); % inclination
S0(3) = DEP_oe(4); % longitute of ascending node
S0(4) = DEP_oe(5); % argument of pericenter
S0(5) = t0;        % departure time
% Lower and upper bounds
lb = [0.95,  0,    0,    0,    t0 - 2*24*3600];
ub = [0.999, pi/2, 2*pi, 2*pi, t0 + 2*24*3600];


% Plot initial guess and check quality. Script can usually find a solution
% even with a bad guess, but a good guess obviously converges faster.
plot_initial_guess(S0)
% Continue or stop execution prompt
kbhit = input('\nType y to continue, anything else to stop: ','s');
if kbhit=='y'
    kbhit = [];
    close all
else
    close all
    cspice_kclear
    return
end 

%%% Targeting - distance and inclination %%%
% Active-set is used as an algorithm for targeting as it can take large
% steps and there are no constraints in this phase.
options = optimoptions("fmincon",...
    "Algorithm","active-set",...
    'Display','iter',...
    "FunctionTolerance",1e-04,...
    "ConstraintTolerance",1e-16,...
    "MaxFunctionEvaluations",1000);
fprintf('\nTargeting...')
target_switch = 1;
[S0_opt, ~] = ...
    fmincon(@cost,S0,[],[],[],[],lb,ub,@constraints,options);
fprintf('Lunar orbit radius: %.2f km\n', dist_moon_min - 1737)
fprintf('Lunar orbit inclination: %.2f rad\n',inclination_clos_app)

%%% Optimization - Delta V %%%
% Interior-point is used here as it's usually faster since it satisfies
% constraints at every iteration and steps can be small in this phase.
% Active-set can be also used; the solutions it finds are better but the
% difference is usually negligible.
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    'Display','iter',...
    "FunctionTolerance",1e-04,...
    "ConstraintTolerance",1e-16,...
    "MaxFunctionEvaluations",1000);
fprintf('Optimizing Delta V...')
target_switch = 2;
[S0_opt, DV_opt] = ...
    fmincon(@cost,S0_opt,[],[],[],[],lb,ub,@constraints,options);


%--------------------------------------------------------------------------
%%% PLOT %%%
% Plot trajectory with optimal parameters
optimal_trajectory(S0_opt)

% Display results
elapsedTime = toc/60;
ToF = (tsol1(ind_dist_min) - S0_opt(5))/(24*3600); 
DepDate = cspice_et2utc(S0_opt(5),'C',0);
fprintf('\nDeparture date: %s\n', DepDate)
fprintf('Total Delta V: %.5f km/s\n', DV_opt)
fprintf('Trans Lunar Injection - Delta V1: %.5f km/s\n', DV1)
fprintf('Lunar Orbit Insertion - Delta V2: %.5f km/s\n', DV2)
fprintf('Time of Flight: %.2f days\n', ToF)
fprintf('Elapsed time: %.2f minutes\n', elapsedTime)
fprintf('Lunar orbit height: %.2f km\n', ...
    dist_moon_min - 1737)
fprintf('Lunar orbit inclination: %.2f deg\n',rad2deg(inclination_clos_app))

% Unload kernels
cspice_kclear
