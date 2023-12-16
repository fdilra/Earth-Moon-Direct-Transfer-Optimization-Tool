close all
clear all

tic

%%% GLOBAL VARIABLES AND KERNEL INITIALIZATION %%%
global GM  DEP_oe timestep rad_arr dist_moon_min target_switch
global tsol1 xsol1 DV1 DV2 ind_dist_min inclination_clos_app inc_arr 

% Load kernels
cspice_furnsh( 'metakr.tm' );


%--------------------------------------------------------------------------
%%% INITIAL GUESS & USER-DEFINED PARAMETERS %%%
% Launch date in UTC
timstr = '2024 31 AUG 12:00:00';  
% Spacecraft mass
mass_sc = 1000; %[kg]

% Departure orbital elements (with respect to J2000/ICRF. Units: [km, rad])
earth_radius = 6371; % [km]
DEP_oe(1) = earth_radius + 250; % radius of parking orbit/departure orbit pericenter
DEP_oe(2) = 0.965; % eccentricity of departure orbit
DEP_oe(3) = 0.3; % inclination of parking orbit/departure orbit
DEP_oe(4) = 1.5; % longitude of ascending node of parking orbit/departure orbit
DEP_oe(5) = 2.6; % argument of pericenter of parking orbit/departure orbit
% Mean anomaly of departure orbit - DO NOT CHANGE
DEP_oe(6) = 0;          

% Lunar orbit radius and inclination
moon_radius = 1737;       % [km]
rad_arr = moon_radius + 100; % [km]
inc_arr = pi/2;           % [rad]

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
S0(1) = DEP_oe(2);
S0(2) = DEP_oe(3);
S0(3) = DEP_oe(4);
S0(4) = DEP_oe(5);
S0(5) = t0;
% Lower and upper bounds
lb = [0.95,  0,    0,    0,    t0 - 2*24*3600];
ub = [0.999, pi/2, 2*pi, 2*pi, t0 + 2*24*3600];
% Options
options = optimoptions("fmincon",...
    "Algorithm","active-set",...
    'Display','iter',...
    "FunctionTolerance",1e-06,...
    "ConstraintTolerance",1e-16,...
    "MaxFunctionEvaluations",1000);


% Plot initial guess and manually check quality. 
% Execution can be stopped here if initial guess is not satisfying. 
% Argument of pericenter of spacecraft should be opposite to the Moon's 
% position at departure for best results. Changing DEP_oe(5) by pi/2
% increments is the easiest way to find a decent initial guess.
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


% Target distance
fprintf('\nTargeting...')
target_switch = 1;
[S0_opt, ~] = ...
    fmincon(@cost,S0,[],[],[],[],lb,ub,@constraints,options);
fprintf('Lunar orbit radius: %.2f km\n', dist_moon_min - 1737)
fprintf('Lunar orbit inclination: %.2f rad\n',inclination_clos_app)
% Optimization
fprintf('Optimizing Delta V...')
target_switch = 2;
[S0_opt, DV_opt] = ...
    fmincon(@cost,S0_opt,[],[],[],[],lb,ub,@constraints,options);


%--------------------------------------------------------------------------
%%% PLOT %%%
% Integrate trajectory with optimal parameters
optimal_trajectory(S0_opt)
% Spacecraft and Moon positions in ECI
x_sc = xsol1(1:ind_dist_min,1:3) - xsol1(1:ind_dist_min,10:12);
x_m  = xsol1(1:ind_dist_min,7:9) - xsol1(1:ind_dist_min,10:12);
plot3(x_sc(:,1),x_sc(:,2),x_sc(:,3),'b')
hold on
plot3(x_m(:,1),x_m(:,2),x_m(:,3),'r')
% Draw Earth
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;
surf(x,y,z, 'FaceColor', '#4DBEEE')
% Draw Moon at lunar insertion
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + x_m(end,1);
y2 = y2 * 1737 + x_m(end,2);
z2 = z2 * 1737 + x_m(end,3);
surf(x2,y2,z2, 'FaceColor', '#808080')
% Options
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
axis equal
hold off

% Print results
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
    norm(x_sc(end,1:3)-x_m(end,1:3)) - 1737)
fprintf('Lunar orbit inclination: %.2f deg\n',rad2deg(inclination_clos_app))

% Unload kernels
cspice_kclear
