function optimal_trajectory(S)

global GM timestep DEP_oe dist_moon_min tsol1 xsol1
global inclination_clos_app DV1 DV2 ind_dist_min

% Trajectory integration
tspan1 = S(5):timestep:S(5)+7*24*3600;
% Transformation of orbital elements to cartesian state vector in ECI
DEP_cart = cspice_conics([DEP_oe(1),S(1:4),DEP_oe(6),S(5),GM(4)]',S(5));
state_earth0 = cspice_spkezr('EARTH', S(5), 'J2000', 'NONE', 'SSB');
SI0 = DEP_cart(1:6) + state_earth0;
% Celestial bodies' state vectors
state_moon = cspice_spkezr('MOON', tspan1, 'J2000', 'NONE', 'SSB')';
state_earth = cspice_spkezr('EARTH', tspan1, 'J2000', 'NONE', 'SSB')';
xmoon = state_moon(:,1:3);
xearth = state_earth(:,1:3);

% Integration
force_fun = force_model(GM);
tsol1 = [];
xsol1 = [];
odeoptions = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tsol1, xsol1] = ode45(force_fun,tspan1,SI0,odeoptions);

% Minimum distance from Moon
dist_moon = zeros(1,length(tsol1));
for j=1:length(tsol1)
    dist_moon(j) = norm(xsol1(j,1:3) - xmoon(j,1:3));
end
[dist_moon_min,ind_dist_min] = min(dist_moon);

% Inclination at closest approach
J2000_to_IAUmoon     = cspice_sxform('J2000','IAU_MOON',tsol1(ind_dist_min));
x_closest_app        = xsol1(ind_dist_min,1:6) - state_moon(ind_dist_min,1:6);
x_clos_app_IAUmoon   = J2000_to_IAUmoon * x_closest_app';
oe_closest_app       = cspice_oscelt(x_clos_app_IAUmoon,tsol1(ind_dist_min),GM(3));
inclination_clos_app = oe_closest_app(3);

%%% COST FUNCTION %%%
DV1 = abs(norm(xsol1(1,4:6) - state_earth(1,4:6)) - sqrt(GM(4)/DEP_oe(1)));
DV2 = norm(xsol1(ind_dist_min,4:6) - ...
    state_moon(ind_dist_min,4:6)) - sqrt(GM(3)/dist_moon_min);

%--------------------------------------------------------------------------
%%% PLOT %%%
% Coordinates wrt ECI %
XSC  = xsol1(:,1:3) - xearth(:,1:3);
XM   = xmoon(:,1:3) - xearth(:,1:3);
% Spacecraft trajectory
plot3(XSC(1:ind_dist_min,1),XSC(1:ind_dist_min,2),XSC(1:ind_dist_min,3),'b')
hold on
% Draw Earth
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;
surf(x,y,z, 'FaceColor', '#4DBEEE')
% Moon orbit
plot3(XM(1:ind_dist_min,1),XM(1:ind_dist_min,2),XM(1:ind_dist_min,3),'r')
% Draw Moon at orbit insertion
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + XM(ind_dist_min,1);
y2 = y2 * 1737 + XM(ind_dist_min,2);
z2 = z2 * 1737 + XM(ind_dist_min,3);
surf(x2,y2,z2, 'FaceColor', '#808080')
% Options
grid on
axis equal
hold off

