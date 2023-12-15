function plot_initial_guess(S)

global GM timestep DEP_oe

% Trajectory integration
tspan1 = S(5):timestep:S(5)+10*24*3600;
% Transformation of orbital elements to cartesian state vector in ECI
DEP_cart = cspice_conics([DEP_oe(1),S(1:4),DEP_oe(6),S(5),GM(4)]',S(5));
% Departure state vector
SI0([1:3,13:15])      = DEP_cart(1:6);
[SI0([4:6,16:18]), ~] = cspice_spkezr('sun',S(5),'j2000','NONE','earth');
[SI0([7:9,19:21]), ~] = cspice_spkezr('moon',S(5),'j2000','NONE','earth');
SI0([10:12,22:24])    = [0,0,0,0,0,0];
% Integration
fbpfun = fourbp(GM);
tsol1 = [];
xsol1 = [];
[tsol1, xsol1] = ode89(fbpfun,tspan1,SI0);

% Minimum distance from Moon
dist_moon = zeros(1,length(tsol1));
for j=1:length(tsol1)
    dist_moon(j) = norm(xsol1(j,1:3) - xsol1(j,7:9));
end
[dist_min_moon,ind_dist_min] = min(dist_moon);
% Display minimum distance from Moon
fprintf('\nMinimum distance from Moon: %.2f km\n', dist_min_moon)


%--------------------------------------------------------------------------
% Coordinates wrt ECI %
xSC  = xsol1(1:end,1:3) - xsol1(1:end,10:12);
xM   = xsol1(1:ind_dist_min,7:9) - xsol1(1:ind_dist_min,10:12);
% Spacecraft trajectory
plot3(xSC(:,1),xSC(:,2),xSC(:,3),'b')
hold on
plot3(xSC(ind_dist_min,1),xSC(ind_dist_min,2),xSC(ind_dist_min,3),'ob')
% Draw Earth
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;
surf(x,y,z, 'FaceColor', '#4DBEEE')
% Moon orbit
plot3(xM(1:120/timestep:end,1),xM(1:120/timestep:end,2),xM(1:120/timestep:end,3),'r')
% Draw Moon at closest approach
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + xM(end,1);
y2 = y2 * 1737 + xM(end,2);
z2 = z2 * 1737 + xM(end,3);
surf(x2,y2,z2, 'FaceColor', '#808080')
% Options
grid on
axis equal
hold off


