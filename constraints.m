function [c,ceq]=constraints(S)

global rad_arr dist_moon_min inclination_clos_app inc_arr target_switch
global rad_tol inc_tol

if target_switch==1
    c = [];
elseif target_switch==2
    c(1) = abs(dist_moon_min - rad_arr) - rad_tol;
    c(2) = abs(inclination_clos_app - inc_arr) - inc_tol;
    c = c';
end



ceq = [];


