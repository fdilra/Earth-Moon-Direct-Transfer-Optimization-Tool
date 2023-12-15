function [c,ceq]=constraints(S)

global rad_arr dist_moon_min inclination_clos_app inc_arr target_switch

if target_switch==1
    c = [];
elseif target_switch==2
    c(1) = abs(dist_moon_min - rad_arr) - 1;
elseif target_switch==3
    c(1) = abs(dist_moon_min - rad_arr) - 1;
    c(2) = abs(inclination_clos_app - inc_arr) - 0.0015;
    c = c';
end



ceq = [];


