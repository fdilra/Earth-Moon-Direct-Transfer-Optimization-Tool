function F = force_model(GM)

F = @(t, u) compute_force(t, u, GM);

function force = compute_force(t, u, GM)
state = cspice_spkezr('SUN', t, 'J2000', 'NONE', 'SSB');
sun_position = state(1:3);

state = cspice_spkezr('MOON', t, 'J2000', 'NONE', 'SSB');
moon_position = state(1:3);

state = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SSB');
earth_position = state(1:3);

norm_sun = norm(sun_position - u(1:3))^3;
norm_moon = norm(moon_position - u(1:3))^3;
norm_earth = norm(earth_position - u(1:3))^3;

force = [
    u(4:6);
    %
    GM(2) * (sun_position - u(1:3)) / norm_sun + ...
    GM(3) * (moon_position - u(1:3)) / norm_moon + ...
    GM(4) * (earth_position - u(1:3)) / norm_earth;
    ];
return
