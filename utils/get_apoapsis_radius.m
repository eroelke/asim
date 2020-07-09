% 
function [r_a, r_p, e] = get_apoapsis_radius(y, mu)

ae = rv2ae(y, mu); %orbital elements
a = ae(1);
e = ae(2);
r_a = a*(1+e);  %apoapsis
r_p = a*(1-e);  %periapsis

end