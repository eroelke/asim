% Function to determine required planet-relative velocity to get a specific
% inertial velocity magnitude.
% Inputs:
% Vel_ini = desired inertial velocity
% 
% Outputs:
% mag = required planet-relative velocity

function mag = find_vel_pp_mag(vel_ini, gamma, omega, r)

A = 1;
B = 2*norm(omega)*norm(r)*cos(gamma);
C = (norm(omega)*norm(r))^2 - vel_ini^2;

p = [A B C];
root = roots(p);

mag = max(root);
end