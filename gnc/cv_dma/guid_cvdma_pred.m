% guid_cvdma_pred.m 
%   -predictor for continuous drag modulation for aerocapture
%   -propagates state to atmospheric exit for a constant reference area
%
% Aeroassist Simulation source code
% Developed by:
%   Space Systems Design Lab, Georgia Tech
%   Colorado Center for Astrodynamics Research, CU Boulder
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   s - NPC state structure
%   p - NPC input structure
%   Aref - current reference area (m2)
%
% Outputs:
%   s - updated state struct with r_ap and delta_ap
%
% Major Revision History:
%   *Created August 2018, E. Roelke
%   *Fixed oe calcs Sep 2018 - E. Roelke

function [ra, delta_ra,pflag] = guid_cvdma_pred(y0, t0, s, p, area_ref)
%#codegen

% Initialize
y = y0; % multi, initial state     y = [r_inertial;v_inertial]
t = t0; % s, initial time
pflag = 0; % nd, termination flag
h = 1/p.traj_rate;    %s, integration time step

K_dens = s.K_dens;  % density corrector
mass = p.mass;
cd = s.cd;
cl = 0;

% for debugging
% i = 1;
% x(:,i) = y0;
% ttemp(i) = t0;
%% Integration (RK4)
while ~pflag
    % use same integration as in simulation (finding error sources): 10/18
    [k1,alt] = npc_eom(t, y, area_ref, K_dens, p, s, cd, cl, mass);
    [k2,~] = npc_eom(t + h/2, y + (h*k1/2), area_ref, K_dens, p, s, cd, cl, mass);
    [k3,~] = npc_eom(t + h/2, y + (h*k2/2), area_ref, K_dens, p, s, cd, cl, mass);
    [k4,~] = npc_eom(t + h, y + k3*h, area_ref, K_dens, p, s, cd, cl, mass);
    y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    t = t + h;
    % save state (for debugging)
%     i = i + 1;
%     x(:,i) = y;  % ensure jettison event(s) test
%     ttemp(i) = t;
%     m(i) = mass;
    
    %% Update dependent states
%     alt = norm(y(1:3)) - p.planet.r_p;
    
    %% Check termination conditions
    if alt > p.planet.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt < p.planet.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
    elseif t > p.t_max
        pflag = 2; % Maximum time exceeded
        break;
    end
end % while ~plag

% calc orbit elements and params - get orbital apoapsis error
r = y(1:3); % m, position vector
v = y(4:6); % m/s, inertial velocity vector
oe = rv2oe([r;v],p.planet.mu);    % get orbital elements from r, v
a = oe(1);
e = oe(2);
ra = a*(1+e);       % apoapsis for given jettison time
delta_ra = ra - (p.ha_tgt + p.planet.r_p);   %apoapsis radius error, m

end % guid_sej_a_pred


