% guid_sej_a_pred.m 
%   Predictor for single-event jettison aerocapture guidance algorithm
%   Propogates the entire trajectory forward including jettison event(s)
%   given jettison time in tj(i)
%
% Aeroassist Simulation source code
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   s - NPC state structure
%   p - NPC input structure
%   j_ind - current jettison index for actual vehicle state
%   dt - small jettison time change (for newton method)
%
% Outputs:
%   r_a - estimated apoapsis radius
%   dr_a - error in estimated apoapsis radius
%
% *Created Sep 2018, E. Roelke
% *Added outside loop for multiple stages, Mar 2019, E. Roelke    
% 
function [r_a,dr_a,pflag] = guid_dej_n_pred(y0, t0, s, p, j_ind, dt, guid)
%#codegen

% Initialize
% y = y0; % multi, initial state     y = [r_inertial;v_inertial]
% t = t0; % s, initial time
% h = p.t_inc_max; % s, (max) increment step size
h = 1/p.traj_rate;    %s, integration time step

K_dens = guid.s.atm.K_dens;

% for debugging
% i = 1;
% x(:,i) = y0;
% ttemp(i) = t0;

y = y0;
t = t0;
jett_flag = logical(false);
pflag = 0;

tj0 = s.tj_next(j_ind);        % jettison time this index
tj0 = tj0 + dt;

%% Integration (RK4)
while ~pflag
    
    if t >= tj0 && ( jett_flag == logical(false) )
        jett_flag = logical(true);
        j_ind = j_ind + 1;
    end
        
    % get area, mass of current stage
    area_ref = p.area_refs(j_ind);
    mass = p.masses(j_ind);
    cd = p.cds(j_ind);
    cl = p.cls(j_ind);
    
    % integrate with rk4
    [k1,alt] = gnc_eom(t, y, area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k2,~] = gnc_eom(t + h/2, y + (h*k1/2), area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k3,~] = gnc_eom(t + h/2, y + (h*k2/2), area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k4,~] = gnc_eom(t + h, y + k3*h, area_ref, K_dens, p, s, cd, cl, mass, guid);
    y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    t = t + h;
    
    %% Check termination conditions
    if alt > guid.p.planet.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt < guid.p.planet.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
%     elseif t > p.t_max
%         pflag = 3; % Maximum time exceeded
%         break;
    end
end % while ~plag

if pflag ~= 2
    r = y(1:3); % m, position vector
    v = y(4:6); % m/s, inertial velocity vector
    ae = rv2ae([r;v], guid.p.planet.mu);    % get orbital elements from r, v
    a = ae(1);
    e = ae(2);     

    r_a = a*(1+e);       % apoapsis for given jettison time
    dr_a = r_a - (p.tgt_ap + guid.p.planet.r_e);   %apoapsis radius error, m
   
    if (e > 1)
        pflag = 3;  % hyperbolic traj
        return;
    end
    
else % surface impact
    r_a = guid.p.planet.r_e; 
    dr_a = r_a - (p.tgt_ap + guid.p.planet.r_e);
    return;
end

end % guid_sej_a_pred
