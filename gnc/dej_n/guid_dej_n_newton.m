% guid_dej_n_newton.m
%   n-stage dicrete event jettison aerocapture numerical predctor-corrector
%   guidance algorithm using Newton Raphson root-finding
% Entry systems Design Lab, University of Colorado Boulder
% Created By: Evan Roelke, Nov 2018
% 

function s = guid_dej_n_newton(y0,t0,s,p,j_ind)
% Initialize
t = t0; % s, current time
flag1 = false;
flag2 = false;
pflag = false;
h = 1/p.traj_rate;    %s, integration time step

tj0 = s.tj_next(j_ind);        % jettison time to try
tj1 = tj0 + h;      % differential jettison time change

jett_flag = logical(false); % predictor jettison flag
K_dens = s.K_dens;

alt1 = norm(y0) - p.planet.r_p;
alt2 = alt1;

while ~pflag
    
    if t >= tj0 && ( jett_flag == logical(false) )
        jett_flag = logical(true);
        j_ind = j_ind + 1;
        y1 = y; % copy pre-jettison state to new var
    end

    % get area, mass of current stage
    area_ref = p.area_refs(j_ind);
    mass = p.masses(j_ind);
    cd = p.cds(j_ind);
    cl = p.cls(j_ind);
    
    % 4th order runge-kutta integrator
    if flag1 == false
        k1 = npc_eom(t, y, area_ref, K_dens, p, s, cd, cl, mass);
        k2 = npc_eom(t + h/2, y + (h*k1/2), area_ref, K_dens, p, s, cd, cl, mass);
        k3 = npc_eom(t + h/2, y + (h*k2/2), area_ref, K_dens, p, s, cd, cl, mass);
        k4 = npc_eom(t + h, y + k3*h, area_ref, K_dens, p, s, cd, cl, mass);
        y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    end
    % get altitude for termination conditions
    alt1 = norm(y(1:3)) - p.planet.r_p;
    
    if t >= tj1 && flag2 == false
        k1 = npc_eom(t, y, area_ref, K_dens, p, s, cd, cl, mass);
        k2 = npc_eom(t + h/2, y + (h*k1/2), area_ref, K_dens, p, s, cd, cl, mass);
        k3 = npc_eom(t + h/2, y + (h*k2/2), area_ref, K_dens, p, s, cd, cl, mass);
        k4 = npc_eom(t + h, y + k3*h, area_ref, K_dens, p, s, cd, cl, mass);
        y1 = y1 + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
        
        alt2 = norm(y1(1:3)) - p.planet.r_p;
    end

    % Check termination conditions
    if alt1 > p.planet.alt_max && flag1 == false
        flag1 = 1; % Maximum altitude exceeded
        r1 = y(1:3); % m, position vector
        v1 = y(4:6); % m/s, inertial velocity vector
    elseif alt1 < p.planet.alt_min  && flag1 == false
        flag1 = 1; % Maximum altitude exceeded
        r1 = y(1:3); % m, position vector
        v1 = y(4:6); % m/s, inertial velocity vector
    end
    
    if alt2 > p.planet.alt_max && flag2 == false
        flag2 = 1; % Maximum altitude exceeded
        r2 = y(1:3); % m, position vector
        v2 = y(4:6); % m/s, inertial velocity vector
    elseif alt2 < p.planet.alt_min  && flag2 == false
        flag2 = 1; % Maximum altitude exceeded
        r2 = y(1:3); % m, position vector
        v2 = y(4:6); % m/s, inertial velocity vector
    end  
    
    if (flag1) && (flag2)
        pflag = true;
        break;
    end
    
    % update time
    t = t + h;
end % while loop

% calc orbit elements and params - get orbital apoapsis errors
oe1 = rv2oe([r1;v1],p.planet.mu);
a = oe1(1);
e = oe1(2);
s.r_ap = a*(1+e);       % apoapsis for given jettison time
s.delta_ap = s.r_ap - (p.tgt_ap + p.planet.r_p);   %curr apoapsis error, m

oe2 = rv2oe([r2;v2],p.planet.mu);
a = oe2(1);
e = oe2(2);
dr2 = a*(1+e) - (p.tgt_ap + p.planet.r_p);  % newton apoapsis error

% newton method
drdt = ( dr2 - s.delta_ap )/h;
update = s.delta_ap/drdt;   % update term, f/f'

% put bounds on t_jett update (in case of hitting surface,
% hyperbolic orbit)
% scale by magnitude of d?
if update >= 30
    update = 30;
elseif update <= -30
    update = -30;
end

% update next jettison time to try
s.tj_next(j_ind) = s.tj_next(j_ind) - update;

end % guid_dej_n_newton






