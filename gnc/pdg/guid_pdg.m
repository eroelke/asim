% guid_pdg.m
%   Pure-Predictive Guidance Algorithm for Aerocapture
%
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Jul 2018
%
% Inputs:
%   i: guidance input structure
%   s: guidance state structure
%   p: guidance parameter structure
% 
function [s] = guid_pdg( i, s, p, guid )
%#codegen
% calculate guidance initiation parameters
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude

if s.A_mag >= p.A_sens_atm %&& s.stage < p.n_jett
    % Estimate density multiplier based on current acceleration
%     if (p.mcflag)
%         s = est_K_dens(i, s, p); % density corrector for monte carlo runs
%     else
%         s.K_dens = 1;
%     end

    % get state
    y0 = [i.R_pci; i.V_inrtl_pci];  % multi, initial state vector
    t0 = i.t;                       % s, current time     
    
    % run predictor with current vehicle state
%     j_ind = s.stage + 1;            % 1-indexed, current jettison stage
%     [ha,hp,zeta] = predict_orbit( y0, t0, j_ind, s, p );

    % run predictor with jettisoned vehicle state
    j_ind = s.stage + 2;            %1-indexed sc jettison phase
    [ha_j,hp_j,zeta_j, s] = predict_orbit( y0, t0, j_ind, s, p, guid );
    
%     dzeta = zeta - zeta_j;  % energy difference     
    
    % determine jettison trigger
    switch (p.trigger)
        case 0 % energy trigger
            % calculate desired orbital energy
            ha_des = (p.ha_tgt + guid.p.planet.r_p) - p.bias(s.stage + 1).tgt_ap;  % desired apoapsis radius
            hp_des = hp_j + guid.p.planet.r_p;  % periapsis radius (min altitude)
            e_des = (ha_des - hp_des)/(ha_des + hp_des);    % desired eccentricity
            sma_des = (ha_des+hp_des)/2; % m, desired semi-major axis
            zeta_des = -guid.p.planet.mu/(2*sma_des);   % desired orbital energy
            
            trig_val = (zeta_j - zeta_des); % difference in orbital energies
            
            % implement an orbital energy tolerance? not very intuitive..
            if trig_val <= 0
                s.jflag(s.stage+1) = logical(true);
            else
                s.jflag(s.stage+1) = logical(false);
            end
        case 1 % altitude trigger
            trig_val = ha_j;     % current jett apoapsis altitude
            % check within tolerance
            if trig_val <= (p.ha_tgt + p.ha_tol) && trig_val >= (p.ha_tgt - p.ha_tol)
                s.jflag(s.stage+1) = logical(true);
            % check if lower than target (haf decreases with time)
            elseif trig_val <= p.ha_tgt
                s.jflag(s.stage+1) = logical(true);
            else
                s.jflag(s.stage+1) = logical(false);
            end
            
        case 2 %velocity trigger
            % to do
            trig_val = norm(y0(4:6));
        case 3 %time trigger
            % to do
            trig_val = t0;
        case 4 %apoapsis altitude trigger
            % to do
            trig_val = ha_j - guid.p.planet.r_p;
        otherwise %energy trigger
            ha_des = p.ha_tgt + guid.p.planet.r_p;  % desired apoapsis radius
            hp_des = hp_j + guid.p.planet.r_p;  % periapsis radius
            e_des = (ha_des - hp_des)/(ha_des + hp_des);    % desired eccentricity
            sma_des = (ha_des+hp_des)/2; % m, desired semi-major axis
            zeta_des = -guid.p.planet.mu/(2*sma_des);   % desired orbital energy
            
            trig_val = (zeta_j - zeta_des); % difference in orbital energies
            
            % implement an orbital energy tolerance? not very intuitive..
            if (trig_val <= 0)
                s.jflag(s.stage+1) = logical(true);
            else
                s.jflag(s.stage+1) = logical(false);
            end
            
    end
    % save values
%     s.ha = ha;
    s.ha_j = ha_j;
%     s.zeta = zeta;
    s.trig_val = trig_val;
else
    s.jflag(s.stage+1) = false;
%     s.ha = nan;
    s.ha_j = nan;
%     s.zeta = nan;
%     s.delta_ap = nan;
    s.trig_val = nan;
end % sens atm check
end % pred_guid_ac



function [ha, hp, zeta, s] = predict_orbit( y0, t0, j_ind, s, p, guid )
% predict the orbit if stage had been jettisoned

% initialize
dt = 1/p.rate;   % current time integration step
flag = false;
y = y0;
t = t0;
% i = 1;
% rmag = norm(y0(1:3));
% rp = rmag;
% rpflag = logical(false);
% get vehicle parameters for current stage
cl = p.cls(j_ind);
cd = p.cds(j_ind);
area = p.area_refs(j_ind);
mass = p.masses(j_ind);
K_dens = guid.s.atm.K_dens;
i = 1;
h = nan(p.t_max/dt,1);
% V = h;
% integrate trajectory to completion
while ~flag   % run until termination condition is met
    
    %% Integrate (Runge-Kutta 4th order)
    
    k1 = gnc_eom(t, y, area, K_dens, p, s, cd, cl, mass, guid);
    k2 = gnc_eom(t + dt/2, y + (k1*dt/2), area, K_dens, p, s, cd, cl, mass, guid);
    k3 = gnc_eom(t + dt/2, y + (k2*dt/2), area, K_dens, p, s, cd, cl, mass, guid);
    k4 = gnc_eom(t + dt, y + k3*dt, area, K_dens, p, s, cd, cl, mass, guid);
    y = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + dt;
    h(i) = norm(y(1:3)) - guid.p.planet.r_p;
%     V(i) = norm(y(4:6));
    
%     vmag = norm(y(4:6));
    % check termination conditions
    alt = norm(y(1:3)) - guid.p.planet.r_p;
    if ((alt < 0) || (alt > guid.p.planet.alt_max) || (t >= p.t_max))
        flag = true;
    end
    i = i+1;

end

% get final orbit
r = y(1:3); % m, position vector
v = y(4:6); % m/s, inertial velocity vector
vmag = norm(v);
rmag = norm(r); % spacecraft radius

% oe = rv2oe([r;v],p.planet.mu);
ae = rv2ae([r;v], guid.p.planet.mu);
a = ae(1);
e = ae(2);
ra = a*(1+e); % apoapsis *radius*
% rp = a*(1-e);   % periapsis radius (inside atm)
hp = min(h);    % periapsis *altitude*
ha = ra - guid.p.planet.r_p; % apoapsis altitude
% rp = min(rmag);

% error - for data output
s.dr_ap_true = ra - (p.ha_tgt + guid.p.planet.r_e);   %apoapsis radius error, m
s.dr_ap = s.dr_ap_true - p.bias(s.stage + 1).tgt_ap; %minus because negative error is below target

zeta = 0.5*vmag^2 - guid.p.planet.mu/rmag;   % energy of orbit with jettison

end % predict_orbit
