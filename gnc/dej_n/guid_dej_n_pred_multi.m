% guid_sej_a_pred.m 
%   Predictor for single-event jettison aerocapture guidance algorithm
%   Propogates the entire trajectory forward including jettison event(s)
%   given jettison time in tj(i)
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
function [r_a,dr_a,pflag] = guid_dej_n_pred_multi(y0, t0, s, p, j_ind, dt, guid)
%#codegen

% Initialize
% y = y0; % multi, initial state     y = [r_inertial;v_inertial]
% t = t0; % s, initial time
pflag = 0; % nd, termination flag
% h = p.t_inc_max; % s, (max) increment step size
h = 1/p.traj_rate;    %s, integration time step
% h = 0.1;

% prevent prediction at previous time
% while tj0 <= t0
%     tj0 = tj0 + 5;
% end

K_dens = guid.s.atm.K_dens;

% for debugging
% i = 1;
% x(:,i) = y0;
% ttemp(i) = t0;

r_a = nan(p.n_jett,1);
dr_a = r_a;

for i = 1:(p.n_jett-(j_ind-1))
% loop through all stages
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
        
        % save jettison point for next stage
        y0 = y;
        t0 = t;
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

    % get altitude
%     alt = norm(y(1:3)) - p.planet.r_e;
    
    % save state (for debugging)
%     i = i + 1;
%     x(:,i) = y;  % ensure jettison event(s) test
%     ttemp(i) = t;
%     m(i) = mass;  

    %% Check expedience conditions
%     if alt <= p.alt_min
%         h = 2; % going to impact, expediate process
%     end
    
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
    oe = rv2ae([r;v], guid.p.planet.mu);    % get orbital elements from r, v
    a = oe(1);
    e = oe(2);     

    r_a(i) = a*(1+e);       % apoapsis for given jettison time
    dr_a(i) = r_a(i) - (p.tgt_ap + guid.p.planet.r_e);   %apoapsis radius error, m
   
    % check hyperbolic traj
    if ((i == 1) && (e > 1))
        pflag = 3;
        return;
    end
    
else
    r_a(1:end) = guid.p.planet.r_e;
    dr_a(1:end) = 0;
    return;
end

end % for loop

% check to see if jettison events worked - for debugging
% r = x(1:3,:);
% v = x(4:6,:);
% for i = 1:size(x,2)
%     vmag(i) = norm(v(:,i));
%     rmag(i) = norm(r(:,i));
% end
% plot(vmag./1000,rmag./1000)
% plot(ttemp,vmag./1000)

end % guid_sej_a_pred



% %% old EOM (had some discrepancy between simulation and this)
% function [ydot] = eom(t, y, area_ref, K_dens, p, s, cd, cl, mass)
% % Equations of motion
% 
% %% Assign state vectors
% R = y(1:3);
% R_mag = norm(R);
% R_u = R/R_mag;
% V_inrtl = y(4:6);
% 
% %% Planet-fixed velocity
% V_pf = V_inrtl - cross(p.planet.omega,R);
% V_pf_mag = norm(V_pf);
% V_pf_u = V_pf/V_pf_mag;
% 
% %% Compute acceleration vector
% % Aerodynamics
% alt = R_mag - p.planet.r_p; % m, geocentric altitude
% dens_model = lin_interp(p.planet.atm_table(:,1), p.planet.atm_table(:,2),alt); % kg/m^3, atmosphere model density
% dens = K_dens * dens_model; % kg/m^3, corrected density
% dynp = 0.5*dens*V_pf_mag*V_pf_mag; % N/m^2, dynamic pressure
% drag_u = -V_pf_u; % nd, unit vector along drag force vector
% % A_drag = (p.cd*dynp*area_ref/p.mass) * drag_u; % m/s^2, drag acceleration vector
% A_drag = (cd*dynp*area_ref/mass) * drag_u; % m/s^2, drag acceleration vector
% % lift_u = [0;0;0]; % kludge for now
% % A_lift = p.cl*dynp*area_ref * lift_u / p.mass; % m/s^2, lift acceleration vector
% % A_aero = A_lift + A_drag; % m/s^2, total aerodynamic acceleration vector
% 
% % Gravity (inverse square)
% A_grav = -(p.planet.mu/R_mag^2) * R_u; % m/s^2, gravity acceleration vector
% 
% % Gravity (j2)
% % g_ii_r2 = p.mu*mass/(R_mag^2);
% % J2_pp = [-mu/R_mag^3*pos_pp(1)*-3/2*J2*(re/pos_pp_mag)^2*(5*pos_pp(3)^2/pos_pp_mag^2-1); ...
% %          0; ...
% %          -mu*pos_pp(3)/pos_pp_mag^3*-3/2*J2*(re/pos_pp_mag)^2*(5*pos_pp(3)^2/pos_pp_mag^2-3)] * mass;
% % J2_ii = L_PI'*J2_pp;
% % gravity_ii = g_ii_r2*(-pos_ii_hat) + J2_ii; % Inertial gravity force [N]
% 
% % Total acceleration vector
% A_total = A_drag + A_grav; % m/s^2
% 
% %% Assign state derivatives
% ydot = [V_inrtl; A_total];
% 
% end % eom