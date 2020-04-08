% guid_sej_e_pred.m 
%   Predictor for single-event jettison entry guidance algorithm
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   p_mode - uint8(1), nd, predictor mode
%   t_jettison - double(1), s, jettison time
%   K_dens - double(1), nd, atmospheric density multiplier
%   v_chute - double(1), m/s, current parachute deploy velocity
%   p - struct(1), multi, SEJ-E parameter structure
%
% Outputs:
%   t - double(1), s, final time
%   y - double(6), multi, final state
%   pflag - int(1), nd, predictor termination flag
%
% Major Revision History:
%   *Created AUG 2012, Z.R.Putnam

function [t,y,pflag] = guid_sej_e_pred(y0, t0, pred_mode, t_jettison, K_dens, v_chute, p)
%#codegen

%% Initialize
y = y0; % multi, initial state
t = t0; % s, initial time
pflag = 0; % nd, termination flag
area_ref = p.area_ref(1); % m^2, set reference area to initial value
mass = p.mass(1); % kg, set mass to initial value
cd = p.cd; % nd, set drag coefficient to initial value
h = p.t_inc_max; % s, step size
chute_deploy_count = 2*round(p.t_inc_max/p.t_inc_min); % use small time steps for 2 large time steps near parachute deploy
count = 0; % nd, initialize counter to zero
% A_sens_mag = double(0); % Sensed acceleration, m/s^2 (initialize to zero)
A_mag = double(0); 
parachute_deploy = false; % set deploy flag to false

% Compute derived quantites based on initial state
% alt = norm(y(1:3)) - p.p_r; % altitude, m (assumes spherical planet)
alt = guid_sej_e_alt(y(1:3),p);
V_pf = y(4:6) - cross(p.omega,y(1:3)); % planet-relative velocity, m/s
V_pf_mag = norm(V_pf); % planet-relative velocity magnitude, m/s

switch pred_mode    
    case 1 % pre-jettison to parachute deploy (no parachute model)
        v_term = v_chute;
        jflag = false;
        
    case 2 % pre-jettison to surface (with parachute model)
        v_term = -1000; % some negative number
        jflag = false;
        
    case 3 % post-jettison to parachute deploy (no parachute modeled)
        v_term = v_chute; % some negative number
        jflag = true; % immediately jettison skirt

    case 4 % post-jettison (with parachute)
        v_term = -1000; % some negative number
        jflag = true; % immediately jettison skirt       
        
    otherwise % duplicate case 1
        v_term = p.v_chute;
        jflag = false;
end


%% Integration (RK4)
while ~pflag

    
    %% Jettison logic    
    if jflag
        area_ref = p.area_ref(2); % m^2
        mass = p.mass(2); % kg
        jflag = false;
    end
    % Set step size
    if (t_jettison > t) && (t_jettison <= t+p.t_inc_max)
        h = t_jettison - t;
        jflag = true;
    else
        h = p.t_inc_max;
    end
    
    if (t_jettison > t) && (V_pf_mag < p.v_j_min)
        jflag = true;
        t_jettison = t;
    end

    % Get total accel
    [k1,dynp] = eom(y, area_ref, mass, cd, K_dens, parachute_deploy, p);
%     A_sens_mag = cd*dynp*area_ref/mass; % m/s^2, sensed acceleration
    A_mag = norm(k1(4:6));
    
    %% Parachute deploy logic
    if (t_jettison < t) ... % post jettison
            && ((p.mode == uint8(2)) || (p.mode == uint8(3))) % ...AND mode 2 or 3

        if (V_pf_mag >= v_chute) ... % chute deploy velocity condition not yet met
                && ((V_pf_mag - A_mag*h) <= v_chute) ... % chute deploy velocity condition will be met in the next time step or so
                && (count == 0)
            h = p.t_inc_min; % switch to smaller step size
            count = count + 1;
        elseif (count ~= 0) && (count <= chute_deploy_count)
            h = p.t_inc_min; % maintain smaller step size
            count = count + 1;
        else
            h = p.t_inc_max; % otherwise maintain maximum step size
        end
        
        if (V_pf_mag <= v_chute)
            parachute_deploy = true;
        end
        
    end
    
    %% Integrate one time step
    [k1,dynp] = eom(y, area_ref, mass, cd, K_dens, parachute_deploy, p);
    [k2,dynp] = eom(y+0.5*h*k1, area_ref, mass, cd, K_dens, parachute_deploy, p);
    [k3,dynp] = eom(y+0.5*h*k2, area_ref, mass, cd, K_dens, parachute_deploy, p);
    [k4,dynp] = eom(y+h*k3, area_ref, mass, cd, K_dens, parachute_deploy, p);
    y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + h; % increment time
    
    %% Update dependent states
%     alt = norm(y(1:3)) - p.p_r;
    alt = guid_sej_e_alt(y(1:3),p);
    V_pf = y(4:6) - cross(p.omega,y(1:3));
    V_pf_mag = norm(V_pf);
    
    %% Check termination conditions
    if V_pf_mag < v_term
        pflag = 4; % Parachute deploy velocity reached
        break;
    elseif alt > p.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt < p.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
    elseif t > p.t_max
        pflag = 3; % Maximum time exceeded
        break;
    end
end

end % guid_sej_a_pred



function [ydot, dynp] = eom(y, area_ref, mass, cd, K_dens, parachute_deploy, p)
% Equations of motion

%% Assign state vectors
R = y(1:3);
R_mag = norm(R);
R_u = R/R_mag;
V_inrtl = y(4:6);

%% Planet-fixed velocity
V_pf = V_inrtl - cross(p.omega,R);
V_pf_mag = norm(V_pf);
V_pf_u = V_pf/V_pf_mag;

%% Compute acceleration vector
% Aerodynamics
% alt = R_mag - p.p_r; % m, geocentric altitude
alt = guid_sej_e_alt(y(1:3),p);
dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3, atmosphere model density
dens = K_dens * dens_model; % kg/m^3, corrected density
dynp = 0.5*dens*V_pf_mag*V_pf_mag; % N/m^2, dynamic pressure
drag_u = -V_pf_u; % nd, unit vector along drag force vector
A_drag = (cd*dynp*area_ref/mass) * drag_u; % m/s^2, drag acceleration vector
if parachute_deploy
    A_drag = A_drag + (p.chute_cd*dynp*p.chute_area_ref/mass) * drag_u;
end
lift_u = [0;0;0]; % kludge for now
A_lift = (p.cl*dynp*area_ref/mass) * lift_u; % m/s^2, lift acceleration vector
A_aero = A_lift + A_drag; % m/s^2, total aerodynamic acceleration vector

% Gravity (inverse square)
% A_grav = -(p.mu/R_mag^2) * R_u; % m/s^2, gravity acceleration vector

% Gravity (j2)
temp_s1 = dot( p.npole, R_u ); % nd
temp_s2 = 1 - 5*temp_s1*temp_s1; % nd
temp_v2 = 2*temp_s1*p.npole; % nd
temp_v3 = temp_s2*R_u; % nd
temp_v2 = temp_v2 + temp_v3; % nd
temp_s1 = p.p_r / R_mag; % nd
temp_s1 = temp_s1*temp_s1; % nd
temp_s1 = 1.5*temp_s1*p.j2; % nd
U_G = R_u + temp_s1*temp_v2;
A_grav = -(p.mu / (R_mag*R_mag)) * U_G; % m/s^2, gravity acceleration vector

% Total acceleration vector
A_total = A_aero + A_grav; % m/s^2

%% Assign state derivatives
ydot = [V_inrtl; A_total];

end % eom