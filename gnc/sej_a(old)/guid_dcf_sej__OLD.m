% guid_dcf_sej.m 
%   Deceleration curve fit aerocapture guidance algorithm for
%       single-event jettison systems
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state structure
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created Jan 2017, M. Werner

function [s] = guid_dcf_sej( i, s, p, init_flag )
%#codegenf

%% First-pass Initialization
if init_flag
    s = init(i, s, p);
end

%% Preliminary calculations
% Vector magnitudes
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude

% % get density factor for atmospheric uncertainties
% s = est_K_dens(i, s, p);

%% Sequencer
% Choose first step according to current phase
switch s.phase
    case 1
        s.next_step = uint8(1);
    case 2
        s.next_step = uint8(2);
    case 3
        s.next_step = uint8(3);
    otherwise
        s.next_step = uint8(4);
end
% Sequence through steps
while (s.next_step ~= uint8(4))
    switch s.next_step
        case 1 % pre-g1 hold
            s = pre_g1( i, s, p );
        case 2 % g2 measurement
            s = get_g2( i, s, p );
        case 3 % exit coast
            s = atm_exit( i, s, p );           
        otherwise % Shouldn't be possible
            s.next_step = uint8(4);
    end % switch s.next_step
end % while

end



function [ s ] = init(i, s, p)
% First-pass initialization function

%% Initialize states
s.jettison(1) = logical(false); % nd, jettison flag
s.jettison(2) = logical(false); % nd, jettison flags
s.phase = uint8(1); % nd, current guidance phase
s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
s.tj_curr = p.tj_ini;   % s
s.cmd_area_ref = p.area_ref(2);

end % init



function [ s ] = pre_g1(i, s, p)
% Pre-entry hold prior to reaching first deceleration value

%% Test for atmospheric entry
% TEST: has the g1 value been reached yet?
if s.A_mag >= p.g1
    % YES, sensible atmosphere has been reached; proceed to g2 phase
    s.g1_time = i.t;    % Time at which g1 was obtained 
    s.g2_time = s.g1_time + p.delT; % Time at which to obtain g2
    s.phase = uint8(2); % g2 phase
    s.next_step = uint8(2); % run g2
else
    % NO, g1 has not been reached; do nothing    
    s.next_step = uint8(3); % Execution complete
end

end % pre_entry



function [ s ] = get_g2(i, s, p)
% Obtains a second deceleration measurement after a set time interval, and
% then uses this measurement with the provided deceleration curve fit to
% obtain an estimate for jettison time.


%% Check for jettison
if s.tj_curr <= i.t
    %% Jettison complete, switch to exit phase
    s.jettison(1) = true; % set flag
    s.phase = uint8(3); % exit phase  
    s.next_step = uint8(4); % exit phase

elseif s.g1_time == 0
    %% g1 time has either not been obtained or is faulty, do nothing
    s.next_step = uint8(4);
    
    % Rounds to two decimal places to more accurately enable comparison
elseif round(s.g2_time * 100)/100 <= round(i.t*100)/100 && s.g2 == 0
    %% Obtain g2 measurement and use it to calculate jettison time

    s.g2 = s.A_mag;
    
    % Calculate jettison time
    s.togo = polyval(p.decel_curve, s.g2);
    s.tj_curr = s.togo + i.t;       % Time-to-go plus current time..not working right?
        
    % Check for jettison this time step (if t_togo + t < t)
    if s.tj_curr <= i.t
        s.jettison(1) = true;       % set flag
        s.phase = uint8(3);         % exit phase  
        s.next_step = uint8(4);     % Execution complete
    else
        s.next_step = uint8(4);     % Execution complete
    end
    
else
    s.next_step = uint8(4); % Execution complete

end

end % get_g2


%% Post-drag-device-jettison hold, performs no actions
function [ s ] = atm_exit(i, s, p)

    % Do nothing
    s.next_step = uint8(4); % execution complete

end % post_jettison




%% Density factor estimator
function [ s ] = est_K_dens(i, s, p)

%     if (p.mode == uint8(2)) && (s.phase == uint8(3))
%         idx = uint8(3);
%     else
%         idx = uint8(1);
%     end

    V_mag = norm(i.V_pf_pci);

    % Estimate density
    dens_est = (2*p.mass*s.A_mag) / (V_mag*V_mag*p.area_ref(1)*p.cd); % kg/m^3, dynam pres
    alt = norm(i.R_pci) - p.p_r; % m
    dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3

    % Density factor
    K_dens = dens_est/dens_model; % nd

    % Limit to min/max values
    K_dens = median([p.K_dens_min, K_dens, p.K_dens_max]); % nd

    % Low-pass filter
    s.K_dens = p.K_dens_gain*K_dens + (1-p.K_dens_gain)*s.K_dens; % nd

end % est_K_dens


function Amag_meas = get_Amag(y, area_ref, K_dens, p, mass)
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
alt = R_mag - p.p_r; % m, geocentric altitude
dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3, atmosphere model density
dens = K_dens * dens_model; % kg/m^3, corrected density
dynp = 0.5*dens*V_pf_mag*V_pf_mag; % N/m^2, dynamic pressure
drag_u = -V_pf_u; % nd, unit vector along drag force vector
% A_drag = (p.cd*dynp*area_ref/p.mass) * drag_u; % m/s^2, drag acceleration vector
A_drag = (p.cd*dynp*area_ref/mass) * drag_u; % m/s^2, drag acceleration vector
% lift_u = [0;0;0]; % kludge for now
% A_lift = p.cl*dynp*area_ref * lift_u / p.mass; % m/s^2, lift acceleration vector
% A_aero = A_lift + A_drag; % m/s^2, total aerodynamic acceleration vector

% Gravity (inverse square)
% A_grav = -(p.mu/R_mag^2) * R_u; % m/s^2, gravity acceleration vector

% Gravity (inverse-square, j2)
% temp_s1 = dot( p.npole, R_u ); % nd
% temp_s2 = 1 - 5*temp_s1*temp_s1; % nd
% temp_v2 = 2*temp_s1*p.npole; % nd
% temp_v3 = temp_s2*R_u; % nd
% temp_v2 = temp_v2 + temp_v3; % nd
% temp_s1 = p.p_r / R_mag; % nd
% temp_s1 = temp_s1*temp_s1; % nd
% temp_s1 = 1.5*temp_s1*p.j2; % nd
% U_G = R_u + temp_s1*temp_v2;
% A_grav = -(p.mu / (R_mag*R_mag)) * U_G; % m/s^2, gravity acceleration vector

% Total acceleration vector
% A_total = A_drag + A_grav; % m/s^2

Amag_meas = norm(A_drag);

%% Assign state derivatives
% ydot = [V_inrtl; A_total];

end % get_Amag