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
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created Jan 2017, M. Werner
%   *Added g3, Jul 2018, E. Roelke

function [s] = guid_dcfj( i, s, p )
%#codegen

% initialize vals
% s.jflag = logical(false);

% get accel data
s.A_mag = norm(i.A_sens_pci);

if (~isnan(s.tj_curr))
    return;
% check if algorithm required
elseif s.A_mag >= p.g1
% else, proceed with measurement
% get first-measurement data (if haven't before)
    if s.g1_time == 0
        s.g1_time = i.t;     % Time at which g1 was obtained (now)
        s.g2_time = s.g1_time + p.dt1; % Time at which to obtain g2
        s.g3_time = s.g2_time + p.dt2; % Time at which to obtain g3
    end

    % Check if within g2 measurement time (within 2 decimals)
    if i.t >= s.g2_time
        s.g2 = s.A_mag;
        % Calculate jettison time with decel curve and g2
        if (p.mcflag)
            s.g_ratio = s.g2/p.g2_nom;
        else
            s.g_ratio = 1;
        end

        s.togo = polyval(p.decel_curve, s.g2/s.g_ratio);
%             s.togo = polyval(p.decel_curve, s.g2);
        s.tj_curr = s.togo + i.t;       % Time-to-go plus current time
%             s.tj_curr = s.togo + i.t - p.dt1;

        % check g3 update - doesn't exist yet
        % if g3 is above curve, atmos is more dense (earlier t_jett)
        % if g3 is below curve, atmos is less dense (later t_jett)

        % Check for jettison this time step (if t_togo + t < t)
        if i.t >= s.tj_curr
            % set jettison
            if s.jflag == false
                s.stage = s.stage + 1;
                s.jflag = logical(true);      % set flag
                s.cmd_area_ref = p.area_ref(2);
            end
        end

    else
        s.jflag = logical(false);     % else, not time to jettison
        return;
    end % g2 measure time
else % no guid algorithm
    s.jflag = logical(false);
end % decel check

end % guid_dcf_sej
    

%% Density factor estimator
% function [ s ] = est_K_dens(i, s, p)
% %     j_ind = s.stage+1;
% %     mass = p.masses(j_ind);
% %     area = p.area_ref
% %     if (p.mode == uint8(2)) && (s.phase == uint8(3))
% %         idx = uint8(3);
% %     else
% %         idx = uint8(1);
% %     end
% 
%     V_mag = norm(i.V_pf_pci);
% 
%     % Estimate density
%     dens_est = (2*p.mass*s.A_mag) / (V_mag*V_mag*p.area_ref(1)*p.cd); % kg/m^3, dynam pres
%     alt = norm(i.R_pci) - p.r_p; % m
%     dens_model = lin_interp(p.atm(:,1), p.atm(:,2),alt); % kg/m^3
% 
%     % Density factor
%     K_dens = dens_est/dens_model; % nd
% 
%     % Limit to min/max values
%     K_dens = median([p.K_bounds(1), K_dens, p.K_bounds(2)]); % nd
% 
%     % Low-pass filter
%     s.K_dens = p.K_gain*K_dens + (1-p.K_gain)*s.K_dens; % nd
% 
% end % est_K_dens
% 
% 
% function Amag_meas = get_Amag(y, area_ref, K_dens, p, mass)
% % Equations of motion
% 
% %% Assign state vectors
% R = y(1:3);
% R_mag = norm(R);
% R_u = R/R_mag;
% V_inrtl = y(4:6);
% 
% %% Planet-fixed velocity
% V_pf = V_inrtl - cross(p.omega,R);
% V_pf_mag = norm(V_pf);
% V_pf_u = V_pf/V_pf_mag;
% 
% %% Compute acceleration vector
% % Aerodynamics
% alt = R_mag - p.r_p; % m, geocentric altitude
% dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3, atmosphere model density
% dens = K_dens * dens_model; % kg/m^3, corrected density
% dynp = 0.5*dens*V_pf_mag*V_pf_mag; % N/m^2, dynamic pressure
% drag_u = -V_pf_u; % nd, unit vector along drag force vector
% % A_drag = (p.cd*dynp*area_ref/p.mass) * drag_u; % m/s^2, drag acceleration vector
% A_drag = (p.cd*dynp*area_ref/mass) * drag_u; % m/s^2, drag acceleration vector
% % lift_u = [0;0;0]; % kludge for now
% % A_lift = p.cl*dynp*area_ref * lift_u / p.mass; % m/s^2, lift acceleration vector
% % A_aero = A_lift + A_drag; % m/s^2, total aerodynamic acceleration vector
% 
% % Gravity (inverse square)
% % A_grav = -(p.mu/R_mag^2) * R_u; % m/s^2, gravity acceleration vector
% 
% % Gravity (inverse-square, j2)
% % temp_s1 = dot( p.npole, R_u ); % nd
% % temp_s2 = 1 - 5*temp_s1*temp_s1; % nd
% % temp_v2 = 2*temp_s1*p.npole; % nd
% % temp_v3 = temp_s2*R_u; % nd
% % temp_v2 = temp_v2 + temp_v3; % nd
% % temp_s1 = p.p_r / R_mag; % nd
% % temp_s1 = temp_s1*temp_s1; % nd
% % temp_s1 = 1.5*temp_s1*p.j2; % nd
% % U_G = R_u + temp_s1*temp_v2;
% % A_grav = -(p.mu / (R_mag*R_mag)) * U_G; % m/s^2, gravity acceleration vector
% 
% % Total acceleration vector
% % A_total = A_drag + A_grav; % m/s^2
% 
% Amag_meas = norm(A_drag);
% 
% %% Assign state derivatives
% % ydot = [V_inrtl; A_total];
% 
% end % get_Amag