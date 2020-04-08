% guid_sej_e.m 
%   Numeric predictor-corrector guidance algorithm for
%       single-event jettison entry systems
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
%   s - struct(1), multi, guidance state struc  ture
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created AUG 2012, Z. R. Putnam

function [s] = guid_sej_e( i, s, p, init_flag )
%#codegen

%% First-pass Initialization
if init_flag
    s = init(i, s, p);
end

%% Preliminary calculations
% Vector magnitudes
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude
s.t = i.t - s.t_ini; % s, guidance time, relative to initialization time


%% Sequencer
% Choose first step according to current phase
switch s.phase
    case 1
        s.next_step = uint8(1);
    case 2
        s.next_step = uint8(2);
    case 3
        s.next_step = uint8(3);
    case 4
        s.next_step = uint8(4);
    otherwise
        s.next_step = uint8(5);
end
% Sequence through steps
while (s.next_step ~= uint8(5))
    switch s.next_step
        
        case 1 % pre-entry hold
            s = pre_entry( i, s, p );
        
        case 2 % numeric predictor-corrector
            s = npc( i, s, p );
        
        case 3 % post-jettison hold
            s = post_jettison( i, s, p );
        
        case 4 % terminal descent
            s = terminal_descent( i, s, p );            
        
        otherwise % Shouldn't be possible
            s.next_step = uint8(5);
    end % switch s.next_step
end % while
   

end % guid_sej_e



function [ s ] = init(i, s, p)
% First-pass initialization function

%% Initialize states
s.jettison = logical(false); % nd, jettison flag
s.phase = uint8(1); % nd, current guidance phase
s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
s.t_inc = p.t_inc_max; % s
s.trust_region = p.trust_region_ini;
s.t_j_curr = p.t_jettison_ini; % s
s.t_j_next = p.t_jettison_ini; % s
s.K_dens = 1; % nd, assume on-board model is correct
s.range_to_tgt = 0; % m
s.flight_range = 0; % m
s.delta_range = 10*p.tol; % m
s.v_chute = p.v_chute;


% time at initialization
s.t_ini = i.t;

% Compute initial target position vector in ECEF
r_tgt_ECEF = zeros(3,1);
temp_s1 = (1.0 - p.flatness)^double(2); % nd
temp_s2 = p.p_r / sqrt( (cos(p.tgt_lat))^double(2) + ...
	temp_s1 * (sin(p.tgt_lat))^double(2) ); % m
temp_s3 = (temp_s2 + 0.0)*cos(p.tgt_lat); % m
r_tgt_ECEF(1) = temp_s3 * cos(p.tgt_lon); % m
r_tgt_ECEF(2) = temp_s3 * sin(p.tgt_lon); % m
r_tgt_ECEF(3) = (temp_s1*temp_s2 + 0.0) * ...
	sin(p.tgt_lat); % m

% Transform Earth spin axis and target position vector to ECI
ECEF_to_ECI = eye(3);
r_tgt_ECI = ECEF_to_ECI * r_tgt_ECEF; % m
s.u_tgt_ini = unit(r_tgt_ECI);


s.u_n_traj = unit(cross(i.R_pci,i.V_inrtl_pci));

% Altitude coefficients
% temp_s1 = (1 - p.flatness)^2;
% s.alt_coeff1 = temp_s1 - 1; % nd
% s.alt_coeff2 = p.p_r*(1 - p.flatness); % ft


end % init



function [ s ] = pre_entry(i, s, p)
% Pre-entry hold prior to reaching sensible atmosphere

%% Test for atmospheric entry
% TEST: has the sensible atmosphere been reached yet?
if s.A_mag > p.A_sens_atm
    % YES, sensible atmosphere has been reached; proceed to NPC phase    
    s.phase = uint8(2); % NPC phase
    s.next_step = uint8(2); % run NPC
else
    % NO, sensible atmosphere has not been reached; do nothing    
    s.next_step = uint8(5); % Execution complete
end

end % pre_entry



function [ s ] = npc(i, s, p)
% Numeric predidctor-corrector targeting algorithm to determine when to
% jettison drag device

%% Check for jettison
if (s.t_j_curr <= s.t) || (s.V_pf_mag <= p.v_j_min)
    
    %% Jettison complete or minimum velocity reached: switch to post-jettison phase
    
    s.jettison = true; % set flag
    s.phase = uint8(3); % exit phase  
    s.next_step = uint8(3); % exit phase

    
elseif p.mode
    
    %% Jettison conditions not reached: execute predictor-corrector targeting
    
    % Estimate density multiplier based on current acceleration
    s = est_K_dens(i, s, p);

    % Select predictor mode
    if (p.mode == uint8(2)) || (p.mode == uint8(3))
        % Modes 2 and 3: model parachute dynamics in predictor; target landing site
        s.pred_mode = uint8(2);
    else
        % Mode 1: do not model parachute dynamics, target parachute deploy point
        s.pred_mode = uint8(1);
    end
    
    % Initialize predictor-corrector
    y0 = [i.R_pci; i.V_inrtl_pci]; % md, initial state vector
    t0 = s.t; %i.t - s.t_ini; % s, initial time = current time - initialization time
    s.iter = uint8(0); % nd, iteration number
    cflag = 0; % nd, corrector flag
    
    while ~cflag
        s.iter = s.iter + uint8(1);
    
        [tf, yf, pflag] = guid_sej_e_pred(y0, t0, s.pred_mode, s.t_j_next, s.K_dens, s.v_chute, p);
        [s, cflag, cstatus] = guid_sej_e_corr(tf, yf, y0, i, s, p);
    
        if s.iter >= p.iter_max
            % ACTION: exit, max iterations exceeded
            cflag = 2; 
        end
    
    end % while cflag
    
    if (s.t_j_curr <= s.t)
    
        %% Jettison now and switch to phase 3
        s.jettison = true; % set flag
        s.phase = uint8(3); % high-beta phase
        s.next_step = uint8(3); % high-beta phase
    end
    
    s.next_step = uint8(5); % Execution complete
    
else
    s.next_step = uint8(5); % Execution complete

end

end % npc



function [ s ] = post_jettison(i, s, p)
% Post-drag-device-jettison hold and parachute deploy

if s.V_pf_mag <= s.v_chute
    s.deploy_decel = true;
    s.phase = uint8(4); % parachute descent
    s.next_step = uint8(4); % parachute decent

elseif p.mode == uint8(3)
    %% Range-based parachute deployment trigger

    s.iter = uint8(0); % nd, iteration number
    range_error = double(zeros(2,1));
    v_chute_dat = p.v_chute_box;
    dummy_tj = 0;
    
    % Estimate density multiplier based on current acceleration
    s = est_K_dens(i, s, p);
    
    if s.V_pf_mag > p.v_chute_box(2) % not in box yet
    
        v_chute_dat = p.v_chute_box;
        local_box = p.v_chute_box;
        
        % Predict to max chute deploy velocity
        s.pred_mode = uint8(3);


        % Initialize predictor-corrector
        y0 = [i.R_pci; i.V_inrtl_pci]; % md, initial state vector
        t0 = s.t;%i.t - s.t_ini; % s, initial time = current time - initialization time
        [t0,y0,pflag] = guid_sej_e_pred(y0,t0,s.pred_mode,dummy_tj,s.K_dens,p.v_chute_box(2),p);

        % Chute deploy at max velocity
        s.pred_mode = uint8(4);
        [tf,yf,pflag] = guid_sej_e_pred(y0,t0,s.pred_mode,dummy_tj,s.K_dens,p.v_chute_box(2),p);
    
    else % in box now

        v_chute_dat = [p.v_chute_box(1); s.V_pf_mag];
        local_box = [p.v_chute_box(1); s.V_pf_mag];
        
        % Initialize predictor-corrector
        y0 = [i.R_pci; i.V_inrtl_pci]; % md, initial state vector
        t0 = s.t;%i.t - s.t_ini; % s, initial time = current time - initialization time
        % Chute deploy at current velocity
        s.pred_mode = uint8(4);
        [tf,yf,pflag] = guid_sej_e_pred(y0,t0,s.pred_mode,dummy_tj,s.K_dens,s.V_pf_mag,p);
    
    end

    % Compute trajectory plane normal vector for downrange calculation
    vr0 = y0(4:6) - cross(p.omega,y0(1:3)); % Initial planet-relative velocity in PCI
    u_n = unit(cross(y0(1:3),vr0)); % Unit normal to trajectory plane
    
    % Compute downrange error
    % Compute target position at termination point
    ang = norm(p.omega)*tf; % angle to rotate target position vector about north pole
    u_tgt = s.u_tgt_ini*cos(ang) + ...
        p.npole*dot(p.npole,s.u_tgt_ini)*(1-cos(ang)) - ...
        cross(s.u_tgt_ini,p.npole)*sin(ang); % unit vector along estimated target position at parachute deploy
    % Compute unit vector along predicted final position vector
    u_R_pred = unit(yf(1:3));
    % Compute in-plane components of final position and target position vectors
    u_tgt_inplane = unit(u_tgt - dot(u_tgt,u_n)*u_n);
    u_R_pred_inplane = unit(u_R_pred - dot(u_R_pred,u_n)*u_n);
    % Compute angle to target (with divide-by-zero protection)
    temp = median([1,-1,dot(u_tgt_inplane,u_R_pred_inplane)]);
    theta = acos(temp);
    % Determine sign of downrange error
    if dot(cross(u_R_pred_inplane,u_tgt_inplane),u_n) < 0
        range_error(2) = theta*p.p_r;
    else
        range_error(2) = -theta*p.p_r;
    end
    
    if range_error(2) >= 0 % overshoot with earliest possible deploy
        s.v_chute = v_chute_dat(2);
        s.delta_range = range_error(2);
        s.ngood = 2;
    else
        % Chute deploy at min velocity
        s.pred_mode = uint8(4);
        [tf,yf,pflag] = guid_sej_e_pred(y0,t0,s.pred_mode,dummy_tj,s.K_dens,p.v_chute_box(1),p);
        
        % Compute downrange error
        % Compute target position at termination point
        ang = norm(p.omega)*tf; % angle to rotate target position vector about north pole
        u_tgt = s.u_tgt_ini*cos(ang) + ...
            p.npole*dot(p.npole,s.u_tgt_ini)*(1-cos(ang)) - ...
            cross(s.u_tgt_ini,p.npole)*sin(ang); % unit vector along estimated target position at parachute deploy
        % Compute unit vector along predicted final position vector
        u_R_pred = unit(yf(1:3));
        % Compute in-plane components of final position and target position vectors
        u_tgt_inplane = unit(u_tgt - dot(u_tgt,u_n)*u_n);
        u_R_pred_inplane = unit(u_R_pred - dot(u_R_pred,u_n)*u_n);
        % Compute angle to target (with divide-by-zero protection)
        temp = median([1,-1,dot(u_tgt_inplane,u_R_pred_inplane)]);
        theta = acos(temp);
        % Determine sign of downrange error
        if dot(cross(u_R_pred_inplane,u_tgt_inplane),u_n) < 0
            range_error(1) = theta*p.p_r;
        else
            range_error(1) = -theta*p.p_r;
        end
        
        if range_error(1) < 0
            s.v_chute = p.v_chute_box(1);
            s.delta_range = range_error(1);
            s.ngood = 1;
        else

            while 1
                s.iter = s.iter + uint8(1);
                s.v_chute = polate( v_chute_dat, range_error, local_box );
                               
                % Run predictor
                s.pred_mode = uint8(4);       
                [tf,yf,pflag] = guid_sej_e_pred(y0,t0,s.pred_mode,dummy_tj,s.K_dens,s.v_chute,p);

                % Compute downrange error
                % Compute target position at termination point
                ang = norm(p.omega)*tf; % angle to rotate target position vector about north pole
                u_tgt = s.u_tgt_ini*cos(ang) + ...
                    p.npole*dot(p.npole,s.u_tgt_ini)*(1-cos(ang)) - ...
                    cross(s.u_tgt_ini,p.npole)*sin(ang); % unit vector along estimated target position at parachute deploy
                % Compute unit vector along predicted final position vector
                u_R_pred = unit(yf(1:3));
                % Compute in-plane components of final position and target position vectors
                u_tgt_inplane = unit(u_tgt - dot(u_tgt,u_n)*u_n);
                u_R_pred_inplane = unit(u_R_pred - dot(u_R_pred,u_n)*u_n);
                % Compute angle to target (with divide-by-zero protection)
                temp = median([1,-1,dot(u_tgt_inplane,u_R_pred_inplane)]);
                theta = acos(temp);
                % Determine sign of downrange error
                if dot(cross(u_R_pred_inplane,u_tgt_inplane),u_n) < 0
                    temp_range_error = theta*p.p_r;
                else
                    temp_range_error = -theta*p.p_r;
                end
                
                if range_error(1) > range_error(2)
                    range_error(1) = temp_range_error; % m, estimated range error 
                    v_chute_dat(1) = s.v_chute;
                else
                    range_error(2) = temp_range_error; % m, estimated range error 
                    v_chute_dat(2) = s.v_chute;
                end                    
                
                if range_error(1) < p.tol
                    s.v_chute = v_chute_dat(1);
                    s.delta_range = range_error(1);
                    s.ngood = 3;
                    break;
                elseif range_error(2) < p.tol
                    s.v_chute = v_chute_dat(2);
                    s.delta_range = range_error(2);        
                    s.ngood = 4;
                    break;
                elseif s.iter >= p.iter_max
                    s.ngood = 5;
                    break;
                end
            end           
        end
    end

    s.next_step = uint8(5); % Execution complete

else
    
    s.next_step = uint8(5); % Execution complete
%     %% Original Parachute deploy logic - velocity trigger
%     if s.V_pf_mag < p.v_chute
%         s.deploy_decel = true;
%     end

end
    
% s.next_step = uint8(4); % Execution complete

end % post_jettison



function [ s ] = terminal_descent(i, s, p)
% Cost on parachutes for now

    % Do nothing
    s.next_step = uint8(5); % Execution complete

end % terminal_descent



function [ s ] = est_K_dens(i, s, p)
% Density factor estimator

if s.jettison
    idx = 2;
else
    idx = 1;
end
    
% Estimate density
dens_est = (2*p.mass(idx)*s.A_mag) / (s.V_pf_mag*s.V_pf_mag*p.area_ref(idx)*p.cd); % kg/m^3
% alt = norm(i.R_pci) - p.p_r; % m
alt = guid_sej_e_alt( i.R_pci, p );
dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3

% Density factor
K_dens = dens_est/dens_model; % nd

% Limit to min/max values
K_dens = median([p.K_dens_min, K_dens, p.K_dens_max]); % nd

% Low-pass filter
s.K_dens = p.K_dens_gain*K_dens + (1-p.K_dens_gain)*s.K_dens; % nd

% Log density estimate
s.dens_est = dens_est;

end % est_K_dens


function [ v_chute ] = polate( v_chute_dat, range_error, v_chute_box )

    m = (v_chute_dat(2)-v_chute_dat(1))/(range_error(2)-range_error(1));

    v_chute = v_chute_dat(1) - m*range_error(1);
    
    v_chute = median([v_chute_box(1) v_chute_box(2) v_chute]);
    
end % polate