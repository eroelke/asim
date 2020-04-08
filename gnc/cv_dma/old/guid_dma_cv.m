% guid_dma_cv.m 
%   Numeric predictor-corrector aerocapture guidance algorithm for
%       continuously-variable drag modulation systems
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
%   *Created DEC 2012, Z.R.Putnam

function [s] = guid_dma_cv( i, s, p, init_flag )
%#codegen

%% First-pass Initialization
if init_flag
    s = init(i, s, p);
end

%% Preliminary calculations
% Vector magnitudes
s.A_mag = norm(i.A_sens_pci); % m/s^2, sensed acceleration magnitude
s.V_pf_mag = norm(i.V_pf_pci); % m/s, planet-relative velocity magnitude

if s.A_mag > 10
    
    s.exit_flag = logical(true);
end


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
        case 1 % pre-entry hold
            s = pre_entry( i, s, p );
        case 2 % numeric predictor-corrector
            s = npc( i, s, p );
        case 3 % exit coast
            s = atm_exit( i, s, p );           
        otherwise % Shouldn't be possible
            s.next_step = uint8(4);
    end % switch s.next_step
end % while
   

end % guid_drag_jettison



function [ s ] = init(i, s, p)
% First-pass initialization function

%% Initialize states
s.phase = uint8(1); % nd, current guidance phase
s.next_step = uint8(1); % nd, next guidnace phase to evaluate on this cycle
s.A_mag = double(0); % m/s^2, sensed acceleration magnitude
s.t_inc = p.t_inc_max; % s
s.ra_next = p.area_ref(1); % s
s.K_dens = 1; % nd, initially assume on-board model is correct
s.r_ap = 0; % m
s.delta_ap = 0; % m
s.trust_region = p.trust_region_ini;

s.cmd_area_ref = p.area_ref(1);

s.ra = [0;0];
s.ae = [0;0];

end % init


%% Pre-entry hold prior to reaching sensible atmosphere
function [ s ] = pre_entry(i, s, p)

%% Test for atmospheric entry
% TEST: has the sensible atmosphere been reached yet?
if s.A_mag > p.A_sens_atm
    % YES, sensible atmosphere has been reached; proceed to NPC phase    
    s.phase = uint8(2); % NPC phase
    s.next_step = uint8(2); % run NPC
else
    % NO, sensible atmosphere has not been reached; do nothing    
    s.next_step = uint8(4); % Execution complete
end

end % pre_entry


%% Numeric predictor-corrector targeting
function [ s ] = npc(i, s, p)


%% Check for atmospheric exit
if s.exit_flag && s.A_mag < p.A_sens_atm

    % atmopsheric exit reached
    s.phase = uint8(3); % exit phase  
    s.next_step = uint8(3); % exit phase

else
    %% Execute predictor-corrector
    
    % Estimate density multiplier based on current acceleration
    s = est_K_dens(i, s, p);

    % Initialize predictor-corrector
    y0 = [i.R_pci; i.V_inrtl_pci]; % multi, initial state vector
    t0 = i.t; % s, initial time
    s.iter = uint8(0); % nd, iteration number
    s.bound = [0;0]; % nd, reset bounding storage
    s.ra = [0;0];
    s.ae = [0;0];
    cflag = 0; % nd, corrector flag   
    
    while ~cflag

        % run predictor to determine final state
        [t, y, pflag] = guid_dma_cv_pred(y0, t0, s.ra_next, s.K_dens, p); 
        
        % run corrector to adjust command area
        [s, cflag, ~] = guid_dma_cv_corr(t, y, pflag, s, p);
    
        if s.iter >= p.iter_max
            % ACTION: exit, max iterations exceeded
            cflag = 2; 
        else
            s.iter = s.iter + uint8(1); % Increment iteration counter
        end
        
    end % while cflag  

    % Assign command
    s.cmd_area_ref = s.ra_next;
    
    
    s.next_step = uint8(4); % Execution complete

end

end % npc




%% Post-drag-device-jettison hold, performs no actions
function [ s ] = atm_exit(i, s, p)

    % Do nothing
    s.next_step = uint8(4); % execution complete

end % post_jettison




%% Density factor estimator
function [ s ] = est_K_dens(i, s, p)

    % Estimate density
    dens_est = (2*p.mass*s.A_mag) / (s.V_pf_mag*s.V_pf_mag*s.cmd_area_ref*p.cd); % kg/m^3
    alt = norm(i.R_pci) - p.p_r; % m
    dens_model = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3

    % Density factor
    K_dens = dens_est/dens_model; % nd

    % Limit to min/max values
    K_dens = median([p.K_dens_min, K_dens, p.K_dens_max]); % nd

    % Low-pass filter
    s.K_dens = p.K_dens_gain*K_dens + (1-p.K_dens_gain)*s.K_dens; % nd

end % est_K_dens
