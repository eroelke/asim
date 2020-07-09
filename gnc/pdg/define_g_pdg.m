% define_g_pdg.m 
%   Define drag modulation aerocapture predictive guidance state structure
%
% Aeroassist Simulation source code (ASIM)
%
% Developed by:
%   CU Boulder
%
% Inputs:
%   none
%   
% Outputs:
%   g_pdg - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created Jun 2018, E. Roelke
% 
function g_pdg = define_g_pdg()
%#codegen
sigmas = struct( ...
    'm', double(0), ... %mass uncertainty percentage
    'cd', double(0), ... %drag coeff 1 sigma uncertinty percentage
    'aref', double(0) ...
);

bias = struct( ...
    'tgt_ap', double(0), ...    %target apoapsis bias
    'tj', double(0) ...     % jettison time bias
);


%% Parameter data structure
p = struct( ...
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
    'ha_tgt', double(0), ... % m, target apoapse altitude
    'ha_tol',double(0), ...  % m, apoapsis altitude tolerance
    'cds', double(zeros(6,1)), ... % nd, drag coefficient (constant across jettison event)
    'cls', double(zeros(6,1)), ... % nd, lift coefficient (constant across jettison event)
    'n_jett',uint8(1), ... % number of jettison events
    'masses', double(zeros(6,1)), ... % kg, spacecraft masses for each stage
    'mcflag', uint8(0), ...     % monte carlo analysis flag (to calc K_dens)
    't_max', double(0), ... % s, predictor maximum integration time
    'tj0', double(zeros(5,1)), ... % s, initial guesses for jettison times
    'rate', double(0), ...      % integration rate
    'trigger', uint8(0), ...    % jettison trigger
    'area_refs', double(zeros(6,1)), ... % m^2, initial (1) and subsequent reference areas
    'bias', repmat(bias,5,1), ...     % biasing structure (per stage)
    'stage_skip', 10, ... %for testing purposes only
    'sigs', sigmas ...   %uncertainty
);

%% State data structure
s = struct( ...
    'jflag', double(zeros(5,1)), ... % nd, jettison flag
    'phase', uint8(1), ... % nd, current guidance phase
    'jett_curr', double(0), ... % current jettison number
    'next_step', uint8(1), ... % nd, next phase to execute on a given cycle
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'cmd_area_ref', double(0), ... % m^2, area_ref command (next area to switch to)
    'K_dens', double(0), ... % nd, current density factor estimate
    'ha_j', double(0), ...  % m, apoapsis estimate w/ jettison
    'trig_val',double(0), ...   % trigger parameter value
    'delta_ap', double(0), ...  % predicted apoapsis error
    'bank', double(0), ...      % bank angle
    'dr_ap', double(0), ...
    'dr_ap_true', double(0), ...
    'stage', double(0) ...  % current stage
); % m, current apoapsis error estimate

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ...     % m/s^2, sensed acceleration magnitude
);

%% Construct combined input data structure
g_pdg = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p ...
);


end
