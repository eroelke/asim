% define_g_dej_n.m 
%   Define n-stage drag modulation aerocapture guidance state structure.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Colorado Center for Astrodynamics Research (CCAR)
%   University of Colorado boulder
%   
% Outputs:
%   g_dej_n - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created August 2018 - E. Roelke

function g_dej_n = define_g_dej_n()
%#codegen
multi_jett = struct( ...
    'pred_all', logical(false), ... %iterate over all stages each guidance step
    'force_jett', logical(false), ...    %force second stage jettison regardless of error change
    'comp_curr', logical(false) ...     % compare s.dr_ap to s.dr_curr AND s.dr_f(j_ind-1)
    );

p_hydra = struct( ...
    'dtj_flag', logical(false), ... %reduce dtj_lim or no
    'dtj_lim', double(20), ...  % max step size on newton method
    'n_dtj_lims', double(5) ...    % number of times bouncing between dtj_lim before its reduced (dtj_flag must be set)
);

%% Parameter data structure
bias = struct( ...
    'tgt_ap', double(0), ...    %target apoapsis bias
    'tj', double(0) ...     % jettison time bias
);

sigmas = struct( ...
    'm', double(0), ... %mass uncertainty percentage
    'cd', double(0), ... %drag coeff 1 sigma uncertinty percentage
    'aref', double(0) ...
);
p = struct( ...
    'iter_max', uint8(1), ... % nd, normal NPC max interal iterations
    'init_iters',uint8(2), ... %initialize each stage with X NPC iterations
    'A_sens_atm', double(0.25), ... % m/s^2, sensible atmosphere decel
    't_init', double(0), ...    %s, simulation time to start guidance
    'tgt_ap', double(0), ... % m, target apoapse altitude (RENAME?)
    'tol_ap', double(0), ... % m, apoapse tolerance
    'cds', zeros(6,1), ... % nd, drag coefficients for each stage
    'cls', zeros(6,1), ... % nd, lift coefficients for each stage
    'masses', double(zeros(6,1)), ... % kg, vehicle masses (n-stage + 1)   
    'tol', double(0), ... % jettison time convergence tolerance (for bisection)
    'tj_ini', double(zeros(5,1)), ... % s, initial guesses for jettison times
    'n_jett', uint8(0), ...     % number of jettison stages feasible
    't_max', double(0), ... %max predictor time
    'area_refs', double(zeros(6,1)), ... % m^2, n+1 stages
    'pred_mode', uint8(0), ...  % 0 reg, 1 saves data (for debugging)
    'npc_mode', uint8(0), ...  % corrector type (0 bisection, 1 newton, 2 multi, 3 ensemble)
    'traj_rate', double(0), ...  % integration rate for NPC (Hz) (bias if != sim rate)
    'hydra', p_hydra, ...      % hybrid NPC struct
    'multi', multi_jett, ...    % multi stage struct settings
    'bias', repmat(bias,5,1), ...     % biasing structure (per stage)
    'sigs', sigmas ...  % 1 sigma uncertainties
);
    

%% State data structure
s_hydra = struct( ...
    'dtj', double(0), ...   % last jettison time update
    'dtj_lim', double(0), ...   %current jettison time update limit
    'dtj_lim_count', double(0) ... % counter of dtj limit update
);

s = struct( ...
    'jettison', false(5,1), ... % nd, jettison flag
    'init_flags', false(5,1), ... %stage inititialization flags
    'stage', uint8(0), ...         % current jettison stage (0 being no jettison yet)
    'phase', uint8(1), ... % nd, current guidance phase
    'next_step', uint8(1), ... % nd, next phase to execute on a given cycle
    'iter', uint8(0), ... % nd, NPC iteration count
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'cmd_area_ref', double(0), ... % m^2, area_ref command (next area to switch to)
    'V_pf_mag', double(0), ... % m/s, planet-relative velocity magnitude
    't_inc', double(0), ... % s, current adjustment increment
    'tj_curr', double(zeros(5,1)), ... % s, current best estimate jettison time
    'tj_next', double(zeros(5,1)), ... % s, next jettison time to try
    'tj', double(zeros(2,1)), ... % s, current (1) and previous (2) jettison times
    'ae', double(zeros(2,1)), ... % s, current (1) and previous (2) jettison times
    'bound', double(zeros(2,1)), ... % s, L and R bounding flags
    'ngood', double(0), ... % s, current (1) and previous (2) jettison times
    'trust_region', double(0), ... % s, current trust region
    'step_size', double(0), ...   % step size for corrector
    'r_ap', double(0), ... % m, current (1) and previous (2) apoapse estimate
    'dr_ap', double(0), ...  % apoapsis radius error
    'dr_ap_true', double(0), ...    %apoapsis radius error (ignoring biasing)
    'dr_f', double(nan(5,1)), ...   % each jettison stage final apoapsis err
    'dr_curr', double(nan), ... % current apoapsis error estimate (for current jettison time)
    'ncalls', double(0), ...    % number of guidance calls
    'bank', double(0), ...       % vehicle bank angle
    'hydra', s_hydra ...    % HYDRA state struct
); % m, current (1) and previous (2) apoapse error estimate for each stage

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ... % m/s^2, sensed acceleration magnitude
);

%% Construct combined input data structure
g_dej_n = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p ...
);





% OLD - changed Oct 2018, E. Roelke
% p = struct( ...
%     'mode', uint8(0), ... % nd, perform targeting tasks or use initial guess for jettison
%     'iter_max', uint8(0), ... % nd, maximum allowable NPC iterations
%     'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
%     'tgt_ap', double(0), ... % m, target apoapse altitude (RENAME?)
%     'tol_ap', double(0), ... % m, apoapse tolerance
%     'cds', zeros(5,1), ... % nd, drag coefficients for each stage
%     'cls', zeros(5,1), ... % nd, lift coefficients for each stage
%     'masses', double(zeros(6,1)), ... % kg, vehicle masses (n-stage + 1)   
%     'K_dens_min', double(0), ... % nd, minimum allowable density factor
%     'K_dens_max', double(0), ... % nd, maximum allowable density factor
%     'K_dens_gain', double(0), ... % nd, density factor low-pass filter gain
%     't_max', double(0), ... % s, predictor maximum integration time
%     'planet', p, ...        % planet model
%     'tol', double(0), ... % jettison time convergence tolerance  
%     'tj_ini', double(zeros(5,1)), ... % s, initial guesses for jettison times
%     't_inc_max', double(0), ... % s, maximum jettison time adjustment increment
%     't_inc_min', double(0), ... % s, minimum jettison time adjustment increment
%     'n_jett', uint8(0), ...     % number of jettison stages feasible
%     'area_refs', double(zeros(6,1)), ... % m^2, n+1 stages
%     'K_flag', uint8(0), ...    % K_density flag (0 for unit density multiplier)
%     'corr_mode', uint8(0) ...      % corrector type (0 bisection, 1 newton)
%     );



end % define_g_seg_a()
