% define_g_dba_n.m
%   input structure for discrete bank-angle switch guidance
% 
% Written by: E. Roelke, Jan 2020
% 
% 
function g_dba = define_g_dba_n()
%% Parameter data structure
p = struct( ...
    'iter_max', uint8(0), ... % nd, maximum iterations (for bisection)
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
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
    'npc_mode', uint8(0), ...  % corrector type (0 bisection, 1 newton, 2 multi, 3 ensemble)
    'traj_rate', double(0), ...  % integration rate for NPC (Hz) (bias if != sim rate)
    'bank_angles', double(zeros(5,1)) ...       % vehicle bank angles to command
);
    

%% State data structure
s = struct( ...
    'switch', false(5,1), ... % nd, switch flag
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
    'dr_f', double(nan(5,1)), ...   % each jettison stage final apoapsis err
    'dr_curr', double(nan), ... % current apoapsis error estimate (for current jettison time)
    'ncalls', double(0), ...    % number of guidance calls
    'bank', double(0) ...  % current bank angle
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
g_dba = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p ...
);

end