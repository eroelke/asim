% define_guid.m 
%   Define guidance data structure
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   none
%   
% Outputs:
%   guid - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created 17 NOV 2011, Z.R. Putnam
%   *3 FEB 2012, Z.R. Putnam, added angle-of-attack commands, DC guidance
%       structure
%   *Modified 9 FEB 2012, M.J. Grant, added DC parameters
%   *AUG 2012, Z.R.Putnam, added single-event jettison entry guidance
%   *Add drag modulation guidance, Sep 2018, E. Roelke
%	*Add density corrector to main guidance struct, Jan 2019, E. Roelke
%	*Add atmospheric estimation params, Feb 2019, E. Roelke
%	*Add drag-modulation guidance mode, Feb 2019, E. Roelke

function guid = define_guid()
%#codegen

%% Command data structure
cmd = struct( ...
    'deploy_decel', logical(false), ...
    'aoa', double(0), ...
    'aoa_rate', double(0), ...
    'bank', double(0), ...
    'bank_rate', double(0), ...
    't_jettison', double(0), ... % jettison time(s), s
    'delta_mass', double(0), ...   % mass to jettison
    'area_ref', double(0) );     % area ref, m^2

%% Parameter data structure
% atmospheric estimation parameter inputs
p_atm = struct( ...
    'Kflag', uint8(0), ...      % density corrector flag (for monte carlo)
    'K_bounds', double([0.1, 2]), ...     % density corrector range
    'K_gain', double(0.1), ...    % density corrector gain
    'mode', uint8(1), ...   % atmospheric estimation mode (default density factor)
    'rss_flag', logical(false), ... %flag to save atm rss error
    'mc_ind', double(0), ...     % true atmosphere index (debugging, plotting purposes only)
    'ens_tol', double(0.01) ...   % ensemble tolerance (% of estimated density)
);

planet = struct( ... 
    'p_ind', uint8(2), ...    % planet index (earth default)
    'alt_max', double(0), ... % m, predictor max altitude
    'alt_min',double(0), ...    %m, predictor min alt
    'r_p', double(0), ...   % planet polar radius
    'r_e', double(0), ...   % planet equatorial radius
    'mu', double(0), ... % m^3/s^2, gravitational parameter
    'j2', double(0), ... % nd, J2 gravity coefficient  
    'npole', double(zeros(3,1)), ... % nd, planetary north pole
    'omega', double(zeros(3,1)), ... % rad/s, planetary rotation rate
    'atm_nom', double(zeros(1000,7)), ...  % expected (nominal) atm table
    'atm_true', double(zeros(1000,7)) ... %true atm table (for debugging/plotting only)
);

% overall paramer input struct
p = struct( ...
    'atm', p_atm, ...       % atmospheric estimation struct
    'planet', planet, ...   % planet input struct
    'aoa_mode', uint8(0), ...   % angle-of-attack modulation mode
    'bank_mode', uint8(0), ...  % bank modulation mode
    'prop_mode', uint8(0), ...  % propulsive guidance mode
    'dm_mode', uint8(0), ...   % drag-modulation mode
    'rate', double(0), ...
    'tgt_lat', double(0), ...
    'tgt_lon', double(0), ...
    'const_aoa', double(0), ...
    'const_bank',double(0), ...
    'const_bank_rate',double(0),...
    'aoa_schedule', double(zeros(100,2)),...
    'aoa_sched_len', double(0),...
    'bank_schedule',double(zeros(100,2)),...
    'bank_sched_len',double(0),...
    'energy_trigger_lambda',double(0),...
    'energy_trigger_h0',double(0),...
    'energy_trigger_v0',double(0),...
    'energy_trigger_hfloor',double(0) ... 
);

%% State data structure
% atmospheric estimation struct
s_atm = struct( ...
    'K_dens', double(1), ...    % density factor from main guidance computer
    'K_true', double(1), ...    % density factor of actual atmosphere
    'rho_est', double(0), ...   % current density estimate
    'atm_ind', double(1), ...   % index of current density history loc
    'atm_hist',double(nan(1000,4)), ...  % atmospheric history (alt, rho, T, P)
    'rho_K', double(0), ...   % density estimate from main guidance comp
    'atm_err', double(0), ...   % std error in the total atmosphere
    'atm_curr', zeros(1000,2), ...  % altitude, rho table of current best atm estimate
    'ind_curr', double(0), ...   %ensemble filter index
    'rss_nom', double(0), ...   %rss error for nominal atm
    'rss_K', double(0), ... % rss error for density corrector
    'rss_ens', double(0) ...    %rss error for ensemble filter
);

% master state struct
s = struct( ...
    'atm', s_atm, ...
    'init_flag',logical(false), ...    
    'sg_ratio', double(0), ...
    'begin_loft',logical(false) ... 
);

%% Algorithm-specific data structures
msl = define_g_msl; % MSL-prime entry guidance
dej_n = define_g_dej_n; % n-stage jettison aerocapture NPC guidance (drag modulation)
% sej_a = define_g_sej_a_old; % 1- and 2-stage jettison npc guidance
sej_e = define_g_sej_e; % Single-event-jettison entry guidance (drag modulation)
dcfj = define_g_dcfj; % Single-event-jettison entry deceleration-curve-fit guidance (drag modulation)
manual = define_g_manual;   % Single-event-jettison with user-input jettison time (drag modulation)
adm = define_g_adm; % Single-event-jettison entry guidance (analytic divert maneuver)
cvdma = define_cv_dma;
gt = define_g_gt; % Gravity turn terminal descent guidance
% dc = define_g_dc; % Dream Chaser numeric predictor-corrector guidance
pd = define_g_pdg; % predictive guidance algorithm for aerocapture

%% Construct combined data structure
guid = struct( ...
    'cmd', cmd, ... % command struct
    'p', p, ...     % param struct
    's', s, ...     % state struct
    'msl', msl, ...     % mars-science lab lift modulation
    'dej_n', dej_n, ... % n-stage discrete jettison NPC guidance
    'sej_e', sej_e, ...
    'dcfj', dcfj, ... % deceleration curve-fit guidance for jettison
    'manual', manual, ...   % manual jettison
    'cvdma', cvdma, ...  % continuous drag modulation
    'adm', adm, ... % analytic guidance
    'gt', gt, ...   % gravity turn
    'pd', pd ...   % predictive aerocapture guidance trigger (n event)
    );

end % define_guid()
