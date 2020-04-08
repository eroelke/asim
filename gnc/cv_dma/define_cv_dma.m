% define_cv_dma.m
%   define input and trajectory structures for continous drag modulation
%   aerocapture
% 
% Entry systems Design Lab (EsDL)
% Colorado Center for Astrodynamics Research (CCAR)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Dec 2018
% 
function g_cvdma = define_cv_dma()

%% Planet model data structure
planet = struct( ... 
    'alt_min', double(0), ... % m, predictor minimum altitude
    'alt_max', double(0), ... % m, predictor maximum altitude 
    'r_p', double(0), ...   % planet polar radius
    'r_e', double(0), ...   % planet equatorial radius
    'mu', double(0), ... % m^3/s^2, gravitational parameter
    'j2', double(0), ... % nd, J2 gravity coefficient  
    'npole', double(zeros(3,1)), ... % nd, planetary north pole
    'omega', double(zeros(3,1)), ... % rad/s, planetary rotation rate
    'atm_table', double(zeros(1000,7)) ...  % nominal atm table
    );

%% Input Parameter Data Structure
p = struct( ... 
    'mode', uint8(0), ...   % guidance mode (0 ref traj, 1 predictor-corrector)
    'decel', double(zeros(1000,2)), ... % [alt decel] table of nom. traj
    'iter_max', uint8(0), ... % nd, maximum iterations
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
    'ha_tgt', double(0), ... % m, target apoapse altitude
    'ha_tol', double(0), ... % m, apoapse tolerance
    'planet', planet, ...        % planet model
    'cd_ini', double(0), ...    % initial drag coefficient
    'rn', double(0), ...    % nose radius, m
    'delta_c_ini', double(0), ...   % intiial cone half-angle (rad)
    'rc_ini', double(0), ...  % initial base radius
    'area_rate',double(0), ...  % max/min area change rate
    'rc_bounds',double(zeros(2,1)), ...   % max/min backshell radii vals (low, high)
    'mass', double(0), ...  % vehicle mass (const)
    'c', double(0), ... % chord length (m)
    'cd_bounds', double(zeros(2,1)), ...    % reasonable cd_bounds
    't_max',double(0), ...   % max predictor time length
    'traj_rate', double(0) ...  % integration rate for NPC (Hz) (bias if != sim rate)
    );

%% Guidance State Structure
s = struct( ... 
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'K_dens', double(0), ...    % density corrector estimate (from outside guidance)
    'Aref', double(0), ...  % reference area (m2)
    'cd',double(0), ... % current drag coefficient
    'delta_c', double(0), ... % current cone half-angle (rad)
    'rc', double(0), ...    % current base radius (m)
    'ra', double(0), ...    % current apoapsis radius estimate (m)
    'delta_ra', double(0), ...  % current apoapsis radius error (m)
    'V_pf_mag', double(0), ...  % planet-relative velocity mag (m/s)
    'iter', uint8(0), ...   % iteration number
    'beta',double(0), ...   % ballistic coefficient (kg/m2)
    'rflag', logical(false), ...  % initialize system
    'sig_cd', double(0), ...    % uncertainty on drag coeff
    'bank', double(0) ...       % vehicle bank angle
    );

%% Guidance Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

g_cvdma = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );



end %define_cv_dma