% define_g_adm.m 
%   Define analytic drag modulation entry guidance state structure.
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
%   g_adm - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created JUL 2015, Z.R.Putnam

function [g_adm] = define_g_adm()
%#codegen

%% Parameter data structure
p = struct( ...
    'mode', uint8(0), ... % nd, 0=use initial guess for jettison (no targeting); 1=nominal; 2=model chute; 3=chute range trigger
    'iter_max', uint8(0), ... % nd, maximum allowable NPC iterations
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
    'tgt_lat', double(0), ... % m, target latitude
    'tgt_lon', double(0), ... % m, target longitude
    'tol', double(0), ... % m, range miss tolerance
    'surface_g', double(0), ... % m/s2, surface gravity acceleration
    'R', double(0), ... % planetary radius
    'H', double(0), ... % m, atmospheric scale height
    'rho_ref', double(0), ... % kg/m3, reference density
    'beta1', double(0), ... % kg/m2, initial ballistic coefficient
    'beta2', double(0), ... % kg/m2, final ballistic coefficient
    'V_j_min', double(0), ... % m/s, minimum skirt jettison velocity
    'V_chute', double(0), ... % m/s, velocity at which to deploy parachute
    'alt_min', double(0), ... % m, predictor minimum altitude
    'alt_max', double(0), ... % m, predictor maximum altitude    
    'K_dens_min', double(0), ... % nd, minimum allowable density factor
    'K_dens_max', double(0), ... % nd, maximum allowable density factor
    'K_dens_gain', double(0), ... % nd, density factor low-pass filter gain
    't_max', double(0), ... % s, predictor maximum integration time
    't_jettison_ini', double(0), ... % s, initial guess for jettison time
    'area_ref', double(zeros(2,1)), ... % m^2, initial (1) and final (2) aerodynamic area reference
    'npole', double(zeros(3,1)), ... % nd, planetary north pole
    'omega', double(zeros(3,1)), ... % rad/s, planetary rotation rate
    'atm_table', double(zeros(100,2)) ); % md, density vs altitude atm model

%% State data structure
s = struct( ...
    'jettison', logical(false), ... % nd, jettison flag
    'deploy_decel', logical(false), ... % nd, parachute deployment flag
    'phase', uint8(1), ... % nd, current guidance phase
    'next_step', uint8(1), ... % nd, next phase to execute on a given cycle
    'u_tgt_ini', double(zeros(3,1)), ... % nd, unit vector along initial target position
    'u_n_traj', double(zeros(3,1)), ... % nd, unit vector perpendictular to trajectory plane
    'gs1', double(0), ... % rad, constant flight-path angle
    'gs2', double(0), ... % rad, constant flight-path angle  
    't_ini', double(0), ... % s, time at guidance initialization
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'K_dens', double(0), ... % nd, current density factor estimate
    'dens_est', double(0), ... % kg/m^3, current density estimate
    'V_pf_mag', double(0), ... % m/s, planet-relative velocity magnitude
    't_inc', double(0), ... % s, current adjustment increment
    't_j_curr',double(0), ... % s, current best estimate jettison time
    't_j_next',double(0), ... % s, next jettison time to try
    't', double(0), ... % s, current time
    'trust_region', double(0), ... % s, current trust region
    'range_to_tgt', double(0), ... % m, current estimated range to target
    'flight_range', double(0), ... % m, current estimated flight range
    'delta_range', double(0), ... % m, current range error estimate
    'V_chute', double(0) ); % m/s, current parachute deploy velocity

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

%% Construct combined input data structure
g_adm = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );


end % define_g_seg_e()
