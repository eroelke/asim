% define_g_sej_e.m 
%   Define drag modulation entry guidance state structure.
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
%   g_sej_a - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created AUG 2012, Z.R.Putnam

function [g_sej_e] = define_g_sej_e()
%#codegen

%% Parameter data structure
p = struct( ...
    'mode', uint8(0), ... % nd, 0=use initial guess for jettison (no targeting); 1=nominal; 2=model chute; 3=chute range trigger
    'iter_max', uint8(0), ... % nd, maximum allowable NPC iterations
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
    'tgt_lat', double(0), ... % m, target latitude
    'tgt_lon', double(0), ... % m, target longitude
    'tol', double(0), ... % m, range miss tolerance
    'cd', double(0), ... % nd, drag coefficient (constant across jettison event)
    'cl', double(0), ... % nd, lift coefficient (constant across jettison event)
    'flatness', double(0), ... % nd, flatness parameter
    'v_j_min', double(0), ... % m/s, minimum skirt jettison velocity
    'v_chute', double(0), ... % m/s, velocity at which to deploy parachute
    'chute_area_ref', double(0), ... % m^2, parachute reference area
    'chute_cd', double(0), ... % nd, parachute drag coefficient
    'alt_min', double(0), ... % m, predictor minimum altitude
    'alt_max', double(0), ... % m, predictor maximum altitude    
    'K_dens_min', double(0), ... % nd, minimum allowable density factor
    'K_dens_max', double(0), ... % nd, maximum allowable density factor
    'K_dens_gain', double(0), ... % nd, density factor low-pass filter gain
    't_max', double(0), ... % s, predictor maximum integration time
    'p_r', double(0), ... % m, planet radius (equatorial)
    'mu', double(0), ... % m^3/s^2, gravitational parameter
    'j2', double(0), ... % nd, J2 gravity coefficient    
    't_jettison_ini', double(0), ... % s, initial guess for jettison time
    't_inc_max', double(0), ... % s, maximum jettison time adjustment increment
    't_inc_min', double(0), ... % s, minimum jettison time adjustment increment
    'trust_region_ini', double(0), ... % s, initial trust region for root-finding
    'mass', double(zeros(2, 1)), ... % kg, vehicle mass
    'area_ref', double(zeros(2,1)), ... % m^2, initial (1) and final (2) aerodynamic area reference
    'npole', double(zeros(3,1)), ... % nd, planetary north pole
    'omega', double(zeros(3,1)), ... % rad/s, planetary rotation rate
    'v_chute_box', double(zeros(2,1)), ... % m/s, parachute deployment box
    'atm_table', double(zeros(100,2)) ); % md, density vs altitude atm model

%% State data structure
s = struct( ...
    'jettison', logical(false), ... % nd, jettison flag
    'deploy_decel', logical(false), ... % nd, parachute deployment flag
    'phase', uint8(1), ... % nd, current guidance phase
    'next_step', uint8(1), ... % nd, next phase to execute on a given cycle
    'pred_mode', uint8(1), ... % nd, predictor mode
    'iter', uint8(0), ... % nd, NPC iteration count
    'u_tgt_ini', double(zeros(3,1)), ... % nd, unit vector along initial target position
    'u_n_traj', double(zeros(3,1)), ... % nd, unit vector perpendictular to trajectory plane
    't_ini', double(0), ... % s, time at guidance initialization
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'K_dens', double(0), ... % nd, current density factor estimate
    'dens_est', double(0), ... % kg/m^3, current density estimate
    'V_pf_mag', double(0), ... % m/s, planet-relative velocity magnitude
    't_inc', double(0), ... % s, current adjustment increment
    't_j_curr',double(0), ... % s, current best estimate jettison time
    't_j_next',double(0), ... % s, next jettison time to try
    't', double(0), ... % s, current time
    'tj', double(zeros(2,1)), ... % s, current (1) and previous (2) jettison times
    're', double(zeros(2,1)), ... % s, current (1) and previous (2) range errors
    'ngood', double(0), ... % s, current (1) and previous (2) jettison times
    'trust_region', double(0), ... % s, current trust region
    'range_to_tgt', double(0), ... % m, current estimated range to target
    'flight_range', double(0), ... % m, current estimated flight range
    'delta_range', double(0), ... % m, current range error estimate
    'v_chute', double(0) ); % m/s, current parachute deploy velocity

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

%% Construct combined input data structure
g_sej_e = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );


end % define_g_seg_e()
