% define_g_dma_cv.m 
%   Define continously-variable drag modulation aerocapture guidance state
%   structure.
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
%   *Created DEC 2012, Z.R.Putnam

function [g_sej_a] = define_g_dma_cv()
%#codegen

%% Parameter data structure
p = struct( ...
    'iter_max', uint8(0), ... % nd, maximum allowable NPC iterations
    'A_sens_atm', double(0), ... % m/s^2, sensible atmosphere boundary
    'tgt_ap', double(0), ... % m, target apoapse altitude (RENAME?)
    'tol_ap', double(0), ... % m, apoapse tolerance
    'cd', double(0), ... % nd, drag coefficient (constant across jettison event)
    'cl', double(0), ... % nd, lift coefficient (constant across jettison event)
    'mass', double(0), ... % kg, vehicle mass
    'alt_min', double(0), ... % m, predictor minimum altitude
    'alt_max', double(0), ... % m, predictor maximum altitude    
    'K_dens_min', double(0), ... % nd, minimum allowable density factor
    'K_dens_max', double(0), ... % nd, maximum allowable density factor
    'K_dens_gain', double(0), ... % nd, density factor low-pass filter gain
    't_max', double(0), ... % s, predictor maximum integration time
    'p_r', double(0), ... % m, planet radius
    'mu', double(0), ... % m^3/s^2, gravitational parameter
    'j2', double(0), ... % nd, J2 gravity coefficient    
    't_inc_max', double(0), ... % s, maximum jettison time adjustment increment
    't_inc_min', double(0), ... % s, minimum jettison time adjustment increment
    'trust_region_ini', double(0), ... % s, initial trust region for root-finding
    'area_ref', double(zeros(2,1)), ... % m^2, minimum and maximum reference areas
    'npole', double(zeros(3,1)), ... % nd, planetary north pole
    'omega', double(zeros(3,1)), ... % rad/s, planetary rotation rate
    'atm_table', double(zeros(100,2)) ); % md, density vs altitude atm model

%% State data structure
s = struct( ...
    'exit_flag', logical(false), ... % nd, flag for atmospheric exit
    'phase', uint8(1), ... % nd, current guidance phase
    'next_step', uint8(1), ... % nd, next phase to execute on a given cycle
    'iter', uint8(0), ... % nd, NPC iteration count
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'cmd_area_ref', double(0), ... % m^2, area_ref command (next area to switch to)
    'K_dens', double(0), ... % nd, current density factor estimate
    'V_pf_mag', double(0), ... % m/s, planet-relative velocity magnitude
    't_inc', double(0), ... % s, current adjustment increment
    'ra_next',double(0), ... % s, next jettison time to try
    'ra', double(zeros(2,1)), ... % s, current (1) and previous (2) jettison times
    'ae', double(zeros(2,1)), ... % s, current (1) and previous (2) jettison times
    'bound', double(zeros(2,1)), ... % s, L and R bounding flags
    'ngood', double(0), ... % s, current (1) and previous (2) jettison times
    'trust_region', double(0), ... % s, current trust region
    'r_ap', double(0), ... % m, current (1) and previous (2) apoapse estimate
    'delta_ap', double(0) ); % m, current (1) and previous (2) apoapse error estimate

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

%% Construct combined input data structure
g_sej_a = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );


end % define_g_dma_cv()
