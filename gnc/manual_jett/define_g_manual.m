% define_g_manual.m 
%   Define drag modulation aerocapture guidance state structure for 
%   user-input jettison time.
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
%   g_manual - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created JAN 2017, M. Werner
%   *Added multi-jettison capabilitiy - Evan Roelke

function [g_manual] = define_g_manual()
%#codegen

%% Parameter data structure
p = struct( ...
    'masses', double(zeros(6,1)), ... % kg, vehicle masses
    'area_refs', double(zeros(6,1)), ... % m^2, reference areas
    'n_jett',uint8(1), ...     % number of jettison stages
    't_jetts', double(zeros(6,1)) ... % s, manual jettison times
    ); 


%% State data structure
s = struct( ...
    'jflag', logical(false(5,1)), ... % nd, jettison flag
    'stage',double(0), ...   % current vehicle stage (0-5)
    'j_ind',double(1) ...   % current jettison index(1-6) to load m,area
    );

%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

%% Construct combined input data structure
g_manual = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );

end % define_g_manual()
