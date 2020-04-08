% define_g_dcf_sej.m 
%   Define drag modulation aerocapture guidance state structure for 
%   deceleration curve fit guidance.
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
%   g_dcf_sej - struct(1), multi, guidance data structure
%
% Major Revision History:
%   *Created JAN 2017, M. Werner
%   *Updated Terms, Oct 2018, E. Roelke

function [g_dcfj] = define_g_dcfj()
%#codegen

%% Parameter data structure
p = struct( ...
    'mass', double(0), ... % kg, vehicle mass
    'mass_jettison', double(0), ... % kg, mass that is jettisoned
    'area_ref', double(zeros(2,1)), ... % m^2, initial (1) and final (2) aerodynamic area reference
    'g1', double(0), ...    % m/s^2, deceleration value at which to initiate guidance
    'dt1', double(0), ... % s, time interval between g1 and g2
    'dt2', double(0), ... % s, time interval between g2 and g3 measurements
    'g2_nom', double(0), ...    % nominal g2 meas. (for decel multiplier)
    'mcflag',double(0), ... % monte carlo sim flag
    'decel_curve', double(zeros(4, 1)) ); % md, coefficients of deceleration linear curve fit


%% State data structure
s = struct( ...
    'stage', double(0), ... % nd, current jettison stage
    'jflag', logical(false), ...    % jettison flag
    'A_mag', double(0), ... % m/s^2, sensed acceleration magnitude
    'cmd_area_ref', double(0), ... % m^2, area_ref command (next area to switch to)
    'g1_time', double(0), ...   % s, time at which g1 was obtained
    'g2_time', double(0), ...   % s, time at which to obtain g2
    'g3_time',double(0), ...    % s, time at which to measure g3
    'g2', double(0), ...    % m/s^2, deceleration value after g1
    'g3', double(0), ...    % m/s2, decel value after g2
    'g_ratio', double(0), ...   % decel multiplier
    'togo', double(0), ...  % s, time-to-go until jettison
    'tj_curr',double(nan) ); % s, current jettison time


%% Input data structure
i = struct( ...
    't', double(0), ... % s, current time
    'R_pci', double(zeros(3,1)), ... % m, current position vector
    'V_inrtl_pci', double(zeros(3,1)), ... % m/s, current inertial velocity vector
    'V_pf_pci', double(zeros(3,1)), ... % m/s, current planet-relative, planet-fixed velocity vector
    'A_sens_pci', double(zeros(3,1)) ); % m/s^2, sensed acceleration magnitude

%% Construct combined input data structure
g_dcfj = struct( ...
    'i', i, ...    
    's', s, ...
    'p', p );

end % define_g_dcf_sej()
