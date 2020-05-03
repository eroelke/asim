% define_veh.m 
%   Define vehicle data structure for use inside trajectory
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
%   veh - struct(1), multi, vehicle data structure
%
% Major Revision History:
%   *Created AUG 2012, Z.R.Putnam
%   *Add cd, E. Roelke

function veh = define_veh()
%#codegen

%% Parameter data structure
p = struct( ...
    'm_jettison', double(0) );

%% State data structure
s = struct( ...
    'not_jettisoned', logical(true), ...
    'compute_trim', uint8(0), ...
    'area_ref', double(0), ...
    'aoa', double(0), ...
    'ssa', double(0), ...
    'bank', double(0), ...
    'cd', double(0), ...
    'mass', double(0), ...
    'prop_remain',logical(true));

%% Construct combined input data structure
veh = struct( ...
    's', s, ...
    'p', p );

end % define_veh()