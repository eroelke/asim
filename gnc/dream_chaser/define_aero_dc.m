% define_aero_dc.m 
%   Define Dream Chaser aerodynamics data structure
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
%   dc - struct(1), multi, Dream Chaser aerodynamics data structure
%
% Major Revision History:
%   *Created 21 FEB 2012, Z.R. Putnam

function dc = define_aero_dc()
%#codegen

% temp = load(['.' filesep 'data' filesep 'aero_dc_v52.mat']);
temp = load('data/aero_data/aero_dc_v52.mat');
dc = temp.dc;

end % define_aero_dc()