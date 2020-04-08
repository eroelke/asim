% build_mex.m 
%   Generates compiled mex file of simulation
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   None
%   
% Outputs:
%   None
%
% Major Revision History:
%   *Created SEP 2011, Z. Putnam
%   *28 SEP 2011, Z. Putnam, functionalized, added ASIM-specific
%       configuration options 
%   *28 MAR 2014, ZR Putnam, added more compiler options, changed main file
%       from "main" to "asim"

function [] = build_mex()

%% Build MEX file

fprintf('Aeroassist Simulation Code Builder');
fprintf('\nBuilding MEX code...\n');

tic; % start timer

% Define input arguments (input structure)
args = { define_in };

% Create custom configuration options file
config = coder.config('MEX');
config.DynamicMemoryAllocation = 'AllVariableSizeArrays'; % Allow variable size arrays
% config.GenCodeOnly = true;

% Build MEX file
codegen -config config -report main1.m -args args;

fprintf('...done. Elapsed time: %6.2fs\n', toc);


end % build_mex()
