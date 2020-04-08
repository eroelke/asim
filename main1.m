% main.m 
%   Top-level main ASIM function; returns output data structure based on
%   input data structure
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - struct, multi, Main input structure; see documentation
%     
% Outputs:
%   dat - struct, multi, Main output structure; see documentation
%
% Major Revision History:
%   *Created 2008, M. Grant
%   *19 SEP 2011, B. Steinfeldt, cleaned extraneous functions
%   *27 SEP 2011, Z. Putnam, changed to structure input/output, 
%      removed inputs, mysave functions

function dat = main1( in ) %#codegen

%% Initialize
[dat, veh] = initialize(in);

%% Run Trajectory
[dat] = trajectory(in,veh,dat);

end % main()