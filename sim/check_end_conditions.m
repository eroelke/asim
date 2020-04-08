% check_end_conditions.m 
%   Checks simulation termination conditions.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - data structure, misc, input data structure
%   dat - data structure, misc, current vehicle trajectory data
%   calcs - data strucutre, misc, current vehicle misc. calculations
%   i - int8(1), nd, current length of the state vector 
%   
% Outputs:
%   flag - int8(1), nd, vehicle drag coefficient at angle of attack
%
% Major Revision History:
%   *Modified from original code
%   *Original: 2008, M. Grant
%   *Revised 20 SEP 2011, B. Steinfeldt, Added modularity to code
%   *28 SEP 2011, Z.R. Putnam, updated input structure variables
%   *7 NOV 2011, Z.R. Putnam, changed inputs to accept current and previous
%       calcs structure instead of dat structure

function [flag] = check_end_conditions(in,calcs_curr,calcs_prev)
%#codegen

% Assign variables
flag = 0;

% Determine which side end trigger is during previous and current time step
switch in.s.term.mode
    case 1 % Altitude
        sign1 = sign(calcs_prev.alt - in.s.term.value);
        sign2 = sign(calcs_curr.alt - in.s.term.value);
    case 2 % Mach
        sign1 = sign(calcs_prev.mach - in.s.term.value);
        sign2 = sign(calcs_curr.mach - in.s.term.value);
    otherwise
        sign1 = nan;
        sign2 = nan;
end

if sign1 + sign2 == 0
    % Value has been crossed, determine direction
    if (in.s.term.direction == 1) && (sign2 == 1)
        flag = 1;
    elseif (in.s.term.direction == -1) && (sign2 == -1)
        flag = 1;
    elseif in.s.term.direction == 0
        flag = 1;
    end
end

% Determine if vehicle mass is below zero
if calcs_curr.mass < 0
    flag = 1;
end

end % check_end_conditions()
