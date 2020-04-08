% termination.m 
%   Checks simulation termination conditions and returns termination flag 
%   and status.
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
%   calcs_curr - struct, misc, current integration time step trajectory 
%       calculations
%   calcs_prev - struct, misc, previous integration time step trajectory 
%       calculations
%   
% Outputs:
%   flag - uint8(1), nd, termination flag
%   cond - uint8(10), nd, vector of flags indicating active conditions
%
% Major Revision History:
%   *Created 9 NOV 2011, Z.R. Putnam from check_end_conditions.m

function [flag,cond] = termination( in, calcs_curr, calcs_prev)
%#codegen

%% Initialize
flag = logical(false);
or_flag = uint8(0);
and_flag = uint8(0);
cond = uint8(zeros(10,1));
or_cond = uint8(zeros(5,1));
and_cond = uint8(zeros(5,1));

%% Check OR conditions
for ii = 1:length(in.s.term.or.var)
    if in.s.term.or.var(ii) == uint8(0) % break out of loop at first unused index
        break;
    else
        or_cond(ii) = check_cond( in.s.term.or.type(ii), ...
                                  in.s.term.or.var(ii), ...
                                  in.s.term.or.value(ii), ...
                                  in.s.term.or.direction(ii), ...
                                  calcs_curr, calcs_prev );
        or_flag = or_flag + or_cond(ii);
    end
end % OR conditions

%% Check AND conditions
if in.s.term.and.var(1) ~= uint8(0)
    and_flag = uint8(1);
    for ii = 1:length(in.s.term.and.var)
        if in.s.term.and.var(ii) == uint8(0) % break out of loop at first unused index
            break;
        else
            and_cond(ii) = check_cond( in.s.term.and.type(ii), ...
                                       in.s.term.and.var(ii), ...
                                       in.s.term.and.value(ii), ...
                                       in.s.term.and.direction(ii), ...
                                       calcs_curr, calcs_prev );
            and_flag = and_flag * and_cond(ii);
        end
    end % AND conditions
end


%% Compute outputs
flag = logical(or_flag + and_flag);
cond = [or_cond; and_cond];


end % termination



function [cond] = check_cond(type,var,value,direction,calcs_curr,calcs_prev)
%% Check if individual condition is satisfied

%% Initialize
cond = uint8(0);

%% Check condition
if type == 0 
    %% Crossing/equality condition
    switch var % Determine which side end trigger is during previous and current time step
        case 1 % altitude
            sign1 = sign(calcs_prev.alt - value);
            sign2 = sign(calcs_curr.alt - value);
        case 2 % mass
            sign1 = sign(calcs_prev.mass - value);
            sign2 = sign(calcs_curr.mass - value);
        case 3 % Mach
            sign1 = sign(calcs_prev.mach - value);
            sign2 = sign(calcs_curr.mach - value);
        case 4 % Planet-fixed velocity magnitude
            sign1 = sign(calcs_prev.vel_pp_mag - value);
            sign2 = sign(calcs_curr.vel_pp_mag - value);
        case 5 % Scaled Energy
            sign1 = sign(calcs_prev.scaled_energy - value);
            sign2 = sign(calcs_curr.scaled_energy - value);
            
        otherwise
            sign1 = NaN;
            sign2 = NaN;
    end % switch

    if (sign1*sign2) <= 0 % Value has been crossed, determine direction and set condition
        if (direction == int8(1)) && (sign2 == 1)
            cond = uint8(1);
        elseif (direction == int8(-1)) && (sign2 == -1)
            cond = uint8(1);
        elseif direction == int8(0)
            cond = uint8(1);
        end
    end 

else 
    %% Inequality condition
    switch var
        case 1 % altitude
            sign1 = sign(calcs_curr.alt - value);
        case 2 % mass
            sign1 = sign(calcs_curr.mass - value);
        case 3 % Mach
            sign1 = sign(calcs_curr.mach - value);
        case 4 % Planet-fixed velocity magnitude
            sign1 = sign(calcs_curr.vel_pp_mag - value);
        otherwise
            sign1 = NaN;
    end % switch                

    if (sign1*type) > 0 % Set condition
        cond = uint8(1);
    end
    
end

end % termination()