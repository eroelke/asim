% define_aero_dc_act.m 
%   Define Dream Chaser actuator deflection data structure
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
%   dc_act - struct(1), multi, Dream Chaser actuator deflection data structure
%
% Major Revision History:
%   *Created 21 FEB 2012, Z.R. Putnam

function dc_act = define_aero_dc_act()
%#codegen

dc_act = struct(... % deg, actuator deflections
    'DLE', double(0), ... 
    'DRE', double(0), ...
    'DUL', double(0), ...
    'DUR', double(0), ...
    'DLL', double(0), ...
    'DLR', double(0), ...
    'DR', double(0), ...
    'GEAR', double(0) );

end % define_aero_dc_act()