% get_dcaero_clcd.m 
%   Returns Dream Chaser lift and drag coefficients based on mach and aoa
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   mach - double(1), nd, current Mach number
%   aoa - double(1), rad, current angle-of-attack
%   basic - struct(1), md, Dream Chaser basic aerodynamic data structure
%   
% Outputs:
%   cl - double(1), nd, vehicle lift coefficient at angle of attack
%   cd - double(1), nd, vehicle drag coefficient at angle of attack
%
% Major Revision History:
%   *Created FEB 2012, Z.R.Putnam

function [cl,cd] = get_dcaero_clcd(mach, aoa, basic)
%#codegen

% Get aoa-polynomial coefficients
cl_p = lin_interp(basic.mach,basic.Cl,mach);
cd_p = lin_interp(basic.mach,basic.Cd,mach);

% Angle-of-attack in degrees
aoa_d = aoa*180/pi;

% Compute lift and drag coefficients
cl = cl_p(1) + cl_p(2)*aoa_d + cl_p(3)*aoa_d^2 + cl_p(4)*aoa_d^3;
cd = cd_p(1) + cd_p(2)*aoa_d + cd_p(3)*aoa_d^2 + cd_p(4)*aoa_d^3;

end % get_dcaero_clcd