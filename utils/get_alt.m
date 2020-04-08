% get_alt.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   pos_pp - double(3), m, position vector in planet-centric frame
%   re - double(1), m, equatorial radius
%   rp - double(1), m, polar radius
%   lat - double(1), rad, latitude
%     
% Outputs:
%   alt - double(1), m, altitude above planetary ellipsoid
%
% Major Revision History:
%   *Created 2012, J. Kelly

function alt = get_alt(pos_pp, re, rp, lat) 
%#codegen

f = (re-rp)/re; % flattening
e2 = 1-(1-f)^2; % ellipticity
s = norm([pos_pp(1) pos_pp(2)]);
N = re / sqrt(1-e2*sin(lat)^2); % radius of curvature in the vertical prime
alt = s*cos(lat) + (pos_pp(3) + e2*N*sin(lat))*sin(lat) - N;

end % get_alt()