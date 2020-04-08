% get_unit_NED.m 
%   Return unit vectors in north, east, and down directions in the PCPF 
%   frame based on current latitude and longitude
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   lat - double(1), rad, latitude
%   lon - double(1), rad, longitude
%     
% Outputs:
%   uN - double(3), nd, unit vector pointing north, PCPF
%   uE - double(3), nd, unit vector pointing east, PCPF
%   uD - double(3), nd, unit vector pointing down, PCPF
%
% Major Revision History:
%   *Created 5 MAY 2012, Z.R.Putnam

function [uN, uE, uD] = get_unit_NED( lat, lon )
%#codegen

    % Compute first in xyz coordinates (z: north pole, x-z plane: contains
    % r, y: completes right-handed set)
    uDxyz = [-cos(lat);0;-sin(lat)];
    uNxyz = [-sin(lat);0;cos(lat)];
    uExyz = [0;1;0];

    % Rotate by longitude to change to PCPF frame
    L3 = [cos(lon), -sin(lon), 0;
          sin(lon),  cos(lon), 0;
                 0,         0, 1];
    uN = L3*uNxyz;
    uE = L3*uExyz;
    uD = L3*uDxyz;

end % get_unit_NED()