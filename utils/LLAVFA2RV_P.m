% LLAVFA2RV_P.m
%   Transformation between planet-relative state scalars and vectors
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   lat - double(1), rad, geodetic latitude
%   lon - double(1), rad, longitude
%   alt - double(1), m, altitude above the ellipsoid
%   gamma_pp - double(1), rad, flight-path angle (negative downward)
%   az - double(1), rad, azimuth angle
%   vel_pp_mag - double(1), m/s, planet-centric velocity magnitude (speed)
%   r_e - double(1), m, equatorial radius
%   r_p - double(1), m, polar radius
%   
% Outputs:
%   pos_pp - double(3), m, position in planet-centric frame (pcpf)
%   vel_pp - double(3), m/s, velocity in planet-centric frame (pcpf)
%
% Major Revision History:
%   *9 OCT 2011, B. Steinfeldt, original creation
%   *14 OCT 2011, B. Steinfeldt, modularity added
%   *13 AUG 2012, Z.R.Putnam, changed from geocentric to geodetic
%   (geocentric still available when r_e == r_p)

function [pos_pp, vel_pp] = LLAVFA2RV_P(lat, lon, alt, gamma_pp, az, vel_pp_mag, r_e, r_p) 
%codegen

% Position vector in PCPF
pos_pp = get_r( lat, lon, alt, r_e, r_p );


% Planet-relative velocity vector in PCPF
[uN,uE,uD] = get_unit_NED(lat,lon); % nd, unit vectors along local NED axes
vD = -vel_pp_mag*sin(gamma_pp); 
vH = vel_pp_mag*cos(gamma_pp);
vE = vH*sin(az);
vN = vH*cos(az);
vel_pp = vN*uN + vE*uE + vD*uD;


end % LLAVFA2RV_P()