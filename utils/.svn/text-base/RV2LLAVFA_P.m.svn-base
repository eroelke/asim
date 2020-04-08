% RV2LLAVFA_P.m
%   Transformation between planet-relative state vectors and scalars
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   pos_pp - double(3), m, position in PCPF frame
%   vel_pp - double(3), m/s, planet-relative velocity in PCPF frame
%   r_e - double(1), m, equatorial radius
%   r_p - double(1), m, polar radius
%   
% Outputs:
%   lat - double(1), rad, geodetic latitude
%   lon - double(1), rad, longitude
%   alt - double(1), m, geodetic altitude
%   gamma_pp - double(1), rad, planet-relative geodetic flight-path angle (negative downward)
%   az - double(1), rad, planet-relative geodetic azimuth angle
%   vel_pp_mag - double(1), m/s, planet-relative velocity magnitude
%
% Major Revision History:
%   *9 OCT 2011, B. Steinfeldt, original creation
%   *14 OCT 2011, B. Steinfeldt, modularity added
%   *13 AUG 2012, Z.R.Putnam, changed from geocentric to geodetic
%   (geocentric still available when r_e == r_p)

function [lat, lon, alt, gamma_pp, az, vel_pp_mag] = RV2LLAVFA_P(pos_pp, vel_pp, r_e, r_p) 
%codegen

% Position scalars
[lat,lon] = get_lat_lon(pos_pp, r_e, r_p); % rad, latitude and longitude
alt = get_alt( pos_pp, r_e, r_p, lat ); % m, altitude


% Velocity scalars
vel_pp_mag = norm(vel_pp,2); % m/s, planet-relative velocity magnitude
[uN,uE,uD] = get_unit_NED(lat,lon); % nd, unit vectors along local NED axes
vN = dot(vel_pp,uN); % m/s, planet-relative velocity magnitude along uN
vE = dot(vel_pp,uE); % m/s, planet-relative velocity magnitude along uE
vD = dot(vel_pp,uD); % m/s, planet-relative velocity magnitude along uD
vH = norm(vE*uE + vN*uN); % m/s, planet-relative velocity magnitude in local horizontal plane
az = atan2(vE,vN); % rad, planet-relative azimuth
gamma_pp = atan2(-vD,vH); % rad, planet-relative flight-path


end % RV2LLAVFA_P()