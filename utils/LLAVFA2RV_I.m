% LLAVFA2RV_I.m
%   Transformation between state vector components and inertial vectors
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
%   omega - double(3), rad/s, angular velocity vector of planet
%   t0 - double(1), s, epoch time
%   t - double(1), s, current time
%   r_e - double(1), m, equatorial radius
%   r_p - double(1), m, polar radius
%   
% Outputs:
%   pos_ii - double(3), m, position in inertial frame
%   vel_ii- double(3), m/s, velocity in inertial frame
%
% Major Revision History:
%   *9 OCT 2011, B. Steinfeldt, original creation
%   *14 OCT 2011, B. Steinfeldt, modularity added

function [pos_ii, vel_ii] = LLAVFA2RV_I(lat, lon, alt, ...
                                        gamma_pp, az, vel_pp_mag, ...
                                        omega, t0, t, r_e, r_p) 
%codegen

%% Get planet-centric state vector
[pos_pp, vel_pp] = LLAVFA2RV_P(lat, lon, alt, gamma_pp, az, vel_pp_mag, r_e, r_p);

%% Convert planet-centric vectors to inertial vectors
pos_ii = PCPF2PCI(pos_pp, omega, t0, t);
vel_ii = PCPF2PCI(vel_pp, omega, t0, t) + cross(omega,pos_ii);

end  % LLAVFA2RV_I()