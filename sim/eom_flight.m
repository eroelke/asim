% eom_flight.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   t - double(1), s, current time
%   y - double(n), misc, current state vector
%                        Currently: [x, y, z, u, v, w, m, downrange,
%                                          crossrange, heat load]'
%   in - data structure, misc, input data structure
%   azi_pp0 - double(1), rad, initial aziumth angle in planet-centric frame
%   veh - struct(1), multi, vehicle data structure
%       s.bank - double(1), rad, current bank angle
%       s.ssa - double(1), rad, sideslip angle
%       s.area_ref - double(1), m^2, current aerodynamic reference area
%       s.compute_trim - uint8(1), nd, trim angle-of-attack calculation flag
%   ctrl - struct(1), multi, control data structure
%   
% Outputs:
%   ydot - double(n), misc, state derivative vector
%
% Major Revision History:
%   *Created 2008, M. Grant
%   *Revised 20 SEP 2011, B. Steinfeldt, Added modularity to code
%   *8 FEB 2012, Z.R.Putnam, added angle-of-attack

function [ydot] = eom_flight( t, y, in, azi_pp0, veh, ctrl )
%#codegen

%% Inputs and Initial Calculations

calcs = traj_calcs( in, t, y, veh, ctrl );


%% Derivative Calculations

% Change in azimuth angle from initial state
d_azi = calcs.azi_pp - azi_pp0;

% Derivatives
ydot1 = zeros(1,10);
ydot1(1:3) = calcs.vel_ii;
ydot1(4:6) = calcs.force_ii/calcs.mass;
ydot1(7) = -norm(calcs.thrust_ii)/(in.c.earth_g*in.v.prop.main.isp);
ydot1(8) = calcs.vel_pp_mag*cos(calcs.gamma_pp)*cos(d_azi);
ydot1(9) = calcs.vel_pp_mag*cos(calcs.gamma_pp)*sin(d_azi);
ydot1(10) = calcs.heat_rate;
ydot = ydot1';

end % eom_flight()
