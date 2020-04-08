% PCPF2PCI.m 
%   Vector transformation between planet-centric and inertial frame
%   
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   vec_pp - double(3), misc, vector in planet-centric frame
%   t0 - double(1), s, epoch time
%   t - double(1), s, time since epoch
%   omega - double(3), rad/s, planet's rotation rate
%   
% Outputs:
%   vec_ii - double(3), misc, vector in inertial frame
%
% Major Revision History:
%   *9 OCT 2011, B. Steinfeldt, original creation
%   *14 OCT 2011, B. Steinfeldt, reduced to solely vector transformations


function [vec_ii] = PCPF2PCI(vec_pp, omega, t0, t) %codegen

% Inertial to planet relative transformation
rot_angle = norm(omega)*(t-t0); % Rotation angle [rad]
[L_PI] = get_LPI(rot_angle);

% Perform the vector transformation
vec_ii = L_PI'*vec_pp; % Position vector planet/planet [m]

return