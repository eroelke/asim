% PCI2PCPF.m 
%   Vector transformation between inertial frame and planet-centric frame
%   
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   vec_ii - double(3), misc, vector in inertial frame
%   t0 - double(1), s, epoch time
%   t - double(1), s, time since epoch
%   omega - double(3), rad/s, planet's rotation rate
%   
% Outputs:
%   vec_pp - double(3), misc, vector in planet-centric frame
%
% Major Revision History:
%   *9 OCT 2011, B. Steinfeldt, original creation
%   *14 OCT 2011, B. Steinfeldt, reduced to solely vector transformations


function [vec_pp] = PCI2PCPF(vec_ii, omega, t0, t) %codegen

% Inertial to planet relative transformation
rot_angle = norm(omega)*(t+t0); % Rotation angle [rad]
[L_PI] = get_LPI(rot_angle);

% Perform the vector transformation
vec_pp = L_PI*vec_ii; % Position vector planet/planet [m]

end % PCI2PCPF()