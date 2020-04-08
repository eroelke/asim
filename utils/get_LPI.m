% get_LPI.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   rot_angle - double(1), rad, angle rotated through since epoch
%   
% Outputs:
%   L_PI - double(3,3), nd, rotation matrix from the inertial frame to the
%                           planet relative frame
%
% Major Revision History:
%   *20 SEP 2011, B. Steinfeldt, original creation

function [L_PI] = get_LPI(rot_angle)
%#codegen

L_PI = [ cos(rot_angle) sin(rot_angle) 0; ...
    -sin(rot_angle) cos(rot_angle) 0; ...
    0              0 1];


end % get_LPI