% cross.m 
%   This function overrides the intrinsic Matlab cross function for speed.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   a - double(3), nd, vector
%   b - double(3), nd, vector
%   
% Outputs:
%   c - double(3), nd, cross product of a and b
%
% Major Revision History:
%   *Created NOV 2011, Z.R. Putnam

function c = cross(a,b)
%#codegen

% Compute cross product using skew-symmetric matrix
c = skew(a)*b;

end % cross()