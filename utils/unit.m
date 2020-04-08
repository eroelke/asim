% unit.m 
%   Compute unit vector
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
%   
% Outputs:
%   u - double(3), nd, unit vector along input, a
%
% Major Revision History:
%   *Created 17 JUL 2012, Z.R.Putnam

function u = unit(a)

    a_mag = sqrt( dot(a,a) );
    u = a ./ a_mag;

end % unit