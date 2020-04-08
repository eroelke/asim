% skew.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   vec - double(3), nd, vector from which the skew symmetric matrix is mad 
%     
% Outputs:
%   n - double(3,3), nd, skew symmetric matrix
%
% Major Revision History:
%   *Created 26 OCT 2008, M. Grant

function n = skew(vec) 
%#codegen

n = [     0  -vec(3)  vec(2); ...
      vec(3)      0  -vec(1); ...
     -vec(2)  vec(1)      0];

end % skew()

