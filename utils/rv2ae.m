% rv2ae.m
%   convert PCI r,v into semi-major axis and eccentricity
%   prevents any errors from acos and whatnot for extreme exit conditions
% Entry systems Design Lab
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Aug 2019
% 
function [ae] = rv2ae(rv,mu)
%#codegen
% extract r,v, from input
r = rv(1:3);
v = rv(4:6);
    
rmag = norm(r);     %magintude of position vector (m)
vmag = norm(v);     %magnitue of velocity vector (m/s)

a = rmag/(2 - (rmag*vmag^2)/mu);        %semimajor axis (vis viva equation)
h = cross(r,v);         %angular momentum
% hmag = norm(h);
e = norm(cross(v,h)/mu - r/rmag); %eccentricity
ae = [a e];
end