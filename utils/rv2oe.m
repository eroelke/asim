%% rv2oe.m
% Given rPCI (m) and vPCI (in m/s),
%  compute orbital elements (oe) in perifocal frame
% 
% Aeroassist Simulation (ASIM)
%   Space Systems Design Laboratory - Georgia Tech
%   Colorado Center for Astrodynamics Research - CU Boulder
% 
% Inputs:
%   rv: 6x1 vector for position, velocity [m; m/s]
%   mu: planet gravitational parameter (m3/s2)
% 
% Outputs:
%   oe: 6x1 vector of orbital elements
%       a: semi-major axis (m)
%       e: eccentricity
%       i: inclination (rad)
%       W: longitude of ascending node (LAN) (rad)
%       w: argument of periapsis (AOP) (rad)
%       f: true anomaly at epoch (rad)
% 
% Created Sep 2015, E. Roelke
% Updated Sep 2018, E. Roelke
% 
function [oe] = rv2oe(rv,mu)
% extract r,v, from input
r = rv(1:3);
v = rv(4:6);
    
rmag = norm(r);     %magintude of position vector (m)
vmag = norm(v);     %magnitue of velocity vector (m/s)

a = rmag/(2 - (rmag*vmag^2)/mu);        %semimajor axis (vis viva equation)
h = cross(r,v);         %angular momentum
hmag = norm(h);
e = cross(v,h)/mu - r/rmag; %eccentricity vector

emag = norm(e);
K = [0 0 1];            %unit vector normal to earth's equatorial plane
n = cross(K,h/hmag);     %nodal vector
nmag = norm(n);
i = acos(h(3)/hmag);    %inclination angle in degrees

W = acos(n(1)/nmag);        %angle of ascending node (rad)
w = acos(dot(n,e)/(nmag*emag));      %argument of periapse (rad)
f = acos(dot(e,r)/(emag*rmag));      %true anomaly (rad)

%test quadrants for W
if n(2) > 0 && W > pi
    W = (2*pi) - W;
elseif n(2) < 0 && W < pi
    W = (2*pi) - W;
end

%test quadrants for w
if e(3) > 0 && w > pi
    w = (2*pi) - w;
elseif e(3) < 0 && w < pi
    w = (2*pi) - w;
end

%test quadrant for f
if dot(r,v) > 0 && f > pi
    %1st/2nd quad
    f = (2*pi) - f;
elseif dot(r,v) < 0 && f < pi
    %3rd/4th quad
    f = (2*pi) - f;
end

e = norm(e);
oe = [a e i W w f]';

end