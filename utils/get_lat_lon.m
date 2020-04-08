% get_lat_lon.m 
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   pos_pp - double(3), m, position vector in planet-centric frame
%   re - double(1), m, equatorial radius
%   rp - double(1), m, polar radius
%   
% Outputs:
%   lat - double(1), rad, planet-detic latitude
%   lon - double(1), rad, longitude
%
% Major Revision History:
%   *20 SEP 2011, B. Steinfeldt, original creation
%   *05 MAR 2012, J. Kelly, converted to planet-detic latitude instead of
%       planet-centric
%   *9 MAY 2012, Z.R. Putnam, removed alt as an output (unused), change
%       long to lon

function [lat,lon] = get_lat_lon(pos_pp,re,rp)
%#codegen

f = (re-rp)/re; % flattening
e2 = 1-(1-f)^2; % ellipticity
pos_pp_mag = norm(pos_pp);
s = norm([pos_pp(1) pos_pp(2)]);

% Calculate initial guesses for reduced latitude (latr) and planet-detic latitude (latd)
latr = atan(pos_pp(3) / ( (1-f)*s ));
latd = atan( ( pos_pp(3) + e2*(1-f)*re*sin(latr)^3/(1-e2) ) / ( s - e2*re*cos(latr)^3 ) );

% Recalculate reduced latitude based on planet-detic latitude
latr2 = atan( (1-f)*sin(latd) / cos(latd) );
diff = latr - latr2;

% Iterate until reduced latitude converges
while diff > 10^-7
    latr = latr2;
    latd = atan( ( pos_pp(3) + e2*(1-f)*re*sin(latr)^3/(1-e2) ) / ( s - e2*re*cos(latr)^3 ) );
    latr2 = atan( (1-f)*sin(latd) / cos(latd) );
    diff = latr - latr2;
end
lat = latd;

% Calculate longitude
lon = asin(pos_pp(2)/sqrt(pos_pp(1)^2+pos_pp(2)^2));

% Quadrant check on longitude
if pos_pp(1) < 0
    lon = pi - lon;
elseif pos_pp(1) >0
    if pos_pp(2) < 0
        lon = 2*pi + lon;
    end
end

end % get_lat_lon