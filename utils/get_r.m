% get_r.m 
%   Returns current position vector based on lat/lon/alt
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   lat - double(1), rad, latitude (geodetic)
%   lon - double(1), rad, longitude
%   alt - double(1), rad, altitude (geodetic)
%   r_e - double(1), m, planet equatorial radius
%   r_p - double(1), m, planet polar radius
%     
% Outputs:
%   r - double(3), m, position vector in PCPF
%
% Major Revision History:
%   *Created 5 MAY 2012, Z.R.Putnam

function r = get_r( lat, lon, alt, r_e, r_p )
%#codegen

    alat = abs(lat);
    
    if alat > 1.5707
        x0 = tan(alat-1.5707)*r_p;
        z0 = r_p;
    else
        A = tan(alat)*(r_p/r_e)^2;
        x0 = 1/sqrt( 1/r_e^2 + (A/r_p)^2 );
        z0 = A*x0;
    end
    
    rxz = [x0+alt*cos(alat); 0; sign(lat)*(z0+alt*sin(alat))];
    L3 = [cos(lon), -sin(lon), 0;
          sin(lon),  cos(lon), 0;
                 0,         0, 1];
    r = L3*rxz;

end % get_r