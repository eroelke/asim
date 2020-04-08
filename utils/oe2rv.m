% AE 6353
% Z. R. Putnam
% 16 SEP 2010

% Homework #2
% Problem 6 A1

% Input - oe
%  oe(1) = semi-major axis (a)
%  oe(2) = eccentricity magnitude (e)
%  oe(3) = inclination (i)
%  oe(4) = argument of periapsis (small omega)
%  oe(5) = longitude of ascending node (capital omega)
%  oe(6) = true anomaly angle (nu)

function [rv]=oe2rv(oe,mu)

sma = oe(1);
ecc = oe(2);
inc = oe(3);
apr = oe(4);
lan = oe(5);
taa = oe(6);

slr = sma * (1 - ecc^2);
r = slr / (1 + ecc*cos(taa));

r_PQW = [r*cos(taa); r*sin(taa); 0];
v_PQW = [-sqrt(mu/slr)*sin(taa); sqrt(mu/slr)*(ecc + cos(taa)); 0];

% form rotation matrix
R3_m_lan = [ cos(-lan) sin(-lan) 0;
            -sin(-lan) cos(-lan) 0;
             0         0         1];
R1_m_inc = [1 0          0;
            0 cos(-inc)  sin(-inc);
            0 -sin(-inc) cos(-inc)];
R3_m_apr = [ cos(-apr) sin(-apr) 0;
            -sin(-apr) cos(-apr) 0;
             0         0         1];

R = R3_m_lan * R1_m_inc * R3_m_apr;
         
r_IJK = R * r_PQW;
v_IJK = R * v_PQW;

rv = [r_IJK; v_IJK];

end