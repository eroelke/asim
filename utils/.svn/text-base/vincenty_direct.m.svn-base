%% Vincenty's Direct Solution of Geodesics on an Ellipsoid
%  Computes geodetic latitude and longitude of a point specified by a
%  geodetic distance and azimuth from a given point
%
% Inputs:
%   phi1: point 1 geodetic latitude, rad
%   lambda1: point 1 longitude, rad
%   s: geodetic distance between point 1 and point 2, length
%   alpha1: azimuth of geodesic, east of north
% Parameters:
%   a: ellipsoid semimajor axis, length
%   b: ellipsoid semiminor axis, length
% Output:
%   phi2: point 2 geodetic latitude, rad
%   lambda2: point 2 longitude, rad
%   alpha2: azimuth of geodesic, east of north
%   e: function status
%     0: Nominal termination
%     1: Algorithm will not converge
%     2: Maximum iterations exceeded
% Functions called:
%   tan, cos, sin, atan, atan2
%
% Reference:
%  Vincenty, T., "Direct and Inverse Solutions of Geodesics on the
%  Ellipsoid with Application of Nested Equations," Survey Review, Vol. 23,
%  No. 176, April 1975.
%
% Revision history:
%   *Created 16 June 2011 by Z.R. Putnam, Space Systems Design Lab, Georgia
%     Institute of Technology
%
function [phi2, lambda2, alpha2, e] = vincenty_direct(phi1, lambda1, ...
    s, alpha1, a, b)
%#codegen

tol = 1e-8; % tolerance on lambda convergence
iter_max = 10; % maximum allowable iterations

f = (a-b)/a; % ellipsoid flattening, nd

cos_alpha1 = cos(alpha1);
sin_alpha1 = sin(alpha1);

tan_U1 = (1-f)*tan(phi1); % reduced latitude, rad
cos_U1 = cos(atan(tan_U1));
sin_U1 = tan_U1*cos_U1;
sigma1 = atan2( tan_U1,cos_alpha1 );

sin_alpha = cos_U1*sin_alpha1;
cos_alpha_sq = (1-sin_alpha)*(1+sin_alpha);
sin_alpha_sq = 1-cos_alpha_sq;
u_sq = cos_alpha_sq * (a^2-b^2)/b^2;

A = 1+ u_sq/16384 * (4096 + u_sq*(-768 + u_sq*(320 - 175*u_sq))); % (3)
B = u_sq/1024 * (256 + u_sq*(-128 + u_sq*(74 - 47*u_sq))); % (4)

sbA = s/(b*A);
sigma = sbA;

dsigma = 2*tol;
e = 0; % assume nominal termination
while (dsigma > tol)
    cos_sigma = cos(sigma);
    sin_sigma = sin(sigma);
    
    cos_2sigma_m = cos(2*sigma1 + sigma);
    
    delta_sigma = B*sin_sigma*(cos_2sigma_m + 0.25*B* ...
        (cos_sigma*(-1 + 2*cos_2sigma_m^2) - (1/6)*B*cos_2sigma_m * ...
        (-3 + 4*sin_sigma^2)*(-3 + 4*cos_2sigma_m^2))); % (6)   
    
    sigma = sbA + delta_sigma;   

    iter = iter+1; % increment iterations
    if iter >= iter_max
        e = 2; % maximum allowable iterations exceeded
        break;
    end
    
end % while (dsigma > tol)

cos_sigma = cos(sigma);
sin_sigma = sin(sigma);
cos_2sigma_m = cos(2*sigma1 + sigma);

phi2 = atan2( sin_U1*cos_sigma + cos_U1*sin_sigma*cos_alpha1, ...
    (1-f)*sqrt(sin_alpha_sq + (sin_U1*sin_alpha - cos_U1*cos_sigma*cos_alpha1)^2) ); % (8)
lambda = atan2( sin_sigma*sin_alpha1, ...
    cos_U1*cos_sigma - sin_U1*sin_sigma*cos_alpha1 ); % (9)

C = (f/16) * cos_alpha_sq * ( 4 + f*(4 - 3*cos_alpha_sq) ); % (10)
L = lambda - (1-C)*f*sin_alpha * ...
    ( sigma + C*sin_sigma*...
    (cos_2sigma_m + C*cos_sigma * (-1 + 2*cos_2sigma_m^2)) ); % (11)

lambda2 = L + lambda1;

alpha2 = atan2( sin_alpha, ...
    -sin_U1*sin_sigma + cos_U1*cos_sigma*cos_alpha1); % (12)

end % vincenty_direct()