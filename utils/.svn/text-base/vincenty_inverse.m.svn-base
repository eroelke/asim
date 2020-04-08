%% Vincenty's Inverse Solution of Geodesics on an Ellipsoid
%  Computes geodetic distance between two points specified by their
%  geodetic latitudes and longitudes
%
% Inputs:
%   phi1: point 1 geodetic latitude, rad
%   lambda1: point 1 longitude, rad
%   phi2: point 2 geodetic latitude, rad
%   lambda2: point 2 longitude, rad
% Parameters:
%   a: ellipsoid semimajor axis, length
%   b: ellipsoid semiminor axis, length
% Output:
%   s: geodetic distance between point 1 and point 2, length
%   alpha1: azimuth of geodesic, east of north
%   alpha2: azimuth of geodesic, east of north
%   e: function status
%     0: Nominal termination
%     1: Algorithm will not converge
%     2: Maximum iterations exceeded
% Functions called:
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
function [s, alpha1, alpha2, e] = vincenty_inverse(phi1, lambda1, ...
    phi2, lambda2, a, b) 
%#codegen

% Initalize parameters
lambda = 0;
sin_lambda = 0;
cos_lambda = 0;
sin_sigma_sq = 0;
sin_sigma = 0;
cos_sigma = 0;
sigma = 0;
sin_alpha = 0;
cos_alpha_sq = 0;
cos_2sigma_m = 0;
C = 0;

tol = 1e-8; % tolerance on lambda convergence
iter_max = 10; % maximum allowable iterations

f = (a-b)/a; % ellipsoid flattening, nd

tan_U1 = (1-f)*tan(phi1); % reduced latitude, rad
cos_U1 = cos(atan(tan_U1));
sin_U1 = tan_U1*cos_U1;

tan_U2 = (1-f)*tan(phi2); % reduced latitude, rad
cos_U2 = cos(atan(tan_U2));
sin_U2 = tan_U2*cos_U2;

L = lambda2-lambda1; % difference in longitude, rad

e = 0; % assume nominal termination


dlambda = 2*tol;
iter = 0;
while (dlambda > tol)

    lambda = L; % (13) initial guess
    sin_lambda = sin(lambda);
    cos_lambda = cos(lambda);   

    sin_sigma_sq = (cos_U2*sin_lambda)^2 + ...
        (cos_U1*sin_U2 - sin_U1*cos_U2*cos_lambda)^2; % (14)
    sin_sigma = sqrt(sin_sigma_sq); % added
    cos_sigma = sin_U1*sin_U2 + cos_U1*cos_U2*cos_lambda; % (15)
    sigma = atan2(sin_sigma,cos_sigma); % (16) - modified
    sin_alpha = cos_U1*cos_U2*sin_lambda / sin_sigma; % (17)
    cos_alpha_sq = 1-sin_alpha^2; % added
    cos_2sigma_m = cos_sigma - 2*sin_U1*sin_U2 / cos_alpha_sq; % (18)

    C = (f/16) * cos_alpha_sq * ( 4 + f*(4 - 3*cos_alpha_sq) ); % (10)
    L = lambda - (1-C)*f*sin_alpha * ...
        ( sigma + C*sin_sigma*...
        (cos_2sigma_m + C*cos_sigma * (-1 + 2*cos_2sigma_m^2)) ); % (11)
    
    if abs(L) > pi
        e = 1; % algorithm will not converge, points are both nearly antipodal
        break;
    end

    dlambda = abs(lambda-L);
    iter = iter+1; % increment iterations

    if iter >= iter_max
        e = 2; % maximum allowable iterations exceeded
        break;
    end
    
end % while (dlambda > tol)

lambda = L;
sin_lambda = sin(lambda);
cos_lambda = cos(lambda);

u_sq = cos_alpha_sq * (a^2-b^2)/b^2;

A = 1+ u_sq/16384 * (4096 + u_sq*(-768 + u_sq*(320 - 175*u_sq))); % (3)
B = u_sq/1024 * (256 + u_sq*(-128 + u_sq*(74 - 47*u_sq))); % (4)

delta_sigma = B*sin_sigma*(cos_2sigma_m + ...
    0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m^2) - (1/6)*B*cos_2sigma_m * ...
    (-3 + 4*sin_sigma_sq)*(-3 + 4*cos_2sigma_m^2))); % (6)

s = b*A*(sigma-delta_sigma); % (19)

alpha1 = atan2( cos_U2*sin_lambda, ...
    cos_U1*sin_U2 - sin_U1*cos_U2*cos_lambda ); % (20)

alpha2 = atan2( cos_U1*sin_lambda, ...
    -sin_U1*cos_U2 - cos_U1*sin_U2*cos_lambda ); % (21)

end % vincenty_inverse()