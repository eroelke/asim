%% PRED GUID Aerocapture Guidance Source Code
% Calculate altitude above reference ellipsoid
function alt = guid_sej_e_alt( R, p )
%#codegen

% Initialize variables
% alt     = double(0);
% temp_s1 = double(0);
% R_mag   = double(0);
% U_R     = double([0;0;0]);

% temp_s1 = (1 - p.flatness)^2;
% % temp_s1 = temp_s1^2;
% p.alt_coeff1 = temp_s1 - 1; % nd
% p.alt_coeff2 = p.p_r*(1 - p.flatness); % ft

temp_s1 = (1 - p.flatness)^2;
alt_coeff1 = temp_s1 - 1; % nd
alt_coeff2 = p.p_r*(1 - p.flatness); % ft


% Compute position vector magnitude and unit vector
R_mag = norm(R,2); % ft
u_R = R/R_mag; % nd

% Compute Earth oblatness correction term
temp_s1 = alt_coeff2 / ...
    sqrt( (1 - dot(u_R,p.npole)^2 )*alt_coeff1 + 1); % ft

% Compute altitude
alt = R_mag - temp_s1; % ft

end % GDE_PREDGUID_pred_calc_alt
