% alt_from_guid.m
%   calculate vehicle guidance from guidance structures, for GNC systems
% 
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% 
% Written By: Evan Roelke, Aug 2019
% 
function alt = alt_from_guid(y, t, p_struct)

r_ii = y(1:3); % Inertial position

% planet params
omega = p_struct.omega;     % planet rotation rate (rad/s)
rp = p_struct.r_p;
re = p_struct.r_e;

% PCI -> PCPF
theta = norm(omega)*t; % Rotation angle [rad]
L_PI = get_LPI(theta);  % DCM for PCI->PCPF

% Inertial and planet centric velocities
r_pp = L_PI*r_ii; % Position vector planet/planet [m]
% v_pp = L_PI*(v_ii - cross(omega,r_ii)); % Velocity vector planet/planet [m/s]

[lat,~] = get_lat_lon(r_pp,re,rp);
alt = get_alt(r_pp,re,rp,lat);


end