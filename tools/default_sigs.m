% default_sigs.m
%   Create uncertainty structure for monte carlo sims
% 
% Inputs: none
% Outputs:
%   sigs: 1sigma uncertainty structure
%       sigs.m: 1sigma on mass (kg)
%       sigs.cl: 1sigma on lift coeff
%       sigs.cd: 1sigma on drag coeff (for all drag coeffs if changing)
%       sigs.aoa: 1sigma on all angle of attacks (rad)
%       sigs.vmag_atm: 1sigma on entry velocity (m/s)
%       sigs.hmag_atm: 1sigma on atmospheric interface (m)
%       sigs.efpa: 1sigma on entry flight path angle (rad)
%       sigs.lat0: 1sigma on entry latitude (rad)
%       sigs.lon0: 1sigma on entry longitude (rad)
%       sigs.az0: 1sigma on entry azimuth (rad)
%       sigs.imu_bias: 1sigma bais on earth G measurements
%       sigs.imu_noise: 1sigma bias on states [x,v,a];
% 

function sigs = default_sigs()

sigs.m = 0.5 / 100 / 3; %.5 3sigma
sigs.cl = 0;
sigs.cd = 5 / 100 / 3;  %5% 3sigma - see set_jettison.m
sigs.aoa = 0;
sigs.vmag_atm = 10 / 3;  %10 m/s 3sigma
sigs.hmag_atm = 150;     %m, 1sigma
sigs.efpa = (0.018 * pi/180) / 3; % correlated MSL
sigs.lat0 = 5e-4;
sigs.lon0 = 5e-4;
sigs.az0 = 5e-4;
sigs.imu_bias = zeros(9,1); % per axis
sigs.imu_noise = zeros(9,1);    % per axis

end