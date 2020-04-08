%% Uranus aerocapture analytical calculations
% JPL Strategic Study
% Sam Albert
% 7.21.19

clear; clc;

in0 = uranus_in;

%% Shared constants
e = exp(1);
g_e = in0.c.earth_g;
k = in0.p.k;
H = in0.p.atm.scale_height;
rho0 = in0.p.atm.dens_ref;

%% Input trajectory properties
v_atm = 25.8915*1e3; % m/s, not including effect of planet rotation

%%%%% With Drag Skirt %%%%%%%%
fpa = -0.142845365914719;
rn = 1/2 * 6;
BC = 30;

% %%%%% Without Drag Skirt %%%%%
% fpa = -0.151160618164708;
% rn = 1/2 * 2;
% BC = 120;


%% Calculations
n_max = v_atm^2 * sin(fpa) / (2 * H * g_e * e);

qs_max = k * (1 / rn)^(1/2) * (BC * sin(fpa) / (3 * H))^(1/2) * ...
    0.6055 * v_atm^3; % result is W/m^2
qa_max = qs_max / 2; % approximate max q on ADEPT as 1/2 max stag p rate
Qs = k * v_atm^2 * (BC/rho0 * (pi * H / (rn * sin(fpa))))^(1/2); % heat load
Qa = Qs / 2; % assume same approximation here as for heat rate





