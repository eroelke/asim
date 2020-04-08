% cvdma_in_in.m 
%   Create input structure for continuously-variable drag-modulation
%   aerocapture (cvdma)
%
% Aeroassist Simulation source code (ASIM)
%
% Developed by:
%   Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
%
% Inputs:
%   h_atm - atmospheric interface altitude (m)
%   planet - string or index for planet model
%     1: Venus
%     2: Earth
%     3: Mars
%     4: Jupiter
%     5: Saturn
%     6: Titan
%     7: Uranus
%     8: Neptune
%   atm_mode - atmospheric mode
%       1 - exponential atm
%       2 - table look up, no winds
%       3 - table look up, winds
%   
% Outputs:
%   in - struct, multi, input structure with parameters for SEJ guidance
%
% Major Revision History:
%   *Created Feb 2019, E. Roelke
% 
function in = cvdma_in(h_atm,planet,atm_mode)
%#codegen

%% Define input structure
in = define_in;

%% Planetary data
if ischar(planet) == true
    models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
    planet = find(strcmpi(planet,models) == 1);
end
if atm_mode > 3 || atm_mode < 1
    fprintf('Warning: Incorrect Atmospheric Mode. Defaulting to Exponential Model\n');
    atm_mode = 1;
end
in.p = get_planetary_data(planet, atm_mode); % Load planet model, GRAM atmosphere + winds

%% Simulation data

% Termination conditions for aerocapture
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination variable (altitude)
in.s.term.or.value(1) = h_atm; % double, m, terminate at atmospheric exit
in.s.term.or.direction(1) = int8(1); % int8, nd, positive crossing (increasing value)

in.s.term.or.type(2) = int8(-1); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = 0.0; % double, m, terminate at atmospheric exit

% Default (planet-relative) entry state
in.s.traj.gamma_pp = -5*pi/180;
in.s.traj.alt = h_atm; % initial altitude, m
in.s.traj.lat = 0*pi/180; % initial latitude, rad
in.s.traj.lon = 0*pi/180; % initial longitude, rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in.s.traj.vel_pp_mag = 10e3;

% Convert scalars to inertial position and velocity vectors
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

% Default simulation settings
in.s.traj.rate = 10; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 750; % double, s, maximum allowable time
in.s.data_rate = 1;     % data storage rate (Hz)
in.v.gnc.g.p.rate = 0.1;  % guidance rate (Hz)
in.v.gnc.c.p.rate = 0.1;  % control rate (Hz)
in.v.gnc.n.p.rate = in.s.traj.rate;     % nav rate (Hz)

%% Default Vehicle data
% Aerodynamics
in.v.aero.mode = uint8(1); % constant aero coeffs
in.v.aero.aoa = 0; % double, rad, angle-of-attack
in.v.aero.cl = 0; % double, nd, lift coefficient
in.v.aero.cd = 1; % double, nd, drag coefficient
in.v.aero.area_ref = 12.56; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = 1; % double, m, nose radius
in.v.mp.m_ini = 1000; % double, kg, vehicle mass

% Guidance
in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.dm_mode = uint8(5);   % cvdma

% cvdma settings
in.v.gnc.g.p_cvdma.mode = uint8(1); %NPC
% in.v.gnc.g.p_cvdma.decel = 0;   %no ref traj
in.v.gnc.g.p_cvdma.iter_max = uint8(2); % nd, max NPC internal iterations
in.v.gnc.g.p_cvdma.A_sens_atm = 0.25; % m/s^2, sensible decel
in.v.gnc.g.p_cvdma.ha_tgt = 1000e3; % m, target apoapsis alt
in.v.gnc.g.p_cvdma.ha_tol = 10e3; % m, apoapsis alt tolerance
in.v.gnc.g.p_cvdma.cd_ini = in.v.aero.cd;   %initial drag coeff
in.v.gnc.g.p_cvdma.rn = in.v.aero.nose_radius;  %nose radius
in.v.gnc.g.p_cvdma.delta_c_ini = 60 * pi/180;    % initial cone half-angle
in.v.gnc.g.p_cvdma.rc_ini = 2;
in.v.gnc.g.p_cvdma.area_rate = pi*(0.25^2); %max/min change rate is 0.25 m backshell
in.v.gnc.g.p_cvdma.rc_bounds = [0.5 2.5];   %max/min backshell radii
in.v.gnc.g.p_cvdma.mass = in.v.mp.m_ini;    %vehicle mass
in.v.gnc.g.p_cvdma.c = in.v.gnc.g.p_cvdma.rc_ini / tan(in.v.gnc.g.p_cvdma.delta_c_ini);
in.v.gnc.g.p_cvdma.t_max = 500;
in.v.gnc.g.p_cvdma.traj_rate = in.s.traj.rate;
in.v.gnc.g.p_cvdma.cd_bounds = [0.25 2]';

% guidance planet model
in.v.gnc.g.p_cvdma.planet.alt_max = h_atm;
in.v.gnc.g.p_cvdma.planet.alt_min = 0.0;
in.v.gnc.g.p_cvdma.planet.r_p = in.p.r_p;
in.v.gnc.g.p_cvdma.planet.r_e = in.p.r_e;
in.v.gnc.g.p_cvdma.planet.mu = in.p.mu; % m^3/s^2
in.v.gnc.g.p_cvdma.planet.j2 = in.p.j2; % nd
in.v.gnc.g.p_cvdma.planet.npole = [0;0;1]; % nd
in.v.gnc.g.p_cvdma.planet.omega = in.p.omega; % rad/s
in.v.gnc.g.p_cvdma.planet.atm_table = in.p.atm.table(:,1:7);    %nominal atmosphere

end % cvdma_in
