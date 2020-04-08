%% ac_earth_MC
%   Runs an end-to-end Monte Carlo simulations using the NPC guidance
%   algorithm running at 5 hz.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
%   Jet Propulsion Laboratory
%   California Institute of Technology
%
% Major Revision History:
%   Jan. 2017 - Created, Michael Werner

%% Prepare workspace
clearvars; close all; clc;

%% Setup - Basic Inputs
% Vehicle properties
R = 0.25;           % Pre-jettison radius (m)
r = 0.1;            % Post-jettison radius (m)
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
m1 = 19.62;          % Pre-jettison mass (kg)
m2 = 14.42;         % Post-jettison mass (kg)
rn = 9/16*r;       % Nose radius (m)
delta = 60;         % Sphere-cone half angle (deg)
cd = 1.53;

% Simulation setup
nruns = 10;       % Number of Monte Carlo runs

% Trajectory propertiesin
ha_tgt = 1760e3;    % Target post-aerocapture apogee altitude (m)
hp_tgt = 180e3;     % Target post-aerocapture perigee altitude (m)

%% Setup - Low-rate NPC Input Structure
in = sej_a_earth_in; % nominal input structure

% Initial state: Trajectory properties
in.s.traj.gamma_pp = -5.04*pi/180;
in.s.traj.alt = 125e3; % initial altitude, m
in.s.traj.lat = 0*pi/180; % initial latitude, rad
in.s.traj.lon = 168.7778*pi/180; % initial longitude, rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in.s.traj.vel_pp_mag = find_vel_pp_mag(10311, in.s.traj.gamma_pp, in.p.omega, (in.p.r_e + in.s.traj.alt)); % initial velocity, m/s

% Convert scalars to inertial position and velocity vectors
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

% Vehicle properties
in.v.mp.m_ini = m1; % kg
in.v.mp.m_jettison = m1-m2;
in.v.aero.area_ref = area(1);
in.v.aero.mode = uint8(1);     % User-specified lift and drag coefficients
in.v.aero.cd = cd;
in.v.aero.cl = 0;
in.v.aero.nose_radius = rn; % Nose radius (m)
in.v.gnc.g.p_sej_a.area_ref = [transpose(area); 1]; % (m^2)

% Sim settings
in.s.traj.rate = 50; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 2000; % double, s, maximum allowable time

% NPC Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 125e3;    % maximum allowable altitude (m)
in.v.gnc.g.p_sej_a.alt_min = min_alt; % m
in.v.gnc.g.p_sej_a.alt_max = max_alt; % m
in.v.gnc.g.p_sej_a.cd = in.v.aero.cd; % nd
in.v.gnc.g.p_sej_a.cl = in.v.aero.cl; % nd
in.v.gnc.g.p_sej_a.mass = in.v.mp.m_ini; % kg
in.v.gnc.g.p_sej_a.mass_jettison = in.v.mp.m_jettison;
in.v.gnc.g.p.rate = 1;  % Guidance execution rate (Hz)

% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Setup - Error Sources
% Error sources for the periapsis lowering burn
err.rp1 = 2;       % error in delivered perigee (km)
err.ra1 = 10;       % error in delivered apogee (km)
err.t = 300;        % error in time to apogee, (seconds)
err.plmdv = 0.017;    % percent of DV magnitude
err.pitch1 = 0;%5;    % burn pitch error (deg)
err.yaw1 = 0;%5;      % burn yaw error (deg)

% Error sources for the atmospheric portion
err.atm_data = load(['atm_earth_gram2007_1000mc.mat']);     % Atmospheric density profiles with dispersion applied
err.v = 0;          % Error in entry velocity (m/s)
err.cn = 0.05;      % Error in normal force coefficient (%)
err.ca = 0.03;      % Error in axial force coefficient (%)

% Error sources for the periapsis raising burn
err.coast = 5;        % coast duration (min)
err.burn = 0;         % burn duration (min)
err.thrust = 0;       % thrust magnitude (%)
err.isp = 0.0106;     % isp (sec)
err.pitch2 = 0;       % burn pitch error (deg)
err.yaw2 = 0;         % burn yaw error (deg)

%% Analysis - Runs Perigee-Lowering MC
out1 = ac_MC_PLM(in, err, nruns);   % Input structure can be any of the 3; only trajectory information is used
% Gathers output
err.fpa = out1.fpa;

%% Analysis - Runs Atmospheric MC for Different Guidance Choices
% Low-rate NPC:
fprintf('\nRunning atmospheric MC. Completion percentage:\n');
out2 = ac_MC_atmospheric(in, err, ha_tgt, nruns);

% Prepare for 3rd MC
err.oe = out2.oe;

%% Analysis - Runs Perigee-Raising MC
out3 = ac_MC_PRM(in, err, nruns, hp_tgt);

% Calculating variables of interest
ra_new = out3.ra_new;
for i = 1:nruns
    % Determines re-entry and escape orbits
    if(ra_new(i) <= 125)
        ra_new(i) = 0;
    elseif(ra_new(i) >= 12000)
        ra_new(i) = 12000;
    end
end

%% Results - Setup
close all;

% Check to save data
% saveDat = true;     % Yes, save data
saveDat = false;    % No, don't save data

% Setup for saving results
fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\AAS GNC\AAS GNC Results';    % Path for home
% fpath = 'C:\Users\mwerner9\Documents\Research\AAS GNC\Figures';    % Path for school
size = 20;
font = 'Times New Roman';

%% Results - Plotting Apoapsis Altitudes
% NPC:
figure(1); 
hist(ra_new, 50); grid on;
set(gca, 'FontSize', size, 'FontName', font);
xlabel({'Final Altitude at Apoapsis, km'; '(a)'}, 'FontSize', size, 'FontName', font);
ylabel('Number of Samples', 'FontSize', size, 'FontName', font);
meanRA = mean(ra_new);
stdRA = std(ra_new);
minRA = min(ra_new);
maxRA = max(ra_new);
fname = 'Earth MC Apoapsis - NPC';
if(saveDat)
    saveas(gca, fullfile(fpath, fname), 'fig');
    saveas(gca, fullfile(fpath, fname), 'tif');
end


% Saving data
if(saveDat)
    save('AC_Earth_MC_Full');
end