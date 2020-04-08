%% ac_mars_auto_jettison
%   Drag-modulated-aerocapture trajectory. A numeric-predictor-corrector 
%   guidance scheme is used as the means of controlling the drag skirt 
%   jettison event.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Major Revision History:
%   Earlier - NPC written, Zach Putnam
%   January 2017 - Developed, Michael Werner

%% Prepare workspace
clearvars; clc;

%% Setup - Basic Inputs
% Vehicle properties
R = 1;           % Pre-jettison radius (m)
r = 0.3;            % Post-jettison radius (m)
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
m1 = 30;          % Pre-jettison mass (kg)
m2 = 20;         % Post-jettison mass (kg)
rn = 9/16*r;       % Nose radius (m)
delta = 70;         % Sphere-cone half angle (deg)
cd = 1.61;          % Drag coefficient

%% Setup - Input Structure
in0 = sej_a_mars_in; % nominal input structure
atm_data = load(['atm_mars_gram2010_1000mc.mat']);     % Atmospheric density profiles with dispersion applied
in0.p.atm.table = atm_data.nom_table;

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = -10.918*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = find_vel_pp_mag(6000, in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt)); % initial velocity, m/s

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = m1; % Initial mass (kg)
in0.v.mp.m_jettison = m1-m2; % Jettisoned mass (kg)
in0.v.aero.area_ref = area(1); % Area before and after jettison (m^2)
in0.v.aero.nose_radius = rn; % Nose radius (m)
in0.v.aero.mode = uint8(1);     % User-specified lift and drag coefficients
in0.v.aero.cd = cd;
in0.v.aero.cl = 0;


% SEJ Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 150e3;    % maximum allowable altitude (m)
in0.v.gnc.g.p_sej_a.alt_min = min_alt; % m
in0.v.gnc.g.p_sej_a.alt_max = max_alt; % m
in0.v.gnc.g.p_sej_a.cd = in0.v.aero.cd;
in0.v.gnc.g.p_sej_a.cl = in0.v.aero.cl;
in0.v.gnc.g.p_sej_a.mass = in0.v.mp.m_ini; % kg
in0.v.gnc.g.p_sej_a.mass_jettison = in0.v.mp.m_jettison;
in0.v.gnc.g.p_sej_a.atm_table = in0.p.atm.table(:, 1:7);

in0.v.gnc.g.p.rate = 5;  % Guidance rate (hz)


% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Inputs
area = [pi*(R)^2, pi*(r)^2];    % m^2
FPA = -10.96*pi/180;    % Rad

in0.v.gnc.g.p_sej_a.area_ref = [transpose(area); 0]; % m^2

ha_tgt = 1760e3;    % m

%% Run trajectories at minimum and maximum allowable FPA
% Case 1: shallowest FPA, minimum beta (shallow side of corridor)
in(1) = in0;
in(1).v.aero.area_ref = area(1);
in(1).s.traj.gamma_pp = FPA;
in(1).v.gnc.g.p_sej_a.tgt_ap = ha_tgt; % m

[in(1).s.traj.r_pci_ini, in(1).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in(1).s.traj.lat, in(1).s.traj.lon, in(1).s.traj.alt, ...
    in(1).s.traj.gamma_pp, in(1).s.traj.az, in(1).s.traj.vel_pp_mag, ...
    in(1).p.omega, in(1).s.traj.t_ini, 0, in(1).p.r_e, in(1).p.r_p);

dat(1) = main1_mex(in(1));

% Compute parameters
beta_1 = in(1).v.mp.m_ini/(in(1).v.aero.area_ref*cd)
beta_2 = (in(1).v.mp.m_ini - in(1).v.mp.m_jettison)/(in(1).v.gnc.g.p_sej_a.area_ref(2)*cd)

if isnan(dat(1).traj.alt(end))
    idx1 = find(isnan(dat(1).traj.alt),1)-1;
else
    idx1 = length(dat(1).traj.alt);
end

rv = [dat(1).traj.pos_ii(idx1,:), dat(1).traj.vel_ii(idx1,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e)/1000

%% Plots
% 
% fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\SciTech\Figures';
size = 16;

% Vel vs Altitude
figure(1); hold on; box on; grid on;
plot(dat.traj.vel_pp_mag/1000, dat.traj.alt/1000,'-k','LineWidth',2);
set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
xlabel({'Mars-Relative Velocity, km/s';'(a)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([3.5 6]);
ylabel('Altitude, km', 'FontSize', size, 'FontName', 'Times New Roman'); ylim([50 150]);
% fname = 'Earth Alt';
% saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');
% 
% % Vel vs Deceleration
figure(2); hold on; box on; grid on;
plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading,'k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Earth-Relative Velocity, km/s';'(b)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([7.5 10]);
% ylabel('Sensed Deceleration, g', 'FontSize', size); ylim([0 4]);
% fname = 'Earth Decel';
% saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');
% 
% % Vel vs Heat Rate
% figure(3); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.heat_rate/10000,'k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Earth-Relative Velocity, km/s'; '(c)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([7.5 10]);
% ylabel('Heat Rate, W/cm^2', 'FontSize', size); ylim([0, 500]);
% fname = 'Earth Heat';
% saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');
% 

% Subplot format:
% figure(1);

% Vel vs Altitude
% subplot(3, 1, 1); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.alt/1000,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Altitude (km)'); ylim([60 130]);
% title('\bfAltitude')

% Vel vs Deceleration
% subplot(3, 1, 2); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Speed Deceleration (g)'); ylim([0 4]);
% title('\bfDeceleration')

% Vel vs Heat Rate
% subplot(3, 1, 3); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.heat_rate/10000,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Heat Rate (W/cm^2)'); ylim([0, 500]);
% title('\bfHeat Rate')