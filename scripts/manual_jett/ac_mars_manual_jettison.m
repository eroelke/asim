%% Prepare workspace
% clear all; close all; clc;
% % close all;
% clear all;

% Turned into function, Dec 2018 E. Roelke


function out = ac_mars_manual_jettison(m,rcs,tj,efpa,v_atm,cd)

%% Setup - Basic Inputs
% Vehicle properties
% R = 1;           % Pre-jettison radius (m)
% r = 0.3;            % Post-jettison radius (m)
R = rcs(1);
r = rcs(2);
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
% m1 = 30;          % Pre-jettison mass (kg)
% m2 = 20;         % Post-jettison mass (kg)
m1 = m(1);
m2 = m(2);
rn = 9/16*r;       % Nose radius (m)
delta = 60;         % Sphere-cone half angle (deg)
% cd = 1.5;          % Drag coefficient

%% Setup - Input Structure
in0 = manual_mars_in; % nominal input structure
atm_data = load(['atm_mars_gram2010_1000mc.mat']);     % Atmospheric density profiles with dispersion applied
in0.p.atm.table = atm_data.nom_table;

% Initial state: Trajectory properties
% in0.s.traj.gamma_pp = -10.918*pi/180;
in0.s.traj.gamma_pp = efpa * pi/180;
in0.s.traj.alt = 125e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = v_atm * 1000;   %m/s

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = m1; % Initial mass (kg)
in0.v.mp.m_jettison = m1-m2; % Jettisoned mass (kg)
in0.v.aero.area_ref = area(1); % Area before and after jettison (m^2)
in0.v.aero.mode = uint8(1);     % User-specified lift and drag coefficients
in0.v.aero.cd = cd;
in0.v.aero.cl = 0;
in0.v.aero.nose_radius = rn; % Nose radius (m)

% Sim settings
in0.s.traj.rate = 100; % double, Hz, integration rate
in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.traj.t_max = 750; % double, s, maximum allowable time

% Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 125e3;    % maximum allowable altitude (m)
in0.v.gnc.g.p_manual.masses(1:2) = [m(1) m(2)];
in0.v.gnc.g.p_manual.t_jetts(1) = tj;  % s
in0.v.gnc.g.p_manual.n_jett = uint8(1);
in0.v.gnc.g.p_manual.area_refs(1:2) = area;

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Run trajectories at minimum and maximum allowable FPA
% Case 1: shallowest FPA, minimum beta (shallow side of corridor)
% in(1) = in0;
% in(1).v.aero.area_ref = area(1);
% in(1).s.traj.gamma_pp = FPA;
% 
% [in(1).s.traj.r_pci_ini, in(1).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
%     in(1).s.traj.lat, in(1).s.traj.lon, in(1).s.traj.alt, ...
%     in(1).s.traj.gamma_pp, in(1).s.traj.az, in(1).s.traj.vel_pp_mag, ...
%     in(1).p.omega, in(1).s.traj.t_ini, 0, in(1).p.r_e, in(1).p.r_p);

out = main1_mex(in0);

% Compute parameters
% beta_1 = in(1).v.mp.m_ini/(in(1).v.aero.area_ref*cd);
% beta_2 = (in(1).v.mp.m_ini - in(1).v.mp.m_jettison)/(in(1).v.gnc.g.p_manual.area_ref(2)*cd);
% out.beta_ratio = beta_2/beta_1;

if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end

rv = [out.traj.pos_ii(idxend,:), out.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
out.haf = round((ra - in0.p.r_e)/1000,5);

xi = -in0.p.mu/oe(1);   % energy = -mu/a
if xi > 0
    fprintf('WARNING: Hyperbolic Exit Trajectory\n')
end
dv = round(out.traj.vel_ii_mag(1) - out.traj.vel_ii_mag(idxend),5);




%% Plots
% 
% fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\SciTech\Figures';
% size = 16;

% % Vel vs Altitude
% figure(1); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.alt/1000,'-k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Mars-Relative Velocity, km/s';'(a)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([3 6]);
% ylabel('Altitude, km', 'FontSize', size, 'FontName', 'Times New Roman'); ylim([50 150]);
% fname = 'Earth Alt';
% saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');
% 
% % Vel vs Deceleration
% figure(2); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading,'k','LineWidth',2);
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

% % Vel vs Altitude
% subplot(3, 1, 1); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.alt/1000,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Altitude (km)'); ylim([60 130]);
% title('\bfAltitude')
% 
% % Vel vs Deceleration
% subplot(3, 1, 2); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Speed Deceleration (g)'); ylim([0 4]);
% title('\bfDeceleration')
% 
% % Vel vs Heat Rate
% subplot(3, 1, 3); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.heat_rate/10000,'k','LineWidth',2);
% xlabel('Planet-relative Velocity (km/s)'); xlim([7.5 10]);
% ylabel('Heat Rate (W/cm^2)'); ylim([0, 500]);
% title('\bfHeat Rate')