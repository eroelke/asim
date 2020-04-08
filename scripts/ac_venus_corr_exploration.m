%% ac_venus_corr_exploration
%   A script to explore the trajectory space for single-event-jettison, 
%   drag-modulated aerocapture mission to Venus. A series of target
%   altitudes and entry velocities are used to calculate entry flight-path
%   angle corridor width. 
%
% Major Revision History:
%   Earlier - Functionality written, Zach Putnam
%   July 2017 - Script developed, Michael Werner
%% Prepare workspace
clc;
clearvars;
close all;

%% Setup - Basic Inputs
% Vehicle properties
R = 0.75;           % Pre-jettison radius (m)
r = 0.175;          % Post-jettison radius (m)
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
m1 = 80;            % Pre-jettison mass (kg)
m2 = 40;            % Post-jettison mass (kg)
rn = 1/2*r;         % Nose radius (m)
delta = 45;         % Sphere-cone half angle (deg)
cd = 1.0;          % Drag Coefficient

%% Setup - Input Structure
in0 = dme_venus_in; % nominal input structure

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = -5.5*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = 10.65e3;

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
in0.s.traj.rate = 50; % double, Hz, integration rate
in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.traj.t_max = 2000; % double, s, maximum allowable time

% Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 150e3;    % maximum allowable altitude (m)

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 150km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Setup - User Inputs

ha_tgt = [400 1000 2000 10000]*1e3;  % m, range of target post-aerocapture apoapse altitudes
v_entry = (10:1:14)*1e3;        % m/s, range of inertial velocities at entry

fpa_guess = [-9 -3.5]*pi/180;    % rad, entry FPA guess for bisection (MUST bracket corridor)

tol_fpa = 0.01*pi/180;          % rad, tolerance on calculation of required FPA
tol_ha = 10e3;                  % m, tolerance on final apoapse altitude

%% Calculations - Initialization
% Set up output data structures
num_ha = length(ha_tgt);        % Number of target apoapses
num_v_entry = length(v_entry);  % Number of entry velocities

fpa_corr = NaN(num_v_entry, num_ha, 2);   % FPA Corridors
corr_width = NaN(num_v_entry, num_ha);    % Corridor widths
ha = fpa_corr;                              % Final altitude for each run

% Prepare generic input structures
% Case 1: Pre-jettison, minimum beta (shallow side of corridor)
in(1) = in0;    
in(1).v.aero.area_ref = area(1);

% Case 2: Post-jettison, maximum beta (steep side of corridor)
in(2) = in0;
in(2).v.aero.area_ref = area(2);
in(2).v.mp.m_ini = m2;

v_entry = 10.65e3;

%% Calculations - Compute Corridor
fprintf('Computing corridors. Progress: \n');
% for i = 1:num_v_entry
    % Re-set initial velocities
    in(1).s.traj.vel_pp_mag = find_vel_pp_mag(v_entry, in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt));
    in(2).s.traj.vel_pp_mag = find_vel_pp_mag(v_entry, in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt));
    
    % Convert scalars to inertial position and velocity vectors
    [in(1).s.traj.r_pci_ini, in(1).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in(1).s.traj.lat, in(1).s.traj.lon, in(1).s.traj.alt, ...
        in(1).s.traj.gamma_pp, in(1).s.traj.az, in(1).s.traj.vel_pp_mag, ...
        in(1).p.omega, in(1).s.traj.t_ini, 0, in(1).p.r_e, in(1).p.r_p);
    
    [in(2).s.traj.r_pci_ini, in(2).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in(2).s.traj.lat, in(2).s.traj.lon, in(2).s.traj.alt, ...
        in(2).s.traj.gamma_pp, in(2).s.traj.az, in(2).s.traj.vel_pp_mag, ...
        in(2).p.omega, in(2).s.traj.t_ini, 0, in(2).p.r_e, in(2).p.r_p);
    
    for j = 1:num_ha
        % Compute FPA required to hit target for pre- and post- jettison
        % configurations
        [fpa_corr(j, 1), ha(j, 1), ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt(j), in(1), tol_fpa, tol_ha);
        [fpa_corr(j, 2), ha(j, 2), ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt(j), in(2), tol_fpa, tol_ha);
        
        % Sets corridor to NaN if final orbits re-enter or fail to capture
        if (any(ha(j, :) <= 200e3) || any(ha(j, :) >= 21500e3))
            fpa_corr(j, :) = NaN;
            ha(j, :) = NaN;
        end
        
        % Compute corridor width
        corr_width(j) = (fpa_corr(j, 1) - fpa_corr(j, 2))*180/pi;
        
        if mod(j,10) == 1
            % Output progress
            fprintf('Run %i\n', j); 
        end
    end
% end
fprintf('\n...done.\n');

%% Post-Processing - Plots
% fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\SciTech\Figures';
% size = 16;
% 
% % Figure 1 - Corridor width for different conditions
% figure(1); hold on; box on; grid on;
% plot(ha_tgt/1e3, corr_width(1, :), 'r', 'LineWidth', 2, 'DisplayName', '10 km/s');
% plot(ha_tgt/1e3, corr_width(2, :), 'g', 'LineWidth', 2, 'DisplayName', '11 km/s');
% plot(ha_tgt/1e3, corr_width(3, :), 'b', 'LineWidth', 2, 'DisplayName', '12 km/s');
% plot(ha_tgt/1e3, corr_width(4, :), '--r', 'LineWidth', 2, 'DisplayName', '13 km/s');
% plot(ha_tgt/1e3, corr_width(5, :), '--g', 'LineWidth', 2, 'DisplayName', '14 km/s');
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% legend('show', 'Location', 'NorthEastOutside');
% xlabel('Post-Aerocapture Altitude at Apoapsis, km', 'FontSize', size, 'FontName', 'Times New Roman'); %xlim([7.5 10]);
% ylabel('Corridor Width, deg', 'FontSize', size, 'FontName', 'Times New Roman'); %ylim([50 130]);
% % fname = 'Venus Corridors';
% % saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');