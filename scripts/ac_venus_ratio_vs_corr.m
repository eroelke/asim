%% ac_venus_ratio_vs_corr
%   A script to explore the trajectory space for single-event-jettison, 
%   drag-modulated aerocapture mission to Venus. A series of ballistic
%   coefficient ratios and entry velocities are used to calculate entry 
%   flight-path angle corridor width. 
%
%   Assumptions: Constant-mass spacecraft and drag skirt. Difference in 
%   ballistic coefficient ratios is caused by increasing drag area.
%% Prepare workspace
clc;
clearvars;
close all;

%% Setup - Basic Inputs
% Vehicle properties
R = 0.75;           % Pre-jettison radius (m)
r = 0.175;          % Post-jettison radius (m)
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
m = 80;            % Pre-jettison mass (kg)
rn = 1/2*r;         % Nose radius (m)
delta = 45;         % Sphere-cone half angle (deg)
cd = 1.0;          % Drag Coefficient

%% Setup - Input Structure
in0 = dme_venus_in; % nominal input structure

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = -8*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = find_vel_pp_mag(10500, in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt)); % initial velocity, m/s

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = m; % Initial mass (kg)
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

ha_tgt = 2000e3;                % km, target post-apoapse altitude

beta_ratio_range = [2 5 10 20];    % nd, range of pre- and post-jettison ballistic coefficient ratios

area_range = area(1)./beta_ratio_range;    % m^2, range of post-jettison areas

v_entry = (10:1:14)*1e3;        % m/s, range of inertial velocities at entry

fpa_guess = [-15 -3.5]*pi/180;    % rad, entry FPA guess for bisection (MUST bracket corridor)

tol_fpa = 0.01*pi/180;          % rad, tolerance on calculation of required FPA
tol_ha = 10e3;                  % m, tolerance on final apoapse altitude

%% Calculations - Initialization
% Set up output data structures
num_ratios = length(beta_ratio_range); % Number of ballistic coefficient ratios
num_v_entry = length(v_entry);      % Number of entry velocities

fpa_corr = NaN(num_v_entry, num_ratios, 2);   % FPA Corridors
corr_width = NaN(num_v_entry, num_ratios);    % Corridor widths
ha = fpa_corr;                              % Final altitude for each run

% Prepare generic input structures
% Case 1: Pre-jettison, minimum beta (shallow side of corridor)
in(1) = in0;
in(1).v.aero.area_ref = area(1);

% Case 2: Post-jettison, maximum beta (steep side of corridor)
in(2) = in0;

%% Calculations - Compute Corridor
fprintf('Computing corridors. Progress: \n');
for i = 1:num_v_entry
    
    % Re-set initial velocities
    in(1).s.traj.vel_pp_mag = find_vel_pp_mag(v_entry(i), in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt));
    in(2).s.traj.vel_pp_mag = find_vel_pp_mag(v_entry(i), in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt));
    
    % Convert scalars to inertial position and velocity vectors
    [in(1).s.traj.r_pci_ini, in(1).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in(1).s.traj.lat, in(1).s.traj.lon, in(1).s.traj.alt, ...
        in(1).s.traj.gamma_pp, in(1).s.traj.az, in(1).s.traj.vel_pp_mag, ...
        in(1).p.omega, in(1).s.traj.t_ini, 0, in(1).p.r_e, in(1).p.r_p);
    
    [in(2).s.traj.r_pci_ini, in(2).s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in(2).s.traj.lat, in(2).s.traj.lon, in(2).s.traj.alt, ...
        in(2).s.traj.gamma_pp, in(2).s.traj.az, in(2).s.traj.vel_pp_mag, ...
        in(2).p.omega, in(2).s.traj.t_ini, 0, in(2).p.r_e, in(2).p.r_p);
    
    for j = 1:num_ratios
        % Set post-jettison area based on selected beta-ratio
        in(2).v.aero.area_ref = area_range(j);
        
        % Compute FPA required to hit target for pre- and post- jettison
        % configurations
        [fpa_corr(i, j, 1), ha(i, j, 1), dat(i, j, 1), ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in(1), tol_fpa, tol_ha);
        [fpa_corr(i, j, 2), ha(i, j, 2), dat(i, j, 2), ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in(2), tol_fpa, tol_ha);
        
        % Sets corridor to NaN if final orbits re-enter or fail to capture
        if (any(ha(i, j, :) <= 200e3) || any(ha(i, j, :) >= 25000e3))
            fpa_corr(i, j, :) = NaN;
            ha(i, j, :) = NaN;
        end
        
        % Compute corridor width
        corr_width(i, j) = (fpa_corr(i, j, 1) - fpa_corr(i, j, 2))*180/pi;
        
        % Output progress
        prog = (((i - 1)*num_ratios + j) / (num_v_entry * num_ratios));
        fprintf('%4.2f percent\n', prog*100); 
    end
end
fprintf('\n...done.\n');

%% Post-Processing - Plots
fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\SciTech\Figures';
size = 16;

% Figure 1 - Corridor width for different conditions
figure(1); hold on; box on; grid on;
plot(beta_ratio_range, corr_width(1, :), 'k-*', 'LineWidth', 1, 'DisplayName', '10 km/s');
plot(beta_ratio_range, corr_width(2, :), 'r-*', 'LineWidth', 1, 'DisplayName', '11 km/s');
plot(beta_ratio_range, corr_width(3, :), 'g-*', 'LineWidth', 1, 'DisplayName', '12 km/s');
plot(beta_ratio_range, corr_width(4, :), 'b-*', 'LineWidth', 1, 'DisplayName', '13 km/s');
plot(beta_ratio_range, corr_width(5, :), 'c-*', 'LineWidth', 1, 'DisplayName', '14 km/s');
set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
legend('show', 'Location', 'NorthEastOutside');
xlabel('\beta_2/\beta_1', 'FontSize', size, 'FontName', 'Times New Roman'); %xlim([7.5 10]);
ylabel('Corridor Width, deg', 'FontSize', size, 'FontName', 'Times New Roman'); %ylim([50 130]);
% fname = 'Venus Corridors';
% saveas(gca, fullfile(fpath, fname), 'fig');
% saveas(gca, fullfile(fpath, fname), 'png');