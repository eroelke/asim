% ac_titan_manual.m
%   manual jettison for aerocapture at titan
% 
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% Colorado Center for Astrodynamics Research
% 
% Created By: Evan Roelke, Dec 2018
% 
function out = ac_titan_manual(m,rcs,tj,efpa,v_atm,cd)

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

min_alt = 0;        % minimum allowable altitude (m)
max_alt = 1270e3;    % maximum allowable altitude (m)

%% simulation settings
in0 = manual_titan_in;
atm = load('./data/atm_titan_gram_10000mc_reg_winds.mat');
in0.p.atm.table = atm.atm_mc_reg_tables.nom_table;

% entry state
in0.s.traj.gamma_pp = efpa * pi/180;
in0.s.traj.alt = max_alt; % initial altitude, m
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
in0.s.traj.t_max = 2000; % double, s, maximum allowable time

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at max altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

% Guidance Settings
in0.v.gnc.g.p_manual.masses(1:2) = [m(1) m(2)];
in0.v.gnc.g.p_manual.t_jetts(1) = tj;  % s
in0.v.gnc.g.p_manual.n_jett = uint8(1);
in0.v.gnc.g.p_manual.area_refs(1:2) = area;


% run nominal trajectory
out = main1_mex(in0);

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


end %ac_titan_manual