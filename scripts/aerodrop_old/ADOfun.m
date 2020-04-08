% ADOfun
%   Uses specified input structure and user-defined parameters to simulate
%   an atmospheric trajectory. Useful for testing quick configuration
%   changes and troubleshooting.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Major Revision History:
%   May 2016 - Created, M. Werner


function [dat,rv,oe,ra,haf,dv,m] = ADOfun(r, Ball, rn, delta, cd, fpa, velinit)


%% Set up input structure

area = pi*(r)^2;
m = Ball * area * cd; % get mass from ballistic coefficient desired

% Staring input structure
% in0 = dme_earth_in; % nominal input structure
in0 = earth_in;

% Initial state: Vehicle properties
in0.v.mp.m_ini = m;
in0.v.aero.area_ref = area;
in0.v.aero.nose_radius = rn;
in0.v.aero.mode = uint8(1);
in0.v.aero.cd = cd;
in0.v.aero.cl = 0; 

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = fpa;
in0.s.traj.alt = 125e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
% in0.s.traj.lon = 168.7778*pi/180; % initial longitude, rad
in0.s.traj.lon = 0*pi/180;
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
% in0.s.traj.vel_pp_mag = find_vel_pp_mag(10311, in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt)); % initial velocity, m/s
in0.s.traj.vel_pp_mag = velinit; % input

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Sim settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 125e3;    % maximum allowable altitude (m)
in0.s.traj.rate = 50; % double, Hz, integration rate
in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.traj.t_max = 2000; % double, s, maximum allowable time

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Run simulation

% dat = main1_mex(in0);
dat = main1(in0);

%% Compute parameters

% Find end of trajectory
if isnan(dat.traj.alt(end))
    idxend = find(isnan(dat.traj.alt),1)-1;
else
    idxend = length(dat.traj.alt);
end


cd = dat(1).traj.cd(idxend);
beta = in0.v.mp.m_ini/(in0.v.aero.area_ref*cd);

rv = [dat.traj.pos_ii(idxend,:), dat.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e)/1000;

dv = dat.traj.vel_ii_mag(1) - dat.traj.vel_ii_mag(idxend);

end
