% run_dej_n_manual.m
%   manual jettison event for given planet
% 
% Written by: Evan Roelke 
% 
% 
function out = run_dej_n_manual(x0,gnc,aero,sim,mc)
format longg

%% Derive values from inputs
stages = gnc.n+1;   % number of stages (including pre-jettison)

areas = pi.*aero.rcs.^2;
max_alt = sim.h_max * 1e3;
min_alt = sim.h_min * 1e3;

%% Base Input Structure
in0 = define_in;

%% Planetary data
if ischar(sim.planet) == true
    models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
    planet = find(strcmpi(sim.planet,models) == 1);
end
if sim.atm_mode > 3 || sim.atm_mode < 1
    fprintf('Warning: Incorrect Atmospheric Mode. Defaulting to Exponential Model\n');
    sim.atm_mode = 1;
end
in0.p = get_planetary_data(planet, sim.atm_mode); % Load planet model, GRAM atmosphere + winds

% Get updated atmospheric profile
[in0.p.atm.table, ~] = parse_atm_data(in0, sim.planet);
in0.p.atm.dh = in0.p.atm.table(2,1) - in0.p.atm.table(1,1);
in0.v.gnc.g.p.planet.atm_nom = in0.p.atm.table(:,1:7);
in0.v.gnc.g.p.planet.atm_true = in0.p.atm.table(:,1:7);

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at top of atmosphere
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

in0.s.traj.t_ini = 0;
in0.s.traj.t_max = sim.t_max;

%% Define Entry State
% in0.s.traj.gamma_pp = -5.5*pi/180;
in0.s.traj.gamma_pp = x0.fpa0*pi/180;
in0.s.traj.alt = max_alt; % initial altitude, m
in0.s.traj.lat = x0.lat0 * pi/180; % initial latitude, rad
in0.s.traj.lon = x0.lon0 * pi/180; % initial longitude, rad
in0.s.traj.az = x0.az0 * pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = x0.v0 * 1e3;    %m/s, entry speed

% Convert PCPF to PCI vectors
if (any(x0.pci) > 0)
    if (size(x0.pci,2) == 6)
        x0.pci = x0.pci';
    end
    in0.s.traj.r_pci_ini = x0.pci(1:3);
    in0.s.traj.v_pci_ini = x0.pci(4:6);
    % 
else
    [in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
        in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
        in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);
end

%% Vehicle properties
in0.v.mp.m_ini = aero.m(1); % Initial mass (kg)
in0.v.mp.m_jettison = 0; % Jettisoned mass (kg) - unused for this guidance algorithm
in0.v.aero.area_ref = areas(1); % initial area difference before and after jettison (m^2)
in0.v.aero.mode = uint8(1); % User-specified, constant lift and drag coefficients
in0.v.aero.cd = aero.cds(1);
in0.v.aero.cl = aero.cls(1);
in0.v.aero.nose_radius = aero.rn; % Nose radius (m)

%% Guidance Settings
in0.v.gnc.g.p.atm.K_gain = 0.2;
in0.v.gnc.g.p.atm.K_bounds = [0.01 2];
in0.v.gnc.g.p.atm.Kflag = uint8(0);     % flag to do density estimate (for monte carlo)
in0.v.gnc.g.p.planet.alt_max = max_alt;
in0.v.gnc.g.p.planet.alt_min = min_alt;
in0.v.gnc.g.p.atm.ens_tol = gnc.p_tol;  % tolerance on ensemble search range (% of density estimate)

% Simulation Rates
in0.s.traj.rate = sim.traj_rate;           % double, Hz, integration rate
in0.v.gnc.g.p.rate = gnc.guid_rate;         % nd, guidance call rate (Hz)
in0.v.gnc.g.p_dej_n.traj_rate = sim.traj_rate;  % guidance integration rate (Hz)
in0.s.data_rate = sim.traj_rate;          % Hz, data storage frequency
in0.v.gnc.n.p.rate = sim.traj_rate;   % navigation rate = traj rate
in0.v.gnc.c.p.rate = sim.traj_rate;    % control rate = traj rate

%% Manual Guidance Settings
in0.v.gnc.g.p.dm_mode = uint8(6);  % manual jettison
% in0.v.gnc.g.p_manual.cds(1:stages) = aero.cds;
in0.v.gnc.g.p_manual.masses(1:stages) = aero.m'; % kg
in0.v.gnc.g.p_manual.area_refs(1:length(areas)) = transpose(areas); % m^2
in0.v.gnc.g.p_manual.n_jett = uint8(gnc.n);
in0.v.gnc.g.p_manual.t_jetts(1:gnc.n) = gnc.tj0;


%% Run Simulation
out = main1_mex(in0);
% out = main1(in0);

% Compute parameters
if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end

beta = nan(length(areas),1);
beta_ratio = nan(length(areas)-1,1);
for i = 1:length(areas)
    beta(i) = aero.m(i)/(areas(i)*aero.cds(i));
    if i > 1
        beta_ratio(i-1) = beta(i)/beta(i-1);
    end
end

rv = [out.traj.pos_ii(idxend,:), out.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e)/1000;

dv = out.traj.vel_ii_mag(1) - out.traj.vel_ii_mag(idxend);

out.haf = round(haf,5);                   %km
out.haf_err = round((haf - gnc.ha_tgt),5);    %km
out.dv = round(dv,5);
out.idxend = idxend;

for i = 1:gnc.n
    tjett = (find(out.g.manual.stage == i,1)-1)/(in0.s.data_rate);
    if (~isempty(tjett))
        out.t_jett(i) = tjett;
        out.idj(i) = find(out.traj.time >= out.t_jett(i),1);
    else
        out.t_jett(i) = nan;
        out.idj(i) = nan;
    end
end
out.betas = round(beta,3);
out.beta_ratios = round(beta_ratio,3);

% TO DO: Monte Carlos

end %run_dej_n_manual.m

