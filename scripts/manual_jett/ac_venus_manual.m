%% Prepare workspace


function dat = ac_venus_manual(fpa,t_jett,traj_rate,ha_tgt,veh)

%% Setup - Basic Inputs
% Vehicle properties
% R = 0.75;           % Pre-jettison radius (m)
% r = 0.175;          % Post-jettison radius (m)

R = veh.r(1);
r = veh.r(2);

area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
% m1 = 80;            % Pre-jettison mass (kg)
% m2 = 40;            % Post-jettison mass (kg)

m1 = veh.m(1);
m2 = veh.m(2);

% rn = 1/2*r;         % Nose radius (m)
rn = veh.rn;
% delta = 45;         % Sphere-cone half angle (deg)
% cd = 1.05;
cd = veh.cd;

%% Setup - Input Structure
in0 = manual_venus_in; % nominal input structure

% temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_JPL.mat');  % home
% temp = load('C:\Users\Evan Roelke\Documents\research\asim\data\atm_venus_mc.mat'); % lab
temp = load('data/atm_data/atm_venus_mc.mat');
in0.p.atm.table = temp.nom_table;

% Initial state: Trajectory properties
% in0.s.traj.gamma_pp = -7.6*pi/180;

in0.s.traj.gamma_pp = fpa*pi/180;

in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = 11e3;

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

% Simulation rates
in0.s.traj.rate = traj_rate; % double, Hz, integration rate
in0.v.gnc.g.p.rate = in0.s.traj.rate; % nd, guidance rate (Hz)
in0.v.gnc.c.p.rate = in0.s.traj.rate;    % control rate = guidance rate
in0.v.gnc.n.p.rate = in0.s.traj.rate;
% in0.s.data_rate = in0.v.gnc.g.p.rate;
in0.s.data_rate = in0.s.traj.rate;

% Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 150e3;    % maximum allowable altitude (m)
in0.v.gnc.g.p_manual.masses(1:2) = [m1 m2];
% in0.v.gnc.g.p_manual.mass_jettison = in0.v.mp.m_jettison;
in0.v.gnc.g.p_manual.area_refs(1:2) = transpose(area);

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
in0.v.gnc.g.p_manual.t_jetts(1) = t_jett;  % s


if (isnan(fpa) == true)
    fpa_tol = 0.05; %deg
%     [fpa_f,~,haf_err,~] = fpa_given_haf(x0.fpa0, gnc.ha_tgt, in0, fpa_tol, gnc.ha_tol);    
    fpa_bounds = [-2, -12];
    [fpa_f, ha, ~, ~, ~] = find_fpa_given_ha(fpa_bounds, ha_tgt, in0, fpa_tol, 10); 
    if (isnan(fpa_f))
        fprintf('WARNING: no solution for EFPA found. Reverting to initial guess\n');
        fpa_f = -5.5;
    end
    
    in0.s.traj.gamma_pp = fpa_f*pi/180;
    % re-do inertial entry state
    [in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
        in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
        in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p); 
end

%% Run Trajectory Simulation
dat = main1_mex(in0);
% dat = main1(in0);

% find last index
if isnan(dat.traj.alt(end))
    idxend = find(isnan(dat.traj.alt),1)-1;
else
    idxend = length(dat.traj.alt);
end

% compute aerocapture parameters
rv = [dat.traj.pos_ii(idxend,:), dat.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e); %m

% dv = dat.traj.vel_ii_mag - dat.traj.vel_ii_mag(idxend);
dat.haf = round(haf/1000,5);    %km
dat.haf_err = round(haf/1000 - ha_tgt,5);   %km
dat.idxend = idxend;
% for i = 1:in0.v.gnc.g.p_manual.n_jett
%     dat.id_jett(i) = find(dat.g.manual.stage == 1,1);
%     dat.t_jett(i) = dat.traj.time(dat.id_jett);
% end



% v = dat.traj.vel_pp_mag;
% h = dat.traj.alt;
% lat = dat.traj.lat;
% lon = dat.traj.lon;
% t = dat.traj.time;
% fpa = dat.traj.gamma_pp;
% save('ref_traj.mat','v','h','lat','lon','t','fpa')

% figure(); hold on
% plot(dat.traj.vel_pp_mag./1000,dat.traj.alt./1000,'LineWidth',2);
% 
% figure(); hold on
% yyaxis left
% plot(dat.traj.time,dat.traj.alt./1000,'r','LineWidth',2)
% yyaxis right
% plot(dat.traj.time,dat.traj.vel_pp_mag./1000,'b','LineWidth',2);

end % ac_venus_manual



%% Plots
% 
% fpath = 'C:\Users\Michael\Documents\School\Research\Year 2\SciTech\Figures';
% size = 16;

% % % Vel vs Altitude
% figure(1); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.alt/1000,'-k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Planet-Relative Velocity, km/s';'(a)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([7 14]);
% ylabel('Altitude, km', 'FontSize', size, 'FontName', 'Times New Roman'); ylim([50 200]);
% % 
% % % Vel vs Deceleration
% figure(2); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.g_loading,'k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Planet-Relative Velocity, km/s';'(b)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([7.5 10]);
% ylabel('Sensed Deceleration, g', 'FontSize', size); ylim([0 4]);
%
% % Vel vs Heat Rate
% figure(3); hold on; box on; grid on;
% plot(dat.traj.vel_pp_mag/1000, dat.traj.heat_rate/10000,'k','LineWidth',2);
% set(gca, 'FontSize', size, 'FontName', 'Times New Roman');
% xlabel({'Planet-Relative Velocity, km/s'; '(c)'}, 'FontSize', size, 'FontName', 'Times New Roman'); xlim([7.5 10]);
% ylabel('Heat Rate, W/cm^2', 'FontSize', size); ylim([0, 500]);
