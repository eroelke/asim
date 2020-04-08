%% ac_earth_dcf_jettison
%   DEFT aerocapture trajectory. A deceleration-curve-fit guidance scheme
%   (similar to Mars Pathfinder) is used as the means of controlling the
%   drag skirt jettison event.
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
% 
% Inputs:
%   fpa_atm: entry flight path angle (deg)
%   trag_rate: integration rate (Hz)
%   ha_tgt: target apoapsis (km)
%   dt: g2 measurment delta t
%   mc: monte carlo struct
%       flag: run monte carlo or not
%       N: number of monte carlo sims
%       sigs: perturbations to initial state
%           [m (kg),cl,cd,aoa,v-atm (m/s),h_atm (m),fpa (rad),lat,lon,az]
% 
%
% Major Revision History:
%   Created Jun 2018 - E. Roelke 
%   Added Monte Carlo: Oct 2018 - E. Roelke

%% Prepare workspace

function out = ac_venus_dcf(fpa_atm,traj_rate,m,rcs,ha_tgt,dt,mc)
%% Setup - Basic Inputs
% Vehicle properties
% b2/b1 = 9.8
% R = 0.75;           % Pre-jettison radius (m)
% r = 0.175;          % Post-jettison radius (m)

m1 = m(1);
m2 = m(2);
r1 = rcs(1);
r2 = rcs(2);

area = [pi*(r1)^2, pi*(r2)^2];    % Drag area (m^2)
% m1 = 80;            % Pre-jettison mass (kg)
% m2 = 40;            % Post-jettison mass (kg)
rn = 1/2*r1;         % Nose radius (m)
% delta = 45;         % Sphere-cone half angle (deg)
cd = 1.05;

% b2/b1 = 4.5
% R = 1.5;
% r = 0.5;
% rn = r/2;
% m1 = 200;
% m2 = 100;
% area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)

b2 = m2/(cd*pi*r2^2);
b1 = m1/(cd*pi*r1^2);
beta_ratio = b2/b1;

%% Setup - Input Structure
in0 = dcfj_venus_in; % nominal input structure
% in0.v.gnc.g.p_dcfj.delT = 20;
% temp = load('./data/atm_venus_JPL.mat');
% in0.p.atm.table = temp.table;

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = fpa_atm*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0 * pi/180; % initial latitude, rad
in0.s.traj.lon = 0 * pi/180; % initial longitude, rad
in0.s.traj.az = 90 * pi/180; % initial azimuth, deg->rad
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

% Sim settings
in0.s.traj.rate = traj_rate; % double, Hz, integration rate (dt = 1/rate)
in0.s.data_rate = traj_rate;
in0.v.gnc.g.p.rate = traj_rate; % guidance rate
in0.v.gnc.c.p.rate = traj_rate; % nav rate
in0.v.gnc.n.p.rate = traj_rate; % control rate

in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.traj.t_max = 500; % double, s, maximum allowable time

% SEJ Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 150e3;    % maximum allowable altitude (m)
in0.v.gnc.g.p_dcfj.mass = in0.v.mp.m_ini; % kg
in0.v.gnc.g.p_dcfj.mass_jettison = in0.v.mp.m_jettison;

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Guidance Parameters

% get polyfit params
if (dt == 10)
    temp = load('./data/dcfj_data/curve_ha_tgt.mat');
    if (ha_tgt == 10000)
        p = temp.p_10k;
    else
        p = temp.p_2k;
    end 
else
    temp = load('./data/dcfj_data/curve_dt.mat');
    if (dt == 15)
        p = temp.p(:,2);
    elseif (dt == 20)
        p = temp.p(:,3);
    end
end

in0.v.gnc.g.p_dcfj.decel_curve = p;
% in0.v.gnc.g.p_dcfj.yfit = yfit;
% in0.v.gnc.g.p_dcfj.ttogo = ttogo;
% in0.v.gnc.g.p_dcfj.g2s = g2;

in0.v.gnc.g.p_dcfj.dt1 = dt;   % time step
in0.v.gnc.g.p_dcfj.dt2 = 0;
in0.v.gnc.g.p_dcfj.g1 = 0.3;    % accel of first g measurement
in0.v.gnc.g.p_dcfj.area_ref = transpose(area); % m^2
in0.v.gnc.g.p_dcfj.mass_jettison = in0.v.mp.m_jettison;
in0.v.gnc.g.p_dcfj.mass = in0.v.mp.m_ini;

%% Run trajectory simulation
in0.v.aero.area_ref = area(1);
% in0.s.traj.gamma_pp = FPA;

% [in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
%     in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
%     in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
%     in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Run nominal case
out = main1_mex(in0);
% fprintf('Running DCFJ in debug mode...\n'); out = main1(in0);

if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end

rv = [out.traj.pos_ii(idxend,:), out.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e);

dv = out.traj.vel_ii_mag(1) - out.traj.vel_ii_mag(idxend);
out.dv = round(dv,5);

idj = find(out.g.dcfj.stage >= 1,1);
tjett = out.traj.time(idj);
if (isempty(tjett))
    out.tjett = nan;
    out.idj = nan;
else
    out.tjett = tjett;
    out.idj = idj;
end

out.haf = round(haf/1000,5);    %km
out.haf_err = round(haf/1000 - ha_tgt,5);   % km
out.idxend = idxend;

tol = 1e-4;
g2_nom = [];
while isempty(g2_nom)
    g2_nom = out.g.dcfj.g2(out.traj.time >= out.g.dcfj.g2_time - tol & ... 
        out.traj.time <= out.g.dcfj.g2_time + tol);
    tol = tol + 1e-4;
end


% check monte carlo
if (mc.flag)
    out_nom = out;
    in0.v.gnc.g.p_dcfj.g2_nom = g2_nom;
    clear out;  % delete nominal trajectory
    out.nom = out_nom;  % nominal trajectory
    % load disperse atmospheric tables
    mc_atm = load('./data/atm_data/atm_venus_mc.mat');
    mc_atm = mc_atm.mc; % 1000 atm tables
    
    % initialize output params
    N = mc.N;
    len = in0.s.traj.t_ini:(1/in0.s.data_rate):in0.s.traj.t_max;
    len = length(len)+1;
    vmag = nan(len,mc.N);
    alt = vmag; fpa = vmag; t = vmag; 
%     rho = vmag; g = vmag;
    t_jett = nan(N,1);
    idj = t_jett;
    idxend = nan(N,1);
    haf = nan(N,1);
    haf_err = haf; dv = haf;
    sigs = mc.sigs;
    flag = double(mc.flag);
    % pre-determine random atmospheres to reduce broadcast overhead
    temp = nan(1000, 7, N);
    for i = 1:N
        rnd = randi(1000);
        temp(:,:,i) = mc_atm(rnd).table;
    end
    
    
    parArg = 6;
    par = parcluster('local');
    par.NumWorkers = 6;
%     parfor (i = 1:N, parArg)
    for i = 1:N
        in_mc = in0;        % copy nominal input struct
        in_mc.v.gnc.g.p_dcfj.mcflag = flag;
%         if mod(i,N/10) == 0
%             fprintf('Monte Carlo %i percent done\n',i*100/N);
%         end
        % apply disperse atmosphere
          % new atmospheric table
          in_mc.p.atm.table = temp(:,:,i);
        
%         figure(1); hold on
%         plot((in_mc.p.atm.table(:,2) - in0.p.atm.table(:,2) )./ ... 
%             in0.p.atm.table(:,2),in_mc.p.atm.table(:,1)./1000)
%         ylim([90 150])
        
        % apply input perturbations (only to traj. guid struct has nominal)
        in_mc = apply_dispersions(in_mc,sigs);
        
%         out_mc = main1(in_mc);    % debugging purposes
        out_mc = main1_mex(in_mc);  % output struct for disperse case
        
        vmag(:,i) = out_mc.traj.vel_pp_mag./1000;   %km/s
        alt(:,i) = out_mc.traj.alt./1000;       %km
        fpa(:,i) = out_mc.traj.gamma_pp.*180/pi;    %deg
        t(:,i) = out_mc.traj.time;                  %s
%         rho(:,i) = out_mc.traj.rho;     % kg/m3
%         g(:,i) = out_mc.traj.g_loading;
        
        jett_ind = (find(out_mc.g.dcfj.stage == 1,1)-1);
        if ~isempty(jett_ind)
            idj(i) = jett_ind;
            t_jett(i) = out_mc.traj.time(idj(i));
        else
            t_jett(i) = nan;    % never jettisoned
            idj(i) = nan;
        end

        % get final index
        if isnan(out_mc.traj.alt(end))
            idxend(i) = find(isnan(out_mc.traj.alt),1)-1;
        else
            idxend(i) = length(out_mc.traj.alt);
        end
        
        % calculate apoapsis
        rv = [out_mc.traj.pos_ii(idxend(i),:), out_mc.traj.vel_ii(idxend(i),:)]';
        oe = rv2oe(rv,in0.p.mu);
        ra = oe(1)*(1+oe(2));
        err = (ra - in0.p.r_e)/1000;
        
        haf(i) = round(err,5);       %km, final apoapsis
        haf_err(i) = round((err - ha_tgt),5);    %km, apoapsis error
        dv(i) = round( out_mc.traj.vel_ii_mag(1) - out_mc.traj.vel_ii_mag(idxend(i)), 5 );
    end % monte carlo loop
    toc
    % store all monte carlo results into out struct
    out.traj.vmag = vmag;
    out.traj.alt = alt;
    out.traj.fpa = fpa;
    out.traj.t = t;
%     out.traj.rho = rho;
%     out.traj.g = g;
    out.haf = haf;
    out.haf_err = haf_err;
    out.dv = dv;
    out.tjett = t_jett;
    out.idj = idj;
    out.idxend = idxend;
    
end % check mc flag

end % ac_venus_dcf