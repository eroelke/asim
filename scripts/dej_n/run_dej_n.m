%% run_dej_n.m
%   -n discrete-event jettison (dej) aerocapture
%
%   Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
%
% Inputs:
%   x0:     Entry Conditions Input Struct
%       fpa0:    entry flight path angle (deg) - [1x1]
%       v0:      entry velocity (km/s)
%       h0:      initial altitude (km)
%       lat:     initial latitude (deg)
%       lon:     initial longitude (deg)
%       az:      initial azimuth (deg)
%   gnc:    Guidance Input Struct
%       mode:       guidance dma mode
%       ha_tgt:     target apoapsis altitude (km) - [1x1]
%       ha_tol:     target apoapsis tolerance (km) - [1x1]
%       n_jett:     number of jettison stages - [1x1]
%       tj0:        initial jettison time(s) guess (sec) - [nx1]
%       dtj_lim:    jettison time update limit
%       npc_mode:   bisection method(0) or newton method(1)
%       guid_rate:  guidance rate (Hz) 
%       atm_mode:   uint8 flag for atmospheric estimation
%       iters:      maximum guidance iterations per call
%       t_init:     guidance initialization time
%       p_tol:      ensemble tolerance
%   aero:   Entry Vehicle Struct
%       rcs:        successive backshell radii (m) - [n+1x1], default [0.75 0.25 0.175 0.1]
%       m:          successive vehicle masses (kg) - [n+1x1], default [80 50 40 30]
%       rn:         nose radius (assumed constant)
%       cds:        drag coeffs per stage
%       cls:        lift coeffs per stage
%   mc: monte carlo input structure
%       mc.flag:    uint8 flag for monte carlo run (and density corrector)
%       mc.N:       number of monte carlo runs
%       mc.sigs:    uncertainty structure [m (kg),cl,cd,aoa,v-atm (m/s),h_atm (m),fpa (rad),lat,lon,az]
%   sim: input simulation structure
%       traj_rate:  trajectory integration rate (Hz)
%       data_rate:  data storage rate (Hz)
%       t_max:      maximum simulation time (s)
%       h_max:      termination above this altitude, also atmospheric interface (km)
%       h_min:      termination below this altitude (km)
%       planet:     planet to model (string or index)
%       atm_mode:   atmospheric mode
%           1:          exponential (default)
%           2:          table look-up, no winds
%           3:          table look-up, winds
%       efpa_flag:      if on, calc optimal nominal efpa before sim
%       parMode:        run in parallel processing (logical)
% 
% Input Notes:
%   n determines number of jettisons even if tj0,m,rcs are more than [nx1]
%   
% Outputs:
%   out:        simulation output data structure
%
% Created Mar 2019, Evan Roelke
% 
function out = run_dej_n(x0,gnc,aero,sim,mc)
format longg

%% Derive values from inputs
stages = gnc.n+1;   % number of stages (including pre-jettison)

areas = pi.*aero.rcs.^2;
max_alt = sim.h_max * 1e3;
min_alt = sim.h_min * 1e3;

%% Base Input Structure
in0 = dej_n_in(gnc.mode, max_alt, sim.planet, sim.atm_mode); % nominal dej-n input structure

% Get updated atmospheric profile
[in0.p.atm.table, mc_atm] = load_atm_data(in0, sim.planet);
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
in0.s.traj.alt = x0.h0 * 1e3; % initial altitude, m
in0.s.traj.lat = x0.lat0 * pi/180; % initial latitude, rad
in0.s.traj.lon = x0.lon0 * pi/180; % initial longitude, rad
in0.s.traj.az = x0.az0 * pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = x0.v0 * 1e3;    %m/s, entry speed

% Convert PCPF to PCI vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

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
in0.v.gnc.g.p.atm.Kflag = false;     % flag to do density estimate (for monte carlo)
in0.v.gnc.g.p.atm.ens_tol = gnc.ecf.p_tol;  % tolerance on ensemble search range (% of density estimate)
in0.v.gnc.g.p.atm.rss_flag = gnc.rss_flag;
in0.v.gnc.g.p.atm.mode = uint8(gnc.atm_mode);   %atmospheric estimation mode
in0.v.gnc.g.p.atm.ecf_mode = uint8(gnc.ecf.mode);
in0.v.gnc.g.p.planet.alt_max = max_alt;
in0.v.gnc.g.p.planet.alt_min = min_alt;
in0.v.gnc.g.p.planet.mode = uint8(sim.atm_mode); %guidance density mode

% Simulation Rates
in0.s.traj.rate = sim.traj_rate;           % double, Hz, integration rate
in0.v.gnc.g.p.rate = gnc.guid_rate;         % nd, guidance call rate (Hz)
in0.v.gnc.g.p_dej_n.traj_rate = sim.traj_rate;  % guidance integration rate (Hz)
in0.s.data_rate = sim.data_rate;          % Hz, data storage frequency
% overall guidance settings
in0.v.gnc.n.p.rate = sim.traj_rate;   % navigation rate = traj rate
in0.v.gnc.c.p.rate = sim.traj_rate;    % control rate = traj rate
in0.v.gnc.g.p.prop_mode = uint8(0);
in0.v.gnc.g.p.dm_mode = uint8(gnc.mode);

switch (gnc.mode)
    case 3 %pdg
        mode = 'pd';
        in0.v.gnc.g.pd.n_jett = uint8(gnc.n);
        in0.v.gnc.g.pd.ha_tgt = gnc.ha_tgt*1000;
        in0.v.gnc.g.pd.ha_tol = gnc.ha_tol*1000;
        in0.v.gnc.g.pd.cds(1:stages) = aero.cds;
        in0.v.gnc.g.pd.cls(1:stages) = aero.cls;
        in0.v.gnc.g.pd.area_refs(1:stages) = transpose(areas); % m^2
        in0.v.gnc.g.pd.masses(1:stages) = aero.m'; % kg
%         in0.v.gnc.g.pd.trigger = uint8(0); %energy trigger
        in0.v.gnc.g.pd.rate = in0.s.traj.rate;    % predictor integration rate
        in0.v.gnc.g.pd.t_max = sim.t_max;
        
        in0.v.gnc.g.pd.sigs.m = 0;
        in0.v.gnc.g.pd.sigs.cd = 0;
        in0.v.gnc.g.pd.sigs.aref = 0;
        
        in0.v.gnc.g.pd.bias = gnc.bias;
        in0.v.gnc.g.pd.stage_skip = gnc.stage_skip;
        
    case 4 % dej-n properties
        mode = 'dej_n';
        in0.v.gnc.g.p_dej_n.cds(1:stages) = aero.cds;
        in0.v.gnc.g.p_dej_n.cls(1:stages) = aero.cls;
        in0.v.gnc.g.p_dej_n.masses(1:stages) = aero.m'; % kg
        in0.v.gnc.g.p_dej_n.area_refs(1:length(areas)) = transpose(areas); % m^2
        in0.v.gnc.g.p_dej_n.t_max = sim.t_max;

        in0.v.gnc.g.p_dej_n.n_jett = uint8(gnc.n);
        in0.v.gnc.g.p_dej_n.iter_max = uint8(gnc.iters); % NPC iteration limit
        in0.v.gnc.g.p_dej_n.tj_ini(1:length(gnc.tj0)) = gnc.tj0;
        in0.v.gnc.g.p_dej_n.tgt_ap = gnc.ha_tgt*1000; % m, target apoapse
        in0.v.gnc.g.p_dej_n.tol_ap = gnc.ha_tol*1000; % m, apoapse tolerance
        in0.v.gnc.g.p_dej_n.npc_mode = uint8(gnc.npc_mode);   % root finding method (1 = newton, 0=bisection)
        in0.v.gnc.g.p_dej_n.pred_mode = uint8(gnc.pred_mode);   % predictor mode (0 normal, 1 verbose)
        in0.v.gnc.g.p_dej_n.t_init = gnc.t_init;

        in0.v.gnc.g.p_dej_n.multi.force_jett = gnc.force_jett;
        in0.v.gnc.g.p_dej_n.multi.comp_curr = gnc.comp_curr;

        in0.v.gnc.g.p_dej_n.hydra.dtj_lim = gnc.dtj_lim;   % max step size for newton method
        in0.v.gnc.g.p_dej_n.hydra.dtj_flag = gnc.hydra_flag;    % hydra dtj reduction on/off

        in0.v.gnc.g.p_dej_n.sigs.m = 0;
        in0.v.gnc.g.p_dej_n.sigs.cd = 0;
        in0.v.gnc.g.p_dej_n.sigs.aref = 0;
        
        in0.v.gnc.g.p_dej_n.bias = gnc.bias;
    case 6 %manual
        mode = 'manual';
        in0.v.gnc.g.p_manual.n_jett = uint8(gnc.n);
        in0.v.gnc.g.p_manual.masses(1:stages) = aero.m'; % kg
        in0.v.gnc.g.p_manual.area_refs(1:length(areas)) = transpose(areas); % m^2
        in0.v.gnc.g.p_manual.t_jetts(1:gnc.n) = gnc.tj0';
        in0.v.gnc.c.p.rate = sim.traj_rate;
end

%% check efpa solving
if (sim.efpa_flag)
    fpa_tol = 0.05; %deg
%     [fpa_f,~,haf_err,~] = fpa_given_haf(x0.fpa0, gnc.ha_tgt, in0, fpa_tol, gnc.ha_tol);    
    fpa_bounds = [x0.fpa0/1.5, x0.fpa0*1.5];
    [fpa_f, ~, ~, ~, ~] = find_fpa_given_ha(fpa_bounds, gnc.ha_tgt, in0, fpa_tol, gnc.ha_tol); 
    if (isnan(fpa_f))
        fprintf('WARNING: no solution for EFPA found. Reverting to initial guess\n');
        fpa_f = x0.fpa0;
    end
    
    in0.s.traj.gamma_pp = fpa_f*pi/180;
    % re-do inertial entry state
    [in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
        in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
        in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p); 
end %solve for efpa


if (gnc.npc_mode == uint8(3))
    % ensemble -> regular for nominal
    in0.v.gnc.g.p_dej_n.npc_mode = uint8(1);
end

if ((~sim.ignore_nominal) && (mc.debug == false || mc.flag == false))
%% Run nominal Trajectory
if (sim.debug)
    fprintf('Running in debug mode...\n'); out = main1(in0);         % debug mode
else
    out = main1_mex(in0);
end

% 
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
    idj = find(out.g.(mode).stage == i,1);
    if ~isempty(idj)
        if (idxend > idj)
            out.idj(i) = idj;
        else
            out.idj(i) = out.idxend;
        end
        out.t_jett(i) = out.traj.time(out.idj(i));
    else
        out.t_jett(i) = nan;
        out.idj(i) = nan;
    end
end
out.betas = round(beta,3);
out.beta_ratios = round(beta_ratio,3);


% if oe(1) < 0
%     fprintf('WARNING: Exiting Orbit is Hyperbolic\n');
% end

else % ignore nominal
    out = struct('efpa', in0.s.traj.gamma_pp*180/pi);
end %check debugging
% check monte carlo

if (mc.flag)
    out_nom = out;
    clear out;  % delete nominal trajectory
    out.nom = out_nom;  % init new output struct
    
    % navigation settings
    if (gnc.nav.rate == 0)
        gnc.nav.rate = 20;
    end
    switch (gnc.nav.mode)
        case 2 %ecrvs
            in0.v.gnc.n.p.mode = uint8(gnc.nav.mode); % Use Markov process/ECRV error model
            in0.v.gnc.n.p.rate = gnc.nav.rate; % Hz
            in0.v.gnc.n.p.seed = uint32(gnc.nav.seed); % nd, Error model seed
            in0.v.gnc.n.p.ecrv.tau = gnc.nav.tau; % s, time constant
            in0.v.gnc.n.p.atm.omega = in0.p.omega; % rad/s, planet angular velocity vector
            in0.v.gnc.n.p.atm.r_e = in0.p.r_e; % m, planet equatorial radius
            in0.v.gnc.n.p.atm.r_p = in0.p.r_p; % m, planet polar radius
            in0.v.gnc.n.p.ecrv.P_SS = gnc.nav.P_SS;
            in0.v.gnc.n.p.ecrv.P_SS(7:9,7:9) = eye(3)*(90*in0.c.earth_g*1e-6 / 3)^2; % m/s^2, 90 micro-g, 3-sigma
        
        case 3  %exponential decay
            in0.v.gnc.n.p.mode = uint8(gnc.nav.mode); % Use Markov process/ECRV error model
            in0.v.gnc.n.p.rate = gnc.nav.rate; % Hz
            in0.v.gnc.n.p.seed = uint32(gnc.nav.seed); % nd, Error model seed
            in0.v.gnc.n.p.atm.omega = in0.p.omega; % rad/s, planet angular velocity vector
            in0.v.gnc.n.p.atm.r_e = in0.p.r_e; % m, planet equatorial radius
            in0.v.gnc.n.p.atm.r_p = in0.p.r_p; % m, planet polar radius
            in0.v.gnc.n.p.n_exp.tau_r = gnc.nav.tau_r;
            in0.v.gnc.n.p.n_exp.tau_v = gnc.nav.tau_v;
            in0.v.gnc.n.p.n_exp.tau_a = gnc.nav.tau_a;
            in0.v.gnc.n.p.n_exp.P0 = gnc.nav.P_SS;
        
        case 5 %exp-decay ecrv
            in0.v.gnc.n.p.mode = uint8(gnc.nav.mode);
            in0.v.gnc.n.p.rate = gnc.nav.rate;
            in0.v.gnc.n.p.seed = uint32(gnc.nav.seed);
            in0.v.gnc.n.p.atm.omega = in0.p.omega; % rad/s, planet angular velocity vector
            in0.v.gnc.n.p.atm.r_e = in0.p.r_e; % m, planet equatorial radius
            in0.v.gnc.n.p.atm.r_p = in0.p.r_p; % m, planet polar radius
            % ecrv data
            in0.v.gnc.n.p.ecrv.tau = gnc.nav.tau_ecrv; % s, ecrv time constant
            % n_exp data
            in0.v.gnc.n.p.n_exp.tau_r = gnc.nav.tau_r;
            in0.v.gnc.n.p.n_exp.tau_v = gnc.nav.tau_v;
            in0.v.gnc.n.p.n_exp.tau_a = gnc.nav.tau_a;
            in0.v.gnc.n.p.n_exp.P0 = gnc.nav.P_SS;
%             in0.v.gnc.n.p.n_exp.ercv0 = ...            
%                 [normrnd(0,gnc.nav.ercv0(1:3),[3,1]); ...      %r_atm uncertainty
%                 normrnd(0,gnc.nav.ercv0(4:6),[3,1]); ...   %v_atm uncertainty
%                 normrnd(0,gnc.nav.ercv0(7:9),[3,1])];
    end
    
    in0.v.gnc.g.p_dej_n.npc_mode = uint8(gnc.npc_mode);
    
    in0.v.gnc.g.p.atm.Kflag = mc.flag;
    in0.v.gnc.g.p.atm.mode = uint8(gnc.atm_mode); % atmospheric capture mode (1-factor, 2-interp, 3-ens filter)
    
    % prellocate data storage variables
    if (~isempty(mc.mcIndex) && mc.N ~= length(mc.mcIndex))
        N = length(mc.mcIndex);
    else
        N = mc.N;
    end
    
    len = in0.s.traj.t_ini:(1/in0.s.data_rate):in0.s.traj.t_max;
    len = length(len) + gnc.n;  % 1 extra per jettison event, and 1 for terminal point
    vmag = nan(len,N);
    alt = vmag; fpa = vmag; t = vmag; 
    rho = vmag;
    t_jett = nan(mc.N,5);
    tjr = t_jett;
    idj = t_jett;
    ha = vmag; ha_err = vmag; 
    rva_err = nan(len,9,N); ercv = rva_err; x_ercv = ercv;
%     ncalls = vmag; 
    tj_curr = vmag; K_dens = vmag; K_true = vmag; rho_est = vmag;
    atm_err = vmag;
    rss_nom = nan(mc.N,1);
    rss_K = vmag; rss_ecf = vmag;
    ind_curr = vmag;
    
    % more preallocation
    haf = nan(N,1);
    haf_err = haf; dv = haf;
    idend = haf;
    
%     tjett = nan;
    sigs = mc.sigs;
    n = gnc.n;
    
    debug = mc.debug;
    ha_tgt = gnc.ha_tgt;
    
    rnd = nan(size(mc_atm,2),1);
    % pre-determine random atmospheres to reduce broadcast overhead
    temp = nan(size(mc_atm,2), 7, N);
    if isempty(mc.mcIndex)
        for i = 1:N
            rnd(i) = randi(size(mc_atm,2));
            temp(:,:,i) = mc_atm(rnd(i)).table; %random disperse atmosphere
    %         temp(:,:,i) = mc_atm(i).table;    %atmospheres in order
        end
    else
        rnd(1:length(mc.mcIndex)) = mc.mcIndex;
        for i = 1:length(mc.mcIndex)
            temp(:,:,i) = mc_atm(rnd(i)).table;
        end
    end

    if (sim.parMode == false || mc.debug == true || ~isempty(mc.mcIndex))
        parArg = 0;
    else
        parArg = sim.nWorkers;
        par = parcluster('local');
        par.NumWorkers = sim.nWorkers;
    end
    
    % nav error values
    nav_mode = gnc.nav.mode;
    if (nav_mode == uint8(5))
        if (isempty(mc.mcIndex))
            len = N;
        else
            len = length(mc.mcIndex);
        end
            
        ecrv0 = zeros(9,len);
        for i = 1:len
            ecrv0(:,i) = [normrnd(0,gnc.nav.ercv0(1:3),[3,1]); ...      %r_atm uncertainty
                normrnd(0,gnc.nav.ercv0(4:6),[3,1]); ...   %v_atm uncertainty
                normrnd(0,gnc.nav.ercv0(7:9),[3,1])];
        end
    end
        
    parfor (i = 1:N, parArg)  % run in parallel processing (need parallel computing toolbox)
%     for i = 1:N
        in_mc = in0;        % copy nominal input struct

        in_mc.p.atm.table = temp(:,:,i);    % new atmospheric table
        in_mc.v.gnc.g.p.planet.atm_true = temp(:,:,i);
%         in_mc.v.gnc.g.p_dej_n.atm.mc_ind = rnd(i);  % atmosphere index
        in_mc.v.gnc.g.p.atm.mc_ind = rnd(i);
        
        % nav error
        switch (nav_mode)
            case 5 %ecrv hybrid
                in_mc.v.gnc.n.p.n_exp.ercv0 = ecrv0(:,i);
        end
        
        % apply input perturbations (only to traj. guid struct has nominal)
        in_mc = apply_dispersions(in_mc,sigs);
        
        % FOR TESTING ONLY
%         in_mc.p.atm.mc = mc_atm;
        
        %%% Run Monte Carlo %%%
        if (debug && parArg == 0)
            fprintf('Running MC in debug mode...\n'); out_mc = main1(in_mc);    % debugging purposes
        else
            out_mc = main1_mex(in_mc);  % output struct for disperse case
        end
        
        % get final index
        if isnan(out_mc.traj.alt(end))
            idxend = find(isnan(out_mc.traj.alt),1)-1;
        else
            idxend = length(out_mc.traj.alt);
        end
        
        vmag(:,i) = out_mc.traj.vel_pp_mag./1000;   %km/s
        alt(:,i) = out_mc.traj.alt./1000;       %km
        fpa(:,i) = out_mc.traj.gamma_pp.*180/pi;    %deg
        t(:,i) = out_mc.traj.time;                  %s
        rho(:,i) = out_mc.traj.rho;     % kg/m3
        rva_err(:,:,i) = out_mc.nav.rva_error;
%         ercv(:,:,i) = out_mc.nav.ercv;
        x_ercv(:,:,i) = out_mc.nav.x_ercv;
        
%         out.traj(i).vmag = out_mc.traj.vel_pp_mag./1000;
%         out.traj(i).alt = out_mc.traj.alt./1000;
%         out.traj(i).fpa = out_mc.traj.gamma_pp;
%         out.traj(i).t = out_mc.traj.time;
%         out.traj(i).rho = out_mc.traj.rho;
        
        for j = 1:n
            tjett = (find(out_mc.g.(mode).stage == j,1)-1)/(in0.s.data_rate);
            if isempty(tjett)
                t_jett(i,j) = nan;
                tjr(i,j) = nan;
                idj(i,j) = nan;
            else
                t_jett(i,j) = tjett;
                tjr(i,j) = tjett / out_mc.traj.time(idxend);
                idjTemp = find(t(:,i) >= t_jett(i,j),1);    % check if jettison final time step
                if isempty(idjTemp)
                    idj(i,j) = nan;
                else
                    idj(i,j) = idjTemp;
                end
            end
        end
       
        ha_err(:,i) = out_mc.g.(mode).dr_ap./1000;  %km
        K_dens(:,i) = out_mc.g.K_dens;  % estimated density variation (nd)
        K_true(:,i) = out_mc.g.K_true;
        
        rho_est(:,i) = out_mc.g.rho_est;    % estimated density (kg/m3)
        switch (mode)
            case 'dej_n'
                ha(:,i) = out_mc.g.dej_n.r_ap./1000;
        %         ncalls(:,i) = out_mc.g.dej_n.ncalls;
                tj_curr(:,i) = out_mc.g.dej_n.tj;
        %         atm_err(:,i) = out_mc.g.dej_n.atm_err;
            case 'pd'
                
            case 'manual'
        end

        % atmospheric estimation stuff
        rss_nom(i) = out_mc.g.rss_nom;
        rss_K(:,i) = out_mc.g.rss_K;
        rss_ecf(:,i) = out_mc.g.rss_ens;
        ind_curr(:,i) = out_mc.g.ind_curr;
        
        % calculate apoapsis
        rv = [out_mc.traj.pos_ii(idxend,:), out_mc.traj.vel_ii(idxend,:)]';
        oe = rv2oe(rv,in0.p.mu);
        ra = oe(1)*(1+oe(2));
        err = (ra - in0.p.r_e)/1000;
        
        haf(i) = round(err,5);       %km, final apoapsis
        haf_err(i) = round((err - ha_tgt),5);    %km, apoapsis error
        dv(i) = round( out_mc.traj.vel_ii_mag(1) - out_mc.traj.vel_ii_mag(idxend), 5 );
        idend(i) = idxend;
        
%         if (gnc.atm_mode == uint8(1))
%             out.atm_model(:,:,i) = out_mc.g.atm_est;
%         else
%             out.atm_model = nan;
%         end
        
        
    end % monte carlo parfor loop
%     toc

    % store all monte carlo results into out struct
    out.traj.vmag = vmag;
    out.traj.alt = alt;
    out.traj.fpa = fpa;
    out.traj.t = t;
    out.traj.rho = rho;
    
    out.g.ha_err = ha_err;
    out.g.K_dens = K_dens;  %density scale factor
    out.g.K_true = K_true;  %true density dispersion factor
    out.g.rho_est = rho_est;    % density estimate
    out.g.mc_inds = rnd;
%     out.g.atm_err = atm_err;    % atmospheric density error
    switch (mode)
        case 'dej_n'
            out.g.ha = ha;
        %     out.g.ncalls = ncalls;
            out.g.tj_curr = tj_curr;

            % atmospheric estimation
            out.g.atm.rss_nom = rss_nom;    % nominal RSS error
            out.g.atm.rss_K = rss_K;        % scale factor RSS error
            out.g.atm.rss_ecf = rss_ecf;    % ensemble filter RSS error
            out.g.atm.ind_curr = ind_curr;  % ecf index
        case 'pd'
            out.g.dr_ap = ha_err;
        case 'manual'
    end
    
    out.haf = haf;
    out.haf_err = haf_err;
    out.dv = dv;
    out.tjett = t_jett;
    out.tjr = tjr;
    out.idj = idj;
    out.idxend = idend;
    
    % nav
    out.g.nav.ecrv0 = ecrv0;
%     out.g.nav.rva_err = rva_err;
%     out.g.nav.ercv = ercv;
%     out.g.nav.x_ercv = x_ercv;
    
    
end % monte carlo check

end % dej_n_earth_auto
