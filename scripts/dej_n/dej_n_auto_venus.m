%% dej_n_auto_venus
%   -Venus drag modulation aerocapture trajectory.
%   -n stage discrete-event jettison (dej)
%
% Developed by:
%   Space Systems Design Lab, Georgia Tech
%   Colorado Center for Astrodynamics Research, CU Boulder
%
% Inputs:
%   fpa_atm:    entry flight path angle (deg) - [1x1]
%   ha_tgt:     target apoapsis altitude (km) - [1x1]
%   n:          number of jettison stages - [1x1]
%   tj0:        initial jettison time(s) guess (sec) - [nx1]
%   rcs:        successive backshell radii (m) - [n+1x1], default [0.75 0.25 0.175 0.1]
%   m:          successive vehicle masses (kg) - [n+1x1], default [80 50 40 30]
%   root_mode:  bisection method(0) or newton method(1)
%   guid_rate:  guidance rate (Hz)
%   mc: monte carlo input structure
%       mc.flag:    flag for monte carlo run (and density corrector)
%       mc.N:       number of monte carlo runs
%       mc.sigs:    [m (kg),cl,cd,aoa,v-atm (m/s),h_atm (m),fpa (rad),lat,lon,az]
%   iters: maximum guidance iterations
%   traj_rate: trajectory integration rate (Hz)
% 
% Input Notes:
%   n determines number of jettisons even if tj0,m,rcs are more than [nx1]
%   
% Outputs:
%   out:        simulation output data structure
%
% Major Revision History:
%   Jul 2018 - Created Function, Evan Roelke

function out = dej_n_auto_venus(fpa_atm, v0, ha_tgt, ha_tol, n, ... 
    tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate)
format longg

%% Define number of stages
stages = n+1;   % number of stages (including no jettison)

%% Check inputs for inconsistencies
% if any(isnan(rcs))==logical(true) || isempty(rcs)==logical(true)
%     R0 = 0.75;
%     rcs = [R0 0.25 0.175 0.1];
% end
% if any(isnan(m))==logical(true) || isempty(m)==logical(true)
%     m0 = 80;
%     m = [m0 50 40 30];    
% end
% 
% if length(tj0) > n
%     fprintf('Warning: tj0 is size [%ix1] instead of [%ix1]\n',(length(tj0)),n);
% elseif length(rcs) > stages
%     fprintf('Warning: rcs is size [%ix1] instead of [%ix1]\n',(length(rcs)),stages);
% elseif length(m) > stages
%     fprintf('Warning: m is size [%ix1] instead of [%ix1]\n',(length(m)),stages);
% end
% 
% if length(tj0) < n
%     error('Input tj0 is size [%ix1] instead of [%ix1]\n',(length(tj0)),n);
% elseif length(rcs) < stages
%     error('Input rcs is size [%ix1] instead of [%ix1]\n',(length(rcs)),stages);
% elseif length(m) < stages
%     error('Input m is size [%ix1] instead of [%ix1]\n',(length(m)),stages);
% end

%% Geometry and Aerodynamics
areas = pi.*rcs.^2;
rn = rcs(2)/2;     % Nose radius (m)
cd = 1.05;
cl = 0;

%% Base Input Structure
in0 = dej_n_venus_in; % nominal input structure

%% Get atmospheric profile
% temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_mc.mat');
% in0.p.atm.table = temp.nom_table;

temp = load('data/atm_data/atm_venus_mc.mat');
% temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_mc.mat');  % home
% temp = load('C:\Users\Evan Roelke\Documents\research\asim\data\atm_venus_mc.mat'); % lab
in0.p.atm.table = temp.nom_table;
% in0.p.atm.table = temp.table;
in0.v.gnc.g.p_dej_n.planet.atm_table = in0.p.atm.table(:,1:7);

%% Define Entry State
% in0.s.traj.gamma_pp = -5.5*pi/180;
in0.s.traj.gamma_pp = fpa_atm*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = v0*1000;

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

%% Vehicle properties
in0.v.mp.m_ini = m(1); % Initial mass (kg)
in0.v.mp.m_jettison = 0; % Jettisoned mass (kg) - unused for this guidance algorithm
in0.v.aero.area_ref = areas(1); % initial area difference before and after jettison (m^2)
in0.v.aero.mode = uint8(1); % User-specified, constant lift and drag coefficients
in0.v.aero.cd = cd;
in0.v.aero.cl = cl;
in0.v.aero.nose_radius = rn; % Nose radius (m)

%% N_DEJ Guidance Settings
max_alt = 150e3;    % maximum allowable altitude (m)
min_alt = 50e3;        % minimum allowable altitude (m)
in0.v.gnc.g.p_dej_n.planet.alt_min = min_alt; % m
in0.v.gnc.g.p_dej_n.planet.alt_max = max_alt; % m
in0.v.gnc.g.p_dej_n.cds(1:stages) = in0.v.aero.cd; % nd, constant cd
in0.v.gnc.g.p_dej_n.cls(1:stages) = in0.v.aero.cl; % nd
in0.v.gnc.g.p_dej_n.masses(1:length(m)) = m'; % kg
% in0.v.gnc.g.p_dej_n.mass_jettison = in0.v.mp.m_jettison;
in0.v.gnc.g.p_dej_n.area_refs(1:length(areas)) = transpose(areas); % m^2

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at top of atmosphere
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

% Simulation Rates
% in0.s.data_rate = 1;          % Hz, data storage frequency
in0.s.traj.rate = traj_rate;           % double, Hz, integration rate
in0.v.gnc.g.p.rate = guid_rate;         % nd, guidance call rate (Hz)
in0.v.gnc.g.p_dej_n.traj_rate = traj_rate;  % guidance integration rate (Hz)
in0.v.gnc.c.p.rate = traj_rate;    % control rate = guidance rate
in0.s.data_rate = 10;
in0.v.gnc.n.p.rate = traj_rate;   % navigation rate = traj rate

in0.v.gnc.g.p_dej_n.n_jett = uint8(n);

% guidance properties
in0.v.gnc.g.p.prop_mode = uint8(0);
in0.v.gnc.g.p.dm_mode = uint8(4);   % dej-n
in0.v.gnc.g.p_dej_n.iter_max = uint8(iters); % NPC iteration limit
% if length(tj0) ~= n
%     error('Number of jettison times not consistent with number of jettison stages')
% else
in0.v.gnc.g.p_dej_n.tj_ini(1:length(tj0)) = tj0;
% end

% in0.v.gnc.g.p_dej_n.K_gain = 0.2;
in0.v.gnc.g.p.K_gain = 0.2;
in0.v.gnc.g.p.K_bounds = [0.1 2];

% update dej properties
in0.v.gnc.g.p_dej_n.tgt_ap = ha_tgt*1000; % m, target apoapse
in0.v.gnc.g.p_dej_n.tol_ap = ha_tol*1000; % m, apoapse tolerance

in0.v.gnc.g.p_dej_n.npc_mode = uint8(root_mode);   % root finding method (1 = newton)
% in0.v.gnc.g.p_dej_n.K_flag = uint8(0);  % density corrector for monte carlo
in0.v.gnc.g.p.Kflag = uint8(0);

% Run nominal Trajectory
out = main1_mex(in0);
% out = main1(in0);         % debug mode
% 
% % Compute parameters
if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end

beta = nan(length(areas),1);
beta_ratio = nan(length(areas)-1,1);
for i = 1:length(areas)
    beta(i) = m(i)/(areas(i)*cd);
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
out.haf_err = round((haf - ha_tgt),5);    %km
out.dv = round(dv,5);
out.idxend = idxend;
for i = 1:n
    tjett = (find(out.g.dej_n.stage == i,1)-1)/(in0.s.data_rate);
    if ~isempty(tjett)
        out.t_jett(i) = tjett;
        out.idj(i) = find(out.traj.time >= out.t_jett(i),1);
    else
        out.t_jett(i) = nan;
        out.idj(i) = nan;
    end
end
out.betas = round(beta,3);
out.beta_ratios = round(beta_ratio,3);


% check monte carlo
if (mc.flag)
    out_nom = out;
    clear out;  % delete nominal trajectory
    out.nom = out_nom;  % nominal trajectory
    
%     in0.v.gnc.g.p_dej_n.K_flag = uint8(mc.flag);    % turn on density corrector
    in0.v.gnc.g.p.Kflag = uint8(mc.flag);
%     in0.v.gnc.g.p.K_gain = in0.v.gnc.g.p_dej_n.K_gain;
%     in0.v.gnc.g.p.K_bounds = in0.v.gnc.g.p_dej_n.K_bounds;
    % load disperse atmospheric tables
    mc_atm = load('./data/atm_data/atm_venus_mc.mat');
    mc_atm = mc_atm.mc; % 1000 atm tables
    

    % preallocate vars
    len = in0.s.traj.t_ini:(1/in0.s.data_rate):in0.s.traj.t_max;
    len = length(len)+1;
    vmag = nan(len,mc.N);
    alt = vmag; fpa = vmag; t = vmag; 
    rho = vmag;
    t_jett = nan(mc.N,5);
    idj = t_jett;
    ha = vmag; ha_err = vmag; ncalls = vmag; 
    tj_curr = vmag; K_dens = vmag;
    haf = nan(mc.N,1);
    haf_err = haf; dv = haf;
%     tjett = nan;
    N = mc.N;
    sigs = mc.sigs;
    % pre-determine random atmospheres to reduce broadcast overhead
    temp = nan(1000, 7, N);
    for i = 1:N
        rnd = randi(1000);
        temp(:,:,i) = mc_atm(rnd).table;
    end
%     tic
    parfor (i = 1:N, 10)
% for i = 1:N
        in_mc = in0;        % copy nominal input struct
%         if mod(i,N/10) == 0
%             fprintf('Monte Carlo %i percent done\n',i*100/N);
%         end
        % apply disperse atmosphere
        
        in_mc.p.atm.table = temp(:,:,i);    % new atmospheric table
        
%         figure(); hold on
%         plot(in0.p.atm.table(:,2),in0.p.atm.table(:,1)./1000)
%         plot(in_mc.p.atm.table(:,2),in_mc.p.atm.table(:,1)./1000)
%         ylim([90 150])
        
        % apply input perturbations (only to traj. guid struct has nominal)
        in_mc = apply_dispersions(in_mc,sigs);
        
%         out_mc = main1(in_mc);    % debugging purposes
        out_mc = main1_mex(in_mc);  % output struct for disperse case
        
        vmag(:,i) = out_mc.traj.vel_pp_mag./1000;   %km/s
        alt(:,i) = out_mc.traj.alt./1000;       %km
        fpa(:,i) = out_mc.traj.gamma_pp.*180/pi;    %deg
        t(:,i) = out_mc.traj.time;                  %s
        rho(:,i) = out_mc.traj.rho;     % kg/m3
        
%         out.traj(i).vmag = out_mc.traj.vel_pp_mag./1000;
%         out.traj(i).alt = out_mc.traj.alt./1000;
%         out.traj(i).fpa = out_mc.traj.gamma_pp;
%         out.traj(i).t = out_mc.traj.time;
%         out.traj(i).rho = out_mc.traj.rho;
        
        for j = 1:n
            tjett = (find(out_mc.g.dej_n.stage == j,1)-1)/(in0.s.data_rate);
            if ~isempty(tjett)
                t_jett(i,j) = tjett;
                idj(i,j) = find(t(:,i) >= t_jett(i,j),1);
            else
                t_jett(i,j) = nan;
                idj(i,j) = nan;
            end
        end
       
        ha(:,i) = out_mc.g.dej_n.r_ap./1000;
        ha_err(:,i) = out_mc.g.dej_n.dr_ap./1000;
        ncalls(:,i) = out_mc.g.dej_n.ncalls;
        tj_curr(:,i) = out_mc.g.dej_n.tj;
%         K_dens(:,i) = out_mc.g.dej_n.K_dens;
        K_dens(:,i) = out_mc.g.K_dens;
        
        % get final index
        if isnan(out_mc.traj.alt(end))
            idxend = find(isnan(out_mc.traj.alt),1)-1;
        else
            idxend = length(out_mc.traj.alt);
        end
        
        % calculate apoapsis
        rv = [out_mc.traj.pos_ii(idxend,:), out_mc.traj.vel_ii(idxend,:)]';
        oe = rv2oe(rv,in0.p.mu);
        ra = oe(1)*(1+oe(2));
        err = (ra - in0.p.r_e)/1000;
        
        haf(i) = round(err,5);       %km, final apoapsis
        haf_err(i) = round((err - ha_tgt),5);    %km, apoapsis error
        dv(i) = round( out_mc.traj.vel_ii_mag(1) - out_mc.traj.vel_ii_mag(idxend), 5 );
    end % monte carlo loop
%     toc
    % store all monte carlo results into out struct
    out.traj.vmag = vmag;
    out.traj.alt = alt;
    out.traj.fpa = fpa;
    out.traj.t = t;
    out.traj.rho = rho;
    out.g.ha = ha;
    out.g.ha_err = ha_err;
    out.g.ncalls = ncalls;
    out.g.tj_curr = tj_curr;
    out.g.K_dens = K_dens;
    out.haf = haf;
    out.haf_err = haf_err;
    out.dv = dv;
    out.tjett = t_jett;
    out.idj = idj;
    
    
end % check mc flag


end % run_ac
