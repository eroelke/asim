% ac_venus_pd
%   predictive guidance for discrete-event drag modulated aerocapture at
%   venus
%
% Inputs:
%   efpa: entry flight path angle (deg)
%   ha_tgt: target apoapsis (km)
%   haf_tol: apoapsis tolerance (km)
%   traj_rate: trajectory integration rate
%   gnc_rate: guidance, navigation, control call rate
%   mc: monte carlo input struct
%       flag: monte carlo sim flag
%       N: number of independent runs
%       sigs: 1sigma standard deviations
%           [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
% 
% Outputs:
%   out: output data struct
% 
function out = ac_venus_pd(efpa,ha_tgt,haf_tol,m,rcs,traj_rate,gnc_rate,trigger,mc)

% initial simulation vals
in0 = pred_guid_venus_in;    % in0put struct for sin0gle jettison event

% set rates
in0.s.traj.rate = traj_rate; % trajectory integration rate (Hz)
in0.s.data_rate = in0.s.traj.rate;
in0.v.gnc.g.p.rate = gnc_rate;   % guidance call rate (Hz)
in0.v.gnc.n.p.rate = traj_rate;   %nav rate
in0.v.gnc.c.p.rate = traj_rate; %control rate

% Define Vals
% m1 = 80;
% m2 = 40;
% r1 = 0.75;
% r2 = 0.175;
m1 = m(1);
m2 = m(2);
r1 = rcs(1);
r2 = rcs(2);
rn = r1/2;
cd = 1.05;
area = pi.*[r1^2 r2^2];

% tgt_ap = 10000; % 10000km altitude target
% ha_tol = 50;    % 50 km tolerance

v_atm = 11e3;
alt_atm = 150e3;
% fpas = -5.4;

% mc.flag = uin0t8(0);
% npc = dej_n_auto_venus(fpa_atm, tgt_ap, ha_tol, 1, 85, [r1 r2], [m1 m2], 1,1,mc,1);
% tj_opt = npc.t_jett;    % 79.92s for this architeture

% Entry State
in0.s.traj.vel_pp_mag = v_atm;
in0.s.traj.alt = alt_atm;
in0.s.traj.gamma_pp = efpa*pi/180;

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = m1; % in0itial mass (kg)
in0.v.aero.area_ref = area(1); % Area before and after jettison (m^2)
in0.v.aero.mode = uint8(1);     % User-specified lift and drag coefficients
in0.v.aero.cd = cd;
in0.v.aero.cl = 0;
in0.v.aero.nose_radius = rn; % Nose radius (m)

% Defin0e guidance parameters
in0.v.gnc.g.pd.ha_tgt = ha_tgt*1e3;
in0.v.gnc.g.pd.ha_tol = haf_tol*1e3;
in0.v.gnc.g.pd.cds(1:2) = cd;
in0.v.gnc.g.pd.area_refs(1:2) = area;
in0.v.gnc.g.pd.masses(1:2) = [m1 m2];
in0.v.gnc.g.pd.trigger = trigger;
in0.v.gnc.g.pd.rate = in0.s.traj.rate;    % predictor integration rate

beta1 = in0.v.mp.m_ini/(in0.v.aero.area_ref*cd);
beta2 = in0.v.gnc.g.pd.masses(2)/(in0.v.gnc.g.pd.area_refs(2)*in0.v.aero.cd);
beta_ratio = beta2/beta1;

% run simulation
% fprintf('Running PDG in debug mode..\n'); out = main1(in0);     % debug mode
out = main1_mex(in0);

% get output parameters
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

% for j = 1:n
out.idj = find(out.g.pd.stage == 1,1);
if ~isempty(out.idj)
    out.tjett = out.traj.time(out.idj);
else
    out.idj = nan;
    out.tjett = nan;
end
% end

out.dv = round(dv,5);    % aerocapture maneuver delta v (m/s)
out.haf = round(haf/1000,5);    %km
out.haf_err = round(haf/1000 - ha_tgt,5);   % km
out.idxend = idxend;

if (mc.flag)    % monte carlo
    out_nom = out;
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
%     t_jett = nan(N,in0.v.gnc.g.pd.n_jett);
    t_jett = nan(N,1);
    idj = t_jett;
    idxend = nan(N,1);
    haf = nan(N,1);
    haf_err = haf; dv = haf;
    sigs = mc.sigs;
    flag = mc.flag;
    % pre-determine random atmospheres to reduce broadcast overhead
    temp = nan(1000, 7, N);
    for i = 1:N
        rnd = randi(1000);
        temp(:,:,i) = mc_atm(rnd).table;
    end
    
%     tic
    parfor (i = 1:N, 6)
%     for i = 1:N
        in_mc = in0;        % copy nominal input struct
        in_mc.v.gnc.g.pd.mcflag = flag;
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
        
%         for j = 1:in0.v.gnc.g.pd.n_jett
            jett_ind = (find(out_mc.g.pd.stage == 1,1)-1);
            if ~isempty(jett_ind)
                idj(i) = jett_ind;
                t_jett(i) = out_mc.traj.time(idj(i));
            else
                t_jett(i) = nan;    % never jettisoned
                idj(i) = nan;
            end
%         end
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
%     toc
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
    
end



end