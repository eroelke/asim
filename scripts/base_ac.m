% E. Roelke, April 2019
function [x0,aero,gnc,sim,mc] = base_ac(isMc)
    x0.fpa0 = -5.9;
    x0.v0 = 12;
    x0.h0 = 125;
    x0.lat0 = 0;
    x0.lon0 = 0;
    x0.az0 = 90;
    
    sim.traj_rate = 100;
    sim.data_rate = 1;
    sim.t_max = 800;
    sim.h_max = x0.h0;
    sim.h_min = 50;
    sim.planet = 'earth';
    sim.atm_mode = uint8(3);   %atm look up with winds
    sim.efpa_flag = false;
    sim.debug = false;
    sim.parMode = false;    % don't run in parallel processing mode
    sim.nWorkers = 10;  %if you do, use 10 workers

    gnc.ha_tgt = 2000;
    gnc.ha_tol = 0.5;
    gnc.n = 1;
    gnc.tj0 = 120;
    gnc.dtj_lim = 30;
    gnc.npc_mode = uint8(1);   %newton method
    gnc.pred_mode = uint8(0);   %regular predictor
    gnc.guid_rate = 0.5;
    gnc.atm_mode = uint8(1);    %density factor
    gnc.iters = uint8(1);
    gnc.t_init = 0;
    gnc.hydra_flag = true;   %hybrid newton-bisection
    gnc.force_jett = false;   % don't force jettison in multi case
    gnc.comp_curr = false;   %compare secondary stage error to only last stage
    gnc.rss_flag = false;   %no atm estimation rss stuff
    gnc.ecf.mode = uint8(0);
    gnc.ecf.p_tol = 0.05;
    
    %navigation
    gnc.nav.mode = 1;   %monte carlo nav mode - perfect nav
    gnc.nav.P_SS = zeros(9,9);
    gnc.nav.rate = sim.traj_rate;
    gnc.nav.seed = 0;
    gnc.nav.tau = 0;

    
    aero.rcs = [0.75 0.175];
    aero.m = [80 40];
    aero.rn = 0.175/2;
    aero.cds = [1.05 1.05];
    aero.cls = [0 0];

    if (isMc)
        mc.flag = true;
        mc.Kflag = true;
        mc.N = 1;
        mc.sigs = default_sigs();
        mc.sigs.cd = .05 * aero.cds(1) / 3;  %5% 3sigma
        mc.debug = false;
    else
        mc.flag = false;
        mc.Kflag = false;
        mc.N = 0;
        mc.sigs = sigs_zero();
        mc.debug = false;
    end
    mc.mcIndex = [];   %force specific monte carlo index/indices
    
    for i = 1:5
        gnc.bias(i,1).tgt_ap = 0;
        gnc.bias(i,1).tj = 0;
    end
end