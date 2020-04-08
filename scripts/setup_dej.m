function [x0,aero,gnc,sim,mc] = setup_dej()
x0.fpa0 = -5.7; % planet-relative EFPA
x0.v0 = 12;     % planet-relative velocity
x0.h0 = 125;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

gnc.ha_tgt = 2000;
gnc.ha_tol = 10;
gnc.n = 1;
gnc.tj0 = 100;
gnc.dtj_lim = 30;
gnc.npc_mode = uint8(1);
gnc.guid_rate = 0.1;
gnc.iters = uint8(1);
gnc.atm_mode = uint8(0);
gnc.t_init = 0;
gnc.p_tol = 0;

aero.rcs = [0.75 0.175];
aero.m = [80 40];
aero.rn = 0.175/2;
aero.cds = [1.05 1.05];
aero.cls = [0 0];

mc.flag = false;
mc.N = 0;
mc.sigs = nan;
mc.debug = false;
mc.sigs = default_sigs();

sim.traj_rate = 250;
sim.data_rate = 1;
sim.t_max = 1000;
sim.h_max = x0.h0;
sim.h_min = 40;
sim.planet = 'earth';
sim.efpa_flag = true;
end