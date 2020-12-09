function [x0,aero,gnc,sim,mc] = base_venus_ac(isMc)
[x0,aero,gnc,sim,mc] = base_ac(isMc);
x0.pci = 0;
x0.fpa0 = -5.4;
x0.v0 = 11;
x0.h0 = 150;
aero.cds = [1.05 1.05];

gnc.tj0 = 100;
gnc.guid_rate = 0.5;
gnc.force_jett = true;
gnc.iters = uint8(1);
gnc.mode = 4;   %dej_n
gnc.stage_skip = 1000;
gnc.atm_mode = 0;   %density factor

sim.planet = 'venus';
sim.atm_mode = 3;
sim.h_max = 150;
sim.h_min = 50;
sim.t_max = 1500;
mc.flag = true;
mc.N = 1000;
sim.data_rate = 1;
sim.ignore_nominal = false;

% beta ratio = 10
aero.m = [150 60]; 
aero.rcs = [1 0.2];
aero.rn = 0.175/2;
aero.cds = [1.05 1.05];

end