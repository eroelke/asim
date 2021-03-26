function [x0,aero,gnc,sim,mc] = base_earth_ac(isMc)
[x0,aero,gnc,sim,mc] = base_venus_ac(isMc);
x0.fpa0 = -5.8;
x0.v0 = 12;
x0.h0 = 125;

sim.planet = 'earth';
sim.h_max = x0.h0;
sim.h_min = 40;
sim.t_max = 1500;

end
