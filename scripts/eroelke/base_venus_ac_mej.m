
function [x0,aero,gnc,sim, mc] = base_venus_ac_mej(b21, isMc)
[x0,aero,gnc,sim, mc] = base_venus_ac(isMc);

b31 = 10;

% entry
x0.v0 = 11;
x0.fpa0 = -5.4;

% guidance
gnc.ha_tgt = 2000;
gnc.ha_tol = 5;
gnc.n = 2;
gnc.tj0 = [100 170];
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);    %MEJ
gnc.dtj_lim = 10;
gnc.force_jett = true;  % unconditional jettison
gnc.comp_curr = false;  % only used for conditional jettison

% mc
mc.debug = false;
mc.flag = logical(isMc);
mc.N = 1000;
    
% aerp
aero.m = [150 60 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * b21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * b31) );

aero.rn = 0.1;
aero.cds = [1.05 1.05 1.05];
aero.cls = [0 0 0];

sim.traj_rate = 100;
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;
if (isMc)
    sim.parMode = true;
end

end % base_venus_ac_mej.m