%% E. Roelke
%
clear;clc

%% corridor width, b2/b1 vs. vel trade study
v0s = 25:1:30;  %km/s
Ni = length(v0s);
% h0 = 1000e3;    %km, atmospheric interface
h0 = (25755603.4663544 - 2.4764e7)/1000;
ratios = 2:2:16;
Nj = length(ratios);
ha_tgt = 430000;    %km
fpa0 = -2.050066770663301e-01 * 180/pi;   %deg
beta1 = 200;
rn = 1;

% make up vals
cd = 1.23;
rcs = [2 0.5];
m1 = beta1 * cd * pi * rcs(1)^2;
m2 = nan(Nj,1);

fpa_guess = [-1 -20];
in = define_in;
in.s.traj.alt = h0;
in.s.traj.gamma_pp = fpa0;
in.s.traj.rate = 100;
in.s.traj.t_max = 5000;

in.p = get_planetary_data(8, 2);

in.v.aero.cd = cd;
in.v.aero.area_ref = pi * rcs(1)^2;
in.v.mp.m_ini = m1;
in.v.aero.nose_radius = rn;
in.v.aero.mode = uint8(1);

tol_fpa = 0.05;
tol_ha = 10;

% basic trajectory
[x0,aero,gnc,sim,mc] = base_jsr();
aero.rcs = rcs(1);
aero.rn = rn;
aero.cds = cd;
aero.cls = 0;
% aero.m = [m1 10*m1*(rcs(2)/rcs(1))^2];
aero.m = [m1];
gnc.guid_rate = 0;
gnc.ha_tgt = ha_tgt;
sim.t_max = 10000;
sim.h_max = 1000;
x0.fpa0 = fpa0;
x0.h0 = h0;
x0.v0 = 31.41749062460320;  %km/s
% x0.v0 = v0s(1);
mc.flag = false;
sim.planet = 'neptune';
sim.atm_mode = 1;

x0.pci = [25755603.4663544, 0, 0,  ... 
          -5423.42400674429, 30948.3904596519, 0];

out = run_dej_n(x0,gnc,aero,sim,mc);

plot(out.traj.vel_pp_mag/1000, out.traj.alt/1000)
keyboard;


for j = 1:Nj
    m2(j) = ratios(j) * m1 * (rcs(2)/rcs(1))^2;
    beta2(j) = m2(j)/(pi * cd * rcs(2)^2);
    
    
    for i = 1:Ni
        in.s.traj.vel_pp_mag = v0s(i);
        
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);
        
        in1 = in;
        in2 = in;
        
        in2.v.mp.m_ini = m2(j);
        in2.v.aero.area_ref = pi * rcs(2)^2;
        
        [fpa_corr(j, i, 1), ~, ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in1, tol_fpa, tol_ha);
        [fpa_corr(j, i, 2), ~, ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in2, tol_fpa, tol_ha);
    end
end




%% Casey asim input
clear; clc;
in = nepcap_asim();
out = main1_mex(in);

if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end
rv = [out.traj.pos_ii(idxend,:), out.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in.p.r_e)/1000;
dv = out.traj.vel_ii_mag(1) - out.traj.vel_ii_mag(idxend);
out.haf = round(haf,5);                   %km
% out.haf_err = round((haf - gnc.ha_tgt),5);    %km
out.dv = round(dv,5);
out.idxend = idxend;
out.xi = oe(1);

% casey comp
dat = load('scripts/eroelke/neptune/asim_out.mat');
dat = dat.aout;
if isnan(dat.traj.alt(end))
    idxend = find(isnan(dat.traj.alt),1)-1;
else
    idxend = length(dat.traj.alt);
end
rv = [dat.traj.pos_ii(idxend,:), dat.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in.p.r_e)/1000;
dv = dat.traj.vel_ii_mag(1) - dat.traj.vel_ii_mag(idxend);
dat.haf = round(haf,5);                   %km
% out.haf_err = round((haf - gnc.ha_tgt),5);    %km
dat.dv = round(dv,5);
dat.idxend = idxend;
dat.xi = oe(1);

figure(); hold on
plot(out.traj.vel_ii_mag./1000, out.traj.alt/1000,'k','LineWidth',2);
plot(dat.traj.vel_ii_mag./1000, dat.traj.alt/1000,'r','LineWidth',2);
legend('My ASIM','Your ASIM','location','ne')


