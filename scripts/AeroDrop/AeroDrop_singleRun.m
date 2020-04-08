%% AeroDrop_main.m
% Sam Albert, EsDL
% Created: 10.4.19

% Simlates flight via run_dej_n()

clear; clc; close all;
tic

%% Find generic CD
% Assume MSL-esque aeroshell, back out MSL CD
m_MSL = 3700; % kg
d_MSL = 4.5; % m, heatshield diameter
BC_MSL = 145; % kg/m^2, ballistic coefficient

A_MSL = pi * (d_MSL)^2; % area of heatshield

CD_MSL = m_MSL / (BC_MSL * A_MSL);

%% Get v0 given vinf
vinf = 0; % km/s, user input
% Set planet params (hardcoded!)
radE = 6378.1363; % km
muE = 3.986004415e5; % km^3/s^2
r_atm = 125; % define atmospheric interface altitude

v_init = sqrt(2*muE/(radE + r_atm) + vinf^2);
v_esc = sqrt(2*muE/(radE + r_atm));

%% Generic setup params
[x0,aero,gnc,sim, mc] = setup_dej();

x0.fpa0 = -5.6; % deg - planet relative! (gamma_pp)
x0.v0 = v_init; % km/s - planet relative! (vel_pp_mag) (v_init set above)

gnc.ha_tgt = 0; % has no effect with gnc off
gnc.n = 0; % 0 jettison events
% gnc.tj0 = [80 120];
gnc.guid_rate = 0;

aero.m = 3700; % MSL-esque mass. Only BC really matters for ballistic
aero.rn = 1.1; 
aero.cds = CD_MSL; % (CD_MSL found above, used as rough guess)
aero.cls = 0;
% Get rc and rn given desired BC (below only works for gnc.n = 0!!!)
BC = 145;
aero.rcs = sqrt(aero.m / (aero.cds * pi * BC));
aero.rn = aero.rcs / 2;

mc.flag = false; % no Monte Carlo (default)

sim.traj_rate = 250; % (default: 250)
sim.planet = 'earth'; % (default: earth)
sim.efpa_flag = false; % manually set efpa

%% Run simulation
out = run_dej_n(x0,gnc,aero,sim,mc);

%% Plots
% Choose plots to include
set.vel_vs_alt = true;
set.alt_vs_time = true;
set.vel_vs_time = true;
set.vel_vs_gs = true;
set.vel_vs_heatRate = true;

AeroDrop_plots(out,set);

%% Analysis
% TODO: set up constraints for apoapsis, heat rate, heat load, g's

% Plot v_esc constraint
figure(1);
plot([v_esc v_esc],[70 125],'r--','LineWidth',1.5);







toc