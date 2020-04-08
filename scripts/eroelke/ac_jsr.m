%% SciTech -> Journal Article Simulations
% Summer 2019
clear; clc;
%% obtain beta ratio
% old vals -> beta ratio 9.18
% m1 = 80; m2 = 40; cd = 1.05; rc1 = 0.75; rc2 = 0.175;
% a1 = pi*rc1^2; a2 = pi*rc2^2;
% b1 = m1/(cd*a1); b2 = m2/(cd*a2);
% % ratios
% rc_ratio = rc2/rc1;
% m_ratio = m2/m1;

% new vals -> beta ratio 10
ratio = 10;
m1 = 150;
r = [1 0.2];
cd = 1.05;
m2 = ratio * m1 * (r(2)/r(1))^2;
b1 = m1/(cd*pi*r(1)^2);
b2 = m2/(cd*pi*r(2)^2);


%% DCF stuff
[x0,aero,gnc,sim,mc] = base_jsr();


%% PDG stuff


%% NPC - bisection


%% NPC - newton
