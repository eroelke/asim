%% Evan Roelke
% Earth Aerocapture Simulations
% 
% clear;clc

%% Basic Lifting Sim
clear; clc;
cd = 1;
beta1 = 80;
b21 = 12;
plotTraj = true;

[x0,aero,gnc,mc,sim] = BaseInputs();
x0.v0 = 13;
x0.fpa0 = -5.6;
mc.flag = false;
x0.az0 = 0;
aero.m = [250 75];
aero.rcs(1) = sqrt(aero.m(1) / (cd * pi * beta1)); %1.231162606;
aero.rcs(2) = sqrt(aero.m(2) / (cd * pi * beta1*b21));
aero.cds(1:2) = [cd cd];
aero.cls(1:2) = [0 0];
sim.efpa_flag = false;
mc.debug = false;
gnc.atm_mode = uint8(1);    %exponential atm
gnc.guid_rate = 1;
gnc.mode = 4;   %dej
gnc.ha_tgt = 42164;
sim.traj_rate = 250;
sim.data_rate = 10;
sim.atm_mode = 1;
sim.t_max = 500;
gnc.tj0 = 100;
gnc.dtj_lim = 10;
sim.h_min = 25;

sim.debug = false;
out = run_dej_n(x0,gnc,aero,sim,mc);
fprintf('tj: %g s\ndh_a: %g km\n',out.t_jett, out.haf_err);

if (plotTraj == true)
figure(); hold on; grid on;
set(gca,'FontSize',14);
plot(out.traj.vel_pp_mag./1000,out.traj.alt./1000,'LineWidth',2)
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')
end

% keyboard

%% Trade Studies EFPA va Vatm - Jan 2019
%{
clear; clc; close all;
betas = [10 25 50 100 150 200 250 500];
% betas = [25 50 100];
Ni = length(betas);

v_atms = linspace(11,16,10);
Nj = length(v_atms);

% vehicle
rc = 1.5;
aref = pi*rc^2;
cd = 1.05;

fpa0 = -5; %initial guess
fpa_tol = 0.025;

ha_tol = 25;        %km
% ha_tgt = 10000;      %km
% ha_tgt = 2000;  % high LEO
% ha_tgt = 35786; % ~GEO
% ha_tgt = 400; %km
% ha_tgt = 200; %km - low LEO

corr0 = [-12 -0.1];

in0 = dme_earth_in;
temp = load('./data/atm_earth_gram2007_1000mc.mat');
in0.p.atm.table = temp.nom_table;
in0.s.traj.alt = 125e3;
in0.s.traj.az = 90*pi/180;
in0.s.traj.gamma_pp = fpa0*pi/180;
in0.v.aero.cd = cd;
in0.v.aero.area_ref = aref;
in0.v.gnc.g.p.prop_mode = uint8(0); % remove guidance
in0.s.traj.rate = 100;
in0.v.gnc.g.p.rate = 0;

% preallocte
m = nan(Ni,1);
fpa_f = nan(Ni,Nj);

for i = 1:Ni
    m(i) = betas(i)*cd*aref;
    for j = 1:length(v_atms)
        in = in0;
        % update input struct
        in.s.traj.vel_pp_mag = v_atms(j)*1000;
        in.v.mp.m_ini = m(i);
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);
        % test sim
        % out = main1_mex(in);
        
%         [fpa_f(i,j),haf,haf_err,~] = fpa_given_haf( ...
%             fpa0, ... %  initial fpa guess, (rad)
%             ha_tgt, ... % target apoapse altitude, (m)
%             in, ... % asim input file, contains plantary parameters, vehicle, and ICs
%             fpa_tol, ... % tolerance on flight path angle to search for (rad)
%             ha_tol); % tolerance on final target apoapse altitude (m)

        % 400km
        [fpa400(i,j), ~,~,~,~] = find_fpa_given_ha(corr0, 400,  ...
            in, fpa_tol, ha_tol);
        
        % 2000km
        [fpa2k(i,j), ~,~,~,~] = find_fpa_given_ha(corr0, 2000,  ...
            in, fpa_tol, ha_tol);
        
        % 10000km
        [fpa10k(i,j), ~,~,~,~] = find_fpa_given_ha(corr0, 10000,  ...
            in, fpa_tol, ha_tol);
        
        % GEO
        [fpaGEO(i,j), ~,~,~,~] = find_fpa_given_ha(corr0, 35786,  ...
            in, fpa_tol, ha_tol);

    end %Nj
end %Ni


for i = 1:Ni
    entries{i} = ['\beta = ' num2str(betas(i)) ' kg/m^2'];
end

% xmin = find(isnan(fpa_f(1,:)) == false,1);

% 400 km
figure(); hold on
set(gca,'FontSize',14)
title(['Earth Apoapsis Altitude Target: ' num2str(400) ' km'],'FontSize',10)
xlabel('Entry Velocity (km/s)')
ylabel('Entry Flight Path Angle (deg)')
plot(v_atms,fpa400,'LineWidth',2)
xlim([min(v_atms) max(v_atms)])
legend(entries,'Location','ne','FontSize',10)
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\efpa_vs_vatm400.fig')
close all

% 2000 km
figure(); hold on
set(gca,'FontSize',14)
title(['Earth Apoapsis Altitude Target: ' num2str(2000) ' km'],'FontSize',10)
xlabel('Entry Velocity (km/s)')
ylabel('Entry Flight Path Angle (deg)')
plot(v_atms,fpa2k,'LineWidth',2)
xlim([min(v_atms) max(v_atms)])
legend(entries,'Location','ne','FontSize',10)
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\efpa_vs_vatm2000.fig')
close all

% 10000 km
figure(); hold on
set(gca,'FontSize',14)
title(['Earth Apoapsis Altitude Target: ' num2str(10000) ' km'],'FontSize',10)
xlabel('Entry Velocity (km/s)')
ylabel('Entry Flight Path Angle (deg)')
plot(v_atms,fpa10k,'LineWidth',2)
xlim([min(v_atms) max(v_atms)])
legend(entries,'Location','ne','FontSize',10)
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\efpa_vs_vatm10000.fig')
close all

% GEO
figure(); hold on
set(gca,'FontSize',14)
title(['Earth Apoapsis Altitude Target: ' num2str(35786) ' km'],'FontSize',10)
xlabel('Entry Velocity (km/s)')
ylabel('Entry Flight Path Angle (deg)')
plot(v_atms,fpaGEO,'LineWidth',2)
xlim([min(v_atms) max(v_atms)])
legend(entries,'Location','ne','FontSize',10)
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\efpa_vs_vatmGEO.fig')
close all
%}


%% Corridor Widths as func of v_atm and beta ratio - January 2019
%{
clear; clc; close all;
ratios = 3:2:15;
Ni = length(ratios);

v_atms = [11 13 15].*1e3;
Nj = length(v_atms);

x0.h0 = 125e3;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;
x0.efpa = -5.4;
cd = 1.05.*ones(2,1);
cl = zeros(2,1);
m1 = 500;
rc = [2 0.5];   % rc ratio 0.25
rn = 0.5.*ones(2,1);
ha_tol = 10;    %km
fpa_tol = 0.025;    %deg
corr0 = [-12 -0.1];

% ha_tgt = 2000;    % high LEO
% ha_tgt = 10000;  % high LEO
% ha_tgt = 35786; % ~GEO
% ha_tgt = 400; %km - low LEO
ha_tgt = 200; %km - low LEO

beta1 = m1/(cd(1)*pi*rc(1)^2);
beta2 = nan(Ni,1);

% preallocate
% fpa_shallow = nan(Ni,Nj);
% fpa_steep = nan(Ni,Nj);
corr_width400 = nan(Ni,Nj);
corr_width2k = nan(Ni,Nj);
corr_width10k = nan(Ni,Nj);
corr_widthGEO = nan(Ni,Nj);

for i = 1:Ni
    beta2(i) = ratios(i)*beta1;
    m2 = beta2(i)*cd(2)*pi*rc(2)^2;
    for j = 1:Nj
        x = x0;
        x.v0 = v_atms(j);
        m = [m1 m2];
        
        % 400 km
        [~,~,~,fpa_corr400] = ...
            calc_fpa_corr('earth',x,cd,cl,m,rn,rc,400,ha_tol,fpa_tol,corr0);
        fpa_shallow = fpa_corr400(1);
        fpa_steep = fpa_corr400(2);
        corr_width400(i,j) = fpa_shallow - fpa_steep;
        
        % 2000 km
        [~,~,~,fpa_corr2k] = ...
            calc_fpa_corr('earth',x,cd,cl,m,rn,rc,2000,ha_tol,fpa_tol,corr0);
        fpa_shallow = fpa_corr2k(1);
        fpa_steep = fpa_corr2k(2);
        corr_width2k(i,j) = fpa_shallow - fpa_steep;
        
        % 10000 km
        [~,~,~,fpa_corr10k] = ...
            calc_fpa_corr('earth',x,cd,cl,m,rn,rc,10000,ha_tol,fpa_tol,corr0);
        fpa_shallow = fpa_corr10k(1);
        fpa_steep = fpa_corr10k(2);
        corr_width10k(i,j) = fpa_shallow - fpa_steep;
        
        % GEO
        [~,~,~,fpa_corrGEO] = ...
            calc_fpa_corr('earth',x,cd,cl,m,rn,rc,35786,ha_tol,fpa_tol,corr0);
        fpa_shallow = fpa_corrGEO(1);
        fpa_steep = fpa_corrGEO(2);
        corr_widthGEO(i,j) = fpa_shallow - fpa_steep;
        
    end %j
end %i

% 400 km
figure(); hold on
set(gca,'FontSize',14)
title(['Corridor Width, h_{a,tgt} = ' num2str(400) ' km'],'FontSize',10)
plot(corr_width400,ratios,'LineWidth',2)
xlabel('Corridor Width (deg)')
ylabel('\beta_2/\beta_1')
legend('V_{atm} = 11','V_{atm} = 13','V_{atm} = 15','location','best')
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\corridor_width\corr_width400.fig')
close all

% 2000 km
figure(); hold on
set(gca,'FontSize',14)
title(['Corridor Width, h_{a,tgt} = ' num2str(2000) ' km'],'FontSize',10)
plot(corr_width2k,ratios,'LineWidth',2)
xlabel('Corridor Width (deg)')
ylabel('\beta_2/\beta_1')
legend('V_{atm} = 11','V_{atm} = 13','V_{atm} = 15','location','best')
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\corridor_width\corr_width2k.fig')
close all

% 10000 km
figure(); hold on
set(gca,'FontSize',14)
title(['Corridor Width, h_{a,tgt} = ' num2str(10000) ' km'],'FontSize',10)
plot(corr_width10k,ratios,'LineWidth',2)
xlabel('Corridor Width (deg)')
ylabel('\beta_2/\beta_1')
legend('V_{atm} = 11','V_{atm} = 13','V_{atm} = 15','location','best')
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\corridor_width\corr_width10k.fig')
close all

% GEO
figure(); hold on
set(gca,'FontSize',14)
title(['Corridor Width, h_{a,tgt} = ' num2str(35786) ' (GEO) km'],'FontSize',10)
plot(corr_widthGEO,ratios,'LineWidth',2)
xlabel('Corridor Width (deg)')
ylabel('\beta_2/\beta_1')
legend('V_{atm} = 11','V_{atm} = 13','V_{atm} = 15','location','best')
saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\corridor_width\corr_widthGEO.fig')
close all

% 
% figure(); hold on
% set(gca,'FontSize',14)
% title(['Corridor Bounds, h_{a,tgt} = ' num2str(ha_tgt/1000) ' km'],'FontSize',10)
% plot(ratios,fpa_shallow,'b','LineWidth',2)
% plot(ratios,fpa_steep,'r','LineWidth',2)
%}

%% Run NPC - tj_opt & haf_err, vs. beta ratios, multiple ha_tgts - Jan 2019
% high fidelity
%{
clear; clc; close all;
v_atm = 12;   %km/s, entry speed
% m = [80 40];
rcs = [0.75 0.175];
rn = 0.5;       % nose radius
cd = 1.05;      % drag coefficient (constant)
% beta1 = 45;
% m1 = beta1*cd*pi*rcs(1)^2;
m1 = 500;
beta1 = m1/(cd*pi*rcs(1)^2);
fpa_atm = -6.5; %deg, efpa

% ha_tgt = 2000;   %km, apoapsis target
ha_tgt = [400 2000 10000];
Nj = length(ha_tgt);
tol_ha = 10;
n = 1;
tj0 = 100;   % initial guess on jettison

ratios = 3:1:15;
beta2 = beta1*ratios;
Ni = length(ratios);

root_mode = uint8(1);
guid_rate = 5;  %guidance rate
mc.flag = uint8(0);
iters = uint8(1);
traj_rate = 1000; % simulation rate;

% preallocate
ldat = nan(Ni,1);
tj_opt400 = ldat;
haf_err400 = ldat;
tj_opt2k = ldat;
haf_err2k = ldat;
tj_opt10k = ldat;
haf_err10k = ldat;

rc2 = rcs(2);

parfor i = 1:Ni
    m2 = beta2(i)*cd*pi*rc2^2;
    m = [m1 m2];
%     for j = 1:Nj

    % 400 km
    out = dej_n_auto_earth(fpa_atm, v_atm, 400, tol_ha, n, ...
        tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
    tj_opt400(i) = out.t_jett/out.traj.time(out.idxend); %s
    haf_err400(i) = out.haf_err;   %km
    
    % 2000 km
    out = dej_n_auto_earth(fpa_atm, v_atm, 2000, tol_ha, n, ...
        tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
    tj_opt2k(i) = out.t_jett/out.traj.time(out.idxend); %s
    haf_err2k(i) = out.haf_err;   %km

    % 10000 km
    out = dej_n_auto_earth(fpa_atm, v_atm, 10000, tol_ha, n, ...
        tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
    tj_opt10k(i) = out.t_jett/out.traj.time(out.idxend); %s
    haf_err10k(i) = out.haf_err;   %km    
        
%     end
end

figure(); hold on
set(gca,'FontSize',14)
title(['Optimal Jettison Time vs. Beta Ratio, v_atm = 12km/s, efpa= ' num2str(fpa_atm) ' deg'],'FontSize',10)
yyaxis left
plot(ratios,tj_opt400,'-o','LineWidth',2)
plot(ratios,tj_opt2k,'--+','LineWidth',2)
plot(ratios,tj_opt10k,':^','LineWidth',2)
ylabel('t_{jettison}/t_{flight}')
yyaxis right
plot(ratios,haf_err400,'-o','LineWidth',2)
plot(ratios,haf_err2k,'--+','LineWidth',2)
plot(ratios,haf_err10k,':^','LineWidth',2)
ylabel('Apoapsis Error (km)')
ylim([-0.5 0.5])
xlim([ratios(1) ratios(end)])
xlabel('\beta_2/\beta_1')
legend('h_{a,tgt} = 400 km', ...
    'h_{a,tgt} = 2000 km', ...
    'h_{a,tgt} = 10000 km', ...
    'h_{a,tgt} = 400 km', ...
    'h_{a,tgt} = 2000 km', ...
    'h_{a,tgt} = 10000 km', ...
    'location','se','FontSize',10)

saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\opt_t_jett\tjOpt_vs_ratio_efpa5.fig')
% close all
%}

%% Run NPC - opt tj & haf_err vs. efpa, multiple ha_tgts
%{
clear; clc; close all
x0.v0 = 12;   %km/s, entry speed
x0.h0 = 125;    %km
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

aero.rcs = [0.75 0.175];
aero.rn = 0.5;       % nose radius
aero.cds = [1.05 1.05];      % drag coefficient (constant)
aero.cls = [0 0];
% beta1 = 45;
% m1 = beta1*cd*pi*rcs(1)^2;
m1 = 80;
beta1 = m1/(aero.cds(1)*pi*aero.rcs(1)^2);
fpa_atm = linspace(-5.8,-6.5,20); %deg, efpa
Ni = length(fpa_atm);

% ha_tgt = 2000;   %km, apoapsis target
ha_tgt400 = 400;
ha_tgt2k = 2000;
ha_tgt10k = 10000;

gnc.ha_tgt = 2000;
gnc.tol_ha = 10;
gnc.n = 1;
gnc.tj0 = 100;   % initial guess on jettison

ratio = 9;
beta2 = ratio*beta1;
m2 = beta2*aero.cds(2)*pi*aero.rcs(2)^2;
aero.m = [m1 m2];
    
gnc.root_mode = uint8(1); % newton method
gnc.guid_rate = 0.2;  %guidance rate
gnc.iters = uint8(1);
gnc.dtj_lim = 20;

sim.traj_rate = 50; % simulation rate;
sim.h_max = 125e3;
sim.h_min = 60e3;
sim.t_max = 750;
sim.data_rate = 1;

mc.flag = uint8(0);

% preallocate
ldat = nan(Ni,1);
tj_opt400 = ldat;
haf_err400 = ldat;
tj_opt2k = ldat;
haf_err2k = ldat;
tj_opt10k = ldat;
haf_err10k = ldat;

for i = 1:Ni
% for i = 1:Ni
    % 400 km
%     out = dej_n_auto_earth(fpa_atm(i), v_atm, ha_tgt400, tol_ha, n, ...
%         tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
%     tj_opt400(i) = out.t_jett/out.traj.time(out.idxend); %s
%     haf_err400(i) = out.haf_err;   %km
    
    % 2000 km
%     out = dej_n_auto_earth(fpa_atm(i), v_atm, ha_tgt2k, tol_ha, n, ...
%         tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
    x0.fpa0 = fpa_atm(i);

    out = dej_n_auto_earth(x0,gnc,aero,mc,sim);
    tj_opt2k(i) = out.t_jett/out.traj.time(out.idxend); %s
    haf_err2k(i) = out.haf_err;   %km

    % 10000 km
%     out = dej_n_auto_earth(fpa_atm(i), v_atm, ha_tgt10k, tol_ha, n, ...
%         tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
%     tj_opt10k(i) = out.t_jett/out.traj.time(out.idxend); %s
%     haf_err10k(i) = out.haf_err;   %km    
        
%     end
end


figure(); hold on
set(gca,'FontSize',14)
title(['Optimal Jettison Time vs. EFPA, h_{a,tgt}; \beta_2/\beta_1= ' num2str(ratio)],'FontSize',10)
yyaxis left
plot(fpa_atm,tj_opt400,'-o','LineWidth',2)
plot(fpa_atm,tj_opt2k,'--+','LineWidth',2)
plot(fpa_atm,tj_opt10k,':^','LineWidth',2)
ylabel('t_{jettison}/t_{flight}')
yyaxis right
plot(fpa_atm,haf_err400,'-o','LineWidth',2)
plot(fpa_atm,haf_err2k,'--+','LineWidth',2)
plot(fpa_atm,haf_err10k,':^','LineWidth',2)
ylabel('Apoapsis Error (km)')
ylim([-0.5 0.5])
xlim([fpa_atm(end) fpa_atm(1)])
xlabel('EFPA (deg)')
legend('h_{a,tgt} = 400 km', ...
    'h_{a,tgt} = 2000 km', ...
    'h_{a,tgt} = 10000 km', ...
    'h_{a,tgt} = 400 km', ...
    'h_{a,tgt} = 2000 km', ...
    'h_{a,tgt} = 10000 km', ...
    'location','se','FontSize',10)

% saveas(gcf,'C:\Users\Evan Roelke\Documents\research\earth_aerocapture\entry_trades\opt_t_jett\tjOpt_vs_efpa_ratio9.fig')
% close all
%}



%% updated atm model NPC sim comparisons - Feb 2019
%{
efpa = -6.0;
v0 = 12;
ha_tgt = 2000;
tol_ha = 10;
n = 1;
tj0 = 120;
rcs = [0.75 0.175];
m = [80 40];
root_mode = uint8(1);   %newton
mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
iters = uint8(1);
atm_mode = uint8(1);    %atm capture on, table size 1000

traj_rate = 100;
guid_rate = 0.5;

tic;
out = dej_n_auto_earth(efpa, v0, ha_tgt, tol_ha, n, ...
    tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate,atm_mode);
dt = toc;

tols = [25 50 100 150];
for i = 1:length(tols)
    num = find( abs(out.haf_err) <= tols(i) );
    pc(i) = length(num)*100/mc.N;
    entries{i} = [num2str(pc(i)) '% within ' num2str(tols(i)) ' km'];
end

figure(); hold on
grid on
set(gca,'FontSize',14)
title(['DEJ-N Earth, AtmCapture, EFPA = ' num2str(efpa) '^{o},h_{a,tgt} = ' ... 
    num2str(ha_tgt) ' km, GNC Rate= ' num2str(guid_rate) ' Hz'],'FontSize',10)
histogram(out.haf_err,'NumBins',50,'FaceColor','r')
xlabel('Apoapsis Error (km)')
ylabel('Occurrences')
text(-150,60,entries)


tjr = nan(mc.N,1);
for i = 1:mc.N
    idend = find( isnan(out.traj.t(:,i)) == 1,1)-1;
    if isempty(idend)
        idend = out.traj.t(end,i);
    end
    tjr(i) = out.tjett(i,1)/out.traj.t(idend,i);
end
    
figure(); hold on
grid on
set(gca,'FontSize',14)
title(['DEJ-N Earth, AtmCapture, EFPA = ' num2str(efpa) '^{o},h_{a,tgt} = ' ... 
    num2str(ha_tgt) ' km, GNC Rate= ' num2str(guid_rate) ' Hz'],'FontSize',10)
plot(tjr,out.haf_err,'*')
xlabel('Jettison Time/Flight Time (s)')
ylabel('Apoapsis Error (km)')
%}

%% updated atm model NPC sim comparisons - Feb 2019
%{
% compare error as func of ha_tgt / jett time
efpa = -5.9;
v0 = 12;
% ha_tgts = [400 800 1000 2000 5000 10000 20000 30000 40000];
ha_tgts = linspace(400,40000,10);
tol_ha = 10;
n = 1;
tj0 = 120;
rcs = [0.75 0.175];
m = [80 40];
cd = 1.05;
root_mode = uint8(1);   %newton

b1 = m(1)/(cd*pi*rcs(1)^2);
b2 = m(2)/(cd*pi*rcs(2)^2);
bratio = b2/b1;

mc.flag = uint8(1);

mc.N = 1000;
% mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
mc.sigs = default_sigs();
iters = uint8(1);

atm_mode = uint8(0);    %atm capture on, table size 1000

traj_rate = 250;
guid_rate = 0.5;

Ni = length(ha_tgts);

x0.fpa0 = -5.9;
x0.v0 = v0;
x0.h0 = 125;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

gnc.ha_tgt = 2000;
gnc.ha_tol = 10;
gnc.n = 1;
gnc.tj0 = 120;
gnc.dtj_lim = 30;
gnc.root_mode = uint8(1);
gnc.guid_rate = 0.5;
gnc.atm_mode = uint8(0);
gnc.iters = uint8(1);

aero.rcs = [0.75 0.175];
aero.m = [80 40];
aero.rn = 0.175/2;
aero.cds = [1.23 1.05];
aero.cls = [0 0];

mc.flag = true;
mc.N = 1;
mc.sigs = nan;
mc.debug = true;

sim.traj_rate = 100;
sim.data_rate = 1;
sim.t_max = 800;
sim.h_max = x0.h0;
sim.h_min = 50;
sim.planet = 'earth';
sim.atm_mode = uint8(3);
sim.efpa_flag = true;

tj_reg = nan(mc.N,Ni);
err_reg = tj_reg;
tjr_reg = tj_reg;

tj_atm = tj_reg;
err_atm = tj_reg;
tjr_atm = tj_reg;

tols = [25 50 100 250 500];
Nj = length(tols);

pc_reg = nan(Ni,Nj);
pc_atm = pc_reg;

for i = 1:Ni
    fprintf('Loop %i, %4.2f%% done\n',i,(i-1)*100/Ni);
%     out_reg = dej_n_auto_earth(efpa, v0, ha_tgts(i), tol_ha, n, ...
%         tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate,uint8(0));
    out_reg = run_dej_n(x0,gnc,aero,sim,mc);
    
    tj_reg(:,i) = out_reg.tjett(:,1);
    err_reg(:,i) = out_reg.haf_err;    
    for j = 1:mc.N
        tjr_reg(j,i) = tj_reg(j,i)/out_reg.traj.t(out_reg.idxend(j),j);
    end
    
    for j = 1:Nj
        num = find( abs(out_reg.haf_err) <= tols(j) );
        pc_reg(i,j) = length(num)*100/mc.N;
    end
    clear out_reg;
    
    out_atm = dej_n_auto_earth(efpa, v0, ha_tgts(i), tol_ha, n, ...
        tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate,atm_mode);
    
    tj_atm(:,i) = out_atm.tjett(:,1);
    err_atm(:,i) = out_atm.haf_err;    
    for j = 1:mc.N
        tjr_atm(j,i) = tj_atm(j,i)/out_atm.traj.t(out_atm.idxend(j),j);
    end
    
    for j = 1:Nj
        num = find( abs(out_atm.haf_err) <= tols(j) );
        pc_atm(i,j) = length(num)*100/mc.N;
    end
    clear out_atm;
    
end

save('C:\Users\Evan Roelke\Documents\research\earth_aerocapture\dej_n\atm_cmp_haTgt.mat', ... 
    'ha_tgts','pc_reg','pc_atm','tols','efpa','rcs','m','guid_rate','tjr_reg','tjr_atm', ... 
    'err_reg','err_atm');

leg_reg = cell(Nj,1);
leg_atm = leg_reg;
for i = 1:Nj
    leg_reg{i} = ['Reg Atm, Tol = ' num2str(tols(i)) ' km'];
    leg_atm{i} = ['Atm Est, Tol = ' num2str(tols(i)) ' km'];
end

header = {['Atm Est. DEJ-N, EARTH, efpa = ' num2str(efpa) ' deg, v0 = ' ... 
    num2str(v0) ' km/s'], ... 
    ['GNC rate = ' num2str(guid_rate) ' Hz, \beta_2/\beta_1 = ' num2str(round(bratio,2))]};
figure(); hold on
grid on
set(gca,'FontSize',14)
title(header,'FontSize',10);
% plot(ha_tgts,err_atm,'LineWidth',2)
plot(ha_tgts./10000,pc_reg(:,1),'-r','LineWidth',1.5)
plot(ha_tgts./10000,pc_reg(:,2),'o-r','LineWidth',1.5)
plot(ha_tgts./10000,pc_reg(:,3),'+-r','LineWidth',1.5)
plot(ha_tgts./10000,pc_reg(:,4),'*-r','LineWidth',1.5)
plot(ha_tgts./10000,pc_reg(:,5),'^-r','LineWidth',1.5)
plot(ha_tgts./10000,pc_atm(:,1),'-b','LineWidth',1.5)
plot(ha_tgts./10000,pc_atm(:,2),'o-b','LineWidth',1.5)
plot(ha_tgts./10000,pc_atm(:,3),'+-b','LineWidth',1.5)
plot(ha_tgts./10000,pc_atm(:,4),'*-b','LineWidth',1.5)
plot(ha_tgts./10000,pc_atm(:,5),'^-b','LineWidth',1.5)
xlabel('Target Apoapsis (10^4 km)')
ylabel('Percent Captured')
legend(vertcat(leg_reg,leg_atm),'location','eastoutside','FontSize',10)


% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% p1 = plot(tjr_reg(:,1),err_reg(:,1),'*r');
% p2 = plot(tjr_atm(:,1),err_atm(:,1),'*b');
% ylim([-250 250])
% legend([p1,p2],'Density Corr.','Atm Estimation','location','nw','FontSize',11)

% figure(); hold on
% histogram(err_reg,'NumBins',5000,'FaceColor','r');
% histogram(err_atm,'NumBins',5000,'FaceColor','b');
% xlim([-1000 1000])
% 
% 
% % post processing
% dat = load('C:\Users\Evan Roelke\Documents\research\earth_aerocapture\dej_n\atm_cmp_efpa.mat');
% ind = 1;
% 
% for ind = 1:length(dat.efpas)
%     efpa = dat.efpas(ind);
%     eReg = dat.err_reg(:,ind);
%     tjrReg = dat.tjr_reg(:,ind);
%     eAtm = dat.err_atm(:,ind);
%     tjrAtm = dat.tjr_atm(:,ind);
% 
%     eRegAvg(ind) = mean(eReg);
%     eRegStd(ind) = std(eReg);
%     eAtmAvg(ind) = mean(eAtm);
%     eAtmStd(ind) = std(eAtm);
% end
% 
% figure(); hold on
% plot(dat.efpas(1:15),3*eRegStd(1:15))
% plot(dat.efpas(1:15),3*eAtmStd(1:15))
% 
% % for j = 1:length(eReg)
% %    rssReg(j) = eReg(j) - eRegAvg;
% %    rssAtm(j) = eAtm(j) - eAtmAvg;
% % end
% 
% n = length(eReg);
% for i = 1:n
%     for j = 1:length(dat.efpas)
%         rssReg(i,j) = ( dat.err_reg(i,j) - mean(dat.err_reg(:,j)) )^2;
%         rssAtm(i,j) = ( dat.err_atm(i,j) - mean(dat.err_atm(:,j)) )^2;
%     end
%     sReg(i) = sqrt( sum(rssReg(i,:))/n ); %1sig value
%     sAtm(i) = sqrt( sum(rssAtm(i,:))/n ); %1sig value
% end
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% p1 = plot(dat.tjr_reg(:,ind),dat.err_reg(:,ind),'*r');
% p2 = plot(dat.tjr_atm(:,ind),dat.err_atm(:,ind),'*b');
% p3 = plot(sReg,'k');
% ylim([-250 250])
% legend([p1,p2],'Density Corr.','Atm Estimation','location','nw','FontSize',11)

%}

%% kalman filtering attempts - Apr 2019
%{
clear;clc;
x0.fpa0 = -5.9;
x0.v0 = 12;
x0.h0 = 125;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

gnc.ha_tgt = 2000;
gnc.ha_tol = 10;
gnc.n = 1;
gnc.tj0 = 120;
gnc.dtj_lim = 30;
gnc.root_mode = uint8(1);
gnc.guid_rate = 0.1;
gnc.atm_mode = uint8(2);
gnc.iters = uint8(1);
gnc.atm_mode = uint8(2);

aero.rcs = [0.75 0.175];
aero.m = [80 40];
aero.rn = 0.175/2;
aero.cds = [1.23 1.05];
aero.cls = [0 0];

mc.flag = true;
mc.N = 1;
mc.sigs = nan;
mc.debug = true;
mc.sigs = default_sigs();

sim.traj_rate = 10;
sim.data_rate = 1;
sim.t_max = 800;
sim.h_max = x0.h0;
sim.h_min = 50;
sim.planet = 'earth';
sim.efpa_flag = false;


out = run_dej_n(x0,gnc,aero,sim,mc);
%}


%% nominal trajectory atmospheric test - Apr 2019
%{
clear;clc;
% get nominal atm.
[x0,aero,gnc,mc,sim] = BaseInputs();
% sim.efpa_flag = false;
% x0.fpa0 = -5.76428769429527;
% sim.efpa_flag = false;
% out = run_dej_n(x0,gnc,aero,sim,mc);
% 
% 
% alt = linspace(min(out.traj.alt),max(out.traj.alt),1000)';
% v = nan(length(alt),1);
% a = v;
% rho = v;
% for i = 1:length(alt)
%     ind = find(out.traj.alt <= alt(i),1);
%     v(i) = out.traj.vel_pp_mag(ind);
% %     a(i) = (out.force_ii(ind)-out.gravity_ii(ind))/out.mass(ind);
%     a(i) = 0;
%     rho(i) = out.traj.rho(ind);
% end
% 
% table = [alt v a rho];
% save('ref.mat','table');
% keyboard;


x0.fpa0 = -5.76428769429527;
mc.flag = true;
sim.efpa_flag = false;
mc.N = 1;
mc.debug = true;
gnc.atm_mode = uint8(3);
gnc.guid_rate = 0.5;
sim.traj_rate = 50;
out = run_dej_n(x0,gnc,aero,sim,mc);
%}

%% visualize atm dispersion
%{
[nom, mc] = parse_atm_data(0, 'earth');
Nk = 20;
avgRho = nan(Nk,1);
stdRho = avgRho;
rho_mc = nan(1000,Nk);
alts = linspace(75,100,Nk).*1000;
for k = 1:Nk
    for i = 1:1000
        rho_mc(i,k) = lin_interp(mc(i).table(:,1),mc(i).table(:,2),alts(k));
    end
    avgRho(k) = mean(rho_mc(:,k));
    stdRho(k) = std(rho_mc(:,k));
end
alts = alts./1000;

% exp. model
p = get_planetary_data( 2, 0 );
h = linspace(50,150,1000).*1000;
rhoExp = nan(1000,1);
for i = 1:1000
    rhoExp(i) = p.atm.dens_ref*exp(-h(i)/p.atm.scale_height);
end
    

figure();
s1 = semilogy(alts.*ones(1000,1),rho_mc,'*b'); hold on
s2 = semilogy(alts, avgRho,'k','LineWidth',2); hold on
s3 = semilogy(alts, avgRho + 3*stdRho,'--k','LineWidth',2); hold on
s4 = semilogy(alts, avgRho - 3*stdRho,'--k','LineWidth',2); hold on
s5 = semilogy(nom(:,1)./1000,nom(:,2),'--g','LineWidth',2); hold on
s6 = semilogy(h./1000,rhoExp,'r','LineWidth',2); hold on
set(gca,'FontSize',14);
xlabel('Altitude (km)')
ylabel('Atmospheric Density (kg/m^3)')
title('\pm3\sigma Atmospheric Density vs. Altitude, Earth')
xlim([alts(1) alts(end)])
legend([s2,s3,s4,s5,s6], ... 
    'Average Density','+3\sigma', ... 
    '-3\sigma','Nominal Model','Exponential Model','location','sw')
%}

%% Atmospheric Modeling Ideas
%{
% [x0,aero,gnc,mc,sim] = BaseInputs();
% x0.fpa0 = -5.76428769429527;
% mc.flag = true;
% sim.efpa_flag = false;
% mc.N = 1000;
% mc.debug = false;
% gnc.atm_mode = uint8(3);
% sim.traj_rate = 50;
% out1 = run_dej_n(x0,gnc,aero,sim,mc);   %sinusoid
% 
% gnc.atm_mode = uint8(4);
% out2 = run_dej_n(x0,gnc,aero,sim,mc);   %only nominal atm
% 
% keyboard;
% 
% figure(); hold on
% histogram(out1.haf_err, 'NumBins',50)
% histogram(out2.haf_err, 'NumBins',50)
% legend('Sinusoid','Nominal','location','ne')


p = get_planetary_data(2,1);
H = p.atm.scale_height;
temp = load('./data/atm_data/atm_earth_gram2016.mat');
nom = temp.nom_table;
mc = temp.mc;

rho_nom = nom(:,2);

h = nan(1000,1);
expAtm = h;
sinAtm = h;
sinAtm2 = h;
for i = 1:1000
    h(i) = nom(i,1);
    expAtm(i) = nom(1,2) * exp(-h(i)/H);
%     sinAtm(i) = expAtm(i)*sin(4*pi*1 + pi*((h(i)/1000)/100));
%     sinAtm(i) = sin(expAtm(i)*(h(i)/1000) / 100);
    sinAtm(i) = expAtm(i)*(1 + (1/10 * sin(pi*h(i)/1000)));
    
    sinAtm2(i) = nom(i,2)*(1 + (1/10 * sin(pi*h(i)/1000)));
    
end

figure()
semilogx(expAtm,h./1000,'k'); hold on
semilogx(sinAtm,h./1000,'r'); hold on
semilogx(nom(:,2),nom(:,1)./1000,'b'); hold on
% semilogx(mc(3).table(:,2),mc(3).table(:,1)./1000,'b'); hold on
ylim([50 140])

% visualize RSS error between K_dens and model
RvarExp = nan(1000,1);
RvarS = RvarExp;
RvarS2 = RvarExp;
RvarMC = RvarExp;
rss = 0;
for j = 1:1000
    RvarExp(j) = expAtm(j)/rho_nom(j); 
    RvarS(j) = sinAtm(j)/rho_nom(j);
    RvarS2(j) = sinAtm2(j)/rho_nom(j);
    RvarMC(j) = mc(2).table(j,2)/rho_nom(j);
end

figure();
semilogx(nom(:,2),nom(:,1)./1000,'k'); hold on
semilogx(mc(400).table(:,2),mc(400).table(:,1)./1000,'r'); hold on
% semilogx(expAtm,h./1000,'b'); hold on
semilogx(sinAtm2,h./1000)


figure(); hold on
% plot(RvarExp,h./1000,'r');
plot(RvarS2,h./1000,'r');
plot(RvarMC,h./1000,'b');
%}

%%%%%%% base setup
function [x0,aero,gnc,mc,sim] = BaseInputs()
[x0,aero,gnc,sim,mc] = base_venus_ac(false);

x0.fpa0 = -5.7;
x0.v0 = 12;
x0.h0 = 125;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

gnc.ha_tgt = 2000;
gnc.ha_tol = 10;
gnc.n = 1;
gnc.tj0 = 100;
gnc.dtj_lim = 30;
gnc.root_mode = uint8(1);
gnc.guid_rate = 0.1;
gnc.iters = uint8(1);
gnc.atm_mode = uint8(0);
gnc.mode = 4;   %dej_n

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
mc.Kflag = true;

sim.traj_rate = 250;
sim.data_rate = 1;
sim.t_max = 1000;
sim.h_max = x0.h0;
sim.h_min = 40;
sim.planet = 'earth';
sim.efpa_flag = true;
sim.atm_mode = 3;
sim.ignore_nominal = false;
end



