%% mars_ac.m
% clear;clc
% 
% % TCW plots
% efpa0 = -6; %initial guess
% v0s = 5:1:10;
% betas = 10:10:250;
% 
% [x0,aero,gnc,sim,mc] = base_jsr();
% sim.planet = 'mars';
% 
% % Ni = length(efpas);
% Nj = length(v0s);
% Nk = length(betas);
% fpa = nan(Nj,Nk);
% 
% cd = 1.2;
% rc = 4.5/2;
% m = nan(Nk,1);
% 
% fpa_guess = [-1 -20];
% ha_tgt = 400;
% tol_fpa = 0.1;  %deg
% tol_ha = 10;    %km
% 
% in = dme_mars_in;
% in.v.aero.area_ref = pi * rc^2;
% in.v.aero.cd = cd;
% in.v.aero.nose_radius = 0.5;
% in.v.aero.mode = uint8(1);
% in.v.gnc.g.p.rate = 0;
% 
% lat = 0;
% lon = 0;
% 
% in.s.traj.alt = 100e3; % atmospheric interface altitude, m
% in.s.traj.lat = lat * pi/180; % initial latitude (72.4deg north), rad
% in.s.traj.lon = lon * pi/180; % initial longitude (338deg west), rad
% % in.s.traj.vel_pp_mag = 6000; % initial velocity, m/s
% % in.s.traj.gamma_pp = -65.4*pi/180; % initial flight-path, rad (Labreton ref)
% in.s.traj.gamma_pp = efpa0 * pi/180;
% in.s.traj.az = 90*pi/180; % initial azimuth, rad
% in.s.traj.t_max = 8000;
% in.s.traj.rate = 100;
% 
% for j = 1:Nj
%     for k = 1:Nk
%             
%         m(k) = betas(k) * pi * cd * rc^2;
%         in.v.mp.m_ini = m(k);
%         in.s.traj.vel_pp_mag = v0s(j)*1000;
% 
%         % Convert scalars to inertial position and velocity vectors
%         [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
%             in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
%             in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
%             in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);
% 
%             fpa(j,k) = find_fpa_given_ha( ...
%                 fpa_guess, ... % 1x2 vector of initial guess on bounds, MUST encompass solution (deg)
%                 ha_tgt, ... % target apoapse altitude (km)
%                 in, ... % asim input file, implicitly contains plantary parameters and ICs
%                 tol_fpa, ... % tolerance on final output FPA (deg)
%                 tol_ha); % tolerance on final target apoapse altitude (km)
% 
% 
% %         [fpa(j,k),~,haf_err,~] = fpa_given_haf( ...
% %             efpa0, ... %  initial fpa guess, (deg)
% %             ha_tgt, ... % target apoapse altitude, km)
% %             in, ... % asim input file, contains plantary parameters, vehicle, and ICs
% %             tol_fpa, ... % tolerance on flight path angle to search for (deg)
% %             tol_ha); % tolerance on final target apoapse altitude km)
%         
%     end
% end
% 
% % save(['C:\Users\Evan Roelke\Documents\research\mars_ac\data\efpa_vs_v0_' ... 
% %     num2str(ha_tgt/1000) 'k.mat']);


%% Monte Carlo
clear;clc
% ha_tgts = [1400 10000];
% [x0,aero,gnc,sim,mc] = base_jsr();

% gnc.ha_tgt = 1400;
gnc.tj0 = 70;
gnc.dtj_lim = 20;
beta1 = 148;
% ratio = 10;
% x0.v0 = 5.7;
% x0.fpa0 = -9;

% fpas = [-11.5 -9.5];
% LtD = 0.25;
% cd = 1.2;
% cl = 1.2 * LtD;

% mc.flag = false;
% mc.debug = true;

% gnc.ha_tgt = 10000;
% haf_err = gnc.ha_tgt;
% while (abs(haf_err) > 5)
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_err = out.haf_err;
%     if abs(haf_err) < 5
%         break;
%     else
%         x0.fpa0 = x0.fpa0 - 0.05;
%     end
% end
% keyboard;

haf_tol = 50;

fpa0s = linspace(-8.5, -12, 20)';
% v0s = linspace(5, 8, 10)';
v0s = 5.7;
LDs = [0.1 0.2 0.3 0.4];
ratios = [5 10 15 20];
% ha_tgts = 1400:200:10000;
ha_tgt = 10000;

Ni = length(fpa0s);
Nj = length(ratios);

tjrL = nan(Ni,Nj);       tjrD = tjrL;
errL = nan(Ni,Nj);       errD = errL;
gLoadL = nan(Ni,Nj);     gLoadD = gLoadL;
heatRateL = nan(Ni,Nj);  heatRateD = heatRateL;
parfor i = 1:Ni
    for j = 1:Nj
        [x0,aero,gnc,sim,mc] = base_jsr();
        sim.planet = 'mars';
        x0.h0 = 100;               
        x0.fpa0 = fpa0s(i);
        x0.v0 = v0s;
        gnc.ha_tgt = ha_tgt;

        sim.h_max = x0.h0;
        sim.h_min = 0;
        mc.sigs.cd = aero.cds(1) * 0.05;    %5%
        mc.sigs.cl = aero.cls(1) * 0.05;
        gnc.npc_mode = uint8(1);
        gnc.guid_rate = 0.1;
        sim.t_max = 2500;
        sim.efpa_flag = false;
        sim.traj_rate = 50;
        sim.parMode = true;

        mc.flag = false;
        mc.N = 0;
        mc.Kflag = false;

        aero.cds = [1.2 1.2];
        aero.rn = 0.2;
        m1 = 5500;
        %%%%%%%%%%%%%%%
        % lift mod
        aero.cls = aero.cds(1) * LDs(j) .* [1 1];
        rc1 = sqrt(m1/(aero.cds(1)*beta1*pi));
        aero.rcs = [rc1 rc1];
        aero.m = [m1 m1];

        out = run_dba_n(x0,gnc,aero,sim,mc);
        errL(i,j) = out.haf_err;
        gLoadL(i,j) = max(out.traj.g_loading);
        heatRateL(i,j) = max(out.traj.heat_rate);

        if (abs(errL(i,j)) < haf_tol)
            tjrL(i,j) = out.t_jett / out.traj.time(out.idxend);
        else
            tjrL(i,j) = nan;
        end

        %%%%%%%%%%%%%%%%%
        % drag mod

        m1 = 5500;
        m2 = 3000;
        rc1 = sqrt(m1/(aero.cds(1)*beta1*pi));
        rc2 = sqrt(m2/(aero.cds(1)*beta1*ratios(j)*pi));
        aero.rcs = [rc1 rc2];
        aero.m = [m1 m2];                
        aero.cls = [0 0];

        out = run_dej_n(x0,gnc,aero,sim,mc);

        errD(i,j) = out.haf_err;
        gLoadD(i,j) = max(out.traj.g_loading);
        heatRateD(i,j) = max(out.traj.heat_rate);

        if (abs(errD(i,j)) < haf_tol)
            tjrD(i,j) = out.t_jett / out.traj.time(out.idxend);
        else
            tjrD(i,j) = nan;
        end
    end
end

keyboard;

% save('C:\Users\Evan Roelke\Documents\research\mars_ac\data\opt_ctrl_v057_1400kmTgt.mat', ... 
%     'errD','errL','tjrD','tjrL','gLoadD','gLoadL','heatRateD','heatRateL', ... 
%     'fpa0s','v0s','LDs','ratios','ha_tgts','Ni','Nj','Nk','Nz');

% save('C:\Users\Evan Roelke\Documents\research\mars_ac\data\opt_ctrl_v057_McCtrl.mat', ... 
%     'errD','errL','tjrD','tjrL','gLoadD','gLoadL','heatRateD','heatRateL', ... 
%     'fpa0s','v0s','LDs','ratios','ha_tgts','Ni','Nj');

% save('C:\Users\Evan Roelke\Documents\research\mars_ac\data\opt_control_events.mat', ... 
%     'errD','errL','tjrD','tjrL','gLoadD','gLoadL','heatRateD','heatRateL', ... 
%     'fpa0s','v0s','LDs','ratios','ha_tgts','Ni','Nj','Nk','Nz');

save('C:\Users\Evan Roelke\Documents\research\mars_ac\data\opt_ctrl_v57_10kTgt.mat', ... 
    'errD','errL','tjrD','tjrL','gLoadD','gLoadL','heatRateD','heatRateL', ... 
    'fpa0s','v0s','LDs','ratios','ha_tgt','Ni','Nj');


% figure(); hold on
% title('1400 km Target')
% set(gca,'FontSize',14)
% plot(fpa0s, tjr, 'LineWidth',2)
% xlabel('EFPA')
% ylabel('Jettison Time Ratio (t_{switch}/t_f)')


% for i = 1:4
%     tjrLift(:,i) = tjrL(:,1,i);
%     tjrDrag(:,i) = tjrD(:,1,i);
% end

% reshape(tjrD, [20,4]);
% tjrL = reshape(tjrL, [20,4]);
% 
% 

figure(); hold on
set(gca,'FontSize',14)
title('1400 km Target, 5.7 km/s','FontSize',10)
plot(fpa0s, tjrL(:,1),'r','Linewidth',2)
plot(fpa0s, tjrL(:,2),'r-.','Linewidth',2)
plot(fpa0s, tjrL(:,3),'r--','Linewidth',2)
plot(fpa0s, tjrL(:,4),'r:','Linewidth',2)
plot(fpa0s, tjrD(:,1),'b','Linewidth',2)
plot(fpa0s, tjrD(:,2),'b-.','Linewidth',2)
plot(fpa0s, tjrD(:,3),'b--','Linewidth',2)
plot(fpa0s, tjrD(:,4),'b:','Linewidth',2)
ylim([0.1 0.4])
ylabel('Control Event Time Ratio, t_*/t_f')
legend( 'L/D = 0.1', ... 
        'L/D = 0.2', ... 
        'L/D = 0.3', ... 
        'L/D = 0.4', ... 
        '\beta_2/\beta_1 = 5', ... 
        '\beta_2/\beta_1 = 10', ... 
        '\beta_2/\beta_1 = 15', ... 
        '\beta_2/\beta_1 = 20', ... 
        'location','nw')
xlabel('Entry Flight Path Angle (degrees)')
%     
% 
% keyboard;

%%%%%
keyboard
clear;clc
% load('C:\Users\Evan Roelke\Documents\research\mars_ac\data\opt_control_events.mat');
tjs = [0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38];
%x(fpa, v0, ctrl, tgt)

% dat1 = nan(Ni,Nz);
% dat2 = dat1;
% for i = 1:20
%     for j = 1:10
%         dat1(i,j) = tjrD(i,3,5,j);
%         dat2(i,j) = tjrL(i,3,5,j);
%         
% %         if (isnan(dat1(i,j))
% %             dat1(i,j) = interp1(dat1(
%     end
% end

a = tjrD;
a = fillmissing(a, 'linear');

b = tjrL;
b = fillmissing(b, 'linear');



figure(); hold on
set(gca,'FontSize',14)
[c1, h1] = contour(fpa0s, ha_tgts, tjrD',tjs,'r','LineWidth',2);
clabel(c1,h1,'Color','r');
[c2, h2] = contour(fpa0s, ha_tgts, b',tjs(3:end),'b','LineWidth',2);
clabel(c2,h2,'Color','b');
xlim([-10.5 fpa0s(1)])
xlabel('Entry Flight Path Angle')
ylabel('Apoapsis Target (km)')
legend('\beta_2/\beta_1 = 10','L/D = 0.25','location','sw')
xlim([-10.1 -8.5])


[c2, h2] = contour(fpa0s, v0s, tjrL(:,:,2,2)', tjs, 'b');
clabel(c2,h2);
% [cs3,h3] = contour(fpa0s, ha_tgts, tof_S2U', dt_contours, 'k');
% clabel(cs3,h3);
xlabel('Entry Flight-Path Angle (degrees)')
ylabel('Planet-Relative Velocity (km/s)')

% beta1
tjrL_57_1400 = reshape(tjrL(:,3,:,3),[20,10]);
tjrD_57_1400 = reshape(tjrD(:,3,:,3), [20,10]);


% figure(); hold on
% plot(fpa0s, tjrL_57_1400)
% plot(fpa0s, tjrD_57_1400)

%% Monte Carlo
clear;clc
ha_tgts = [1400 10000];
[x0,aero,gnc,sim,mc] = base_jsr();
sim.planet = 'mars';
x0.h0 = 100;
gnc.ha_tgt = 1400;
gnc.tj0 = 70;
gnc.dtj_lim = 20;
beta1 = 148;
ratio = 10;
x0.v0 = 5.7;
x0.fpa0 = -11;

% fpas = [-11.5 -9.5];
% LtD = 0.25;
cd = 1.2;
cl = 0;

aero.cds = [cd cd];
aero.cls = [cl cl];
aero.rn = 0.2;
% rc1 = aero.rn * 2.5;
m1 = 5500;
m2 = 3000;
rc1 = sqrt(m1/(aero.cds(1)*beta1*pi));
rc2 = sqrt(m2/(aero.cds(1)*beta1*ratio*pi));
aero.rcs = [rc1 rc2];
aero.m = [m1 m2];
sim.h_max = x0.h0;
sim.h_min = 0;
mc.sigs.cd = aero.cds(1) * 0.05;    %5%
gnc.npc_mode = uint8(1);
sim.t_max = 2500;
sim.efpa_flag = false;
sim.traj_rate = 250;
sim.parMode = true;

% mc.flag = false;
% mc.debug = true;

% gnc.ha_tgt = 10000;
% haf_err = gnc.ha_tgt;
% while (abs(haf_err) > 5)
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_err = out.haf_err;
%     if abs(haf_err) < 5
%         break;
%     else
%         x0.fpa0 = x0.fpa0 - 0.05;
%     end
% end
% keyboard;

out1 = run_dej_n(x0,gnc,aero,sim,mc);


for i = 1:length(ha_tgts)
    fprintf('Case %i\n', i);
    gnc.ha_tgt = ha_tgts(i);
%     x0.fpa0 = fpas(i);
    out(i) = run_dej_n(x0,gnc,aero,sim,mc);
end

keyboard;

out1 = out(1);
out2 = out(2);

% dat = load('C:\Users\Evan Roelke\Documents\research\mars_ac\data\mcHa_10deg.mat');
% out1 = dat.out1;
% out2 = dat.out2;

tols = 0:1:100; %error percent
Ni = length(tols);
pc1 = nan(Ni,1);
pc2 = pc1;

for i = 1:Ni  
    err1 = abs(out1.haf_err / ha_tgts(1)) * 100;
    err2 = abs(out2.haf_err / ha_tgts(2)) * 100;

    n1 = find(err1 <= tols(i));
    n2 = find(err2 <= tols(i));
    pc1(i) = length(n1) * 100/1000;
    pc2(i) = length(n2) * 100/1000;
end

figure(); hold on
grid on;
set(gca,'FontSize',14)
plot(tols, pc1, 'r','LineWidth',2)
plot(tols, pc2, 'b','LineWidth',2)
legend('1400 km Target','10000 km Target','Location','nw')
xlabel('Apoapsis Percent Error')
ylabel('Percent Captured')


