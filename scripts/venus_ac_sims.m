%% Venus aerocapture simulations
% Evan Roelke
% %  EsDL
% 
clear; clc;
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% ENTRY TRADES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% entry trade studies - feb 2019
%{
% clear; clc; close all;
% betas = [25 50 100 200 500];
% Ni = length(betas);
% 
% v_atms = linspace(10,15,10);
% Nj = length(v_atms);
% 
% % vehicle
% [x0,aero,gnc,sim,mc] = base_jsr();
% aref = pi*aero.rcs(1)^2;
% cd = aero.cds(1);
% b1 = aero.m(1)/(cd*aref);
% b2 = aero.m(2)/(cd*pi*aero.rcs(2)^2);
% 
% 
% fpa0 = -5.4; %initial guess
% 
% 
% fpa_tol = 0.025;
% 
% ha_tol = 25;        %km
% % ha_tgt = 10000;      %km
% ha_tgt = 2000;  % high LEO
% % ha_tgt = 35786; % ~GEO
% % ha_tgt = 400; %km
% % ha_tgt = 200; %km - low LEO
% 
% corr0 = [-10 -2];
% 
% in0 = dme_venus_in;
% temp = load('./data/atm_data/atm_venus_mc.mat');
% in0.p.atm.table = temp.nom_table;
% in0.s.traj.alt = x0.h0 * 1e3;
% in0.s.traj.az = 90*pi/180;
% in0.s.traj.gamma_pp = fpa0*pi/180;
% in0.v.aero.cd = cd;
% in0.v.aero.area_ref = aref;
% in0.v.gnc.g.p.dm_mode = uint8(0); % remove guidance
% in0.s.traj.rate = 100;
% in0.v.gnc.g.p.rate = 0;
% in0.v.mp.m_ini = aero.m(1);
% in0.v.aero.nose_radius = aero.rn;
% 
% % preallocte
% m = nan(Ni,1);
% fpa_f = nan(Ni,Nj);
% 
% for i = 1:Ni
%     m(i) = betas(i)*cd*aref;
%     for j = 1:length(v_atms)
%         in = in0;
%         % update input struct
%         in.s.traj.vel_pp_mag = v_atms(j)*1000;
%         in.v.mp.m_ini = m(i);
%         [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
%             in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
%             in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
%             in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);
%         % test sim
%         out = main1_mex(in);
%         
% %         [fpa_f(i,j),haf,haf_err,~] = fpa_given_haf( ...
% %             fpa0, ... %  initial fpa guess, (rad)
% %             ha_tgt, ... % target apoapse altitude, (m)
% %             in, ... % asim input file, contains plantary parameters, vehicle, and ICs
% %             fpa_tol, ... % tolerance on flight path angle to search for (rad)
% %             ha_tol); % tolerance on final target apoapse altitude (m)
% 
%         [fpa_f(i,j), ~,~,~,~] = find_fpa_given_ha(corr0, ha_tgt,  ...
%             in, fpa_tol, ha_tol);
% 
%     end %Nj
% end %Ni
% 
% 
% for i = 1:Ni
%     entries{i} = ['\beta = ' num2str(betas(i)) ' kg/m^2'];
% end
% 
% 
% figure(); hold on
% set(gca,'FontSize',14)
% title(['Venus Apoapsis Altitude Target: ' num2str(ha_tgt) ' km'],'FontSize',10)
% xlabel('Entry Velocity (km/s)')
% ylabel('Entry Flight Path Angle (deg)')
% plot(v_atms,fpa_f,'LineWidth',2)
% xlim([min(v_atms) max(v_atms)])
% legend(entries,'Location','ne','FontSize',10)
% 
% keyboard;
%}

%% beta ratio as a function of entry corridor,a couple different entry velocities
%{
clearvars -except save_path; clc;
v_atm = [10,11,12];
Ni = length(v_atm);
ha_tgt = [2000 10000];
Nj = length(ha_tgt);

m1 = 80;
% m1 = 150;
Nk = 10;
ratios = linspace(2,20,Nk);

rc1 = 0.75;
rc2 = 0.175;
aero.rc = [rc1 rc2];
aero.cd = 1.05;
% rc1 = 1;
% rc2 = 0.2;
m = nan(2,Nk);
% fpa_corr = nan(Ni,Nj,Nk); %[v h b]
corr1 = nan(Ni,Nk);
corr2 = corr1;
run = 0;
for i = 1:Ni %v_atms
%     for j = 1:Nj % ha_tgts
    for k = 1:Nk %beta ratios

        m(:,k) = [m1, m1*ratios(k)*((rc2/rc1)^2)];
        aero.m = m(:,k);

        percent = (run/(Ni*Nk))*100;
        fprintf('Run %i/%i, %i/%i, %i/%i, %2.1f percent.\n',i,Ni,Nj,k,Nk,percent);

        % 2k target
        corr1(i,k) = ac_venus_find_corr(aero,v_atm(i),ha_tgt(1));
        % 10k target
        corr2(i,k) = ac_venus_find_corr(aero,v_atm(i),ha_tgt(2));
        
        run = run+1;
    end
%     end
end

save([save_path 'entry_trades/beta_vs_corr_width.mat'],'corr1','corr2','v_atm','ha_tgt','m','ratios','aero');


% dat = load('beta_vs_corr_width.mat');

keyboard;

% corr1 = nan(Ni,Nk);
% corr2 = corr1;
% for i = 1:Ni %v
%     for k = 1:Nk %b
%     corr1(i,k) = fpa_corr(i,1,k); %2k
%     corr2(i,k) = fpa_corr(i,2,k); %10k
%     end
% end

entries = cell(Ni,1);
for i = 1:Ni
    entries{i} = ['v_{0} = ' num2str(v_atm(i)) ' km/s'];
end

figure(); hold on
set(gca,'FontSize',14);
title('2000 km target')
plot(corr1, ratios,'LineWidth',2)
legend(entries,'location','nw')
xlabel('Corridor Width (degrees)')
ylabel('\beta_{2}/\beta_{1}')

figure(); hold on
set(gca,'FontSize',14);
title('10000 km target');
plot(corr2, ratios,'LineWidth',2)
legend(entries,'location','nw')
xlabel('Corridor Width (degrees)')
ylabel('\beta_{2}/\beta_{1}')

figure(); hold on
set(gca,'FontSize',14);
title('Target Apoapsis Difference');
plot(corr2-corr1, ratios,'LineWidth',2)
legend(entries,'location','ne')
xlabel('\Delta Corridor Width (degrees)')
ylabel('\beta_{2}/\beta_{1}')
%}

%% optimal jettison time(s)
%{
clear;clc;
fpa_bounds = [-6.3, -7.2];
% ha_tgt = 2000;
ratios = 4:4:20;
% ratios = [10];
Nj = length(ratios);
m2 = nan(Nj,1);

efpas = linspace(fpa_bounds(1),fpa_bounds(2),20)';
Ni = length(efpas);
tj = nan(Ni,Nj); haf_err = tj; tjr = tj;
tj0 = 120;
tol = 5;    % stop within x km
veh.m = [150 60];
veh.r = [1 0.2];
veh.rn = 0.175/2;
veh.cd = 1.05;

haf_tol = 100;
[x0,aero,gnc,sim,mc] = base_venus_ac;
% gnc.ha_tgt = 400;
mc.flag = false;
sim.t_max = 2000;
sim.traj_rate = 1000;
sim.data_rate = 1000;
sim.h_min = 35;
x0.v0 = 15;
gnc.iters = uint8(5);
gnc.npc_mode = uint8(0);
gnc.guid_rate = 1;
gnc.ha_tgt = 105500;

entries = cell(Nj,1);
% for k = 1:Nk %tgts
for j = 1:Nj %betas
    fprintf('Outer Loop %i\n',j);
%     m2(j) = ratios(j) * veh.m(1) * (veh.r(2)/veh.r(1))^2;
%     veh.m(2) = m2(j);

    m2(j) = ratios(j) * aero.m(1) * (aero.rcs(2)/aero.rcs(1))^2;
    aero.m(2) = m2(j);
    
    for i = 1:Ni %fpa
%         tj_init = tj0 - 5*(i-1);
%         tj_init = tj0;
%         [tj(i,j),tjr(i,j),haf_err(i,j),~] = find_opt_tj(ha_tgt,efpas(i),tj_init,tol,veh);
        
        x0.fpa0 = efpas(i);
        out = run_dej_n(x0,gnc,aero,sim,mc);
        tj(i,j) = out.t_jett(1);
        tjr(i,j) = tj(i,j) / out.traj.time(out.idxend);
        haf_err(i,j) = out.haf_err;
        
        if (tj(i,j) < 0 || abs(haf_err(i,j)) > haf_tol)
            tj(i,j) = nan;
            tjr(i,j) = nan;
            haf_err(i,j) = nan;
        end
    end
end

% keyboard;
save(['../venus_ac/entry_trades/opt_t_jett/v' num2str(x0.v0) '_'  ... 
    num2str(gnc.ha_tgt/1000) 'k.mat'], ... 
    'tjr','tj','efpas','haf_err','ratios');
% end
keyboard;
% for j = 1:Nj
%     entries{j} = ['\beta_2/\beta_1 = ' num2str(ratios(j))];
% end

% figure(); hold on
% set(gca,'FontSize',14)
% % yyaxis left
% plot(efpas, tjr,'LineWidth',2)
% ylabel('Jettison Time (s)')
% % yyaxis right
% plot(efpas, haf_err,'LineWidth',2)
% ylabel('Apoapsis Error (km)')
% ylim([-tol-2 tol+2])
% legend(entries,'location','nw')
% xlabel('Entry Flight-Path Angle (deg)')
% ylim([0.1 0.6])

% keyboard;

%%%%% post process
p = '..\venus_ac\entry_trades\opt_t_jett\';

% dat = load([p 'v11_0.4k.mat']);
dat = load([p 'v15_0.4k.mat']);

[err,tj,tjr] = brush_tj_opts(dat.haf_err,dat.tj, ... 
    dat.tjr,dat.ratios,'linear',50);

Nj = length(dat.ratios);
betas = cell(Nj,1);
for j = 1:Nj
    betas{j} = ['\beta_2/\beta_1 = ' num2str(dat.ratios(j))];
end

t = '15 km/s, 400 km tgt';

%tjr
figure(); hold on
set(gca,'FontSize',14)
grid on;
title(t);
plot(dat.efpas, tjr,'LineWidth',2)
legend(betas, 'location','nw')
xlabel('EFPA (deg)')
ylabel('t/t_f')

%tj
figure(); hold on
set(gca,'FontSize',14)
grid on;
title(t);
plot(dat.efpas, tj,'LineWidth',2)
legend(betas, 'location','nw')
xlabel('EFPA (deg)')
ylabel('Time (s)')

keyboard;
%}

%% beta ratio targeting specific apopapsis as a function of entry FPA
%{
% 2/28/18
% clear; clc
% ha_tgts = [400 1000 2000].*1e3;
% 
% % fpa corridors (10.64543 km/s v_atm)
%     % 400 km: -5.2016 to -5.6428 deg
%     % 1000 km: -5.1953 to -5.6371 deg
%     % 2000 km: -5.1847 to -5.6278 deg
% 
%     % with venus GRAM atomsphere
%     % 400 km: -5.1945 to -5.6265
%     % 1000 km: -5.1884 to -5.6199
%     % 2000 km: -5.1784 to -5.6095
% 
% Ni = length(ha_tgts);
% Nj = 65;
% Nk = 75;
% 
% 
% v_atm = 10.64543e3; %m/s
% h_atm = 150e3;  %m
% 
% fpa_atm = linspace(-5.17, -5.65, Nj)';
% ratios = linspace(2.5, 20, Nk)';
% 
% rcs = [0.75 0.175];
% m1 = 80;
% 
% m = nan(Nk,2);
% run = 0;
% for i = 1:Ni        %ha_tgts
%     for j = 1:Nj        %fpas
%         x_atm = [h_atm,v_atm,fpa_atm(j),0,0,90];
%         
%         for k = 1:Nk        %betas
%            percent = (run/(Ni*Nj*Nk))*100;
%            fprintf('Run %i/%i, %i/%i, %i/%i. %2.1f percent.\n',i,Ni,j,Nj,k,Nk,percent);            
%             
%             m(k,:) = [m1, m1*ratios(k)*( rcs(2)/rcs(1) )^2];
%             
%             [out(i,j,k), ~, ~, beta_ratio(i,j,k)] = run_aerocapture(...
%                 [0 1 0], 'venus', 3, ...
%                 rcs, [1.05 1.05], [0 0], [0 0], [45 45], 0.5.*rcs, m(k,:), 0,  ...    %spacecraft geometry
%                 x_atm, 0, ...                   %entry state vector and uncertanties
%                 0, ha_tgts(i), 1e3, 0.01, ...
%                 2, 0, 0, 0);
%             
%             run = run + 1;
%         end
%     end
% end
% 
% 
% save('beta_vs_ha_trajs_course.mat')
% 
% 
% 
% dat = load('beta_vs_ha_trajs_course.mat');
% 
% for i = 1:Ni
%     for j = 1:Nj
%         for k = 1:Nk
%             if isnan(out(i,j,k).nom.orbit.haf) ~= 1
%                 haf(i,j,k) = out(i,j,k).nom.orbit.haf;
%             else
%                 haf(i,j,k) = nan;
%         end
%     end
% end
%}


%% probabilty of capture vs. beta ratio
%{
clear; clc
% ha_tgt = 2000e3;
ha_tgts = [400 2000 10000 35000].*1000;

ha_tol = 10e3;  %10 km
% haf_corr = [ha_tgt - ha_tol, ha_tgt + ha_tol]./1000;
v_atm = 10.65e3;
fpa_atm = mean([-5.1945, -5.6265]);

x_atm = [150e3, v_atm, fpa_atm, 0,0,90];
rcs = [0.75 0.175];
ratios = [2 4 6 8 10 12 14 16 18 20];
m1 = 80;

sigs = [0.25, 0.03, 0.03, 0, 10, 50, 0.005, 0.005, 0.005, 0.005];

iter = 7;

Ni = length(ratios);
m2 = nan(Ni,1);
% percent = nan(Ni,1);

N_mc = 1000;

for j = 1:length(ha_tgts)

for i = 1:Ni
    fprintf('Iter Run %i\n',i);
    
    m2(i) = m1*ratios(i)*(rcs(2)/rcs(1))^2;
    
    % run sim
    out = ....
        venus_sej_auto(m2(i), fpa_atm,v_atm,ha_tgts(j), ...
        0,iter,N_mc,sigs,ha_tol);
    
    captured = find(out.haf >=  ha_tgts(j) - ha_tol & ...
        out.haf <= ha_tgts(j) + ha_tol);
    if ~isempty(captured)
        percent(j,i) = length(captured)/length(out.haf) * 100;
    else
        percent(j,i) = 0;
    end
    
end

end

save 'venus_ac/capture_prob_betas.mat'

figure(); hold on
grid on
set(gca,'FontSize',14)
plot(ratios,percent,'LineWidth',2)
xlabel('\beta_{2}/\beta_{1}')
ylabel('Capture Rate (%)')
legend('Target Apoapsis: 400 km', ... 
    'Target Apoapsis: 2000 km', ... 
    'Target Apoapsis: 10000 km', ... 
    'Target Apoapsis: 35000 km', ... 
    'location','nw');
xlim([min(ratios) max(ratios)])
%}

%{
figure(); hold on
plot(ratios,haf./1000,'LineWidth',2)
plot(ratios,haf_corr(1).*ones(1,Ni),'--k','LineWidth',1.5)
plot(ratios,haf_corr(2).*ones(1,Ni),'--k','LineWidth',1.5)
xlabel('\beta_{2}/\beta_{1}')
ylabel('Final r_{a} (km)')

figure(); hold on
for i = 1:Ni
    plot(dat(i).traj.vel_pp_mag./1000,dat(i).traj.alt./1000, 'LineWidth', 1.5) 
end
xlabel('Velocity (m/s)')
ylabel('Altitude (km)');
%}

%% run probability of capture vs. iter_max (2000 km+ target) - 4/28/18
%       note that ha_tol doesn't appear to affect results...
%{
clear;clc;

% iters = [5 7 9 11 13 15]; % iter_max vals
iters = 5:5:30;
% iters = 5;

ha_tgt = 2000e3;
ha_tol = 10e3;  %10 km orbit tolerance
v_atm = 10.65e3; %m/s
h_atm = 150e3;  %m
fpa_atm = mean([-5.2016, -5.6428]); %mean of fpa corr
% m = [80 40];
m1 = 72.5;
m2 = 37.8;
rcs = [0.75 0.1];
rns = 0.5.*rcs;

N_mc = 1;

% state vectors
% x_atm = [h_atm,v_atm,fpa_atm,0,0,90];
sigs = [0.25, 0.03, 0.03, 0, 10, 50, 0.005, 0.005, 0.005, 0.005];

percent = nan(length(iters),1);

for i = 1:length(iters)
    
%     if mod(i,5) == 0
        fprintf('Iter Run %i\n',i);
%     end
    
    % run sim
    out = ...
        ac_venus_auto_jettison(m2, fpa_atm,v_atm,ha_tgt, ... 
        0,iters(i),N_mc,sigs,ha_tol);
    
    haf(i) = out.haf;
    
%     % calculate capture percentage
%     captured = find(out(i).haf >=  ha_tgt - ha_tol & ... 
%         out(i).haf <= ha_tgt + ha_tol);
%     if ~isempty(captured)
%         percent(i) = length(captured)/length(out(i).haf) * 100;
%     else
%         percent(i) = 0;
%     end
    
end

% keyboard;
% figure(); hold on
% hist(out(6).haf,100);
% plot((ha_tgt - ha_tol).*ones(length(out(6).haf),1),... 
%     linspace(0,max(max(out(6).haf)),length(out(6).haf)))
% plot((ha_tgt + ha_tol).*ones(length(out(6).haf),1),linspace(0,max(max(out(6).haf,length(out(6).haf)))))

figure(); hold on
grid on
set(gca,'FontSize',14)
% plot(iters,percent)
plot(iters,haf./1000,'LineWidth',2)
legend('Target 2000 km','location','se')
xlabel('Max Iterations')
ylabel('Final Apoapsis (km)')

% figure(); hold on
% grid on
% hist(out.haf(i)./1000,100)
% 
% iter = 1;
% for i = 1:1000
%     if out.haf(i) > (ha_tgt - 2*ha_tol) && out.haf(i) < (ha_tgt + 2*ha_tol)
%         x(iter) = out.haf(i)/1000;
%         iter = iter + 1;
%     end
% end
% 
% figure(); hold on
% set(gca,'FontSize',14)
% hist(x,100);
% xlabel('Apoapsis Altitude (km)')
% ylabel('Frequency')
% legend('10km Tol = 24.4% Capture','location','nw')
%}

%% probability of capture vs. haf tolerance - 4/28/18
%{
clear;clc
tols = [0 20 50 100 200 500 1000].*1000;

ha_tgt = 2000e3;

ha_tgts = [400 2000 10000 35000].*1000;

% ha_tol = 50e3;  %10 km orbit tolerance
v_atm = 10.65e3; %m/s
h_atm = 150e3;  %m
fpa_atm = mean([-5.2016, -5.6428]); %mean of fpa corr
% m = [80 40];
m1 = 80;
m2 = 40;
rcs = [0.75 0.175];
rns = 0.5.*rcs;

N_mc = 1000;

% state vectors
% x_atm = [h_atm,v_atm,fpa_atm,0,0,90];
sigs = [0.25, 0.03, 0.03, 0, 10, 50, 0.005, 0.005, 0.005, 0.005];

iter = 5;

for j = 1:length(ha_tgts)

for i = 1:length(tols)
    
    fprintf('Iter Run %i\n',i);
    
    % run sim
    out = ....
        venus_sej_auto(m2, fpa_atm,v_atm,ha_tgts(j), ...
        0,iter,N_mc,sigs,tols(i));
    
    captured = find(out.haf >=  ha_tgts(j) - tols(i) & ...
        out.haf <= ha_tgts(j) + tols(i));
    if ~isempty(captured)
        percent(j,i) = length(captured)/length(out.haf) * 100;
    else
        percent(j,i) = 0;
    end
end

end


save 'venus_ac/capture_prob_tols.mat'
% % % a = load('capture_prob.mat');


figure(); hold on
grid on
set(gca,'FontSize',14)
plot(tols./1000,percent,'LineWidth',2)
xlabel('Apoapsis Tolerance (km)')
ylabel('Capture Percentage')
legend( ['Target r_{a}: ' num2str(ha_tgts(1)./1000) ' km'], ... 
    ['Target r_{a}: ' num2str(ha_tgts(2)./1000) ' km'], ... 
    ['Target r_{a}: ' num2str(ha_tgts(3)./1000) ' km'], ... 
    ['Target r_{a}: ' num2str(ha_tgts(4)./1000) ' km'], ... 
    'location','nw')

%}

%% differences in atmospheres - 4/28/18
%{
temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_JPL.mat');
atm_table = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_mc.mat');
orig = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus.mat');

%
figure(); hold on
plot(atm_table.nom_table(:,2:4),atm_table.nom_table(:,1)./1000,'LineWidth',2)
legend('Dens','Pres','Temp','location','ne')
ylim([90 150])

figure(); hold on
plot(atm_table.nom_table(:,2) - temp.table(:,2),temp.table(:,1)./1000)

figure(); hold on
plot(temp.table(:,2:4),temp.table(:,1)./1000,'LineWidth',2)
legend('Dens','Pres','Temp','location','ne')
ylim([90 150])

figure(); hold on
plot(orig.table(:,2:4),orig.table(:,1)./1000,'LineWidth',2)
legend('Dens','Pres','Temp','location','ne')
ylim([90 150])

% plot alt vs. dens
figure(); hold on
plot(temp.table(:,2),temp.table(:,1)./1000,'b','LineWidth',2)
plot(atm_table.nom_table(:,2),atm_table.nom_table(:,1)./1000,'--k','LineWidth',2)
plot(orig.table(:,2),orig.table(:,1)./1000,'-.r','LineWidth',2)
legend('JPL','VenusGRAM','Original')
xlabel('Density')
ylabel('Altitude')

figure(); hold on
plot(temp.table(:,4)-orig.table(:,4),temp.table(:,1)./1000,'b','LineWidth',2)
plot(atm_table.nom_table(:,4)-orig.table(:,4),atm_table.nom_table(:,1)./1000,'--k','LineWidth',2)
% plot(orig.table(:,2),orig.table(:,1)./1000,'-.r','LineWidth',2)
legend('JPL','VenusGRAM','location','ne')
xlabel('Temperature Difference')
ylabel('Altitude')
ylim([90 200])

% plot alt vs. pressure
figure(); hold on
plot(temp.table(:,1)./1000,temp.table(:,3),'b','LineWidth',2)
plot(atm_table.nom_table(:,1)./1000,atm_table.nom_table(:,3),'--k','LineWidth',2)
plot(orig.table(:,1)./1000,orig.table(:,3),'-.r','LineWidth',2)
legend('JPL','VenusGRAM','Original')

figure(); hold on
plot(temp.table(:,1)./1000,temp.table(:,4),'b','LineWidth',2)
plot(atm_table.nom_table(:,1)./1000,atm_table.nom_table.table(:,4),'--k','LineWidth',2)
legend('JPL','VenusGRAM')

figure(); hold on
plot(temp.table(:,1)./1000,temp.table(:,5),'b','LineWidth',2)
plot(atm_table.nom_table(:,1)./1000,atm_table.nom_table(:,5),'--k','LineWidth',2)
legend('JPL','VenusGRAM')

figure(); hold on
plot(temp.table(:,1)./1000,temp.table(:,6),'b','LineWidth',2)
plot(atm_table.nom_table(:,1)./1000,atm_table.nom_table(:,6),'--k','LineWidth',2)
legend('JPL','VenusGRAM')


figure(); hold on
plot(orig.table(:,2:7),orig.table(:,1)./1000,'LineWidth',2)
legend('Density','Pressure','Temp', 'Wind E','Wind N','Wind U','location','best')

%}

%% probability of capture vs. corrector time steps (Jun 2018)
%{
clear; clc
ha_tol = 10e3;
ha_tgts = [400 2000 10000 35000].*1000;

v_atm = 10.65e3; %m/s
h_atm = 150e3;  %m
fpa_atm = -5.5; %mean of fpa corr
% m = [80 40];
m1 = 80;
m2 = 40;
rcs = [0.75 0.175];
rns = 0.5.*rcs(1);

N_mc = 1000;

% state vectors
% x_atm = [h_atm,v_atm,fpa_atm,0,0,90];
sigs = [0.25, 0.03, 0.03, 0, 10, 50, 0.005, 0.005, 0.005, 0.005];

iter = 5;

t_incs = [0.01 2];

t_steps = [ 1    2; ...
            2    4; ...
            4    8; ...
            8    10; ... 
            10   20; ...
            20   40];

debug = 1;
        
for j = 1:length(ha_tgts)

for i = 1:size(t_steps,1)
    
    fprintf('Iter Run %i\n',i);
    
    % run sim
    out = ...
        venus_sej_auto(m2, fpa_atm,v_atm,ha_tgts(j), ...
        0,iter,N_mc,sigs,ha_tol,t_incs,t_steps(:,i),debug);
    

    
    capt = find(out.haf >=  ha_tgts(j) - ha_tol & ...
        out.haf <= ha_tgts(j) + ha_tol);
    if ~isempty(capt)
        pc(j,i) = length(capt)/length(out.haf) * 100;
    else
        pc(j,i) = 0;
    end
  
end

end

save 'venus_ac/data/capture_tSteps.mat'
%}

%% basic NPC case
%{
clear;clc

N_mc = 1;
ha_tgt = 2000e3;
% Ni = length(ha_tgts);

ha_tol = 1e3;  % orbit tolerance

v_atm = 10.65e3; %m/s
h_atm = 150e3;  %m
fpa_atm = -5.4;
% m = [80 40];
m1 = 72.5;
m2 = 37.8;
rcs = [0.75 0.1];
rns = 0.5.*rcs;
sigs = [0.25, 0.03, 0.03, 0, 10, 50, 0.005, 0.005, 0.005, 0.005];

iter = 5;

dflag = uint8(0);  % debug mode

t_incs = [0.01 2];  % npc time increments
t_steps = [10 20];
% for i = 1:Ni
%     fprintf('Run %i\n',i);
    
    % run sim
    out = ....
        venus_sej_auto(m2,fpa_atm,v_atm,ha_tgt, ...
        0,iter,N_mc,sigs,ha_tol,t_incs, t_steps, 0.1, dflag);
    
% end

captured = find(out.haf >=  ha_tgt - ha_tol & ...
    out.haf <= ha_tgt + ha_tol);

if ~isempty(captured)
    pc = length(captured)/length(out.haf) * 100;
else
    pc = 0;
end

err = (out.haf - ha_tgt)./1000;
fprintf('Apoapse Error: %4.2f km\n',err)

figure(); hold on
plot(out.traj.vel_pp_mag./1000,out.traj.alt./1000,'LineWidth',2)
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')

% figure(); hold on
% set(gca,'FontSize',14)
% histogram(out.haf./1000,N_mc/2);
% % xlim( [ha_tgt/1000 - 10, ha_tgt/1000 + 10] )
% legend('Target Apoapsis: 35000 km','location','nw')


% figure(); hold on
% for i = 1:N_mc
%     plot(out.traj(i).vel,out.traj(i).alt)
% end
% 


% find time of peak accel


% figure(); hold on
% plot(out.traj.time,out.traj.fpa_atm.*180/pi)
% plot(out.traj.time,out.traj.g)
% legend('FPA','G-Loading','location','se')
%}

%% opt jettison time vs. error (Sep 2018)
%{
clear;clc
% temp = load('C:\Users\Evan Roelke\Documents\research\asim\venus_ac\ac_guid_algorithms\DCF\t_exact_dt10.mat');
% t0 = temp.t_jett;
% m = [150 60];
% rcs = [1 0.2];

m = [80 40];
rcs =[0.75 0.175];

ra_tgt = [2e3 10e3];
ra_tol = 50;
n = 1;  % single jett
% fpas = temp.fpas;
Ni = 20;
fpas = linspace(-5.2,-5.6,Ni);
guid_rate = 0.5;    %guidacne rate, Hz
mc.flag = false;
mc.sigs = 0;

traj_rate = 250;

t_2k = nan(Ni,1);
t_10k=t_2k;t_2k = t_2k;t_10k = t_2k;
err_2k=t_2k;err_10k=t_2k;err_2k=t_2k;err_10k=t_2k;
err1_n1 = t_2k; err2_n1=t_2k;
for i = 1:Ni
    if mod(i,Ni/2)==0
        fprintf('50 percent done\n')
    end
    t0 = 150-(5*i); % initial jettison time guess
    n1 = dej_n_auto_venus(fpas(i),ra_tgt(1),n,t0(i),rcs,m,1,guid_rate,mc,1,traj_rate);
    n2 = dej_n_auto_venus(fpas(i),ra_tgt(2),n,t0(i),rcs,m,1,guid_rate,mc,1,traj_rate);
    
    
    t_2k(i) = n1.t_jett;
    t_10k(i) = n2.t_jett;
    err_2k(i) = n1.haf_err;
    err_10k(i) = n2.haf_err;
%     err1_n1(i) = n1.g.dej_n.delta_ap(n1.idj)/1000;    % still bad - newton method diverged
%     b1 = dej_n_auto_venus(fpas(i),ra_tgt(1),ra_tol,n,t0(i),rcs,m,0,guid_rate);
%     b2 = dej_n_auto_venus(fpas(i),ra_tgt(2),ra_tol,n,t0(i),rcs,m,0,guid_rate);
    
%     t_2k(i) = b1.t_jett;
%     t_10k(i) = b2.t_jett;
%     err_2k(i) = b1.haf_err;
%     err_10k(i) = b2.haf_err;
end

% diff = err1_n1 - err2_n1;
% save('venus_ac/ac_guid_algorithms/NPC/nominal_errs.mat','t_2k','t_10k','err_2k','err_10k','fpas')

figure(); hold on
grid on
set(gca,'FontSize',12)
title('Nominal NPC Errors for 2000 and 10000 km r_{a}')
yyaxis left
plot(fpas,err_2k,'b-o','LineWidth',1.5)
plot(fpas,err_10k,'b--o','LineWidth',1.5)
ylim([-100 150])
ylabel('Apoapsis Error (km)')
yyaxis right
plot(fpas,t_2k,'r-o','LineWidth',1.5)
plot(fpas,t_10k,'r--o','LineWidth',1.5)
ylabel('Jettison Time (s)')
xlabel('Entry FPA (deg)')
legend('r_{a}=2000km','r_{a}=10000km','r_{a}=2000km','r_{a}=10000km','location','nw')
%}

%% comparison of bisection method and Newton method (Sep 2018)
%{
% single jettison comparison
clear;clc
m = [80 40];
rcs = [0.75 0.175];
ra_tgt = 10e3;
ra_tol = 25;
n = 1;  % single jett
guid_rate = 0.1;
tj0 = 90;  %initial guess, sec
mc.flag = uint8(0);
ob = dej_n_auto_venus(-5.4,ra_tgt,ra_tol,n,tj0,rcs,m,0,guid_rate,mc);
on = dej_n_auto_venus(-5.4,ra_tgt,ra_tol,n,tj0,rcs,m,1,guid_rate,mc);
% plot(on.traj.vel_pp_mag./1000,on.traj.alt./1000);


% % check difference between guidance err and sim err
% figure(); hold on
% title('Comparison of Newton Method and Bisection Method for NPC')
% plot(on.traj.time,on.g.dej_n.delta_ap./1000,'b','LineWidth',1.5)
% plot(on.traj.time,on.haf_err.*ones(length(on.traj.time),1),'--k','LineWidth',1.5)
% xlim([0 150])
% ylim([-10 30])
% legend(['Diff = ' num2str(on.haf_err - on.g.dej_n.delta_ap(on.idj)/1000) ' km'])

ymin = -1;
ymax = 2;
x = on.g.dej_n.ncalls(on.idj)+1;
figure(); hold on
set(gca,'FontSize',14)
plot(ob.g.dej_n.ncalls,ob.g.dej_n.delta_ap/1000 / 1e4,'-.r','LineWidth',2);
plot(on.g.dej_n.ncalls,on.g.dej_n.delta_ap/1000 / 1e4,'b','LineWidth',2);
plot([x x],[ymin ymax],'--k','LineWidth',1.5)
xlabel(['Number of Guidance Iterations, (' num2str(guid_rate) ' Hz)'])
ylabel('Target Apoapsis Error (10^{4} km)')
% ax1 = gca;
% ax1.XColor = 'r';
% ax1.YColor = 'r';
% ax1_pos = ax1.Position;
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
% xlim([0 round(on.traj.time(on.idxend),0)])
% xlabel('Time (sec)')
ylim([ymin,ymax])
legend(['Bisection Method: \Deltar_{a}=' num2str(abs(ob.g.dej_n.delta_ap(ob.idxend)/1000)) ' km'], ... 
    ['Newton Method: \Deltar_{a}=' num2str(abs(on.g.dej_n.delta_ap(on.idxend)/1000)) ' km'])



% plot(on.g.dej_n.ncalls,(on.g.dej_n.dr_opt - ob.g.dej_n.dr_opt)/1000,'b');

% 
% ylim1 = min(ob.g.dej_n.delta_ap)/1000;
% ylim2 = max(ob.g.dej_n.delta_ap)/1000;
% ylim = ([ylim1 ylim2]);
% figure(); hold on
% plot(ob.traj.time,ob.g.dej_n.delta_ap./1000)
% plot([ob.traj.time(ob.idj) ob.traj.time(ob.idj)],[ylim1 ob.g.dej_n.delta_ap(ob.idj)/1000],'k','LineWidth',2)
% 
% 
% figure(); hold on
% plot(on.traj.vel_pp_mag./1000,on.traj.alt./1000,'b')
% plot(ob.traj.vel_pp_mag./1000,ob.traj.alt./1000,'r')


% observe error as func of jettison time


% figure(); hold on
% yyaxis right
% plot(out.traj.time,out.g.dej_n.stage)
% yyaxis left
% plot(out.traj.time,out.g.dej_n.dr_opt./1000)
% % plot(out.traj.time,out.g.dej_n.delta_ap./1000)
% ylim([-100 100])



b2k = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\NPC\1stage\data\mcBisection_fpas_hafs_5hz.mat')
tols = [100 500];
num = find(abs(b2k.haf_b_2k(:,2)-2000) <= tols(1));
pc = length(num)*100/1000
num = find(abs(b2k.haf_b_2k(:,2)-2000) <= tols(2));
pc = length(num)*100/1000

n2k = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\NPC\1stage\data\mcNewton_fpas_hafs_5hz.mat');
num = find(abs(n2k.haf_n_2k(:,2)-2000) <= tols(1));
pc = length(num)*100/1000
num = find(abs(n2k.haf_n_2k(:,2)-2000) <= tols(2));
pc = length(num)*100/1000

histogram(n2k.haf_n_2k(:,2)-2000,'NumBins',50);
xlim([-750 750])

%}

%% multiple jettison events (Sep-Oct 2018)
%{
clear;clc
m = [80 45 30];
rcs = [0.75 0.2 0.1];
ra_tgt = 10e3;
ra_tol = 25;
n = 2;  % single jett
tj0 = [70 120];  %initial guess, sec
guid_rate = 0.2;

on = dej_n_auto_venus(-5.4,ra_tgt,ra_tol,n,tj0,rcs,m,0,guid_rate);

figure(); hold on
y_min = -10000;
% y_max = 1000;
% ylim([y_min y_max])
xlim([0, on.traj.time(on.idxend)])
plot(on.traj.time,on.g.dej_n.delta_ap./1000,'LineWidth',2)
plot([on.traj.time(on.idj(1)), on.traj.time(on.idj(1))], ... 
    [y_min 10000],'--k','LineWidth',1.5)
plot([on.traj.time(on.idj(2)), on.traj.time(on.idj(2))], ... 
    [y_min 10000],'--k','LineWidth',1.5)
% text(5,-5,'Jettison 1')

% plot(on.traj.time, on.g.dej_n.dr_opt)
%}



%% DCF monte carlo vs. fpa - Nov 2018
%{
ha_tgt = 2000;
% m = [80 40];
% rcs = [0.75 0.175];
m = [150 60];
rcs = [1 0.2];

dt = 10;
fpas = [-5.4 -5.6];

Ni = length(fpas);

mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = default_sigs();
tols = [10 100 500 1000];
Nj = length(tols);

pc = nan(Ni,Nj);
traj_rate = 250;
for i = 1:Ni
    out(i) = ac_venus_dcf(fpas(i), traj_rate, m, rcs, ha_tgt, dt, mc);
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        pc(i,j) = length(num)*100/mc.N;
    end
end

keyboard;

haf_errs = [out(1).haf_err, out(2).haf_err];
hafs = [out(1).haf, out(2).haf];
dvs = [out(1).dv, out(2).dv];

% save results
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';

save([save_path '/dej_n/dcfj/data/mc/mc_fpa2k.mat'], ...
    'tols','pc','fpas','haf_errs','dvs','hafs')

load([save_path '/dej_n/dcfj/data/mc/mc_fpa2k.mat']);

% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Hist DCF, h_{a,tgt} = 10000 km,Traj Rate 100 Hz','FontSize',10)
histogram(hafs(:,1)./1000,'NumBins',50,'FaceColor','r')
histogram(hafs(:,2)./1000,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Altitude (10^3 km)')
ylabel('Occurence')
legend('EFPA = -5.4^{o}','EFPA = -5.6^{o}','location','ne')

% save figure
% saveas(gca, ... 
%     'C:\Users\Evan\Documents\Academics\Research\Venus_ac\DCF\mc_fpa_3_5.fig')
%}

%% DCF monte carlo vs. ha tgt - Nov 2018
%{
ha_tgt = [2000 10000];
m = [80 40];
rcs = [0.75 0.175];

dt = 10;
% fpas = [-5.3 -5.5];
efpa = -5.4;

Ni = length(ha_tgt);

mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
tols = [10 100 500 1000];
Nj = length(tols);

traj_rate = 500;
pc = nan(Ni,Nj);

for i = 1:Ni
    out(i) = ac_venus_dcf(efpa, traj_rate, m, rcs, ha_tgt(i), dt, mc);
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        pc(i,j) = length(num)*100/mc.N;
    end
end
keyboard;

% save results
save('C:\Users\Evan\Documents\Academics\Research\Venus_ac\DCF\percent_ha.mat', ...
    'tols','pc','fpas')

load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\DCFJ\data\dcf_mc_ha.mat');


% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Hist DCF, EFPA = -5.4^{o}, Traj Rate 100 Hz','FontSize',10)
histogram(out(1).haf,'NumBins',50,'FaceColor','r')
histogram(out(2).haf,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Altitude (km)')
ylabel('Occurence')
legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km','location','ne')
% save figure
saveas(gca, ... 
    'C:\Users\Evan\Documents\Academics\Research\Venus_ac\DCF\mc_ha.fig')
%}


%% PDG nominal errors vs. fpa - fidelity checks
%{
clear; clc; close all;
fpas = round(linspace(-5.4,-5.8,19),2)';
m = [150 60];
rcs = [1 0.2];
ha_tgts = [2000 10000];
trigger = uint8(0);
Ni = length(fpas);
Nj = length(ha_tgts);
traj_rate = 100;
gnc_rate = 10;   % gnc call rate, hz
t2g = traj_rate/gnc_rate;

ha_tol = 100;          % lower tol and check performance
tj = nan(Ni,Nj);
err = tj;
haf = tj; dv = tj;
mc.flag = uint8(0);
for j = 1:Nj % ha tgts
    fprintf('Target %i\n',j);
    for i  = 1:Ni % fpas
        out = ac_venus_pd(fpas(i),ha_tgts(j),ha_tol,m,rcs,traj_rate,gnc_rate,trigger,mc);
        tj(i,j) = out.tjett;    %s
        haf(i,j) = out.haf; %km
        dv(i,j) = out.dv;   %m/s
        err(i,j) = out.haf_err; %km
    end
end

figure(); hold on
grid on
set(gca,'FontSize',14)
title('Nominal Apoapsis Errors, \beta_2/\beta_1 = 9.18, GNC Rate = 5 Hz, Traj/GNC Ratio = 100')
plot(fpas,err,'-o','LineWidth',2)
xlabel('EFPA (deg)')
ylabel('Apoapsis Error (km)')
legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km', 'location', 'best')
ylim([-250 100])
% saveas(gca, ... 
%     'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\pdg_err_ha_low_fidel.fig')
%}
%% PDG monte carlo analysis as func of FPA and tgt AP - Energy Trigger
%{
% energy trigger, fpa performance
clear; clc; close all
fpas = [-5.4 -5.6];
% fpas = -5.8;
m = [150 60];
rcs = [1 0.2];
ha_tgt = [2 10].*1e3;
Nk = length(ha_tgt);
Ni = length(fpas);
traj_rate = 100;
gnc_rate = 0.1;
ha_tol = 10;
mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = default_sigs();
tols = [10 50 100 250 500 1000];
Nj = length(tols);
pc_2k = nan(Ni,Nj);
pc_10k = pc_2k;
trigger = uint8(0); % energy trigger

haf_2k = nan(mc.N,Ni);
haf_10k = haf_2k;
dv_2k = haf_2k;
dv_10k = haf_2k;
for i = 1:Ni % fpas
    fprintf('PDG FPA Sim %i\n',i);
    out1 = ac_venus_pd(fpas(i),ha_tgt(1),ha_tol,m,rcs, ... 
        traj_rate,gnc_rate,trigger,mc);
%     out2 = ac_venus_pd(fpas(i),ha_tgt(2),ha_tol,m,rcs, ... 
%         traj_rate,gnc_rate,trigger,mc);
    for j = 1:Nj
        num = find(abs(out1.haf_err) <= tols(j));
        pc_2k(i,j) = length(num)*100/mc.N;
%         num = find(abs(out2.haf_err) <= tols(j));
%         pc_10k(i,j) = length(num)*100/mc.N;
    end
    
    haf_2k(:,i) = out1.haf;
    dv_2k(:,i) = out1.dv;
%     
%     haf_10k(:,i) = out2.haf;
%     dv_10k(:,i) = out2.dv;
end

keyboard;
% save results
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';
save([save_path 'dej_n\PDG\data\mc_fpaBoth_0.1hz.mat'], ...
    'tols','pc_2k','haf_2k','dv_2k','fpas')
keyboard;


a = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\PDG\data\mc_fpa54_56_5hz.mat');
b = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\PDG\data\mc_fpa58_5hz.mat');
% histogram of results - diff fpa, 2k target
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo, h_{a,tgt} = 2000 km,GNC Rate=5hz Hz,Traj Rate 100 Hz')
h3 = histogram(b.haf,'NumBins',50,'FaceColor','g');
h2 = histogram(a.haf_2k(:,2),'NumBins',50,'FaceColor','b');
h1 = histogram(a.haf_2k(:,1),'NumBins',50,'FaceColor','r');
xlabel('Apoapsis Altitude (km)')
ylabel('Occurence')
legend([h1,h2,h3],'EFPA = -5.4^{o}', ... 
    'EFPA = -5.6^{o}', ... 
    'EFPA = -5.8^{o}', ... 
    'location','ne')

% histogram of results - diff fpa, 10k target
figure(); hold on; grid on
set(gca,'FontSize',14)
title(['Monte Carlo, h_{a,tgt} = 10000 km,GNC Rate ' num2str(gnc_rate) ' Hz,Traj Rate 100 Hz'])
histogram(haf_10k(:,2)./1000,'NumBins',50,'FaceColor','r')
histogram(haf_10k(:,1)./1000,'NumBins',20,'FaceColor','b')
xlabel('Apoapsis Altitude (10^3 km)')
ylabel('Occurence')
legend('EFPA = -5.6^{o}','EFPA = -5.4^{o}','location','nw')

% histogram - diff ha_tgt, same fpa (-5.4)
figure(); hold on; grid on
set(gca,'FontSize',14)
title(['Monte Carlo Histogram, EFPA=-5.4, GNC Rate ' num2str(gnc_rate) ' Hz,Traj Rate 100 Hz'])
histogram(haf_2k(:,1) - 2000,'NumBins',50,'FaceColor','r')
histogram(haf_10k(:,1) - 10000,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('2000 km Target','10000 km Target','location','nw')
% histogram - diff ha_tgt, same fpa (-5.6)
figure(); hold on; grid on
set(gca,'FontSize',14)
title(['Monte Carlo Histogram, EFPA=-5.6, GNC Rate ' num2str(gnc_rate) ' Hz,Traj Rate 100 Hz'])
histogram(haf_10k(:,2) - 10000,'NumBins',50,'FaceColor','b')
histogram(haf_2k(:,2) - 2000,'NumBins',50,'FaceColor','r')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('10000 km Target','2000 km Target','location','ne')

% save figure
% saveas(gca, ... 
%     'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\monte_carlo\pdg_mc_fpa_2_4.fig')
%}

%% PDG monte carlo as func of ha_tgt (Oct 2018)
%{
clear; clc; close all
efpa = -5.4;
ha_tgt = [2000];
Ni = length(ha_tgt);
traj_rate = 100;
gnc_rate = 5;
haf_tol = 10;
mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = default_sigs();
tols = [10 50 100 250 500];
Nj = length(tols);
percent = nan(Ni,Nj);
haf_med = nan(Ni,1);
haf_std = nan(Ni,1);
trigger = uint8(0);
m = [150 60];
rcs = [1 0.2];
for i = 1:Ni % ha
    fprintf('Starting PDG Case %i\n',i);
    out = ac_venus_pd(efpa,ha_tgt(i),haf_tol,m,rcs,traj_rate,gnc_rate,trigger,mc);
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        percent(i,j) = length(num)*100/mc.N;
    end
    haf_med(i) = median(out(i).haf);
    haf_std(i) = std(out(i).haf);
end

keyboard;
% save capture results
% save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\PDG\monte_carlo\percent_ha.mat', ... 
%     'tols','percent','haf_med','haf_std')

% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, EFPA = -5.4 Hz, GNC Rate = 5 Hz')
h1 = histogram(haf_10k(:,1) -10000,'NumBins',50,'FaceColor','b');
h2 = histogram(haf_2k(:,1) - 2000,'NumBins',50,'FaceColor','r');
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend([h2, h1],'h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km','location','nw')

figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, EFPA = -5.6 Hz, GNC Rate = 5 Hz')
h1 = histogram(haf_10k(:,2) -10000,'NumBins',50,'FaceColor','b');
h2 = histogram(haf_2k(:,2) - 2000,'NumBins',50,'FaceColor','r');
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend([h2, h1],'h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km','location','nw')

% save figure
saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\monte_carlo\pdg_mc_ha.fig')
%}
%% PDG error vs. GNC rate
%{
m = [150 60];
rcs = [1 0.2];
ha_tgt = 2000;
ha_tol = 1;
trigger = uint8(0);
efpa = -5.4;
mc.flag = false;

a = 0.1;
b = 5000;
n = 25;
gnc_rates = nan(n,1);
growth = (b/a)^(1/(n - 1));
for i = 1:n
    gnc_rates(i) = a*growth^(i-1);
end
Ni = length(gnc_rates);
traj_rate0 = 1000;

haf_err = nan(Ni,1);
dv = haf_err;
tjr = haf_err;
for i = 1:Ni
    if gnc_rates(i) > traj_rate0
        traj_rate = round(gnc_rates(i),0);
    else
        traj_rate = traj_rate0;
    end
    out = ac_venus_pd(efpa,ha_tgt,ha_tol,m,rcs,traj_rate,gnc_rates(i),trigger,mc);
    haf_err(i) = out.haf_err;
    dv(i) = out.dv;
    tjr(i) = out.tjett / out.traj.time(out.idxend);
    
%     keyboard
end

save([save_path 'dej_n/PDG/data/err2k_vs_gncRate_fpa5.4.mat'], ...
    'dv','haf_err','tjr','efpa','ha_tgt');

ind = 14;
fig = figure(); hold on
set(gca,'FontSize',14)
title('2k tgt, -5.4deg, PDG, Venus, \beta_2/\beta_1=10')
yyaxis left
plot(gnc_rates(1:ind),haf_err(1:ind),'LineWidth',2)
ylabel('Nominal Apoapsis Altitude Error (km)')
ylim([-200 5])
yyaxis right
plot(gnc_rates(1:ind),tjr(1:ind)-tjr(ind),'LineWidth',2)
ylabel('\Delta t_{jett}/t_f')
ylim([-0.025 0])
xlabel('GNC Rate (Hz)')
xlim([0 gnc_rates(ind)])
%}

%%%%%%%%%%%%%%%
% NPC STUFF - SEJ
%%%%%%%%%%%%%%%
%% Basic SEJ event
%{
clear;clc;
[x0,aero,gnc,sim,mc] = base_jsr();
b21 = 10;
aero.m = [150 60];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / (aero.m(1) * b21) );

gnc.guid_rate = 5;
gnc.npc_mode = uint8(1);
gnc.tj0 = 90;
gnc.ha_tgt = 10000;
sim.traj_rate = 250;
mc.N = 1000;
gnc.dtj_lim = 10;
gnc.hynap = true;
sim.parMode = true;
sim.nWorkers = 6;
x0.fpa0 = -5.6;
fprintf('SEJ Run 1\n');
out = run_dej_n(x0,gnc,aero,sim,mc);
% 
keyboard;

save(['C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\HYDRA\1stage\data\mc_n_' ... 
    num2str(gnc.ha_tgt/1000) 'k_' ... 
    num2str(gnc.guid_rate) 'hz_' ... 
    num2str(abs(x0.fpa0)*10) 'deg.mat']);
% 
% figure();
% histogram(out.haf_err,'NumBins',50);
% xlim([-2000 2000])
% 
keyboard;
%}




%% High Fidelity NPC newton/bisection method monte carlo analysis as func of ha_tgt, fpa
% guidance rate of 5 Hz
%{
clear;clc;close all;
ha_tgt = [2 10].*1000;
fpas = [-5.4 -5.6];
[x0,aero,gnc,sim,mc] = base_jsr();
gnc.guid_rate = 0.5;
sim.parMode = true;
sim.t_max = 2000;


Ni = length(fpas);

mc.N = 1000;

% preallcate
haf_b_2k = nan(mc.N,Ni); haf_b_10k = haf_b_2k;
haf_n_2k = haf_b_2k; haf_n_10k = haf_b_2k;
dv_b_2k = haf_b_2k; dv_n_2k = haf_b_2k;
dv_b_10k = haf_b_2k; dv_n_10k = haf_b_2k;

% mc.debug = true;

for i = 1:Ni
    x0.fpa0 = fpas(i);
    
    gnc.ha_tgt = 2000;
    fprintf('FPA %4.2f, Target Ap %4.2f\n',fpas(i), gnc.ha_tgt);
%     gnc.npc_mode = uint8(0);
%     out_b_2k = run_dej_n(x0,gnc,aero,sim,mc);
    gnc.npc_mode = uint8(1);
    out_n_2k = run_dej_n(x0,gnc,aero,sim,mc);
    
    gnc.ha_tgt = 10000;
    fprintf('FPA %4.2f, Target Ap %4.2f\n',fpas(i), gnc.ha_tgt);
%     gnc.npc_mode = uint8(0);
%     out_b_10k = run_dej_n(x0,gnc,aero,sim,mc);
    gnc.npc_mode = uint8(1);
    out_n_10k = run_dej_n(x0,gnc,aero,sim,mc);    
    
%     haf_b_2k(:,i) = out_b_2k.haf;
%     haf_b_10k(:,i) = out_b_10k.haf;
    haf_n_2k(:,i) = out_n_2k.haf;
    haf_n_10k(:,i) = out_n_10k.haf;   
    
%     dv_b_2k(:,i) = out_b_2k.dv;
%     dv_b_10k(:,i) = out_b_10k.dv;
    dv_n_2k(:,i) = out_n_2k.dv;
    dv_n_10k(:,i) = out_n_10k.dv;
    
%     for j = 1:Nj
%         num = find(abs(out(i).haf_err) <= tols(j));
%         percent(i,j) = length(num)*100/mc.N;
%     end
%     haf_med(i) = median(out(i).haf);
%     haf_std(i) = std(out(i).haf);
end
keyboard;

% toc
% save capture results
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\HYDRA\';
save([save_path '1stage\data\mcNewton_fpas_hafs_' num2str(gnc.guid_rate) 'hz.mat'], ... 
    'haf_n_2k','haf_n_10k','dv_n_2k','dv_n_10k','fpas','ha_tgt','x0','gnc');

% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title(['NPC, EFPA=-5.4^{o},GNC rate = ' num2str(gnc.guid_rate) 'Hz, ha_tgt=' num2str(gnc.ha_tgt)])
histogram(haf_b_2k(:,1),'NumBins',50,'FaceColor','r')
histogram(haf_n_2k(:,1),'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Altitude (km)')
ylabel('Occurence')
legend('Bisection Method','Newton Method','location','nw')
% % save figure
% % saveas(gca, ... 
% %     'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\NPC\monte_carlo\npc_mc_ha_5Hz.fig')


n = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\NPC\1stage\data\mcNewton_fpas_hafs_0.5hz.mat');
b = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\NPC\1stage\data\mcBisection_fpas_hafs_0.5hz.mat');

figure(); hold on
histogram(n.haf_n_2k(:,1)/1000,'NumBins',50,'FaceColor','b')
histogram(b.haf_b_2k(:,1)/1000,'NumBins',50,'FaceColor','r')
xlim([0 4])


%% GNC Rate NPC monte carlo as func of iter_max - compare methods (Nov 2018)
%{
% guidance rate of 0.1 Hz
clear;clc;close all;
m = [80 40];
fpa_atm = -5.4;
rcs = [0.75 0.175];
ha_tgt = 10000;
ha_tol = 25;
n = 1;  % single jett
% root_mode = 1;  % newton method
tj0 = 85;  %initial guess, sec
sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
% % % [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
mc.flag = uint8(1);
mc.sigs = sigs;
mc.N = 1000;
iters = [1 3 5];    % num iterations per guidance call
traj_rate = 100;   %1000 Hz
tols = [10 50 100 250 500];

Nj = length(tols);
Ni = length(iters);
guid_rate = 0.1;  % guidance rate

percentn = nan(Ni,Nj);
haf_med_n = nan(Ni,1);
haf_std_n = nan(Ni,1);
percentb = nan(Ni,Nj);
haf_med_b = nan(Ni,1);
haf_std_b = nan(Ni,1);
tic
for i = 1:length(iters)
    fprintf('NPC Simulation %i\n',i);
    out_n(i) = dej_n_auto_venus(fpa_atm, ha_tgt, ha_tol, ... 
        n, tj0, rcs, m, 1, guid_rate,mc,iters(i),traj_rate);
    out_b(i) = dej_n_auto_venus(fpa_atm, ha_tgt, ha_tol, ... 
        n, tj0, rcs, m, 0, guid_rate,mc,iters(i),traj_rate);
    for j = 1:Nj
        num = find(abs(out_n(i).haf_err) <= tols(j));
        percentn(i,j) = length(num)*100/mc.N;
        num = find(abs(out_b(i).haf_err) <= tols(j));
        percentb(i,j) = length(num)*100/mc.N;
    end
    haf_med_n(i) = median(out_n(i).haf);
    haf_std_n(i) = std(out_n(i).haf);
    haf_med_b(i) = median(out_b(i).haf);
    haf_std_b(i) = std(out_b(i).haf);
end
toc
% save capture results
save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\NPC\monte_carlo\percent_iter_0_1Hz.mat', ... 
    'tols','percentn','percentb','haf_med_n','haf_med_b','haf_std_n','haf_std_b')

% histogram of results - newton
figure(); hold on; grid on
set(gca,'FontSize',14)
title({'Monte Carlo Histogram, Newton Method' ... 
    'EFPA = -5.4^{o},h_{a,tgt} = 10000 km, GNC Rate=0.1 Hz'},'FontSize',10)
histogram(out_n(1).haf_err,'NumBins',50,'FaceColor','r')  %1 iter
histogram(out_n(2).haf_err,'NumBins',50,'FaceColor','b')  %3 iters
histogram(out_n(3).haf_err,'NumBins',50,'FaceColor','g')  %5 iters
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('Iter_{max} = 1','Iter_{max} = 3','Iter_{max} = 5','location','nw')
% save figure
saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\NPC\monte_carlo\npc_n_mc_iter_0_1Hz.fig')

% histogram of results - bisection
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, Bisection Method, EFPA = -5.4^{o},h_{a,tgt} = 10000 km,GNC Rate=0.1 Hz')
histogram(out_b(1).haf_err,'NumBins',50,'FaceColor','r')  %2k tgt
histogram(out_b(2).haf_err,'NumBins',50,'FaceColor','b')  %10k tgt
histogram(out_b(3).haf_err,'NumBins',50,'FaceColor','g')  %30k tgt
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('Iter_{max} = 1','Iter_{max} = 3','Iter_{max} = 5','location','ne')
% save figure
saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\NPC\monte_carlo\npc_b_mc_iter0_1Hz.fig')
%}

%% High GNC Rate NPC monte carlo as func of iter_max, compare methods
%{
% guidance rate of 5 Hz
clear;clc;close all;
m = [80 40];
fpa_atm = -5.4;
rcs = [0.75 0.175];
ha_tgt = 10000;
ha_tol = 25;
n = 1;  % single jett
% root_mode = 1;  % newton method
tj0 = 85;  %initial guess, sec
sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
% % % [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
mc.flag = uint8(1);
mc.sigs = sigs;
mc.N = 1000;
iters = [1 3 5];    % num iterations
traj_rate = 100;   %1000 Hz
tols = [10 50 100 250 500];

Nj = length(tols);
Ni = length(iters);
guid_rate = 5;

percentn = nan(Ni,Nj);
haf_med_n = nan(Ni,1);
haf_std_n = nan(Ni,1);
percentb = nan(Ni,Nj);
haf_med_b = nan(Ni,1);
haf_std_b = nan(Ni,1);
tic
for i = 1:length(iters)
    fprintf('NPC Simulation %i\n',i);
    out_n(i) = dej_n_auto_venus(fpa_atm, ha_tgt, ha_tol, ... 
        n, tj0, rcs, m, 1, guid_rate,mc,iters(i),traj_rate);
    out_b(i) = dej_n_auto_venus(fpa_atm, ha_tgt, ha_tol, ... 
        n, tj0, rcs, m, 0, guid_rate,mc,iters(i),traj_rate);
    for j = 1:Nj
        num = find(abs(out_n(i).haf_err) <= tols(j));
        percentn(i,j) = length(num)*100/mc.N;
        num = find(abs(out_b(i).haf_err) <= tols(j));
        percentb(i,j) = length(num)*100/mc.N;
    end
    haf_med_n(i) = median(out_n(i).haf);
    haf_std_n(i) = std(out_n(i).haf);
    haf_med_b(i) = median(out_b(i).haf);
    haf_std_b(i) = std(out_b(i).haf);
end
toc
% save capture results
save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\NPC\monte_carlo\percent_iter_5Hz.mat', ... 
    'tols','percentn''percentb','haf_med_n','haf_med_b','haf_std_n','haf_std_b')

% histogram of results - newton
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, Newton Method, EFPA = -5.4^{o},h_{tgt} = 10k km,GNC Rate=5 Hz')
histogram(out_n(1).haf_err,'NumBins',50,'FaceColor','r')  %2k tgt
histogram(out_n(2).haf_err,'NumBins',50,'FaceColor','b')  %10k tgt
histogram(out_n(3).haf_err,'NumBins',50,'FaceColor','g')  %30k tgt
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('Iter_{max} = 1','Iter_{max} = 3','Iter_{max} = 5','location','nw')
% save figure
saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\NPC\monte_carlo\npc_n_mc_iter_5Hz.fig')

% histogram of results - bisection
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, Bisection Method, EFPA = -5.4^{o}, h_{a,tgt} = 10000 km,GNC Rate=5 Hz')
histogram(out_b(1).haf_err,'NumBins',50,'FaceColor','r')  %2k tgt
histogram(out_b(2).haf_err,'NumBins',50,'FaceColor','b')  %10k tgt
histogram(out_b(3).haf_err,'NumBins',50,'FaceColor','g')  %30k tgt
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('Iter_{max} = 1','Iter_{max} = 3','Iter_{max} = 5','location','nw')
% save figure
saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\NPC\monte_carlo\npc_b_mc_iter5Hz.fig')

%}

%%%%% to do (11/11/18)
%%%%% check if removing the <= error check improves newton method
%%% line 123 in guid_dej_n
%%%%% one simulation with K_dens every time step instead of func call?


%% probability of capture as func of beta ratio for each func (Nov 2018)
%{
% ratios = [6 7 9 11];
ratio = 10;
m1 = 150;
rc_ratio = .175/.75;
rc1 = 1;
% rc2 = rc1*rc_ratio;
rc2 = 0.2;
rcs = [rc1 rc2];
m2 = ratio.*(m1*rcs(2)^2/rcs(1)^2);
cd = 1.05;
v0 = 11;
m = [m1 m2];

Ni = length(m2);    %outer loop num

tols = (10:10:1000)';
Nj = length(tols);
pcNPCn = nan(Ni,Nj);
pcPDG = pcNPCn;
pcDCF = pcNPCn;
pcNPCb = pcNPCn;

sigs = default_sigs();

% sim settings
efpas = [-5.4 -5.6];
Ni = length(efpas);

for i = 1:2
    b(i) = m(i)/(cd * pi * rcs(i)^2);
end

ha_tgt = 10000;
ha_tol = 5;
nj = 1;
tj0 = 90;
mc.sigs = sigs;
mc.flag = uint8(1);
trigger = uint8(0);
dt = 10;
mc.N = 800;

gnc_rate = 0.5;
traj_rate = 100;
iters = 1;

pc_n2k = nan(Ni,Nj);
pc_n10k = pc_n2k;
pc_b2k = pc_n2k;
pc_b10k = pc_n2k;

haf_n_med2k = nan(Ni,1);
haf_n_std2k = haf_n_med2k; 
haf_n_med10k = haf_n_med2k; 
haf_n_std10k = haf_n_med2k; 
dv_n_med2k = haf_n_med2k; dv_n_med10k = haf_n_med2k;
dv_n_std2k = haf_n_med2k; dv_n_std10k = haf_n_med2k;
haf_b_med2k = nan(Ni,1);
haf_b_std2k = haf_b_med2k; 
haf_b_med10k = haf_b_med2k; 
haf_b_std10k = haf_b_med2k; 
dv_b_med2k = haf_b_med2k; dv_b_med10k = haf_b_med2k;
dv_b_std2k = haf_b_med2k; dv_b_std10k = haf_b_med2k;
for i = 1:Ni
    % display progress update
    c = clock;
    fprintf('Beta Ratio %i: %4.0fh:%4.0fm:%4.0fs\n',i,c(4),c(5),c(6));
    
    out2k = dej_n_auto_venus(efpas(i), v0, 2000, ha_tol, ... 
        nj, tj0, rcs, m, 1, gnc_rate,mc,iters,traj_rate);
    out10k = dej_n_auto_venus(efpas(i), v0, 10000, ha_tol, ... 
        nj, tj0, rcs, m, 1, gnc_rate,mc,iters,traj_rate);
%     out_npcb = dej_n_auto_venus(efpa, ha_tgt, ha_tol, ... 
%         n, tj0, rcs, m(i,:), 0, gnc_rate,mc,iters,traj_rate);
%     out_pdg = ac_venus_pd(efpa,ha_tgt,ha_tol,m(i,:),rcs, ... 
%         traj_rate,gnc_rate,trigger,mc);
%     out_dcf = ac_venus_dcf(efpa,traj_rate,m(i,:),rcs, ... 
%         ha_tgt,dt,mc);

    b2k = dej_n_auto_venus(efpas(i), 2000, ha_tol, ... 
        nj, tj0+10, rcs, m, 0, gnc_rate,mc,iters,traj_rate);
    b10k = dej_n_auto_venus(efpas(i), 10000, ha_tol, ... 
        nj, tj0, rcs, m, 0, gnc_rate,mc,iters,traj_rate);
    
    % get capture percentages for this beta ratio
    for j = 1:Nj
       num = find(abs(out2k.haf_err) <= tols(j));
       pc_n2k(j) = length(num)*100/mc.N;
       num = find(abs(out10k.haf_err) <= tols(j));
       pc_n10k(j) = length(num)*100/mc.N;
       
       
       num = find(abs(b2k.haf_err) <= tols(j));
       pc_b2k(j) = length(num)*100/mc.N;
       num = find(abs(b10k.haf_err) <= tols(j));
       pc_b10k(j) = length(num)*100/mc.N;
       
%        nNPCb = find(abs(out_npcb.haf_err) <= tols(j));
%        pcNPCb(i,j) = length(nNPCb)*100/mc.N;
       
%        nPDG = find(abs(out_pdg.haf_err) <= tols(j));
%        pcPDG(i,j) = length(nPDG)*100/mc.N;
       
%        nDCF = find(abs(out_dcf.haf_err) <= tols(j));
%        pcDCF(i,j) = length(nDCF)*100/mc.N;
        
    end
    
%     haf_med_DCF(i) = median(out_dcf.haf);
%     haf_std_DCF(i) = std(out_dcf.haf);
    
    haf_n_med2k = median(out2k.haf);
    haf_n_std2k = std(out2k.haf);
    haf_n_med10k = median(out10k.haf);
    haf_n_std10k = std(out10k.haf);
    haf_b_med2k = median(b2k.haf);
    haf_b_std2k = std(b2k.haf);
    haf_b_med10k = median(b10k.haf);
    haf_b_std10k = std(b10k.haf);
    
    dv_n_med2k = median(out2k.dv);
    dv_n_med10k = median(out10k.dv);
    dv_n_std2k = std(out2k.dv);
    dv_n_std10k = std(out10k.dv);
    dv_b_med2k = median(b2k.dv);
    dv_b_med10k = median(b10k.dv);
    dv_b_std2k = std(b2k.dv);
    dv_b_std10k = std(b10k.dv);
    
    save([save_path 'dej_n\NPC\1stage\data\mcBoth_' num2str(abs(efpas(i))) ...
        'deg_comp_' num2str(gnc_rate) 'hz.mat'], ...
     'haf_n_med2k','haf_n_med10k','haf_n_std2k','haf_n_std10k', ... 
     'dv_n_med2k','dv_n_med10k','dv_n_std2k','dv_n_std10k', ... 
     'haf_b_med2k','haf_b_med10k','haf_b_std2k','haf_b_std10k', ... 
     'dv_b_med2k','dv_b_med10k','dv_b_std2k','dv_b_std10k', ...
     'pc_n2k','pc_n10k','pc_b2k','pc_b10k', ...
     'm','efpas','gnc_rate','rcs');
    
%     haf_med_NPCb(i) = median(out_npcb.haf);
%     haf_std_NPCb(i) = std(out_npcb.haf);

%     haf_med_PDG(i) = median(out_pdg.haf);
%     haf_std_PDG(i) = std(out_pdg.haf);
end

% save DCF capture results
% save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\DCF\mc_percent_beta', ... 
%     'ratios','pcDCF','efpa','gnc_rate','m1','rcs','haf_med_DCF','haf_std_DCF')

% save PDG capture results
% save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\PDG\mc_percent_beta', ... 
%     'ratios','pcPDG','efpa','gnc_rate','m1','rcs','haf_med_PDG','haf_std_PDG')

% % save NPC capture results
% save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\NPC\mc_percent_beta', ... 
%     'ratios','pcNPCn','efpa','gnc_rate','m1','rcs','haf_med_NPCn','haf_std_NPCn', ... 
%     'haf_std_NPCb','haf_med_NPCb')

% save([save_path 'dej_n\NPC\1stage\data\mc_bothTgts_bothEfpa_newton_beta10.mat'], ...
%     'haf_med2k','haf_med10k','haf_std2k','haf_std10k', ... 
%     'dv_med2k','dv_med10k','dv_std2k','dv_std10k', ... 
%     'out2k','out10k')
keyboard;

figure(); hold on
grid on;
set(gca,'FontSize',14)
title('Percent vs. Beta Ratio, Traj Rate 100 Hz, GNC Rate 5 Hz, 10k Target, Newton Method','FontSize',10)
plot(tols,pcNPCn,'LineWidth',2)
xlabel('Apoapsis Altitude Error Tolerance (km)')
ylabel('Percent Captured')
legend('\beta_2/\beta_1 = 6','\beta_2/\beta_1 = 7', ... 
    '\beta_2/\beta_1 = 9','\beta_2/\beta_1 = 11', ... 
    'location','se')
%}

%}

%% NPC check decoupled K_dens - Dec 2018
%{
clear;clc;close all;
m = [80 40];
fpa_atm = -5.4;
% fpa_atm = -5.2;
rcs = [0.75 0.175];
% ha_tgt = [2 10].*1e3;    % ha targets
ha_tgt = 2000;
ha_tol = 25;
n = 1;  % single jett
root_mode = 1;  % newton method
tj0 = 85;  %initial guess, sec
sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
% % % [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
mc.flag = uint8(1);
mc.sigs = sigs;
mc.N = 1;
iters = 1;      % single iteration
traj_rate = 100;   %integration rate
tols = [10 50 100 250 500];

Ni = length(ha_tgt);
Nj = length(tols);
guid_rate = 0.5;

percent = nan(Ni,Nj);
haf_med = nan(Ni,1);
haf_std = nan(Ni,1);

out = dej_n_auto_venus(fpa_atm, ha_tgt, ha_tol, ... 
    n, tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
for j = 1:Nj
    num = find(abs(out.haf_err) <= tols(j));
    percent(j) = length(num)*100/mc.N;
end
haf_med = median(out.haf);
haf_std = std(out.haf);


min_ind = find(out.traj.alt == min(out.traj.alt));
Rvar = nan(min_ind,1);
rss = 0;
for i = 1:min_ind
    Rvar(i) = out.traj.rho(i)/out.nom.traj.rho(i);  % actual variation
    diff(i) = out.g.K_dens(i) - Rvar(i);
    rss = rss + ( diff(i)^2 );
end

figure(); hold on
grid on
set(gca,'FontSize',14)
title(['Density Corrector Performance, GNC Rate = ' num2str(guid_rate) ... 
    ', RSS = ' num2str(round(rss,3)) ... 
    ', Kgain = 0.2'],'FontSize',10)
plot(Rvar, out.traj.alt(1:min_ind),'k','LineWidth',2);
plot(out.g.K_dens(1:min_ind), out.traj.alt(1:min_ind),'b','LineWidth',2);
ylim([100 122]);
xlabel('Density Variation')
ylabel('Altitude (km)')
legend('Actual Variation','Variation Estimate','location','se')


% large K_dens spike at jettison time step!!
%}
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 1-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clearvars -except save_path; clc; close all;
% ha_tgt = [.4 2 10 105.5].*1000;
ha_tgt = 2000;
% fpas = [-6.5 -6.55 -6.6 -6.65];
fpas = -5.4;
[x0,aero,gnc,sim,mc] = base_venus_ac();
gnc.guid_rate = 0.5;
gnc.tj0 = 80;
gnc.hydra_flag = true;
gnc.dtj_lim = 10;
sim.parMode = true;
sim.nWorkers = 10;
sim.t_max = 2000;
sim.h_min = 25;
sim.traj_rate = 50;
sim.data_rate = 10;

x0.v0 = 11;

Ni = length(fpas);
Nj = length(ha_tgt);

mc.N = 500;

% mc.debug = true;

for i = 1:Ni %fpas
    for j = 1:Nj %ha tgts
        x0.fpa0 = fpas(i);
        gnc.ha_tgt = ha_tgt(j);
        
        fprintf('FPA %4.2f, Target Ap %4.2f\n',fpas(i), gnc.ha_tgt);
        gnc.npc_mode = uint8(0);
        out_b = run_dej_n(x0,gnc,aero,sim,mc);
        gnc.npc_mode = uint8(1);
        out_h = run_dej_n(x0,gnc,aero,sim,mc);
        
%         if (gnc.ha_tgt < 1000)
%             save([save_path 'dej_n\HYDRA\1stage\data\' num2str(floor(gnc.ha_tgt)) ...
%                 '\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg.mat'])
%         else
%             save([save_path 'dej_n\HYDRA\1stage\data\' num2str(floor(gnc.ha_tgt/1000)) ...
%                 'k\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg.mat'])
%         end
        
    end %j
end %i

out_stats(out_b)
out_stats(out_h)
fprintf('Finished 1stage sim. '); print_current_time();
keyboard;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 2-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regular 2stage
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac();

b31 = 10;
b21s = [3 5 7 9];

x0.v0 = 11;
efpas = [-5.4 -5.6];
Nk = length(efpas);

gnc.ha_tgt = 10000;
gnc.ha_tol = 10;
gnc.n = 2;
tj0s = [90 170;70 140];
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);

for j = 1:length(b21s)
b21 = b21s(j);
aero.m = [150 60 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * b21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * b31) );

if (aero.rcs(2) > aero.rcs(1)) || (aero.rcs(3) > aero.rcs(2))
    error('Impossible Geometry')
end

b = nan(3,1);   %betas
br = nan(2,1);  %beta ratios
for i = 1:3
    b(i) = aero.m(i) / (aero.cds(1) * pi * aero.rcs(i)^2);
end
for i = 1:2
    br(i) = b(i+1)/b(i);
end
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br31 = b(3)/b(1);

aero.rn = 0.1;
aero.cds = [1.05 1.05 1.05];
aero.cls = [0 0 0];

sim.traj_rate = 200;
sim.planet = 'venus';
sim.efpa_flag = false;
    
mc.debug = false;

mc.flag = true;
sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage, b2/b1 = %2.0f, ha* = %4.0f\n',br21, gnc.ha_tgt);
gnc.force_jett = true;
gnc.comp_curr = false;

for k = 1:Nk
    x0.fpa0 = efpas(k);
    gnc.tj0 = tj0s(k,:);
    out = run_dej_n(x0,gnc,aero,sim,mc);
    
    betas = ['b_ijs' num2str(b21) '_' num2str(b31)];

    if (gnc.ha_tgt < 1000)
        save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt)) ... 
            '\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '.mat'])
    else
        save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt/1000)) ... 
            'k\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '.mat'])
    end
end %fpas
end %b21s
fprintf('Finished 2stage sims. '); print_current_time();
keyboard;


% % % % post process
dat = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\v11_5.4deg_b21=3_b31=10.mat');
out = dat.out;
haf1 = nan(1000,1); haf2 = haf1;
for i = 1:1000
    if (~isnan(out.tjett(i,2)))
        haf2(i) = out.haf_err(i);
        haf1(i) = nan;
    else
        haf1(i) = out.haf_err(i);
        haf2(i) = nan;
    end
end
figure(); hold on
% plot(out.tjett(:,1),haf1,'*b','LineWidth',1.5)
% plot(out.tjett(:,2),haf2,'*r','LineWidth',1.5)
h2=histogram(haf2,'NumBins',50,'FaceColor','b');
h1=histogram(haf1,'NumBins',50,'FaceColor','r');
width = 20;
% h1.Normalization = 'probability';
% h2.Normalization = 'probability';
h1.BinWidth = width;
h2.BinWidth = width;
legend([h1 h2],'Jettison 1','Jettison 2')
%}

%% FJ & CC
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac();

b31 = 10;
b21s = [3 5 7 9];

x0.v0 = 15;
x0.fpa0 = -6.5;

gnc.ha_tgt = 400;
gnc.ha_tol = 5;
gnc.n = 2;
gnc.tj0 = [80 120];
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);

for j = 1:length(b21s)
b21 = b21s(j);
    
% aero.rcs = [1 0.45 0.15];
aero.m = [150 60 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * b21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * b31) );

if (aero.rcs(2) > aero.rcs(1)) || (aero.rcs(3) > aero.rcs(2))
    error('Impossible Geometry')
end

b = nan(3,1);   %betas
br = nan(2,1);  %beta ratios
for i = 1:3
    b(i) = aero.m(i) / (aero.cds(1) * pi * aero.rcs(i)^2);
end
for i = 1:2
    br(i) = b(i+1)/b(i);
end
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br31 = b(3)/b(1);

aero.rn = 0.1;
aero.cds = [1.05 1.05 1.05];
aero.cls = [0 0 0];

mc.flag = true;

sim.traj_rate = 200;
sim.planet = 'venus';
sim.efpa_flag = false;
    
mc.debug = false;
mc.flag = true;
sim.parMode = true;
gnc.dtj_lim = 10;

% x0.fpa0 = fpas;
fprintf('2stage, b2/b1 = %2.0f, ha* = %4.0f\n',br21, gnc.ha_tgt);
gnc.force_jett = false;
gnc.comp_curr = false;
out1 = run_dej_n(x0,gnc,aero,sim,mc);

fprintf('2stage, b2/b1 = %2.0f, ha* = %4.0f, FJ\n',br21, gnc.ha_tgt);
gnc.force_jett = true;
gnc.comp_curr = false;
out2 = run_dej_n(x0,gnc,aero,sim,mc);

fprintf('2stage, b2/b1 = %2.0f, ha* = %4.0f, CC\n',br21, gnc.ha_tgt);
gnc.force_jett = false;
gnc.comp_curr = true;
out3 = run_dej_n(x0,gnc,aero,sim,mc);

fprintf('2stage, b2/b1 = %2.0f, ha* = %4.0f, FJ & CC\n',br21, gnc.ha_tgt);
gnc.force_jett = true;
gnc.comp_curr = true;
out4 = run_dej_n(x0,gnc,aero,sim,mc);

betas = ['b21=' num2str(br21) '_b31=' num2str(b31)];

% keyboard;

if (gnc.ha_tgt < 1000)
    save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt)) ... 
        '\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '_FJ_CC.mat'])
else
    save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt/1000)) ... 
        'k\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '_FJ_CC.mat'])
end

end

fprintf('Finished 2stage sims. '); print_current_time();
keyboard;
%}

%% 2stage Post Process
%{
% dat = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\hydra\mc2k_fpa5_4_0.5Hz, b31=10_b21=3_b32=3.3333.mat');
dat = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\rss_data\v11_5.4deg_b21=3_b31=10.mat');

err = dat.out.haf_err;
tj = dat.out.tjett;
err1 = nan(1000,1);
err2 = err1;
for i = 1:1000
%     err1(ind1) = err(i);
%     tj1(ind1) = tj(i,1);
%     ind1 = ind1 + 1;
%     
%     if (~isnan(tj(i,2)))
%         err2(ind2) = err(i);
%         tj2(ind2) = tj(i,2);
%         ind2 = ind2 + 1;
%     end

    if (isnan(tj(i,2)))
        err1(i) = err(i);
%         ind1 = ind1 + 1;
    else
        err1(i) = nan;
        err2(i) = err(i);
    end
end
ind = find(err1 == min(err(err1 < 0))); %max err

for i = 1:size(dat.out.traj.rho(:,ind),1)
    dens_err(i) = (log(dat.out.traj.rho(i,ind)) - log(dat.out.g.rho_est(i,ind)))/log(dat.out.traj.rho(i,ind));
end

tjr1 = dat.out.traj.t(dat.out.idj(ind,1),ind)/dat.out.traj.t(dat.out.idxend(ind),ind);

min_err = -500;
max_err = 1000;

figure(); hold on
grid on
set(gca,'FontSize',14)
yyaxis left
plot(dat.out.traj.t(:,ind)/dat.out.traj.t(dat.out.idxend(ind),ind), dat.out.g.ha_err(:,ind),'LineWidth',2)
ylabel('Estimated Apoapsis Error (km)')
ylim([min_err max_err])
yyaxis right
plot(dat.out.traj.t(:,ind)/dat.out.traj.t(dat.out.idxend(ind),ind), dens_err * 100,'LineWidth',2)
ylim([-2 3])
xlabel('Normalized Time, t/t_f')
ylabel('% Error of Log(\rho) (%)')
plot([tjr1 tjr1],[min_err max_err],'k--','LineWidth',2)
legend('Estimated Error (km)','Density Error','Jettison 1')

% dens RSS err vs. estimated err at jett2
dens_err1 = nan(1000,1);
dens_err2 = dens_err1;
est_err = dens_err1;
true_err = dens_err1;
dens_diff = dens_err1;
rss_K = dens_err1;
ind_tot = length(find(~isnan(dat.out.tjett(:,2))));
ind1 = 0; ind2 = 0; ind3 = 0; ind4 = 0;
for i = 1:1000
    if (~isnan(tj(i,2)))
        ind = dat.out.idj(i,2);
        if (dat.out.traj.t(dat.out.idxend(i),i) < tj(i,2))
            ind = dat.out.idxend(i);
        end
        dens_err1(i) = (log(dat.out.traj.rho(dat.out.idj(i,1),i)) - ... 
            log(dat.out.g.rho_est(dat.out.idj(i,1),i)))/log(dat.out.traj.rho(dat.out.idj(i,1),i));
        dens_err2(i) = (log(dat.out.traj.rho(ind,i)) - ... 
            log(dat.out.g.rho_est(ind,i)))/log(dat.out.traj.rho(ind,i));
        est_err(i) = dat.out.g.ha_err(ind,i);
        true_err(i) = err(i);
        rss_K(i) = dat.out.g.atm.rss_K(ind,i);
        
        if (rss_K(i) > 1 && est_err(i) < 0)
            ind1 = ind1 + 1; %bottom right
        elseif (rss_K(i) > 1 && est_err(i) > 0)
            ind2 = ind2 + 1; %top right
        elseif (rss_K(i) < 1 && est_err(i) > 0)
            ind3 = ind3 + 1; %top left
        elseif (rss_K(i) < 1 && est_err(i) < 0)
            ind4 = ind4 + 1; %bottom left
        end
            
        
        % density difference between tj2 and tj1
        dens_diff(i) = dens_err2(i) - dens_err1(i);
    end
end

ymin = -200;
ymax = 500;
% xmin = floor(min(rss_K));
% xmax = ceil(max(rss_K));
xmin = 0.5;
xmax = 1.5;
figure(); hold on
title(['betas: ' dat.betas ', ha*=' num2str(dat.gnc.ha_tgt)]);
grid on
set(gca,'FontSize',14)
plot(rss_K, est_err,'*b','LineWidth',1.5)
% plot(dens_err, true_err,'*r','LineWidth',1.5)
ylim([ymin ymax])
xlim([xmin xmax])
plot([1 1], [ymin ymax],'k--','LineWidth',1.5)
plot([xmin xmax], [1 1],'k--','LineWidth',1.5)
xlabel('Density RSS Error (kg/m^3)')
ylabel('\Delta h_a at Jettison 2 (km)')
% bottom right
dim = [0.7 0.5 0.6 0.3];
str = {[num2str(100 * ind1 / ind_tot) '%']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% top right
dim = [0.7 0.1 0.6 0.2];
str = {[num2str(100 * ind2 / ind_tot) '%']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% top left
dim = [0.2 0.5 0.3 0.3];
str = {[num2str(100 * ind3 / ind_tot) '%']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% bottom left
dim = [0.2 0.1 0.3 0.2];
str = {[num2str(100 * ind4 / ind_tot) '%']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 3-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_jsr();
br41 = 10;
br21s = [3 5 7 9];
br31s = [5 7 9 9.5];

% br21s = [5];
% br31s = [7];

fpas = [-5.4 -5.6];
ha_tgts = 2000;

% br21s = [9];
% br31s = [9.5];

% fpas = [-5.4 -5.6];
% 
for jj = 1:length(ha_tgts)
for kk = 1:length(fpas)
for k = 1:length(br31s)
%     
%     br31 = br31s(k);
% 
% % ratios = [4 3 3];
% % m1 = 200;

br21 = br21s(k);
br31 = br31s(k);

aero.m = [150 60 40 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * br21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * br31) );
aero.rcs(4) = sqrt( (aero.m(4) * aero.rcs(1)^2) / ...
    (aero.m(1) * br41) );

if (aero.rcs(4) > aero.rcs(3)) || (aero.rcs(3) > aero.rcs(2)) ... 
 || (aero.rcs(2) > aero.rcs(1))
    error('Impossible Geometry')
end
% 
aero.rn = 0.75/2;
aero.cds = [1.05 1.05 1.05 1.05];
aero.cls = [0 0 0 0];
% 
% gnc.ha_tgt = 10000;
gnc.ha_tgt = ha_tgts(jj);
gnc.n = 3;
gnc.tj0 = [90 170 230];
gnc.guid_rate = 0.5;

sim.t_max = 1500;
% 
for i = 1:4
    b(i) = aero.m(i)/( aero.cds(1) * pi * aero.rcs(i)^2 );
end
br41 = b(4)/b(1);
br31 = b(3)/b(1);
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br43 = b(4)/b(3);
br42 = b(4)/b(2);
% 
mc.flag = true;
mc.N = 1000;
mc.sigs = default_sigs();
mc.sigs.cd = (0.05 * aero.cds(1) / 3);
% mc.Kflag = uint8(0);

sim.traj_rate = 200;
sim.planet = 'venus';
sim.efpa_flag = false;
gnc.dtj_lim = 10;

gnc.iters = uint8(1);
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2); %multi stage
gnc.comp_curr = false;

% mc.debug = true;

sim.parMode = true;
sim.nWorkers = 10;

t = datetime();
fprintf('Starting 3stage Monte Carlo %i, %i/%i/%i, %2.0f:%2.0f:%2.0f\n', ... 
    k,t.Month, t.Day, t.Year, t.Hour, t.Minute, t.Second);
x0.fpa0 = fpas(kk);
fprintf(['Starting 3stage, betas = ' ... 
    num2str(br21) ',' num2str(br31) ',' num2str(br41)]); print_current_time();
out = run_dej_n(x0,gnc,aero,sim,mc);
fprintf(['Finished 3stage sim, betas = ' ... 
    num2str(br21) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

mc.N = 1000;
for i = 1:4
    b(i) = aero.m(i)/(1.05*pi*aero.rcs(i)^2);
end
br21 = b(2)/b(1);
br31 = b(3)/b(1);
br41 = b(4)/b(1);
betas = ['b41=' num2str(br41) '..b31=' ...
    num2str(br31) '..b21=' num2str(br21)];

for j = 1:1000
%     ttf0(j) = out.tjett(j,1)./out.traj.t(out.idxend(j),j);
    ttf1(j) = out.tjett(j,1)./out.traj.t(out.idxend(j),j);
    ttf2(j) = out.tjett(j,2)./out.traj.t(out.idxend(j),j);
    ttf3(j) = out.tjett(j,3)./out.traj.t(out.idxend(j),j);
end

haf_err2 = nan(mc.N,1);
haf_err1 = haf_err2;
haf_err3 = haf_err2;
n3 = 0; % num 3 stage jettisons
n2 = 0; % num 2 stage jettisons
n1 = 0; % num 1 stage jettisons

ha_est1 = haf_err2;
ha_est2 = haf_err2;
ha_est3 = haf_err2;

for i = 1:mc.N
    if ~isnan(out.tjett(i,3))
        haf_err3(i) = out.haf_err(i);
        n3 = n3 + 1;
    elseif ~isnan(out.tjett(i,2))
        haf_err2(i) = out.haf_err(i);
        n2 = n2 + 1;
    else
        haf_err1(i) = out.haf_err(i);
        n1 = n1 + 1;
    end
    
    ind = out.idj(i,1);
    if ~isnan(ind)
        ha_est1(i) = out.g.ha_err(ind,i); %-1 b/c GNC run before store_dat.m
    else
        ha_est1(i) = nan;
    end
    
    ind = out.idj(i,2);
    if ~(isnan(ind))
        ha_est2(i) = out.g.ha_err(ind,i);
    else
        ha_est2(i) = nan;
    end    
    
    ind = out.idj(i,3);
    if ~(isnan(ind))
        ha_est3(i) = out.g.ha_err(ind,i);
    else
        ha_est3(i) = nan;
    end     
    
end

% keyboard;

if (gnc.ha_tgt < 1000)
    save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt)) ... 
        '\rss_data\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '_FJ.mat'])
else
    save([save_path 'dej_n\HYDRA\2stage\data\' num2str(floor(gnc.ha_tgt/1000)) ... 
        'k\rss_data\v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_' betas '_FJ.mat'])
end



end %k betas
end %kk fpas
end %jj ha tgts
keyboard;
% 
% 


header = ['\beta_n/\beta_1=10,beta_ratios=' ... 
    num2str(out.nom.beta_ratios(1)) ' & ' ... 
    num2str(out.nom.beta_ratios(2)) ' & ' ...
    num2str(out.nom.beta_ratios(3))];

% need to redo? need all time indices to match
% figure(); hold on
% set(gca,'FontSize',14)
% title(['beta ratios ' num2str(out.nom.beta_ratios(1)) ' & ' ... 
%     num2str(out.nom.beta_ratios(2)) ', h_{a,tgt} = 2000 km, ' ... 
%     'efpa* = ' num2str(out.nom.traj.gamma_pp(1)*180/pi) ... 
%     ', v0 = ' num2str(x0.v0)],'FontSize',10)
% % histogram(ttf0,'NumBins',50,'FaceColor','r')
% histogram(ttf1,'NumBins',50,'FaceColor','b')
% histogram(ttf2,'NumBins',100,'FaceColor','g')
% histogram(ttf3,'NumBins',100,'FaceColor','c')
% legend('3 Stage, Jettison 1', ... 
%     '3 Stage, Jettison 2', ... 
%     '3 Stage, Jettison 3', ...
%     'FontSize',12)
% xlabel('t_{jett}/t_{f}')
% ylabel('Occurrence')
% % xlim([75 250])
% saveas(gcf, [save_path 'dej_n/3stage/mc' ...
%     num2str(gnc.ha_tgt/1000) 'k_tj_' betas '.fig']);



% figure(); hold on
% set(gca,'FontSize',14)
% title([header 'h_{a,tgt} = 2000 km, ' ... 
%     'efpa* = ' num2str(out.nom.traj.gamma_pp(1)*180/pi) ... 
%     ', v0 = ' num2str(x0.v0)],'FontSize',10)
% histogram(out.haf_err,'NumBins',1000,'FaceColor','r')
% % histogram(out2.haf_err,'NumBins',1000,'FaceColor','b')
% histogram(haf_err1-gnc.ha_tgt,'NumBins',round(n1),'FaceColor','b')
% histogram(haf_err2-gnc.ha_tgt,'NumBins',round(n2),'FaceColor','g')
% histogram(haf_err3-gnc.ha_tgt,'NumBins',round(n3),'FaceColor','r')
% legend(['3 Stage, Single Jettison Event: n = ' num2str(n1*100/mc.N) '%'], ... 
%     ['3 Stage, Two Jettison Events: n = ' num2str(n2*100/mc.N) '%'], ... 
%     ['3 Stage, 3 Jettison Events: n = ' num2str(n3*100/mc.N) '%'], ... 
%     'FontSize',10)
% xlabel('Apoapsis (km)')
% ylabel('Occurrence')
% xlim([-200 400])
% saveas(gcf, [save_path 'dej_n/3stage/mc' ...
%     num2str(gnc.ha_tgt/1000) 'k_haf_' betas '.fig']);

% xlim([-1000 1000])

ind = 1; ind2 = 1; ind3 = 1;
for i = 1:1000
    if ~(isnan(haf_err1(i)))
        ha1(ind) = haf_err1(i);
        ind = ind+1;
    end
    if ~(isnan(haf_err2(i)))
        ha2(ind2) = haf_err2(i);
        ind2 = ind2+1;
    end
    if ~isnan(haf_err3(i))
        ha3(ind3) = haf_err3(i);
        ind3 = ind3 + 1;
    end
end

figure(); hold on
set(gca,'FontSize',14)
title(['beta ratios ' betas ', h_{a,tgt} = 2000 km, ' ... 
    'efpa* = ' num2str(out.nom.traj.gamma_pp(1)*180/pi) ... 
    ', v0 = ' num2str(x0.v0)],'FontSize',10)
% histogram(out1.haf_err,'NumBins',1000,'FaceColor','r')
% histogram(out2.haf_err,'NumBins',1000,'FaceColor','b')
w = 20;
h3=histogram(ha3,'NumBins',50,'FaceColor','r');
h3.BinWidth = w; h3.FaceAlpha = 0.7;
h2=histogram(ha2,'NumBins',50,'FaceColor','g');
h2.BinWidth = w; h2.FaceAlpha = 0.7;
h1=histogram(ha1,'NumBins',50,'FaceColor','b');
h1.BinWidth = w; h1.FaceAlpha = 0.7;
legend([h1 h2 h3],['3 Stage, Single Jettison Event: n = ' num2str(n1*100/mc.N) '%'], ... 
    ['3 Stage, Two Jettison Events: n = ' num2str(n2*100/mc.N) '%'], ... 
    ['3 Stage, Three Jettison Events: n = ' num2str(n3*100/mc.N) '%'], ... 
    'FontSize',10)
xlabel('Apoapsis Error (km)')
ylabel('Occurrence')
xlim([-200 400])


figure(); hold on
set(gca,'FontSize',14)
title(['beta ratios ' num2str(out.nom.beta_ratios(1)) ' & ' ... 
    num2str(out.nom.beta_ratios(2)) ', h_{a,tgt} = ' num2str(gnc.ha_tgt) ' km, ' ... 
    'efpa* = ' num2str(out.nom.traj.gamma_pp(1)*180/pi) ... 
    ', v0 = ' num2str(x0.v0)],'FontSize',10)
plot(ttf1, haf_err1,'b*')
plot(ttf2, haf_err2,'g*')
plot(ttf3, haf_err3,'r*')
ylim([-500 500])
xlabel('t_{jett}/t_{f}')
ylabel('Apoapsis Altitude Error (km)')
legend(['Single Jettison, n = ' num2str(n1*100/mc.N) '%'], ... 
    ['Two Jettison Events: n = ' num2str(n2*100/mc.N) '%'], ... 
    ['Three Jettison Events: n = ' num2str(n3*100/mc.N) '%'], ... 
    'FontSize',10, 'location','south')


figure(); hold on
set(gca,'FontSize',14)
title(['beta ratios ' betas ', h_{a,tgt} = ' num2str(gnc.ha_tgt) ' km, ' ... 
    'efpa* = ' num2str(out.nom.traj.gamma_pp(1)*180/pi) ... 
    ', v0 = ' num2str(x0.v0)],'FontSize',10)
plot(ttf1, ha_est1,'b*');
plot(ttf2, ha_est2,'g*');
plot(ttf3, ha_est3, 'r*');
xlabel('t_{jett}/t_{f}')
ylabel('Estimated Apoapsis Error (km)')
legend('1^{st} Jettison Event', ... 
    '2^{nd} Jettison Event', ... 
    '3^{rd} Jettison Event', ... 
    'FontSize',10, 'location','north')
ylim([-500 500])
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% N-Stage Comparions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clear;clc;
p = 'C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\HYDRA\';
p2 = [p '2stage\data\'];
p3 = [p '3stage\data\'];
%%%%%%%%% efpa = -5.4, 2k
d1 = load([p '1stage\data\mcNewton_fpas_hafs_0.5hz.mat']);
d2_1 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=3_b32=3.3333.mat']);
d2_2 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=7_b32=1.4286.mat']);
d2_4 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=9_b32=1.1111.mat']);
d3_1 = load([p3 'mc_2k_fpa5.4_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_2k_fpa5.4_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_2k_fpa5.4_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_2k_fpa5.4_b41=10..b31=9.5..b21=9.mat']);

% d2_52 = load([p '2stage\data\mc2k_fpaBoth_5Hz, b31=10_b21=5_b32=2.mat']);
% d2_25 = load([p '2stage\data\mc2k_fpaBoth_5Hz, b31=10_b21=2_b32=5.mat']);
% d3_43 = load([p '3stage\data\mc_2k_fpa5.4_b41=10..b31=4..b21=3.mat']);
% d3_63 = load([p '3stage\data\mc_2k_fpa5.4_b41=10..b31=6..b21=3.mat']);
% d3_83 = load([p '3stage\data\mc_2k_fpa5.4_b41=10..b31=8..b21=3.mat']);

d1_haf_errs = d1.haf_n_2k(:,1) - 2000;

% names = {'d1','d2_25','d2_52','d3_43','d3_63','d3_83'};
tols = (0:10:500)';
Ni = length(tols);
pc1 = nan(Ni,1);
pc2_1 = pc1; pc2_2 = pc1; pc2_3 = pc1; pc2_4 = pc1;
pc3_1 = pc1; pc3_2 = pc1; pc3_3 = pc1; pc3_4 = pc1;

mcN = 1000;
for i = 1:Ni
    num = find(abs(d1_haf_errs) <= tols(i));
    pc1(i) = length(num)/mcN * 100;
    
    num = find(abs(d2_1.out.haf_err) <= tols(i));
    pc2_1(i) = length(num)/mcN * 100;
    num = find(abs(d2_2.out.haf_err) <= tols(i));
    pc2_2(i) = length(num)/mcN * 100;
    num = find(abs(d2_3.out.haf_err) <= tols(i));
    pc2_3(i) = length(num)/mcN * 100;
    num = find(abs(d2_4.out.haf_err) <= tols(i));
    pc2_4(i) = length(num)/mcN * 100;
    
    num = find(abs(d3_1.out.haf_err) <= tols(i));
    pc3_1(i) = length(num)/mcN * 100;
    num = find(abs(d3_2.out.haf_err) <= tols(i));
    pc3_2(i) = length(num)/mcN * 100;
    num = find(abs(d3_3.out.haf_err) <= tols(i));
    pc3_3(i) = length(num)/mcN * 100;
    num = find(abs(d3_4.out.haf_err) <= tols(i));
    pc3_4(i) = length(num)/mcN * 100;
end

figure(); hold on
set(gca,'FontSize',14)
title('2k tgt, -5.4deg, 0.5hz, venus')
plot(tols, pc1,'k--','LineWidth',1.5)
plot(tols, pc2_1,'LineWidth',1.5)
plot(tols, pc2_2,'LineWidth',1.5)
plot(tols, pc2_3,'LineWidth',1.5)
plot(tols, pc2_4,'LineWidth',1.5)
plot(tols, pc3_1,'-.','LineWidth',1.5)
plot(tols, pc3_2,'-.','LineWidth',1.5)
plot(tols, pc3_3,'-.','LineWidth',1.5)
plot(tols, pc3_4,'-.','LineWidth',1.5)
legend( ... 
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','se', 'FontSize',10);
xlabel('Apoapsis Error Tolerance (km)')
ylabel('Capture Rate (%)')
xlim([0 200])


%%%%%%%%% efpa = -5.6, 2k
clearvars -except p p2 p3 d1 pc1
d2_1 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=3_b32=3.33.mat']);
d2_2 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=7_b32=1.43.mat']);
d2_4 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=9_b32=1.11.mat']);
d3_1 = load([p3 'mc_2k_fpa5.6_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_2k_fpa5.6_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_2k_fpa5.6_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_2k_fpa5.6_b41=10..b31=9.5..b21=9.mat']);

d1 = load([p '1stage\data\mc_n_2k_0.5hz_56deg.mat']);

% d1_haf_errs = d1.haf_n_2k(:,2) - 2000;
d1_haf_errs = d1.out.haf_err;

tols = (0:10:500)';
Ni = length(tols);
pc1 = nan(Ni,1);
pc2_1 = pc1; pc2_2 = pc1; pc2_3 = pc1; pc2_4 = pc1;
pc3_1 = pc1; pc3_2 = pc1; pc3_3 = pc1; pc3_4 = pc1;
mcN = 1000;
for i = 1:Ni
    num = find(abs(d1_haf_errs) <= tols(i));
    pc1(i) = length(num)/mcN * 100;
    
    num = find(abs(d2_1.out.haf_err) <= tols(i));
    pc2_1(i) = length(num)/mcN * 100;
    num = find(abs(d2_2.out.haf_err) <= tols(i));
    pc2_2(i) = length(num)/mcN * 100;
    num = find(abs(d2_3.out.haf_err) <= tols(i));
    pc2_3(i) = length(num)/mcN * 100;
    num = find(abs(d2_4.out.haf_err) <= tols(i));
    pc2_4(i) = length(num)/mcN * 100;
    
    num = find(abs(d3_1.out.haf_err) <= tols(i));
    pc3_1(i) = length(num)/mcN * 100;
    num = find(abs(d3_2.out.haf_err) <= tols(i));
    pc3_2(i) = length(num)/mcN * 100;
    num = find(abs(d3_3.out.haf_err) <= tols(i));
    pc3_3(i) = length(num)/mcN * 100;
    num = find(abs(d3_4.out.haf_err) <= tols(i));
    pc3_4(i) = length(num)/mcN * 100;
end

figure(); hold on
set(gca,'FontSize',14)
title('2k tgt, -5.6deg, 0.5hz, venus')
plot(tols, pc1,'k--','LineWidth',1.5)
plot(tols, pc2_1,'LineWidth',1.5)
plot(tols, pc2_2,'LineWidth',1.5)
plot(tols, pc2_3,'LineWidth',1.5)
plot(tols, pc2_4,'LineWidth',1.5)
plot(tols, pc3_1,'-.','LineWidth',1.5)
plot(tols, pc3_2,'-.','LineWidth',1.5)
plot(tols, pc3_3,'-.','LineWidth',1.5)
plot(tols, pc3_4,'-.','LineWidth',1.5)
legend( ... 
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','se', 'FontSize',10);
xlabel('Apoapsis Error Tolerance (km)')
ylabel('Capture Rate (%)')
% xlim([0 200])


%%%%%%%%% efpa = -5.4, 10k
clearvars -except p p2 p3 d1 pc1
d2_1 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=3_b32=3.33.mat']);
d2_2 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=7_b32=1.43.mat']);
d2_4 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=9_b32=1.11.mat']);
d3_1 = load([p3 'mc_10k_fpa5.4_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_10k_fpa5.4_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_10k_fpa5.4_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_10k_fpa5.4_b41=10..b31=9.5..b21=9.mat']);

% d1_haf_errs = d1.haf_n_10k(:,1) - 10000;

d1 = load([p '1stage\data\mc_n_10k_0.5hz_54deg.mat']);
d1_haf_errs = d1.out.haf_err;

tols = (0:10:500)';
Ni = length(tols);
pc1 = nan(Ni,1);
pc2_1 = pc1; pc2_2 = pc1; pc2_3 = pc1; pc2_4 = pc1;
pc3_1 = pc1; pc3_2 = pc1; pc3_3 = pc1; pc3_4 = pc1;
mcN = 1000;
for i = 1:Ni
    num = find(abs(d1_haf_errs) <= tols(i));
    pc1(i) = length(num)/mcN * 100;
    
    num = find(abs(d2_1.out.haf_err) <= tols(i));
    pc2_1(i) = length(num)/mcN * 100;
    num = find(abs(d2_2.out.haf_err) <= tols(i));
    pc2_2(i) = length(num)/mcN * 100;
    num = find(abs(d2_3.out.haf_err) <= tols(i));
    pc2_3(i) = length(num)/mcN * 100;
    num = find(abs(d2_4.out.haf_err) <= tols(i));
    pc2_4(i) = length(num)/mcN * 100;
    
    num = find(abs(d3_1.out.haf_err) <= tols(i));
    pc3_1(i) = length(num)/mcN * 100;
    num = find(abs(d3_2.out.haf_err) <= tols(i));
    pc3_2(i) = length(num)/mcN * 100;
    num = find(abs(d3_3.out.haf_err) <= tols(i));
    pc3_3(i) = length(num)/mcN * 100;
    num = find(abs(d3_4.out.haf_err) <= tols(i));
    pc3_4(i) = length(num)/mcN * 100;
end

figure(); hold on
set(gca,'FontSize',14)
title('10k tgt, -5.4deg, 0.5hz, venus')
plot(tols, pc1,'k--','LineWidth',1.5)
plot(tols, pc2_1,'LineWidth',1.5)
plot(tols, pc2_2,'LineWidth',1.5)
plot(tols, pc2_3,'LineWidth',1.5)
plot(tols, pc2_4,'LineWidth',1.5)
plot(tols, pc3_1,'-.','LineWidth',1.5)
plot(tols, pc3_2,'-.','LineWidth',1.5)
plot(tols, pc3_3,'-.','LineWidth',1.5)
plot(tols, pc3_4,'-.','LineWidth',1.5)
legend( ... 
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','se', 'FontSize',10);
xlabel('Apoapsis Error Tolerance (km)')
ylabel('Capture Rate (%)')
xlim([0 200])

%%%%%%%%% efpa = -5.6, 10k
clearvars -except p p2 p3 d1 pc1
d2_1 = load([p2 'mc10k_fpa5.6_0.5Hz__b31=10_b21=3_b32=3.33.mat']);
d2_2 = load([p2 'mc10k_fpa5.6_0.5Hz__b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc10k_fpa5.6_0.5Hz__b31=10_b21=7_b32=1.43.mat']);
d2_4 = load([p2 'mc10k_fpa5.6_0.5Hz__b31=10_b21=9_b32=1.11.mat']);
d3_1 = load([p3 'mc_10k_fpa5.6_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_10k_fpa5.6_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_10k_fpa5.6_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_10k_fpa5.6_b41=10..b31=9.5..b21=9.mat']);

d1_haf_errs = d1.haf_n_10k(:,2) - 10000;

d1 = load([p '1stage\data\mc_n_10k_0.5hz_56deg.mat']);
d1_haf_errs = d1.out.haf_err;

tols = (0:10:500)';
Ni = length(tols);
pc2_1 = pc1; pc2_2 = pc1; pc2_3 = pc1; pc2_4 = pc1;
pc3_1 = pc1; pc3_2 = pc1; pc3_3 = pc1; pc3_4 = pc1;
mcN = 1000;
for i = 1:Ni
    num = find(abs(d1_haf_errs) <= tols(i));
    pc1(i) = length(num)/mcN * 100;
    
    num = find(abs(d2_1.out.haf_err) <= tols(i));
    pc2_1(i) = length(num)/mcN * 100;
    num = find(abs(d2_2.out.haf_err) <= tols(i));
    pc2_2(i) = length(num)/mcN * 100;
    num = find(abs(d2_3.out.haf_err) <= tols(i));
    pc2_3(i) = length(num)/mcN * 100;
    num = find(abs(d2_4.out.haf_err) <= tols(i));
    pc2_4(i) = length(num)/mcN * 100;
    
    num = find(abs(d3_1.out.haf_err) <= tols(i));
    pc3_1(i) = length(num)/mcN * 100;
    num = find(abs(d3_2.out.haf_err) <= tols(i));
    pc3_2(i) = length(num)/mcN * 100;
    num = find(abs(d3_3.out.haf_err) <= tols(i));
    pc3_3(i) = length(num)/mcN * 100;
    num = find(abs(d3_4.out.haf_err) <= tols(i));
    pc3_4(i) = length(num)/mcN * 100;
end

figure(); hold on
set(gca,'FontSize',14)
title('10k tgt, -5.6deg, 0.5hz, venus')
plot(tols, pc1,'k--','LineWidth',1.5)
plot(tols, pc2_1,'LineWidth',1.5)
plot(tols, pc2_2,'LineWidth',1.5)
plot(tols, pc2_3,'LineWidth',1.5)
plot(tols, pc2_4,'LineWidth',1.5)
plot(tols, pc3_1,'-.','LineWidth',1.5)
plot(tols, pc3_2,'-.','LineWidth',1.5)
plot(tols, pc3_3,'-.','LineWidth',1.5)
plot(tols, pc3_4,'-.','LineWidth',1.5)
legend( ... 
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','nw', 'FontSize',10);
xlabel('Apoapsis Error Tolerance (km)')
ylabel('Capture Rate (%)')
% xlim([0 200])

%%%%%%%%%%%%
% histograms
clearvars -except p p2 p3

% 2stage, 2k, -5.4deg
d1 = load([p '1stage\data\mcNewton_fpas_hafs_0.5hz.mat']);
d2_1 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=3_b32=3.3333.mat']);
d2_2 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=7_b32=1.4286.mat']);
d2_4 = load([p2 'mc2k_fpa5_4_0.5Hz, b31=10_b21=9_b32=1.1111.mat']);
d3_1 = load([p3 'mc_2k_fpa5.4_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_2k_fpa5.4_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_2k_fpa5.4_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_2k_fpa5.4_b41=10..b31=9.5..b21=9.mat']);

figure(); hold on
set(gca,'FontSize',14)
title('2stage, 2k tgt, -5.4deg, 0.5hz, venus')
w = 15;
alpha = 0.6;
h0 = histogram(d1.haf_n_2k,'FaceColor','k');
h0.BinWidth = w; h0.FaceAlpha = alpha;
h1 = histogram(d2_1.out.haf_err + 2000,'FaceColor','r');
h1.BinWidth = w; h1.FaceAlpha = alpha;
h2 = histogram(d2_2.out.haf_err + 2000,'FaceColor','b');
h2.BinWidth = w; h2.FaceAlpha = alpha;
h3 = histogram(d2_3.out.haf_err + 2000,'FaceColor','g');
h3.BinWidth = w; h3.FaceAlpha = alpha;
h4 = histogram(d2_4.out.haf_err + 2000,'FaceColor','m');
h4.BinWidth = w; h4.FaceAlpha = alpha;

xlabel('Apoapsis Altitude (km)')
ylabel('Occurrence')
legend([h0 h1 h2 h3 h4], ...
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    'location','ne','FontSize',10)
xlim([-400 600] + 2000)
% 3stage, 2k, -5.4deg
figure(); hold on
set(gca,'FontSize',14)
title('2stage, 2k tgt, -5.4deg, 0.5hz, venus')
w = 15;
alpha = 0.6;
h0 = histogram(d1.haf_n_2k,'FaceColor','k');
h0.BinWidth = w; h0.FaceAlpha = alpha;
h1 = histogram(d3_1.out.haf_err + 2000,'FaceColor','r');
h1.BinWidth = w; h1.FaceAlpha = alpha;
h2 = histogram(d3_2.out.haf_err + 2000,'FaceColor','b');
h2.BinWidth = w; h2.FaceAlpha = alpha;
h3 = histogram(d3_3.out.haf_err + 2000,'FaceColor','g');
h3.BinWidth = w; h3.FaceAlpha = alpha;
h4 = histogram(d3_4.out.haf_err + 2000,'FaceColor','m');
h4.BinWidth = w; h4.FaceAlpha = alpha;

xlabel('Apoapsis Altitude (km)')
ylabel('Occurrence')
legend([h0 h1 h2 h3 h4], ...
    '1 Stage', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','ne','FontSize',10)
xlim([-400 600] + 2000)

%%%% 2stage, 10k, -5.4deg
clearvars -except d1 p p2 p3
d2_1 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=3_b32=3.33.mat']);
d2_2 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=7_b32=1.43.mat']);
d2_4 = load([p2 'mc10k_fpa5.4_0.5Hz__b31=10_b21=9_b32=1.11.mat']);
d3_1 = load([p3 'mc_10k_fpa5.4_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_10k_fpa5.4_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_10k_fpa5.4_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_10k_fpa5.4_b41=10..b31=9.5..b21=9.mat']);

figure(); hold on
set(gca,'FontSize',14)
title('2stage, 10k tgt, -5.4deg, 0.5hz, venus')
w = 0.15;
alpha = 0.6;
h0 = histogram(d1.haf_n_10k/1000,'FaceColor','k');
h0.BinWidth = w; h0.FaceAlpha = alpha;
h1 = histogram((d2_1.out.haf_err + 10000)/1000,'FaceColor','r');
h1.BinWidth = w; h1.FaceAlpha = alpha;
h2 = histogram((d2_2.out.haf_err + 10000)/1000,'FaceColor','b');
h2.BinWidth = w; h2.FaceAlpha = alpha;
h3 = histogram((d2_3.out.haf_err + 10000)/1000,'FaceColor','g');
h3.BinWidth = w; h3.FaceAlpha = alpha;
h4 = histogram((d2_4.out.haf_err + 10000)/1000,'FaceColor','m');
h4.BinWidth = w; h4.FaceAlpha = alpha;
xlabel('Apoapsis Altitude (10^3 km)')
ylabel('Occurrence')
legend([h0 h1 h2 h3 h4], ...
    '1 Stage', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{3,3.33\}', ... 
    '2 Stage, \beta_{i+1}/\beta_i = \{5,2\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{7,1.42\}', ...
    '2 Stage, \beta_{i+1}/\beta_i = \{9,1.11\}', ...
    'location','ne','FontSize',10)
xlim([-1 4]+10)

%3 stage, 10k, -5.4deg
figure(); hold on
set(gca,'FontSize',14)
title('3stage, 10k tgt, -5.4deg, 0.5hz, venus')
w = 0.12;
alpha = 0.6;
h0 = histogram(d1.haf_n_10k/1000,'FaceColor','k');
h0.BinWidth = w; h0.FaceAlpha = alpha;
h1 = histogram((d3_1.out.haf_err + 10000)/1000,'FaceColor','r');
h1.BinWidth = w; h1.FaceAlpha = alpha;
h2 = histogram((d3_2.out.haf_err + 10000)/1000,'FaceColor','b');
h2.BinWidth = w; h2.FaceAlpha = alpha;
h3 = histogram((d3_3.out.haf_err + 10000)/1000,'FaceColor','g');
h3.BinWidth = w; h3.FaceAlpha = alpha;
h4 = histogram((d3_4.out.haf_err + 10000)/1000,'FaceColor','m');
h4.BinWidth = w; h4.FaceAlpha = alpha;
xlabel('Apoapsis Altitude (10^3 km)')
ylabel('Occurrence')
legend([h0 h1 h2 h3 h4], ...
    '1 Stage', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{3,1.67,2\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{5,1.4,1.43\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{7,1.28,1.11\}', ...
    '3 Stage, \beta_{i+1}/\beta_i = \{9,1.05,1.05\}', ...
    'location','ne','FontSize',10)
xlim([-1 4]+10)

%%%% 2stage, 2k, -5.6deg
clearvars -except d1 p p2 p3
d2_1 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=3_b32=3.33.mat']);
d2_2 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=5_b32=2.mat']);
d2_3 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=7_b32=1.43.mat']);
d2_4 = load([p2 'mc2k_fpa5.6_0.5Hz__b31=10_b21=9_b32=1.11.mat']);
d3_1 = load([p3 'mc_2k_fpa5.6_b41=10..b31=5..b21=3.mat']);
d3_2 = load([p3 'mc_2k_fpa5.6_b41=10..b31=7..b21=5.mat']);
d3_3 = load([p3 'mc_2k_fpa5.6_b41=10..b31=9..b21=7.mat']);
d3_4 = load([p3 'mc_2k_fpa5.6_b41=10..b31=9.5..b21=9.mat']);

figure(); hold on
set(gca,'FontSize',14)
title('2stage, 2k tgt, -5.4deg, 0.5hz, venus')
w = 15;
alpha = 0.6;
h0 = histogram(d1.haf_n_2k,'FaceColor','k');
h0.BinWidth = w; h0.FaceAlpha = alpha;
h1 = histogram(d2_1.out.haf_err + 2000,'FaceColor','r');
h1.BinWidth = w; h1.FaceAlpha = alpha;
h2 = histogram(d2_2.out.haf_err + 2000,'FaceColor','b');
h2.BinWidth = w; h2.FaceAlpha = alpha;
h3 = histogram(d2_3.out.haf_err + 2000,'FaceColor','g');
h3.BinWidth = w; h3.FaceAlpha = alpha;
h4 = histogram(d2_4.out.haf_err + 2000,'FaceColor','m');
h4.BinWidth = w; h4.FaceAlpha = alpha;

%%%% 2stage, 10k, -5.6deg
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONTINUOUS DRAG MODULATION (CVDMA) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cvdma monte carlos - Mar/April 2019
%{
clear;clc; close all
ha_tgt = 2000;  %km
ha_tol = 10;    %km   
drc = 0.5;
area_rate = pi*drc^2;   %max/min area rate change
rc_bounds = [0.5 0.75];   %min/max backshell radii
rc_bounds_2 = [1.4 1.5]./2;
efpa = -5.1;
v_atm = 10.65;  %km/s
m = 68.2;
% cd0 = 1.05;
rc0 = 0.5;
delta_c0 = 50;
rn = rc0/2;
chord = rc0/(tan(delta_c0*pi/180));

betaDat = getBetas(rn,rc_bounds,delta_c0,rc0,m, chord);
cd0 = betaDat.cd0;

betaDat2 = getBetas(rn,rc_bounds_2,delta_c0,rc0,m, chord);

[x0,aero,gnc,sim,mc] = setup_cvdma();
x0.v0 = v_atm;
x0.efpa = efpa;
sim.efpa_flag = uint8(1); % calc fpa from ha_tgt   

aero.m = m;
aero.cd0 = betaDat.cd0;
aero.rc0 = rc0;
aero.rn = rn;
aero.delta_c0 = delta_c0;

gnc.guid_rate = 2;
gnc.rc_bounds = rc_bounds;
gnc.ha_tgt = ha_tgt;
gnc.area_rate = area_rate;

mc.flag = uint8(1);
mc.sigs = default_sigs();
mc.sigs.cd = (0.05 * cd0 / 3);  %5% 3sigma
mc.N = 1000;

% sim.traj_rate = 10;

t = datetime();
fprintf('Starting First Monte Carlo %i/%i/%i, %2.0f:%2.0f:%2.0f\n', ... 
    t.Month, t.Day, t.Year, t.Hour, t.Minute, t.Second);
% control authority betaDat
% out1 = cvdma_venus(x0,aero,gnc,sim,mc);

out1 = run_cvdma(x0,aero,gnc,sim,mc);

t = datetime();
fprintf('Starting Second Monte Carlo %i/%i/%i, %2.0f:%2.0f:%2.0f\n', ... 
    t.Month, t.Day, t.Year, t.Hour, t.Minute, t.Second);
% control authority betaDat2
gnc.rc_bounds = rc_bounds_2;
out2 = cvdma_venus(x0,aero,gnc,sim,mc);

tols = [25 50 100 250 500 1000];
pc1 = nan(length(tols),1);
pc2 = pc1;
entries1 = cell(length(tols),1);
entries2 = entries1;
for i = 1:length(tols)
    num1 = find(abs(out1.haf_err) <= tols(i));
    num2 = find(abs(out2.haf_err) <= tols(i));
    
    pc1(i) = length(num1)*100/mc.N;
    pc2(i) = length(num2)*100/mc.N;
    
    entries1{i} = [num2str(pc1(i)) '% within ' num2str(tols(i)) ' km'];
    entries2{i} = [num2str(pc2(i)) '% within ' num2str(tols(i)) ' km'];
end





% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% yyaxis left
% plot(out1.traj.t,out1.g.cvdma.rc,'LineWidth',2)
% ylabel('Backshell Radius (m)')
% yyaxis right
% plot(out1.traj.t,out1.g.ha_err./1000,'LineWidth',2)
% ylabel('Apoapsis Error (km)')
% xlabel('Time (sec)')


% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% title(['Venus CVDMA, ha_{tgt} = ' num2str(ha_tgt) ' km'],'FontSize',10)
% yyaxis left
% plot(out.traj.time,out.g.cvdma.delta_ra./1000,'LineWidth',2)
% xlabel('Time (sec)')
% ylabel('Apoapsis Error (km)')
% yyaxis right
% % figure(); hold on
% % grid on
% % set(gca,'FontSize',14)
% plot(out.traj.time,out.g.cvdma.rc,'LineWidth',2)
% xlabel('Time (sec)')
% ylabel('Backshell Radius (m)')




% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% title(['CVDMA, EFPA = -5.1^{o},h_{a,tgt} = ' ... 
%     num2str(ha_tgt) ' km, GNC Rate= ' num2str(guid_rate) ' Hz'],'FontSize',10)
% plot(out.traj.vmag(:,10),out.traj.alt(:,10),'LineWidth',1.5)
% plot(out.traj.vmag(:,76),out.traj.alt(:,76),'LineWidth',1.5)
% plot(out.traj.vmag(:,254),out.traj.alt(:,254),'LineWidth',1.5)
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')

header = {['CVDMA, EFPA* = ' num2str(out1.nom.traj.gamma_pp(1)*180/pi) '^{o}, h_{a,tgt} = ' ... 
    num2str(gnc.ha_tgt) ' km, GNC Rate= ' num2str(gnc.guid_rate) ' Hz, rc_{bounds} = [' ... 
    num2str(rc_bounds(1)) ' ' num2str(rc_bounds(2)) '] m'], ... 
    ['\beta_{max}/\beta_{min} = ' num2str(betaDat.beta_ratio)]};

figure(); hold on
grid on
set(gca,'FontSize',14)
title(header,'FontSize',10)
histogram(out1.haf_err,'NumBins',50,'FaceColor','r')
xlabel('Apoapsis Error (km)')
ylabel('Occurrences')
text(-500,100,entries1)



header2 = {['CVDMA, EFPA* = ' num2str(out2.nom.traj.gamma_pp(1)*180/pi) '^{o}, h_{a,tgt} = ' ... 
    num2str(gnc.ha_tgt) ' km, GNC Rate= ' num2str(gnc.guid_rate) ' Hz, rc_{bounds} = [' ... 
    num2str(rc_bounds_2(1)) ' ' num2str(rc_bounds_2(2)) '] m'], ... 
    ['\beta_{max}/\beta_{min} = ' num2str(betaDat2.beta_ratio)]};

figure(); hold on
grid on
set(gca,'FontSize',14)
title(header2,'FontSize',10)
histogram(out2.haf_err,'NumBins',50,'FaceColor','r')
xlabel('Apoapsis Error (km)')
ylabel('Occurrences')
text(-400,70,entries2)



%}





%% Param sweep of rc_bounds (min -> 0.75)
%{
clear; clc
rc1 = linspace(0.05,0.7,10);
Ni = length(rc1);
rc_max = 0.75;

rc0 = 0.5;

ha_tgt = 2000;  %km
ha_tol = 10;    %km   

drc = 0.5;
area_rate = pi*drc^2;   %max/min area rate change
efpa = -5.1;
v_atm = 11;  %km/s
m = 80;
% cd0 = 1.05;
delta_c0 = 75;
rn = rc_max/2;

chord = rc0/tan(delta_c0*pi/180);

[x0,aero,gnc,sim,mc] = setup_cvdma();
x0.v0 = v_atm;
x0.efpa = efpa;
sim.efpa_flag = uint8(1); % calc fpa from ha_tgt   

aero.m = m;
aero.rc0 = rc0;
aero.rn = rn;
aero.delta_c0 = delta_c0;

gnc.guid_rate = 2;
gnc.ha_tgt = ha_tgt;
gnc.area_rate = area_rate;

mc.flag = uint8(1);
mc.sigs = default_sigs();

mc.N = 1000;
sim.traj_rate = 250;
sim.data_rate = 1;

tols = [10 50 100 250 500 1000];
Nj = length(tols);
pc = nan(Ni,Nj);

haf_err = nan(mc.N,Ni);
beta_ratios = nan(Ni,1);
for i = 1:Ni
    fprintf('Loop %i\n',i);
    rc_bounds = [rc1(i) rc_max];
    
    betaDat = getBetas(rn,rc_bounds,delta_c0,rc0,m, chord);

    aero.cd0 = betaDat.cd0;
    mc.sigs.cd = (0.05 * betaDat.cd0 / 3);  %5% 3sigma
    gnc.rc_bounds = rc_bounds;

    out = run_cvdma(x0,aero,gnc,sim,mc);
    
    haf_err(:,i) = out.haf_err;
    beta_ratios(i) = betaDat.beta_ratio;
    
    for j = 1:Nj
        num = find(abs(out.haf_err) <= tols(j));
        pc(i,j) = length(num)*100/mc.N;
    end
    
end %Ni

fpa0 = out.nom.traj.gamma_pp(1)*180/pi;
save('../venus_aerocapture/cvdma/catchRate.mat','beta_ratios','haf_err','fpa0')


keyboard;

for i = 1:Ni
    entry{i} = ['\beta_{max}/\beta_{min} = ' num2str(beta_ratios(i))];
end
entry = flip(entry);

figure(); hold on
title('CVDMA Venus - Capture Rate for Beta Ratio')
set(gca,'FontSize',14)
plot(tols,pc(2:7,:)','LineWidth',2)
legend(entry(2:7),'location','nw','FontSize',10)
ylabel('Percent Captured')
xlabel('Tolerance (km)')

%}


%% base trajectory
% fpa = nan;
% t_jett = 200000;
% traj_rate = 50;
% ha_tgt = 2000;
% dat = ac_venus_manual(fpa,t_jett,traj_rate,ha_tgt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dat] = getBetas(rn,rcs,qc0,rc0,m,chord)
% Inputs:
%   rn: nose radius (assumed constant), m
%   rcs: backshell radii bounds, m
%   qc0: initial spherecone half angle, deg
%   rc0: initial backshell radius, m
%   m: mass (assumed constant), kg

qc0 = qc0 * pi/180;

% function handles
% c0 = rc0/(tan(qc0));    % initial chord length, m
delta_c = @(rc) atan(rc/chord);
cd = @(rn,rc,qc) (1-(sin(qc)^4))*(rn/rc)^2 + ...
    (2*(sin(qc)^2))*(1-(rn/rc)^2 * (cos(qc)^2));

% fix rcs if needed
if rcs(2) < rcs(1)
    flip(rcs);
end

% minimum values
qc_min = delta_c(rcs(1));
dat.cd_min = cd(rn,rcs(1),qc_min);
dat.beta_max = m/(dat.cd_min*pi*rcs(1)^2);

% maximum values
qc_max = delta_c(rcs(2));
dat.cd_max = cd(rn,rcs(2),qc_max);
dat.beta_min = m/(dat.cd_max*pi*rcs(2)^2);

dat.beta_ratio = dat.beta_max/dat.beta_min;

% initial vals
dat.cd0 = cd(rn,rc0,qc0);
dat.beta0 = m/(dat.cd0*pi*rc0^2);

end



% base inputs for cvmda
function [x0,aero,gnc,sim,mc] = setup_cvdma()
x0.v0 = 12;
x0.efpa = -5;
x0.h0 = 150;
x0.lat0 = 0;
x0.lon0 = 0;
x0.az0 = 90;

aero.m = 80;
aero.cd0 = 1.05;
aero.rc0 = 0.75;
aero.rn = aero.rc0/2;
aero.delta_c0 = 45;

gnc.guid_rate = 0.5;
gnc.guid_mode = uint8(1);
gnc.rc_bounds = [0.5 1.5];
gnc.ha_tgt = 2000;
gnc.ha_tol = 10;
gnc.area_rate = pi*0.15^2;
gnc.iter_max = uint8(1);
gnc.cd_bounds = [0.25 2];

sim.traj_rate = 250;
sim.data_rate = 1;
sim.planet = 'venus';
sim.atm_mode = 3;
sim.t0 = 0;
sim.t_max = 1000;
sim.alt_max = 150;
sim.alt_min = 0;
sim.efpa_flag = uint8(0);

mc.flag = uint8(1);
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
mc.N = 1000;

end %



