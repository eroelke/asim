clear;clc
% testing
% fpa_corr = [-5.2 -5.6];    %deg
% 
% % Ni = 50;
% fpas = linspace(fpa_corr(1),fpa_corr(2),Ni);  % fpas to investigate


%% Time Difference
temp = load('.\venus_ac\ac_guid_algorithms\NPC\nominal_errs.mat');
fpas = temp.fpas;
Ni = length(fpas);
dt = [10,15,20];
Nj = length(dt);
ha_tgt = 2000;
haf = nan(Ni,Nj);
err = haf;
tj = haf;
for j = 1:Nj
    for i = 1:Ni
       out = ac_venus_dcf( fpas(i), 100, ha_tgt, dt(j) );

       haf(i,j) = out.haf;
       err(i,j) = haf(i,j) - ha_tgt;
       tj(i,j) = out.tjett;
    end
end

% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% % yyaxis left
% plot(fpas,err,'-o','LineWidth',2)
% xlabel('FPA_{atm} (deg)')
% ylabel('Apoapsis Error (km)')
% title('Apoapsis Target: 2000 km')
% legend('\Delta t = 10s','\Delta t = 15s','\Delta t = 20s', ... 
%     'location','ne')
% ylim([-1000 1200])
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% title('Nominal Jettison Time Comparison')
% % yyaxis right
% plot(fpas,tj-temp.t_2k,'-o','LineWidth',2);
% % plot(fpas,temp.t_2k,'--k','LineWidth',2);
% % legend('DCF Algorithm','NPC Algorithm')
% legend('\Delta t = 10s','\Delta t = 15s','\Delta t = 20s', ... 
%     'location','sw')
% xlabel('EFPA (deg)')
% ylabel('\Delta t_{jett} (s)')



% Apoapsis Target Difference
% temp = load('.\venus_ac\ac_guid_algorithms\NPC\nominal_errs.mat');
clear;clc;
temp = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\DCFJ\data\opt_t_jett_beta10.mat');
fpas = temp.fpas;
Ni = length(fpas);
ha_tgts = [2 10].*1e3;
Nj = length(ha_tgts);
dt = 10;
mc.flag = false;
haf = nan(Ni,Nj);
err = haf;
tj = haf;

rcs = [1 0.2];
m = [150 60];

for j = 1:Nj
    for i = 1:Ni
        out = ac_venus_dcf(fpas(i),500,m,rcs,ha_tgts(j),dt,mc);
        
        haf(i,j) = out.haf;
        err(i,j) = haf(i,j) - ha_tgts(j);
        tj(i,j) = out.tjett;
        
    end
end

figure(); hold on
grid on
set(gca,'FontSize',14)
% yyaxis left
plot(fpas(4:end),err(4:end,:),'-o','LineWidth',2)
xlabel('FPA_{atm} (deg)')
ylabel('Apoapsis Error (km)')
% title('Apoapsis Target: 2000 km')
legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km', ... 
    'location','se')
ylim([-1000 1200])

figure(); hold on
grid on
set(gca,'FontSize',14)
title('Nominal Jettison Time Comparison, \Delta t = 10s')
% yyaxis right
plot(fpas,tj(:,1)-temp.t_2k,'-o','LineWidth',2);
plot(fpas,tj(:,2)-temp.t_10k,'-o','LineWidth',2);
% plot(fpas,temp.t_2k,'--k','LineWidth',2);
% legend('DCF Algorithm','NPC Algorithm')
legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km', ... 
    'location','se')
xlabel('EFPA (deg)')
ylabel('\Delta t_{jett} (s)')
ylim([-3 3])


%% monte carlo for two different fpas
clear;clc;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
% % [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
mc.N = 1000;

fpas = [-5.3 -5.5];
% fpas = -5.4;
Ni = length(fpas);
ha_tgts = 2000;
% Ni = length(ha_tgts);

mc.flag = 1;
% flags = [0 1];
% Ni = length(flags);
% tols = [10 50 100 250 500];
% Nj = length(tols);
% percent = nan(Ni,Nj);
traj_rate = 100;
% 
% mc1 = mc;
% mc1.flag = 0;
% mc2 = mc;
% mc2.flag = 1;

for i = 1:Ni
%     fprintf('Monte Carlo %i\n',i)
    out(i) = ac_venus_dcf(fpas(i),traj_rate,ha_tgts,10,mc);
%     out2 = ac_venus_dcf(fpas,traj_rate,ha_tgts,10,mc2);
    % calculate capture percentage based on apoapsis tolerance
%     toc
%     for j = 1:Nj
%         num = find(abs(out(i).haf_err) <= tols(j));
%         percent(i,j) = length(num)*100/mc.N;
%     end
end
% figure(); hold on
% plot(out.nom.traj.vel_pp_mag./1000,out.nom.traj.alt./1000,'--k')
% for i = 1:mc.N
% plot(out.traj.vmag(:,i),out.traj.alt(:,i))
% end

figure(); hold on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, h_{a,tgt} = 2000 km')
histogram(out(1).haf_err,'NumBins',50,'FaceColor','r')
histogram(out(2).haf_err,'NumBins',50,'FaceColor','b')
% histogram(out2.haf_err,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('EFPA = -5.3^{o}','EFPA = -5.5^{o}','location','nw')

saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\DCF\monte carlo\dcf_mc_fpa.fig')

%% monte carlo for two different ha tgts
clear;clc;
close all;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';
% % [m,cl,cd,aoa,vmag_atm,hmag_atm,fpa,lat,lon,az]
mc.N = 1000;

% fpas = [-5.3 -5.5];
fpas = -5.4;
% Ni = length(fpas);
ha_tgts = [2000 10000];
Ni = length(ha_tgts);

mc.flag = 1;
% flags = [0 1];
% Ni = length(flags);
% tols = [10 50 100 250 500];
% Nj = length(tols);
% percent = nan(Ni,Nj);
traj_rate = 100;
% 
% mc1 = mc;
% mc1.flag = 0;
% mc2 = mc;
% mc2.flag = 1;

for i = 1:Ni
%     fprintf('Monte Carlo %i\n',i)
    out(i) = ac_venus_dcf(fpas,traj_rate,ha_tgts(i),10,mc);
%     out2 = ac_venus_dcf(fpas,traj_rate,ha_tgts,10,mc2);
    % calculate capture percentage based on apoapsis tolerance
%     toc
%     for j = 1:Nj
%         num = find(abs(out(i).haf_err) <= tols(j));
%         percent(i,j) = length(num)*100/mc.N;
%     end
end
% figure(); hold on
% plot(out(1).nom.traj.vel_pp_mag./1000,out(1).nom.traj.alt./1000,'--k')
% for i = 1:mc.N
% plot(out(1).traj.vmag(:,i),out(1).traj.alt(:,i))
% end

for i = 1:Ni
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        percent(i,j) = length(num)*100/mc.N;
    end
end
figure(); hold on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, h_{a,tgt} = 2000 km')
histogram(out(1).haf_err,'NumBins',50,'FaceColor','r')
histogram(out(2).haf_err,'NumBins',50,'FaceColor','b')
% histogram(out2.haf_err,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km','location','nw')

save -v7.3 'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\DCF\dcf_mc_ha.mat'

saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\DCF\monte carlo\dcf_mc_ha.fig')


