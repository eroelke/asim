%% nominal energy trigger
% clear; clc; close all;
% fpas = linspace(-5.6,-5.2,20)';
% % fpas = -5.4;
% ha_tgt = [2000 10000];
% % ha_tgt = 10000;
% Ni = length(fpas);
% Nj = length(ha_tgt);
% 
% traj_rate = 100;
% gnc_rate = 10;
% haf_tol = 100;
% 
% tj = nan(Ni,Nj);
% err = tj;
% 
% mc.flag = uint8(0);
% trigger = 'energy';
% 
% for j = 1:Nj % ha_tgt
%     fprintf('ha_tgt %i\n',j);
%     for i = 1:Ni % fpas
%         out = ac_venus_pd(fpas(i),ha_tgt(j),haf_tol,traj_rate,gnc_rate,trigger,mc);
%         tj(i,j) = out.tjett;    %s
%         err(i,j) = out.haf_err; %km
%     end
% end
% 
% 
% figure(); hold on
% grid on
% set(gca,'FontSize',14)
% title('Nominal Apoapsis Errors, \beta_2/\beta_1 = 9.18')
% plot(fpas,err,'LineWidth',2)
% xlabel('EFPA (deg)')
% ylabel('Apoapsis Error (km)')
% legend('h_{a,tgt} = 2000 km','h_{a,tgt} = 10000 km', 'location', 'best')
% ylim([-150 60])
% 
% saveas(gca, ... 
%     'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\pdg_err_ha.fig')


%% nominal errors vs. trigger value
clear; clc; close all;
fpas = linspace(-5.6,-5.2,20)';
% fpas = [-5.3 -5.5];
ha_tgt = 2000;
triggers = [uint8(0) uint8(1)];

Ni = length(fpas);
Nj = length(triggers);

traj_rate = 100;
gnc_rate = 5;
haf_tol = 100;

tj = nan(Ni,Nj);
err = tj;

mc.flag = uint8(0);

% for j = 1:Nj % ha_tgt
%     fprintf('ha_tgt %i\n',j);
for j = 1:Nj % triggers
    fprintf('trigger %i\n',j);
    for i  = 1:Ni % fpas
        out = ac_venus_pd(fpas(i),ha_tgt,haf_tol,traj_rate,gnc_rate,triggers(j),mc);
        tj(i,j) = out.tjett;    %s
        err(i,j) = out.haf_err; %km
    end
end
% end

figure(); hold on
grid on
set(gca,'FontSize',14)
title('Nominal Apoapsis Errors, \beta_2/\beta_1 = 9.18, h_{a,tgt} = 10000 km')
plot(fpas,err,'LineWidth',2)
xlabel('EFPA (deg)')
ylabel('Apoapsis Error (km)')
legend('Energy Trigger','Altitude Trigger, h_{a,tol} = 100 km', 'location', 'best')
ylim([-100 100])

% saveas(gca, ... 
%     'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\pdg_err_trig.fig')

% keyboard
%% monte carlo analysis, fpa
% energy trigger, fpa performance
clear; clc; close all
fpas = [-5.3 -5.5];
ha_tgt = 2000;
Ni = length(fpas);
% Nj = length(ha_tgt);

traj_rate = 100;
gnc_rate = 5;
haf_tol = 10;

% tj = nan(Ni,Nj);
% err = tj;

mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';

tols = [10 50 100 250 500];
Nj = length(tols);
percent = nan(Ni,Nj);

trigger = 'energy';

for i = 1:Ni % fpas
    out(i) = ac_venus_pd(fpas(i),ha_tgt,haf_tol,traj_rate,gnc_rate,trigger,mc);
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        percent(i,j) = length(num)*100/mc.N;
    end
end

save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\PDG\monte_carlo\percent_fpa.mat', ...
    'tols','percent')

% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, h_{a,tgt} = 2000 km')
histogram(out(1).haf_err,'NumBins',50,'FaceColor','r')
histogram(out(2).haf_err,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('EFPA = -5.3^{o}','EFPA = -5.5^{o}','location','nw')

saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\monte_carlo\pdg_mc_fpa.fig')


%% monte carlo, ha_tgt performance, energy trigger
clear; clc; close all
fpas = -5.4;
ha_tgt = [2000 10000];
Ni = length(ha_tgt);

traj_rate = 100;
gnc_rate = 5;
haf_tol = 10;

% tj = nan(Ni,Nj);
% err = tj;

mc.flag = uint8(1);
mc.N = 1000;
mc.sigs = [0.25, 0, 0.003, 0, 10, 50, 5e-4, 5e-4, 5e-4, 5e-4]';

tols = [10 50 100 250 500];
Nj = length(tols);
percent = nan(Ni,Nj);

haf_med = nan(Ni,1);
haf_std = nan(Ni,1);

trigger = 'energy';

for i = 1:Ni % ha
    out(i) = ac_venus_pd(fpas,ha_tgt(i),haf_tol,traj_rate,gnc_rate,trigger,mc);
    for j = 1:Nj
        num = find(abs(out(i).haf_err) <= tols(j));
        percent(i,j) = length(num)*100/mc.N;
    end
    haf_med(i) = median(out(i).haf);
    haf_std(i) = std(out(i).haf);
end

save('C:\Users\Evan Roelke\Documents\research\venus_aerocapture\data\PDG\monte_carlo\percent_ha.mat', ... 
    'tols','percent','haf_med','haf_std')

% histogram of results
figure(); hold on; grid on
set(gca,'FontSize',14)
title('Monte Carlo Histogram, h_{a,tgt} = 2000 km')
histogram(out(1).haf_err,'NumBins',50,'FaceColor','r')
histogram(out(2).haf_err,'NumBins',50,'FaceColor','b')
xlabel('Apoapsis Error (km)')
ylabel('Occurence')
legend('EFPA = -5.3^{o}','EFPA = -5.5^{o}','location','nw')

saveas(gca, ... 
    'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\PDG\monte_carlo\pdg_mc_fpa.fig')

