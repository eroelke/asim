%% Atmoospheric Modeling Sims
% Evan Roelke
clear;clc

%% get venus atms
% atm = load('./data/atm_data/atm_venus_mc.mat');
% mcs = nan(1000,1001);
% mcs(:,1) = atm.nom_table(:,1); %alt
% for i = 1:1000
%     mcs(:,i+1) = atm.mc(i).table(:,2);
% end
% save('./data/atm_data/atm_venus_all.mat','mcs');

%% get earth atms
% atm = load('./data/atm_data/atm_earth_gram2016.mat');
% mcs = nan(1000,1001);
% mcs(:,1) = atm.nom_table(:,1); %alt
% for i = 1:1000
%     mcs(:,i+1) = atm.mc(i).table(:,2);
% end
% save('./data/atm_data/atm_earth_all.mat','mcs');

%% run test case

% [x0,aero,gnc,sim,mc] = base_ac();
% gnc.atm_mode = uint8(0);
% mc.sigs = default_sigs();
% mc.debug = false;
% mc.flag = true;
% sim.traj_rate = 100;
% gnc.guid_rate = 0.0001;
% sim.efpa_flag = true;
% out = run_dej_n(x0,gnc,aero,sim,mc);
% 
% figure(); hold on
% plot(out.g.rho_est);
% plot(out.traj.rho,'k')

% save('ac_test.mat','out')
% keyboard


%% test case for correlation test
%{
% out = load('ac_test.mat');
% out = out.out;
% 
% profile = [out.traj.alt, out.traj.vmag, out.traj.fpa, out.traj.t, out.traj.rho];
% nominal = [out.nom.traj.alt./1000, out.nom.traj.vel_pp_mag./1000, ... 
%     out.nom.traj.gamma_pp .* 180/pi, out.nom.traj.time, out.nom.traj.rho];
% 
% % test data
% ind = out.idxend;
% prof = profile(1:ind,:);
% nom = nominal(1:ind,:);
% 
% % prop forward
% for i = 1:ind
%     alt = prof(i,1);
%     rho_nom = nom(i,end);
%     rho_true = prof(i,end);
%     
%     [L1,L3,sig]  = getAtmL(alt,'rho');
%     f(i) = 0.05 * sin(2*pi*alt);
%     
%     rL(i) = (prof(i,1)+6370)/(1000*L3);
%     erL(i) = exp(rL(i));
%     
%     rho_est(i) = rho_nom*f(i);
%     
%     if (alt >= 65 && alt <= 85)
%         wave(i) = cos((2*pi*alt/20));
%     else
%         wave(i) = 1;
%     end
%     rho_wave(i) = rho_nom*wave(i);
%     
%     L(i) = L3/sqrt(L3^3 * cos(prof(i,3))^2);
%     
%     rho_both(i) = (f(i) + wave(i))*rho_nom;
%     
%     
%     err_nom(i) = (rho_nom - rho_true)/rho_true;
%     err_est(i) = (rho_est(i) - rho_true)/rho_true;
%     err_wave(i) = (rho_wave(i) - rho_true)/rho_true;
%     err_both(i) = (rho_both(i) - rho_true)/rho_true;
%     
%     rho_K(i) = rho_nom*out.g.K_dens(i);
%     err_K(i) = (rho_K(i) - rho_true)/rho_true;
% end
% 
% figure(); hold on
% set(gca,'FontSize',14)
% title('Density vs. Time')
% plot(prof(:,4)./prof(end,4), log(prof(:,end)),'k','LineWidth',1.5)
% plot(nom(1:out.nom.idxend,4)./nom(out.nom.idxend,4), log(nom(1:out.nom.idxend,end)),'LineWidth',1.5)
% plot(prof(:,4)./prof(end,4), log(rho_K),'LineWidth',1.5)
% xlabel('Time Ratio, t_0/t_f')
% ylabel('Ln \rho')
% legend('Truth','Nominal Profile','Density Corrector')
% xlim([0 1])
% ylim([-18.5 -9.5])
% 
% figure(); hold on
% plot(err_nom,'k');
% plot(err_K)
% plot(err_both);
%}

%% 3-Sigma Error vs. Time for various estimation techniques
%{
dat = load('ac_test.mat');
dat = dat.out;
nom = dat.nom;
dat.nom = [];

for i = 1:dat.idxend
    sig(i) = ((dat.g.K_dens(i)*nom.traj.rho(i) - nom.traj.rho(i)))/nom.traj.rho(i);
    sig2(i) = (dat.traj.rho(i) - nom.traj.rho(i))/nom.traj.rho(i);
    sig3(i) = dat.g.dej_n.atm_err(i);
end


figure(); hold on
set(gca,'FontSize',14)
title('3sig error vs. Time For Each Estimation Method. Sample Monte Carlo Run','FontSize',10)
plot(dat.traj.t(1:dat.idxend)/dat.traj.t(dat.idxend),sig*3,'b','LineWidth',1.5)
plot(dat.traj.t(1:dat.idxend)/dat.traj.t(dat.idxend),sig3*3,'r','LineWidth',1.5)
plot(dat.traj.t(1:dat.idxend)/dat.traj.t(dat.idxend),sig2*3,'k','LineWidth',1.5)
legend('Density Factor','No Estimation','location','best')
xlabel('t/t_{max}')
ylabel('3-\sigma Error (%)')
ylim([-10 10])
%}

%% Nonlinear regression attempts
%{
[x0,aero,gnc,sim,mc] = base_ac();
gnc.atm_mode = uint8(5);
mc.sigs = default_sigs();
mc.debug = true;
mc.flag = true;
sim.traj_rate = 100;
gnc.guid_rate = 0.1;
out = run_dej_n(x0,gnc,aero,sim,mc);

%}


%% error bounds as func of time
% kinda not working :(
%{
temp = load('atm_indices_vs_time.mat');
atm_inds = temp.guid.dej_n.s.atm_inds;
atm_true = temp.dat.traj.rho;
t = temp.dat.traj.time;
atm_j = temp.guid.dej_n.s.atm_j;
alt = temp.dat.traj.alt;
idend = find(isnan(alt)==true,1)-1;
atm_est = temp.dat.g.atm_est;

atm = load('data/atm_data/atm_earth_gram2016.mat');
mc = atm.mc;
nom = atm.nom_table;

% figure(); hold on
% atms = nan(atm_j,1000,1000);
% k = 1;
% loop through columns of atm_inds
% for i = 1:atm_j
%     prof = atm_inds(:,i);
%     % loop through rows of column
%     for j = 1:1000
%         if (prof(j) ~= 0)
%            % load this index's mc density profile
%            %atms(i,:,k) = mc(prof(j)).table(:,2);
%            plot(nom(i,1)./1000,mc(j).table(i,2),'*');
%            k = k + 1;
%         end
%         
%     end
% end


rhos = nan(idend,1000);

% figure(); hold on
for i = 1:idend
    alt_ind = find(atm_est(:,1) <= alt(i),1,'first'); %index of atm_inds
    mc_inds = atm_inds(:,alt_ind); % mc atms possible at this alt
    
    for j = 1:1000
        if (mc_inds(j) ~= 0)
            rhos(i,j) = interp1(mc(j).table(:,1),mc(j).table(:,2),alt(i));
        end
    end
end

figure(); hold on
plot(t(1:idend),rhos,'b*');
plot(t(1:idend),atm_true(1:idend),'k','LineWidth',1.5);


figure();
semilogy(alt./1000,atm_true,'k'); hold on
semilogy(atm_est(:,1)./1000,atm_est(:,2))
%}


%% atmospheres as func of solar flux / geomagnetic index
%{
path = './data/atm_data/earth_variations/';

temp1 = load([path 'solarLow_geoLow' '.mat']);
temp2 = load([path 'solarLow_geoMid' '.mat']);
temp3 = load([path 'solarLow_geoHigh' '.mat']);
temp4 = load([path 'solarMid_geoLow' '.mat']);
temp5 = load([path 'solarMid_geoMid' '.mat']);
temp6 = load([path 'solarMid_geoHigh' '.mat']);
temp7 = load([path 'solarHigh_geoLow' '.mat']);
temp8 = load([path 'solarHigh_geoMid' '.mat']);
temp9 = load([path 'solarHigh_geoHigh' '.mat']);


entries = {'F_{10} = 70, ap = 0', ...
    'F_{10} = 70, ap = 0', ...
    'F_{10} = 70, ap = 0', ...
    'F_{10} = 140, ap = 50', ...
    'F_{10} = 140, ap = 50', ...
    'F_{10} = 140, ap = 50', ...
    'F_{10} = 210, ap = 200', ...
    'F_{10} = 210, ap = 200', ...
    'F_{10} = 210, ap = 200'};
    
entries2 = {'F_{10} = 70', ...
    'F_{10} = 140', ...
    'F_{10} = 210'};


figure(); hold on
set(gca,'FontSize',14)
yyaxis left
ylabel('Ln \rho [Ln(kg^3/m^3)]')
plot(temp1.nom_table(:,1)./1000,log(temp1.nom_table(:,2)),'LineWidth',2);
% plot(temp2.nom_table(:,1)./1000,log(temp2.nom_table(:,2)),'LineWidth',2)
% plot(temp3.nom_table(:,1)./1000,log(temp3.nom_table(:,2)),'LineWidth',2)
plot(temp4.nom_table(:,1)./1000,log(temp4.nom_table(:,2)),'LineWidth',2);
% plot(temp5.nom_table(:,1)./1000,log(temp5.nom_table(:,2)),'LineWidth',2)
% plot(temp6.nom_table(:,1)./1000,log(temp6.nom_table(:,2)),'LineWidth',2)
plot(temp7.nom_table(:,1)./1000,log(temp7.nom_table(:,2)),'LineWidth',2);
% plot(temp8.nom_table(:,1)./1000,log(temp8.nom_table(:,2)),'LineWidth',2)
plot(temp9.nom_table(:,1)./1000,log(temp9.nom_table(:,2)),'k','LineWidth',2)
% for i = 1:1000
%     plot(temp1.mc(i).table(:,1)./1000,log(temp1.mc(i).table(:,2)),'b');
% end
yyaxis right
ylabel('Temperature (K)')
plot(temp1.nom_table(:,1)./1000,temp1.nom_table(:,4),'LineWidth',2);
plot(temp4.nom_table(:,1)./1000,temp4.nom_table(:,4),'LineWidth',2);
plot(temp7.nom_table(:,1)./1000,temp7.nom_table(:,4),'LineWidth',2);
xlim([90 200])
xlabel('Altitude [km]')
legend(entries2);


% figure(); hold on
% plot(temp1.nom_table(:,1)./1000,temp1.nom_table(:,4),'LineWidth',2)
% plot(temp2.nom_table(:,1)./1000,temp2.nom_table(:,4),'LineWidth',2)
% plot(temp3.nom_table(:,1)./1000,temp3.nom_table(:,4),'LineWidth',2)
% plot(temp4.nom_table(:,1)./1000,temp4.nom_table(:,4),'LineWidth',2)
% plot(temp5.nom_table(:,1)./1000,temp5.nom_table(:,4),'LineWidth',2)
% plot(temp6.nom_table(:,1)./1000,temp6.nom_table(:,4),'LineWidth',2)
% plot(temp7.nom_table(:,1)./1000,temp7.nom_table(:,4),'LineWidth',2)
% plot(temp8.nom_table(:,1)./1000,temp8.nom_table(:,4),'LineWidth',2)   
% plot(temp9.nom_table(:,1)./1000,temp9.nom_table(:,4),'LineWidth',2)
% xlim([70 200])
% xlabel('Altitude (km)')
% ylabel('Temperature (K)')
% legend(entries);


temp1 = load([path 'solarLow_geoLow' '.mat']);
temp2 = load([path 'case1_2' '.mat']);

figure(); hold on
plot(temp1.nom_table(:,1),log(temp1.nom_table(:,2)));
plot(temp2.nom_table(:,1),log(temp2.nom_table(:,2)));

plot(temp1.nom_table(:,1),temp1.nom_table(:,2)-temp2.nom_table(:,2))

%}


%% ensemble filter base case - July 2019
%{
[x0,aero,gnc,sim,mc] = atm_sim();
mc.N = 1;
sim.traj_rate = 100;
mc.debug = true;
sim.parMode = false;
gnc.p_tol = 0.01;
gnc.t_init = 60;
fprintf('Running Ensemble Filter Case\n');
out = run_dej_n(x0,gnc,aero,sim,mc);

keyboard;
% % 
% % 
% 
% tj = out.traj.t(out.idj(1))/out.traj.t(out.idxend);
% 
% figure(); hold on
% set(gca,'FontSize',14)
% % plot(out.traj.t,out.g.rss_nom(1:length(out.traj.t)),'k');
% plot(out.traj.t./out.traj.t(out.idxend), out.g.rss_ens,'b','LineWidth',2);
% plot(out.traj.t./out.traj.t(out.idxend), out.g.rss_K,'r','LineWidth',2);
% plot([tj tj], [0 3.5e-7],'k--','LineWidth',1.5)
% xlabel('t/t_{f}');
% ylabel('RSS_{\rho} Error')
% legend('\rho_{ensemble}','\rho_{factor}','FontSize',12)

% save('ensemble_test.mat', 'out');
% 
%}

%% ensemble filter tolerance trade study
%{
% apoapsis error vs. 
% [x0,aero,gnc,sim,mc] = base_ac(true);
% gnc.t_init = 40;
% gnc.atm_mode = uint8(1);
% gnc.npc_mode = uint8(3);   % ensemble kalman
% 
% gnc.ha_tgt = 2000;
% sim.efpa_flag = true;
% sim.parMode = true;
% 
% tols = round(linspace(0.01, 0.2, 9),2)';
% % tols = linspace(0.01,0.1,5);
% Ni = length(tols);
% 
% aero.m = [100 54.4];
% aero.cds = [1.05 1.05];
% aero.rcs = [0.75 0.175];
% 
% ratios = 6:2:10;
% Nj = length(ratios);
% m2 = nan(Nj,1);
% 
% sim.traj_rate = 100;
% gnc.guid_rate = 0.5;
% mc.N = 500;
% 
% pc = nan(Nj,Ni);
% b1 = aero.m(1)/(aero.cds(1) * pi * aero.rcs(1)^2);
% m2 = ratios * (aero.m(1) * aero.rcs(2)^2)/(aero.rcs(1)^2);
% 
% haf_errs1 = nan(mc.N, Ni);
% haf_errs2 = haf_errs1;
% haf_errs3 = haf_errs1;
% haf_errs4 = haf_errs1;

% mc.debug = true;    %for now

% % ratio 1 - 6
% fprintf('Starting Beta ratio set 1...\n')
% aero.m(2) = m2(1);
% for i = 1:Ni
%     gnc.p_tol = tols(i);
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_errs1(:,i) = out.haf_err;
% end
% save('ensemble_beta6.mat','haf_errs1');

% ratio 2 - 8

% t = datetime;
% fprintf('Starting Beta ratio set 2, %i/%i, %i:%i\n', t.Month,t.Day,t.Hour,t.Minute);
% aero.m(2) = m2(2);
% for i = 1:Ni
%     gnc.p_tol = tols(i);
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_errs2(:,i) = out.haf_err;
% end
% save('ensemble_beta8.mat','haf_errs2');
% 
% % ratio 3 - 10
% t = datetime;
% fprintf('Starting Beta ratio set 3, %i/%i, %i:%i\n', t.Month,t.Day,t.Hour,t.Minute);
% aero.m(2) = m2(3);
% for i = 1:Ni
%     gnc.p_tol = tols(i);
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_errs3(:,i) = out.haf_err;
% end
% save('ensemble_beta10.mat','haf_errs3');


% % ratio 4
% fprintf('Starting Beta ratio set 4...\n')
% aero.m(2) = m2(4);
% for i = 1:Ni
%     gnc.p_tol = tols(i);
%     out = run_dej_n(x0,gnc,aero,sim,mc);
%     haf_errs4(:,i) = out.haf_err;
% end

% % % post-processing
% r6 = load('ensemble_beta6.mat');
% r8 = load('ensemble_beta8.mat');
% r10 = load('ensemble_beta10.mat');
% 
% 
% entries = {};
% for j = 1:Nj
%     entries{j} = ['\beta_2/\beta_1 = ' num2str(ratios(j))];
% end
% % 
% tol1 = 50;
% tol2 = 250;
% pc1_1 = nan(Ni,1);
% pc1_2 = pc1_1;
% pc2_1 = pc1_1;
% pc2_2 = pc1_1;
% pc3_1 = pc1_1;
% pc3_2 = pc1_1;
% for i = 1:Ni
%     % ratio 1
%     num = find(abs(r6.haf_errs1(:,i)) <= tol1);
%     pc1_1(i) = length(num)/mc.N * 100;
%     num = find(abs(r6.haf_errs1(:,i)) <= tol2);
%     pc1_2(i) = length(num)/mc.N * 100;
%     % ratio 2
%     num = find(abs(r8.haf_errs2(:,i)) <= tol1);
%     pc2_1(i) = length(num)/mc.N * 100;
%     num = find(abs(r8.haf_errs2(:,i)) <= tol2);
%     pc2_2(i) = length(num)/mc.N * 100;
%     % ratio 3
%     num = find(abs(r10.haf_errs3(:,i)) <= tol1);
%     pc3_1(i) = length(num)/mc.N * 100;
%     num = find(abs(r10.haf_errs3(:,i)) <= tol2);
%     pc3_2(i) = length(num)/mc.N * 100;
% end
% 
% % keyboard
% figure(); hold on
% set(gca,'FontSize',14);
% xlabel('Ensemble Tolerance (%)')
% ylabel('Capture Rate (%)')
% plot(tols,pc1_1,'r','LineWidth',2)
% plot(tols,pc1_2,'r--','LineWidth',2)
% plot(tols,pc2_1,'b','LineWidth',2)
% plot(tols,pc2_2,'b--','LineWidth',2)
% plot(tols,pc3_1,'g','LineWidth',2)
% plot(tols,pc3_2,'g--','LineWidth',2)
% legend(['\beta_2/\beta_1 = 6, h_{a} Tol = ' num2str(tol1) ' km'], ...
%     ['\beta_2/\beta_1 = 6, h_{a} Tol = ' num2str(tol2) ' km'], ...
%     ['\beta_2/\beta_1 = 8, h_{a} Tol = ' num2str(tol1) ' km'], ...
%     ['\beta_2/\beta_1 = 8, h_{a} Tol = ' num2str(tol2) ' km'], ...
%     ['\beta_2/\beta_1 = 10, h_{a} Tol = ' num2str(tol1) ' km'], ...
%     ['\beta_2/\beta_1 = 10, h_{a} Tol = ' num2str(tol2) ' km'], ...
%     'FontSize',12,'location','se');
% xlim([tols(1) tols(end)])

% p1=plot(tols, haf_errs1,'LineWidth',2);
% p2=plot(tols, haf_errs2, 'b','LineWidth',2);
% p3=plot(tols, haf_errs3, 'g','LineWidth',2);
% p4=plot(tols, haf_errs4, 'm','LineWidth',2);
% legend([p1,p2,p3,p4],entries)



% total results - all search tolerances/betas
% figure(); hold on
% title('2k ha tgt, 12 km/s, -5.9deg, all search tolerances')
% set(gca,'FontSize',14)
% histogram(r6.haf_errs1,'NumBins',50,'FaceColor','b')
% histogram(r8.haf_errs2,'NumBins',50,'FaceColor','r')
% histogram(r10.haf_errs3,'NumBins',50,'FaceColor','g')
% xlim([-1000 1000])
% xlabel('Apoapsis Error (km)')
% ylabel('Occurrence')
% legend(entries,'location','ne')

% one beta ratio, multiple search tolerances
% figure(); hold on
% title('2k ha tgt, 12 km/s, -5.9deg, search tolerances for single beta')
% set(gca,'FontSize',14)
% histogram(r10.haf_errs3(:,1),'NumBins',50,'FaceColor','b')
% histogram(r10.haf_errs3(:,3),'NumBins',50,'FaceColor','r')
% histogram(r10.haf_errs3(:,5),'NumBins',50,'FaceColor','g')
% xlim([-500 1500])
% xlabel('Apoapsis Error (km)')
% ylabel('Occurrence')
% legend('Search Tolerance = 1%','Search Tolerance = 5%','Search Tolerance = 10%', ... 
%     'location','ne','FontSize',12)
%}

%% Comparison Simulation - Earth, July 2019
%{
% density factor
[x0,aero,gnc,sim,mc] = atm_sim();
sim.traj_rate = 100;
gnc.atm_mode = uint8(0);
gnc.npc_mode = uint8(1);
mc.N = 500;
t = datetime();
fprintf('Factor Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save('rho_factor.mat','out');
clear out;
% density interpolator
gnc.atm_mode = uint8(1);
gnc.npc_mode = uint8(1);
t = datetime();
fprintf('Interp Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save('rho_interp.mat','out');
clear out;
% ensemble filter
gnc.t_init = 30;
% % % % 1 percent
gnc.p_tol = 0.01;
t = datetime();
fprintf('Ensemble 1%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save('ensemble_1percent.mat','out');
clear out;
% % % % 5 percent
gnc.p_tol = 0.05;
t = datetime();
fprintf('Ensemble 5%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save('ensemble_5percent.mat','out');
clear out;
% % % % 10 percent
gnc.p_tol = 0.1;
t = datetime();
fprintf('Ensemble 10%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save('ensemble_10percent.mat','out');
clear out;

% figure(); hold on
% plot(out3.traj.t(1:out3.idxend)./out3.traj.t(out3.idxend), ... 
%     out3.g.rss_K(1:out3.idxend),'b','LineWidth',2)
% plot(out3.traj.t(1:out3.idxend)./out3.traj.t(out3.idxend), ...
%     out3.g.rss_ens(1:out3.idxend),'r--','LineWidth',2)

% % % % post processing
factor = load('rho_factor.mat');
interp = load('rho_interp.mat');
ens1 = load('ensemble_1percent.mat');
ens5 = load('ensemble_5percent.mat');
ens10 = load('ensemble_10percent.mat');
% 
a = ens1.out;
% 
clear tj j;
j = 1;
for i = 1:50
    if ~isnan(a.tjett(i,1))
        tj(j) = a.tjett(i,1) / a.traj.t(a.idxend(i),i);
        j = j+1;
    end
end

for i = 1:1000
    haf_err1(i) = normrnd(median(factor.out.haf_err),std(factor.out.haf_err));
    haf_err2(i) = normrnd(median(interp.out.haf_err),std(interp.out.haf_err));
    haf_err3(i) = normrnd(median(ens1.out.haf_err),std(ens1.out.haf_err));
    haf_err4(i) = normrnd(median(ens5.out.haf_err),std(ens5.out.haf_err));
end
    
figure(); hold on
set(gca,'FontSize',14)
histogram(ens1.out.haf_err,'NumBins',100,'FaceColor','r')
histogram(interp.out.haf_err,'NumBins',100,'FaceColor','b')
histogram(factor.out.haf_err,'NumBins',100,'FaceColor','g')
legend('Ensemble Filter, 1% Tol','Density Interpolator','Density Factor','location','ne','FontSize',12)
xlabel('Apoapsis Altitude Error (km)')
ylabel('Occurrence')
xlim([-1000 1000])
% histogram(haf_err4,'NumBins',50)
%}

%% Comparison Simulation - Venus, July 2019
% density factor
[x0,aero,gnc,sim,mc] = atm_sim();
sim.planet = 'venus';
sim.h_max = 150;
sim.efpa_flag = false;
sim.data_rate = 10;

x0.h0 = 150;
x0.fpa0 = -5.4;
x0.v0 = 11;

aero.rcs = [1 0.2];
aero.m = [150 60];

for i = 1:2
    b(i) = aero.m(i)/(aero.cds(i) * pi * aero.rcs(i)^2);
end

gnc.tj0 = 100;
gnc.ha_tol = 1;

gnc.atm_mode = uint8(0);
gnc.npc_mode = uint8(1);
mc.N = 1000;

p = 'C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\HYDRA\1stage\atm_est\';
t = datetime();
fprintf('Factor Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save([p 'rho_factor.mat'],'out');
clear out;
% density interpolator
gnc.atm_mode = uint8(1);
gnc.npc_mode = uint8(1);
t = datetime();
fprintf('Interp Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save([p 'rho_interp.mat'],'out');
clear out;
% ensemble filter
gnc.t_init = 30;
% % % % 1 percent
gnc.p_tol = 0.01;
t = datetime();
fprintf('Ensemble 1%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save([p 'ensemble_1percent.mat'],'out');
clear out;
% % % % 5 percent
gnc.p_tol = 0.05;
t = datetime();
fprintf('Ensemble 5%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save([p 'ensemble_5percent.mat'],'out');
clear out;
% % % % 10 percent
gnc.p_tol = 0.1;
t = datetime();
fprintf('Ensemble 10%% Sim. %i/%i %i:%i\n',t.Month,t.Day,t.Hour,t.Minute);
out = run_dej_n(x0,gnc,aero,sim,mc);
save([p 'ensemble_10percent.mat'],'out');
clear out;












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x0,aero,gnc,sim,mc] = atm_sim()
[x0,aero,gnc,sim,mc] = base_ac(true);
gnc.t_init = 40;
gnc.atm_mode = uint8(1);
gnc.npc_mode = uint8(3);   % ensemble kalman

gnc.ha_tgt = 2000;
sim.efpa_flag = true;
sim.parMode = true;
sim.nWorkers = 6;

aero.m = [100 54.4];        % beta ratio 10
aero.cds = [1.05 1.05];
aero.rcs = [0.75 0.175];

sim.traj_rate = 250;
gnc.guid_rate = 0.5;
gnc.p_tol = 0.05;
mc.N = 1000;
end

