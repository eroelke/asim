%% venus_ac_atm.m
%   atmospheric estimation code/sims for venus aerocapture
% 
clear;clc;
p = '..\venus_ac\dej_n\atmospheric_estimation\';

%% Monte Carlo nav error scaling post-process
% ${

d100 = load([p 'data/mc/mc500_54deg_navScale_100.mat']);
d80 = load([p 'data/mc/mc500_54deg_navScale_80.mat']);
d60 = load([p 'data/mc/mc500_54deg_navScale_60.mat']);
d40 = load([p 'data/mc/mc500_54deg_navScale_40.mat']);

scales = d100.scales(1:4);
labels = {'out_k','out_di','out_ecf','out_hyb'};
names = {'DSF','DI','ECF','ECFH'};
percentile = 99.9;
Ni = length(labels);
for i = 1:Ni
   mus(:,i) =  [nanmean(d100.(labels{i}).haf_err), ... 
       nanmean(d80.(labels{i}).haf_err), ... 
    nanmean(d60.(labels{i}).haf_err), ... 
    nanmean(d40.(labels{i}).haf_err)]';
    
    stds(:,i) = [nanstd(d100.(labels{i}).haf_err(find(abs(d100.(labels{i}).haf_err) <= prctile(d100.(labels{i}).haf_err, percentile)))), ... 
       nanstd(d80.(labels{i}).haf_err(find(abs(d80.(labels{i}).haf_err) <= prctile(d80.(labels{i}).haf_err, percentile)))), ... 
    nanstd(d60.(labels{i}).haf_err(find(abs(d60.(labels{i}).haf_err) <= prctile(d60.(labels{i}).haf_err, percentile)))), ... 
    nanstd(d40.(labels{i}).haf_err(find(abs(d40.(labels{i}).haf_err) <= prctile(d40.(labels{i}).haf_err, percentile))))]';

    iqrs(:,i) = [iqr(d100.(labels{i}).haf_err), ... 
       iqr(d80.(labels{i}).haf_err), ... 
    iqr(d60.(labels{i}).haf_err), ... 
    iqr(d40.(labels{i}).haf_err)]';
end

colors = get_plot_colors(true);


figure(); hold on; grid on
set(gca,'FontSize',16)
for i = 1:Ni
    plot(scales * 100, stds(:,i),'color',colors(:,i), 'linewidth',2);
end
set ( gca, 'xDir', 'reverse' )
legend(names,'location','nw')
ylabel('Apoapsis Error \sigma (km)');
xlabel('Navigation Error Scale (%)')

%}

%% Monte Carlo: nav error scaling
% %{
clear;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 500;
% mc.mcIndex = 1;
% mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1000;
mc.sigs = default_sigs();
gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
gnc = exp_decay_ecrv(gnc);
ecrv0 = gnc.nav.ercv0;
P0 = gnc.nav.P_SS;
gnc.ecf.p_tol = 0.02;
% mc.debug = true;

scales = flip(0:0.2:1);

x0.fpa0 = -5.4;

for i = 1:length(scales)
gnc.nav.ercv0 = ecrv0 .* scales(i);
gnc.nav.P_SS = P0 .* scales(i);

gnc.atm_mode = uint8(1);    %dsf
out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(2);    %di
out_di = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    % ecf
gnc.ecf.mode = 1;	% no scaling
out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

gnc.atm_mode = uint8(4);    %ecfh
gnc.ecf.mode = 1;   % no scaling
out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% % pert on
% gnc.ecf.pert = true;
% 
% gnc.atm_mode = uint8(3);    % ecf
% gnc.ecf.mode = 1;	% no scaling
% out_ecf_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator
% 
% % hybrid free, no scaling
% gnc.atm_mode = uint8(4);
% gnc.ecf.mode = 1;   % no scaling
% out_hyb_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

save(['../venus_ac/dej_n/atmospheric_estimation/data/mc/mc' ... 
    num2str(mc.N) '_' num2str(abs(x0.fpa0*10)) 'deg_navScale_' ... 
    num2str(scales(i)*100) '.mat'])

end


%}



%% Monte Carlo Post Process
%{
if (~exist('d', 'var'))
%     d = load([p 'data/mc/mc2500_54deg.mat']);
%     d = load([p 'data/mc/mc1k_55deg.mat']);
    
%     d = load([p 'data/mc/mc2500_2k_54deg_guid05_navON_new.mat']);
%     d = d.d;

    d = load([p '\data\mc\old\mc2500_2k_54deg_guid05_navON.mat']);

end

errs = get_default_atm_errs(d);
% % histogram
% labels = default_atm_labels();
% haf_errs = [err_k err_di err_ecf err_ecfh err_ecf_p err_ecfh_p];
% histogram_plot(200, haf_errs, labels, 2000, [-200 200])

% % err vs tj
% figure(); hold on
% set(gca,'FontSize',14)
% plot(d.out_k.tjr(:,1), abs(err_k),'b+')
% plot(d.out_di.tjr(:,1), abs(err_di),'r^')
% plot(d.out_ecf.tjr(:,1), abs(err_ecf),'g*')
% % plot(d.out_hyb.tjr(:,1), abs(err_ecfh),'m*')
% 
% ylim([0 500])
% xlabel('t/t_f')
% ylabel('|\Delta r_a| (km)')
% legend('DSF','DI','ECF')


% % capture rate
haTgt = 2000;
tols = 0:0.1:100;
N = d.mc.N;
pc = nan(length(tols),size(errs,2));
for i = 1:length(tols)
    for j = 1:size(errs,2)
        pc(i,j) = 100 * length(find((100*abs(errs(:,j))/haTgt) <= tols(i))) / N;
    end
end

colors = nan(3,8);
colors(:,1) = [0;0;0];
colors(:,2:8) = get_plot_colors();

labels = default_atm_labels();

figure(); hold on
grid on
set(gca,'FontSize',16)
title(['v' num2str(d.x0.v0) ', h_a^* = ' num2str(haTgt) ', EFPA = ' num2str(d.x0.fpa0)])
for j = 1:size(errs,2)
    plot(tols, pc(:,j),'color',colors(:,j),'LineWidth',2.5)
end
xlim([0 10])
legend(labels);
xlabel('Tolerance (%)')
ylabel('Capture Rate')


% % ecf inds
get_ecf_accuracies(d);

%}

%% Post Process- Box Plots
%{

d = load([p 'data/mc/mc2500_2k_54deg_guid05_navON.mat']);

errs = [d.out_k.haf_err, d.out_di.haf_err, d.out_ecf.haf_err, d.out_hyb.haf_err, d.out_ecf_p.haf_err, d.out_hyb_p.haf_err];


g1 = repmat({'DSF'},size(errs,1),1);
g2 = repmat({'DI'},size(errs,1),1);
g3 = repmat({'ECF'},size(errs,1),1);
g4 = repmat({'ECFH'},size(errs,1),1);
g5 = repmat({'ECF-P'},size(errs,1),1);
g6 = repmat({'ECFH-P'},size(errs,1),1);
g = [g1;g2;g3;g4;g5;g6];

boxplot(errs,g);
grid on
set(gca,'FontSize',14)
ylabel('Apoapsis Error (km)')

%}

%% Post Process 'Nominal' Monte Carlo - dispersed GRAM index vs. err
%{
d = load([p 'data/entry_trades/2k_54deg_gramIdxTrade.mat']);
% d = load([p 'data/entry_trades/2k_55deg_gramIdxTrade.mat']);

inds = 1:1:1000;

haf_errs = get_default_atm_errs(d);
labels = default_atm_labels();

histogram_plot(250, haf_errs, labels, 2000, [-100 100])

% figure(); hold on; grid on
% title(['EFPA Trade Study, GRAM Index ' num2str(d.mc.mcIndex)])
% set(gca,'FontSize',14)
% plot(1:1000, err_k,'k--*','LineWidth',2)
% plot(1:1000, err_di,'r*','LineWidth',2)
% plot(1:1000, err_ecf,'b*','LineWidth',2)
% plot(1:1000, err_ecfh,'g*','LineWidth',2)
% plot(1:1000, err_ecf_p,'b*','LineWidth',2)
% plot(1:1000, err_ecfh_p,'g*','LineWidth',2)
% legend('DSF','DI','ECF','ECFH','location','ne')
% xlabel('Perturbed GRAM Index')
% ylabel('Apoapsis Error (km)')
% xlim([-5.7, -5.35])
% ylim([-200 500])


%}

%% 'Nominal' Monte Carlo - dispersed GRAM index vs. error
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1000;
mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1000;
mc.sigs = zero_sigs();
% gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
% gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;

sim.traj_rate = 100;
sim.data_rate = 1;
gnc.guid_rate = 0.5;

ha_tgts = 2000;

% efpas = [-5.35];
efpas = -5.4;

Ni = length(efpas);
Nj = length(ha_tgts);
% Nk = length(inds);

err_k = nan(Ni,Nj); err_di = err_k; err_ecf = err_k; err_ecfh = err_k; err_ecf_p = err_k; err_ecfh_p = err_k;
% tj_k = nan(Ni,Nj); tj_di = tj_k; tj_ecf = tj_k; tj_ecfh = tj_k;
tjr_k = nan(Ni,Nj); tjr_di = tjr_k; tjr_ecf = tjr_k; tjr_ecfh = tjr_k; tjr_ecf_p = tjr_k; tjr_ecfh_p = tjr_k;
% dv_k = nan(Ni,Nj); dv_di = dv_k; dv_ecf = dv_k; dv_hyb = dv_k;
dvCirc_k = nan(Ni,Nj); dvCirc_di = dvCirc_k; dvCirc_ecf = dvCirc_k; dvCirc_ecfh = dvCirc_k; dvCirc_ecf_p = dvCirc_k; dvCirc_ecfh_p = dvCirc_k;
for i = 1:Ni % efpa
    fprintf('efpa %4.1f deg\n',efpas(i));
    x0.fpa0 = efpas(i);
    for j = 1:Nj %tgt
        gnc.ha_tgt = ha_tgts(j);
        % dsf
        gnc.atm_mode = uint8(1);    %dsf
        out_k = run_dej_n(x0,gnc,aero,sim,mc);
        err_k = out_k.haf_err;
        tjr_k = out_k.tjr(:,1);
        dvCirc_k = out_k.dv;
        
        % di
        gnc.atm_mode = uint8(2);    %di
        out_di = run_dej_n(x0,gnc,aero,sim,mc);
        err_di = out_di.haf_err;
        tjr_di = out_di.tjr(:,1);
        dvCirc_di = out_di.dv;
        
        % ecf free, no scaling
        gnc.atm_mode = uint8(3);    % ecf
        gnc.ecf.mode = 1;	% no scaling
        out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        err_ecf = out_ecf.haf_err;
        tjr_ecf = out_ecf.tjr(:,1);
        dvCirc_ecf = out_ecf.dv;
        
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_ecfh = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        err_ecfh = out_ecfh.haf_err;
        tjr_ecfh = out_ecfh.tjr(:,1);
        dvCirc_ecfh = out_ecfh.dv;
        
        gnc.ecf.pert = true;
        
        % ecf free, no scaling, pert
        gnc.atm_mode = uint8(3);    % ecf
        gnc.ecf.mode = 1;	% no scaling
        out_ecf_p = run_dej_n(x0,gnc,aero,sim,mc);        
        err_ecf_p = out_ecf_p.haf_err;
        tjr_ecf_p = out_ecf_p.tjr(:,1);
        dvCirc_ecf_p = out_ecf_p.dv;
        
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_ecfh_p = run_dej_n(x0,gnc,aero,sim,mc);     
        err_ecfh_p = out_ecfh_p.haf_err;
        tjr_ecfh_p = out_ecfh_p.tjr(:,1);
        dvCirc_ecfh_p = out_ecfh_p.dv;
        
%         keyboard;
        save(['../venus_ac/dej_n/atmospheric_estimation/data/entry_trades/' ...
            num2str(gnc.ha_tgt/1000) 'k_' num2str(abs(x0.fpa0*10)) 'deg_gramIdxTrade.mat']);

    end %j - ha tgt
end %i - efpa

keyboard
%}

%% Monte Carlo: pert ON/OFF, nav errors
%{
clear;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 2500;
% mc.mcIndex = 1;
% mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1000;
mc.sigs = default_sigs();
gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;
% mc.debug = true;

x0.fpa0 = -5.4;

gnc.atm_mode = uint8(1);    %dsf
out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(2);    %di
out_di = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    % ecf
gnc.ecf.mode = 1;	% no scaling
out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

gnc.atm_mode = uint8(4);    %ecfh
gnc.ecf.mode = 1;   % no scaling
out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% pert on
gnc.ecf.pert = true;

gnc.atm_mode = uint8(3);    % ecf
gnc.ecf.mode = 1;	% no scaling
out_ecf_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% hybrid free, no scaling
gnc.atm_mode = uint8(4);
gnc.ecf.mode = 1;   % no scaling
out_hyb_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

save(['../venus_ac/dej_n/atmospheric_estimation/data/mc/mc' ... 
    num2str(mc.N) '_' num2str(abs(x0.fpa0*10)) 'deg.mat'])

%}


%%  ECF performance vs. efpa - 'nominal', perturbed profile
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1;

low_ind = 12;   %rss 19 (non-pert is ~5)
mid_ind = 500;  %rss ~44
high_ind = 998;   %rss ~97 (non-pert is ~20)

inds = [mid_ind, high_ind];
% inds = mid_ind;

% mc.mcIndex = 500;
% mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1200;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
% gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;

sim.traj_rate = 100;
sim.data_rate = 100;
gnc.guid_rate = 1;

ha_tgts = 2000;
efpas = linspace(-5.7,-5.2, 25);
Ni = length(efpas);
Nj = length(ha_tgts);
Nk = length(inds);

for k = 1:Nk %GRAM ind
mc.mcIndex = inds(k);
    
err_k = nan(Ni,Nj); err_di = err_k; err_ecf = err_k; err_hyb = err_k;
tj_k = nan(Ni,Nj); tj_di = tj_k; tj_ecf = tj_k; tj_hyb = tj_k;
tjr_k = nan(Ni,Nj); tjr_di = tjr_k; tjr_ecf = tjr_k; tjr_hyb = tjr_k;
dv_k = nan(Ni,Nj); dv_di = dv_k; dv_ecf = dv_k; dv_hyb = dv_k;
dvCirc_k = nan(Ni,Nj); dvCirc_di = dvCirc_k; dvCirc_ecf = dvCirc_k; dvCirc_hyb = dvCirc_k;
for j = 1:Nj % ha tgt
    fprintf('ha tgt %4.0f\n',ha_tgts(j));
    gnc.ha_tgt = ha_tgts(j);
    for i = 1:Ni %efpa
        x0.fpa0 = efpas(i);        
        % dsf
        gnc.atm_mode = uint8(1);    %dsf
        out_k = run_dej_n(x0,gnc,aero,sim,mc);
        [err_k(i,j), tj_k(i,j), tjr_k(i,j), dv_k(i,j), dvCirc_k(i,j)] = ...
            get_ac_results(out_k);
        
        % di
        gnc.atm_mode = uint8(2);    %di
        out_di = run_dej_n(x0,gnc,aero,sim,mc);
        [err_di(i,j), tj_di(i,j), tjr_di(i,j), dv_di(i,j), dvCirc_di(i,j)] = ...
            get_ac_results(out_di);
        
        % ecf free, no scaling
        gnc.atm_mode = uint8(3);    % ecf
        gnc.ecf.mode = 1;	% no scaling
        gnc.ecf.pert = true;  %perturbed profiles
        out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_ecf(i,j), tj_ecf(i,j), tjr_ecf(i,j), dv_ecf(i,j), dvCirc_ecf(i,j)] = ...
            get_ac_results(out_ecf);
        
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_hyb(i,j), tj_hyb(i,j), tjr_hyb(i,j), dv_hyb(i,j), dvCirc_hyb(i,j)] = ...
            get_ac_results(out_hyb);
        
    end %i
end %j

save(['../venus_ac/dej_n/atmospheric_estimation/data/entry_trades/' ...
    'pert_gramIdx_' num2str(inds(k)) '.mat']);

end %k

keyboard
%}

%% add simulated variations to density profile
%{
if ~exist('atm','var')
    atm = load('./data/atm_data/atm_venus_mc.mat');
    mc = atm.mc;
    nom = atm.nom_table;
end
alt = nom(:,1)./1000;
N = 1000;

K = nan(N,length(nom(:,2)));
for i = 1:N
    dens = mc(i).table(:,2);
    K(:,i) = dens ./ nom(:,2);
end

% save('./data/atm_data/atm_venus_K.mat','K');

indices = 1:1000;
Amax = 0.1;
max = 150;
min = 80;
Knew = nan(1000,1000);
pert = Knew;
rss = nan(1000,1);
A = nan(1000,1); shift = A; period = A;
for i = 1:1000 %profile
    A(i) = Amax*rand(1); % random amplitude
    shift(i) = 2*pi*rand(1);
    period(i) = (max - min) * randi([1 2]);
    temp = 0;
    for j = 1:1000  %alt
        pert(j,i) = A(i)*sin(2*pi*alt(j)/period(i) + shift(i));
        Knew(j,i) = K(j,i) + pert(j,i);
        temp = temp + (Knew(j,i) - nom(j,2))^2;
    end
    rss(i) = sqrt(temp);
%     Knew(:,i) = K(:,i) + pert(:,i);
end

% v = period ./ 70;
% lab = 'Period';
% 
% figure(); hold on; grid on
% set(gca,'FontSize',14)
% xlabel('Perturbed GRAM Index')
% plot(indices, v, '*')
% ylabel(lab);
% yline(mean(v),'linewidth',1.5);
% yline(mean(v) + 3*std(v),'k--','linewidth',1.5);
% yline(mean(v) - 3*std(v),'k--','linewidth',1.5);

ind = find(A >= 0.095);
ind = ind(randi([1 length(ind)]));

figure(); hold on
set(gca,'FontSize',14)
plot(pert, alt,'Color',[0 0 0] + 0.4)
plot(pert(:,ind),alt,'r','linewidth',1.5)
ylim([min max]);
xlabel('K_{pert}')
ylabel('Altitude (km)')


figure(); hold on; grid on
title(['Amplitude=' num2str(round(A(ind),2)) ',period=2\pi/' num2str(period(ind)) ... 
    ' km, shift=' num2str(round(shift(ind)/pi,2)) '\pi'])
set(gca,'FontSize',16)
plot(K(:,ind),alt,'k','linewidth',1.5)
plot(Knew(:,ind),alt,'b','linewidth',1.5);
ylim([80 150]);
ylabel('Altitude (km)')
xlabel('Density Variation')
legend('GRAM','Perturbed GRAM','location','ne')

% perturb densities
mc_new = mc;
for i = 1:1000 %profile
    mc_new(i) = mc(i);
    for j = 1:1000 %alt
        rho(j,i) = mc(i).table(j,2);
        mc_new(i).table(j,2) = Knew(i,j) * mc(i).table(j,2);
        rho_new(j,i) = mc_new(i).table(j,2);
    end
end

% figure();
% title(['Amplitude=' num2str(round(A(ind),2)) ',period=2\pi/' num2str(period(ind)) ... 
%     ' km, shift=' num2str(round(shift(ind)/pi,2)) '\pi'])
% % set(gca,'FontSize',16)
% semilogy(alt,rho(:,ind),'b','linewidth',1); hold on
% semilogy(alt,rho_new(:,ind),'k','linewidth',1);
% ylim([80 150]);
% ylabel('Altitude (km)')
% xlabel('Density Variation')
% legend('GRAM','Perturbed GRAM','location','ne')

% save new profiles
nom_table = atm.nom_table;
mc = mc_new;
save('./data/atm_data/atm_venus_mc_pert.mat','nom_table','mc');

% ECF data array
atm = load('./data/atm_data/atm_venus_mc_pert.mat');
alt = atm.nom_table(:,1);
mc = atm.mc;
mcs = nan(1000,1001);
mcs(:,1) = alt;
for i = 1:1000 %prof
    mcs(:,i+1) = mc(i).table(:,2);
end
% 
save('./data/atm_data/atm_venus_all_pert.mat','mcs');

% keyboard
%}


%% performance vs. ha tgt - 'nominal'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1;

low_ind = 16;   %rss ~5
mid_ind = 500;  %rss ~11
high_ind = 516;   %rss ~20


inds = [low_ind high_ind];
% inds = mid_ind;

% mc.mcIndex = 500;
% mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1000;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
% gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;

sim.traj_rate = 100;
sim.data_rate = 100;
gnc.guid_rate = 1;


ha_tgts = logspace(3,5,30); ha_tgts = ha_tgts * .4;

% efpas = [-5.4, -5.6];
efpas = -5.5;

Ni = length(efpas);
Nj = length(ha_tgts);
Nk = length(inds);

for k = 1:Nk %GRAM ind
mc.mcIndex = inds(k);
    
err_k = nan(Ni,Nj); err_di = err_k; err_ecf = err_k; err_hyb = err_k;
tj_k = nan(Ni,Nj); tj_di = tj_k; tj_ecf = tj_k; tj_hyb = tj_k;
tjr_k = nan(Ni,Nj); tjr_di = tjr_k; tjr_ecf = tjr_k; tjr_hyb = tjr_k;
dv_k = nan(Ni,Nj); dv_di = dv_k; dv_ecf = dv_k; dv_hyb = dv_k;
dvCirc_k = nan(Ni,Nj); dvCirc_di = dvCirc_k; dvCirc_ecf = dvCirc_k; dvCirc_hyb = dvCirc_k;
for i = 1:Ni % efpa
    fprintf('efpa %4.1f deg\n',efpas(i));
    x0.fpa0 = efpas(i);
    for j = 1:Nj %tgt
        gnc.ha_tgt = ha_tgts(j);
        % dsf
        gnc.atm_mode = uint8(1);    %dsf
        out_k = run_dej_n(x0,gnc,aero,sim,mc);
        [err_k(i,j), tj_k(i,j), tjr_k(i,j), dv_k(i,j), dvCirc_k(i,j)] = ...
            get_ac_results(out_k);
        
        % di
        gnc.atm_mode = uint8(2);    %di
        out_di = run_dej_n(x0,gnc,aero,sim,mc);
        [err_di(i,j), tj_di(i,j), tjr_di(i,j), dv_di(i,j), dvCirc_di(i,j)] = ...
            get_ac_results(out_di);
        
        % ecf free, no scaling
        gnc.atm_mode = uint8(3);    % ecf
        gnc.ecf.mode = 1;	% no scaling
        out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_ecf(i,j), tj_ecf(i,j), tjr_ecf(i,j), dv_ecf(i,j), dvCirc_ecf(i,j)] = ...
            get_ac_results(out_ecf);
        
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_hyb(i,j), tj_hyb(i,j), tjr_hyb(i,j), dv_hyb(i,j), dvCirc_hyb(i,j)] = ...
            get_ac_results(out_hyb);
        
    end %j - ha tgt
end %i - efpa

save(['../venus_ac/dej_n/atmospheric_estimation/data/entry_trades/' ...
    num2str(abs(x0.fpa0*10)) 'deg_haTrade_gramIdx_' num2str(inds(k)) '.mat']);

end %k
%}


%% nominal performance vs. haTGt post process
%{
d = load([p 'data\entry_trades\haTrade_gramIdx_516.mat']);

ind = 1;    %-5.4deg

figure(); hold on; grid on
title(['EFPA Trade Study, GRAM Index ' num2str(d.mc.mcIndex)])
set(gca,'FontSize',14)
plot(d.ha_tgts, d.err_k(ind,:),'k--','LineWidth',2)
plot(d.ha_tgts, d.err_di(ind,:),'r','LineWidth',2)
plot(d.ha_tgts, d.err_ecf(ind,:)','b','LineWidth',2)
plot(d.ha_tgts, d.err_hyb(ind,:)','g','LineWidth',2)
legend('DSF','DI','ECF','ECFH','location','ne')
xlabel('Apoapsis Altitude Target (km)')
ylabel('Apoapsis Error (km)')
% xlim([-5.7, -5.35])
ylim([-200 500])


%}


%%  performance vs. efpa - 'nominal'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1;

low_ind = 16;   %rss ~5
mid_ind = 500;  %rss ~11
high_ind = 516;   %rss ~20

inds = [low_ind high_ind];
% inds = mid_ind;

% mc.mcIndex = 500;
% mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 1000;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
% gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;

sim.traj_rate = 100;
sim.data_rate = 100;
gnc.guid_rate = 1;

ha_tgts = 2000;
efpas = linspace(-5.8,-5.2, 30);
Ni = length(efpas);
Nj = length(ha_tgts);
Nk = length(inds);

for k = 1:Nk %GRAM ind
mc.mcIndex = inds(k);
    
err_k = nan(Ni,Nj); err_di = err_k; err_ecf = err_k; err_hyb = err_k;
tj_k = nan(Ni,Nj); tj_di = tj_k; tj_ecf = tj_k; tj_hyb = tj_k;
tjr_k = nan(Ni,Nj); tjr_di = tjr_k; tjr_ecf = tjr_k; tjr_hyb = tjr_k;
dv_k = nan(Ni,Nj); dv_di = dv_k; dv_ecf = dv_k; dv_hyb = dv_k;
dvCirc_k = nan(Ni,Nj); dvCirc_di = dvCirc_k; dvCirc_ecf = dvCirc_k; dvCirc_hyb = dvCirc_k;
for j = 1:Nj % ha tgt
    fprintf('ha tgt %4.0f\n',ha_tgts(j));
    gnc.ha_tgt = ha_tgts(j);
    for i = 1:Ni %efpa
        x0.fpa0 = efpas(i);        
        % dsf
        gnc.atm_mode = uint8(1);    %dsf
        out_k = run_dej_n(x0,gnc,aero,sim,mc);
        [err_k(i,j), tj_k(i,j), tjr_k(i,j), dv_k(i,j), dvCirc_k(i,j)] = ...
            get_ac_results(out_k);
        
        % di
        gnc.atm_mode = uint8(2);    %di
        out_di = run_dej_n(x0,gnc,aero,sim,mc);
        [err_di(i,j), tj_di(i,j), tjr_di(i,j), dv_di(i,j), dvCirc_di(i,j)] = ...
            get_ac_results(out_di);
        
        % ecf free, no scaling
        gnc.atm_mode = uint8(3);    % ecf
        gnc.ecf.mode = 1;	% no scaling
        out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_ecf(i,j), tj_ecf(i,j), tjr_ecf(i,j), dv_ecf(i,j), dvCirc_ecf(i,j)] = ...
            get_ac_results(out_ecf);
        
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        [err_hyb(i,j), tj_hyb(i,j), tjr_hyb(i,j), dv_hyb(i,j), dvCirc_hyb(i,j)] = ...
            get_ac_results(out_hyb);
        
    end %i
end %j

save(['../venus_ac/dej_n/atmospheric_estimation/data/entry_trades/' ...
    'gramIdx_' num2str(inds(k)) '.mat']);

end %k
%}

%% nominal performance vs. efpa post process
%{
% d = load([p 'data\entry_trades\gramIdx_516.mat']);

% d = load([p 'data\entry_trades\pert_gramIdx_500.mat']);
d = load([p 'data\entry_trades\pert_gramIdx_12.mat']);

ind = 1;    %2000 km tgt

figure(); hold on; grid on
title(['EFPA Trade Study, GRAM Index ' num2str(d.mc.mcIndex)])
set(gca,'FontSize',14)
plot(d.efpas, d.err_k(:,ind),'k--','LineWidth',2)
plot(d.efpas, d.err_di(:,ind),'r','LineWidth',2)
plot(d.efpas, d.err_ecf(:,ind)','b','LineWidth',2)
plot(d.efpas, d.err_hyb(:,ind)','g','LineWidth',2)
legend('DSF','DI','ECF','ECFH','location','ne')
xlabel('EFPA (deg)')
ylabel('Apoapsis Error (km)')
xlim([-5.7, -5.35])
ylim([-200 500])


%}

%% preallocate ecrv errors
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 100;
mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 600;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
gnc = exp_decay_ecrv(gnc, mc.N, sim, false);    % identical initial nav errors


%}

%% investigate GRAM indices
%{
temp = load(['./data/atm_data/atm_venus_mc.mat']);
nom = temp.nom_table;
mc = temp.mc;

inds = [427];

alt = nom(:,1)./1000;
K = nan(length(inds),length(nom(:,2)));
for i = 1:length(inds)
    dens = mc(inds(i)).table(:,2);
    for j = 1:length(dens)
        K(i,j) = dens(j) / nom(j,2);
    end
end


%}

%% Check potential of hybrid method - cheating mode, MC
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 100;
% mc.mcIndex = 1;
mc.ordered = true;
sim.ignore_nominal = true;
sim.t_max = 600;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
tvec = (0 : 1/sim.traj_rate : sim.t_max)'; % Trajectory time vector, s
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);    % identical initial nav errors
% gnc = exp_decay_ecrv(gnc, mc.N, sim, true);
gnc = exp_decay_ecrv(gnc);
gnc.ecf.p_tol = 0.02;
% mc.debug = true;

% dsf
gnc.atm_mode = uint8(1);    %dsf
out_k = run_dej_n(x0,gnc,aero,sim,mc);

% di
gnc.atm_mode = uint8(2);    %di
out_di = run_dej_n(x0,gnc,aero,sim,mc);

% ecf free, no scaling
gnc.atm_mode = uint8(3);    % ecf
gnc.ecf.mode = 1;	% no scaling
out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% hybrid free, no scaling
gnc.atm_mode = uint8(4);
gnc.ecf.mode = 1;   % no scaling
out_hyb = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

lab{1} = ['DSF, 1\sigma=' num2str(nanstd(out_k.haf_err))];
lab{2} = ['DI, 1\sigma=' num2str(nanstd(out_di.haf_err))];
lab{3} = ['ECF, 1\sigma=' num2str(nanstd(out_ecf.haf_err))];
lab{4} = ['ECF Hybrid, 1\sigma=' num2str(nanstd(out_hyb.haf_err))];

histogram_plot(mc.N / 5, [out_k.haf_err, out_di.haf_err, ... 
    out_ecf.haf_err, out_hyb.haf_err], lab, 2000);

% for i = 1:mc.N
%     vSig(i) = norm(gnc.nav.ecrv0(4:6,i));
% end
% 
% % plot(out_ecf.g.mc_inds, vSig,'*')
% 
% figure(); hold on; grid on
% set(gca,'FontSize',14)
% p1=plot(out_ecf.g.mc_inds, out_ecf.haf_err,'b*');
% p0=plot(out_di.g.mc_inds, out_di.haf_err,'r*');
% 
% 
% figure(); hold on; grid on
% set(gca,'FontSize',14)
% p1=plot(out_ecf.tjr(:,1),out_ecf.haf_err,'b*');
% p0=plot(out_di.tjr(:,1),out_di.haf_err,'r*');
% legend([p0 p1],'DI','ECF','location','se')
% xlabel('t_{jett}/t_f')
% ylabel('Apoapsis Error (km)')

% save('../venus_ac/dej_n/atmospheric_estimation/data/ecf/test_mc.mat')

%}

%% Uniqueness of GRAM
%{
atm = load('./data/atm_data/atm_venus_mc.mat');
% atm = load('./data/atm_data/atm_venus_mc_pert.mat');

alt = atm.nom_table(:,1);
nom = atm.nom_table(:,2);
% disp = nan(1000,1000);
rss_err = nan(1000,1);
indices = 1:1:1000;
for i = 1:1000 %profile
    disp = atm.mc(i).table(:,2);

    % compute rss error to nom
    rss = 0;
    for j = 1:1000 %alt
        err = (disp(j) - nom(j))^2;
        rss = rss + err;
    end
    rss_err(i) = sqrt(rss);
end

min_ind = find(rss_err == min(rss_err));
max_ind = find(rss_err == max(rss_err));

save('../venus_ac/atmosphere/gram_rss_errors.mat');


figure(); hold on; grid on
title('VenusGRAM RSS Errors')
set(gca,'FontSize',14)
% histogram(rss_err,'NumBins',50)
plot(indices, rss_err,'*')
% plot(100, rss_err(100),'k*','LineWidth',2)
% plot(369,rss_err(369),'r*','LineWidth',2)
% plot(68,rss_err(68),'b*','LineWidth',2)
% plot(227,rss_err(227),'g*','LineWidth',2)
% plot(780,rss_err(780),'y*','LineWidth',2)
% plot(798,rss_err(798),'c*','LineWidth',2)
xlabel('GRAM Index')
ylabel('RSS Error')
% legend('Range','Truth-100','369','68','227','780','798')

%}

%% DI/ECF hybrid - mc
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 1;
mc.mcIndex = randi([1 1000]);
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = true;
sim.parMode = 0;
% gnc.pred_mode = 1;  % verbose
% gnc = exp_decay_ecrv(gnc);
% 
% gnc.atm_mode = uint8(1);    %dsf
% out_k = run_dej_n(x0,gnc,aero,sim,mc);
% 
% gnc.atm_mode = uint8(2);    %di
% out_di = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    %hybrid
gnc.ecf.p_tol = 0.1;
out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator


% save('../venus_ac/dej_n/atmospheric_estimation/data/hybrid/mc_test.mat');

for i = 1:mc.N
    if ~isnan(out_ecf.idj(i,1))
        ind_j(i) = out_ecf.g.atm.ind_curr(out_ecf.idj(i,1),i);
    else
        ind_j(i) = out_ecf.g.atm.ind_curr(out_ecf.idxend(i),i);
    end
end

% figure(); hold on; grid on
% histogram(ind_j, 'NumBins',200)
% xline(100);
% 
% figure(); hold on; grid on
% set(gca,'FontSize',14)
% plot(ind_j, out_ecf.haf_err,'+','LineWidth',1.5)
% xlabel('GRAM Index at Jettison')
% ylabel('Apoapsis Error (km)')
% 
% figure(); hold on; grid on
% set(gca,'FontSize',14)
% plot(out_ecf.traj.t, out_ecf.g.atm.ind_curr,'LineWidth',1.2)
% p1=plot(out_ecf.tjett(:,1), ind_j,'k*')
% yline(100);
% xlabel('Time (s)')
% ylabel('GRAM Index')
% legend([p1],'Jettison')
%}

%% DI/ECF hybrid - hybrid sigma v vs. performance
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 100;
mc.N = 1;
sim.ignore_nominal = true;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% gnc.pred_mode = 1;  % verbose
gnc.atm_mode = uint8(4);    %hybrid
gnc.ecf.p_tol = 0.02;

sigs = linspace(0, 1, 100)'; % percent

out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator
err(1) = out.haf_err;
t(:,1) = out.traj.t;
inds(:,1) = out.g.atm.ind_curr;

% add error
gnc = exp_decay_ecrv(gnc);
for i = 2:length(sigs)
gnc0 = gnc;
gnc0.nav.P_SS = gnc0.nav.P_SS .* sigs(i);
out = run_dej_n(x0,gnc0,aero,sim,mc);    %perfect nav, interpolator

err(i) = out.haf_err;
t(:,i) = out.traj.t;
inds(:,i) = out.g.atm.ind_curr;

end

% ECF index over time
figure(); hold on; grid on
plot(t, inds,'LineWidth',1.5)
yline(100,'k--','LineWidth',1.5)
xlabel('Time (s)')
ylabel('GRAM Index')
title(['GRAM Index ' num2str(mc.mcIndex)]);

% error vs. sigs
figure(); hold on; grid on
set(gca,'FontSize',14)
plot(sigs * 100, err,'+','LineWidth',2)
xlabel('\sigma scale, %')
ylabel('Apoapsis Error (km)')
title(['GRAM Index ' num2str(mc.mcIndex)]);

%}

%% DI/ECF hybrid - mc on single index
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 100;
mc.N = 1;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = true;
gnc = exp_decay_ecrv(gnc);
gnc.atm_mode = uint8(4);    %hybrid
gnc.ecf.p_tol = 0.001;
% mc.debug = true;

N = 100;
err = nan(N,1);

for i = 1:N
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator
err(i) = out.haf_err;
t(:,i) = out.traj.t;
inds(:,i) = out.g.atm.ind_curr;

ecrv0(:,i) = out.g.nav.ecrv0;
end

% save('../venus_ac/dej_n/atmospheric_estimation/data/hybrid/gram_idx100_mc.mat')
for i = 1:N
%     idx0 = find(~isnan(inds(:,i)),1);
    ind_f(i) = inds(find(isnan(t(:,i)),1)-1,i);
end

% figure(); hold on; grid on
% set(gca,'FontSize',14)
% histogram(ind_f,'NumBins',100);
% xline(mc.mcIndex,'k--','Color',[0 0 0] + 0.5)
% 
figure(); hold on; grid on
title(['Tol: ' num2str(100*gnc.ecf.p_tol) '%, Success Rate: ' num2str(length(find(ind_f == 100))) ' /' ... 
    num2str(N) ', (' ... 
    num2str(100 * length(find(ind_f == 100)) / N) '%), 1\sigma=' ... 
    num2str(round(std(err),1)) ' km'])
set(gca,'FontSize',14)
p0=plot(ind_f, err,'+','LineWidth',2);
yline(0)
xlabel('Final Ecf Index')
ylabel('Apoapsis Error (km)')
xline(mc.mcIndex,'k--','Color',[0 0 0] + 0.5)


%try and get entry uncertainty relationship to final index
for i = 1:N
    rSig(i) = norm(ecrv0(1:3,i));
    vSig(i) = norm(ecrv0(4:6,i));
    aSig(i) = norm(ecrv0(7:9,i));
    sig(i) = norm(ecrv0(:,i));
end

figure(); hold on; grid on
set(gca,'FontSize',14)
xline(100,'k--','LineWidth',1)
plot(ind_f, vSig,'+','LineWidth',2)
xlabel('Final GRAM Index')
ylabel('|\sigma_v| (m/s)')

%}

%% DI/ECF hybrid - nominal
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);

% mc.mcIndex = 1;
mc.ordered = false;
% mc.mcIndex = 454;   %427    % problem children
mc.N = 1;
mc.mcIndex = 1;
sim.ignore_nominal = true;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
sim.traj_rate = 100;
% sim.data_rate = 10;
gnc.guid_rate = 0.5;
gnc.tj0 = 120;
x0.fpa0 = -5.4;
gnc.pred_mode = 1;  % verbose
% sim.atm_mode = uint8(2);  %no winds
sim.atm_mode = uint8(3);    %winds
% sim.atm_mode = uint8(4);    %perfect

gnc = exp_decay_ecrv(gnc);

% gnc.atm_mode = uint8(1);    %k
% out_k = run_dej_n(x0,gnc,aero,sim,mc);
% 
% gnc.atm_mode = uint8(2);    %DI
% out_di = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    %ecf
% gnc.ecf.p_tol = 0.1;
gnc.ecf.mode = 0x81;    %give index, no scaling
% gnc.pred_mode = 1;
% mc.debug = true;
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% figure(); hold on
% yyaxis left
% plot(out.traj.t,out.g.atm.ind_rss)
% yyaxis right
% plot(out.traj.t,out.g.atm.ind_curr)


% figure(); hold on; grid on
% plot(out_k.traj.t,out_k.g.K_dens,'Linewidth',2)
% plot(out_di.traj.t,out_di.g.K_model,'Linewidth',2)
% plot(out_ecf.traj.t,out_ecf.g.K_model,'Linewidth',2)
% plot(out_ecf.traj.t,out_ecf.g.K_true,'k--','LineWidth',2)
% xlabel('Time')
% ylabel('rho_{ECF} Error')
% legend('DSF','DI','ECF Model','Truth','location','ne')

% for i = 1:out.idxend
%     if (out.g.rho_model(i) == 0)
%         err(i) = 0;
%     else
%         err(i) = 100 * ((out.g.rho_model(i) - out.traj.rho(i))) / out.traj.rho(i);
%     end
% end

% figure(); hold on; grid on
% plot(out.traj.t(1:out.idxend),err,'Linewidth',2)
% % plot(out.traj.t,out.traj.rho,'k--','LineWidth',2)
% xlabel('Time')
% ylabel('rho_{ECF} Error (%)')
% legend('ECF Model','Truth','location','ne')

% gram index 'nominal' MC performance
% figure(); hold on; grid on
% plot(1:1000,out.haf_err,'*')

%}

%% ensemble correlation filter - 'nominal'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = randi([1 100]);
mc.N = 1;
sim.ignore_nominal = true;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
gnc.pred_mode = 1;  % verbose

% out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    %ecf
% gnc.ecf.pert = true;
gnc.ecf.p_tol = 0.02;
gnc.ecf.mode = 1;
% mc.debug = true;
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% calc rss errors of index up to min altitude
atm = load('./data/atm_data/atm_venus_mc.mat');
min_h_ind = find(out.traj.alt == min(out.traj.alt));
alt = atm.nom_table(:,1);
nom = flip(atm.nom_table(:,2));
mc_prof = flip(atm.mc(mc.mcIndex).table(:,2));
% disp = nan(1000,1000);
rss_err_nom = zeros(1000,1);
rss_err_mc = rss_err_nom;
indices = 1000;
for i = 1:1000 %profile
    disp = flip(atm.mc(i).table(:,2));

    % compute rss error to nom
    temp1 = 0;
    temp2 = 0;
    for j = 1:min_h_ind %alt
        err_nom = (disp(j) - nom(j))^2;
        err_mc = (disp(j) - mc_prof(j))^2;
        temp1 = temp1 + err_nom;
        temp2 = temp2 + err_mc;
    end
    rss_err_nom(i) = sqrt(temp1);
    rss_err_mc(i) = sqrt(temp2);
end

% get rss error of ECF gram index each time step
rss = zeros(length(out.traj.t),1);
rss2 = rss;
for i = 1:length(out.traj.t)
    if isnan(out.g.atm.ind_curr(i))
        rss(i) = 0;
        rss2(i) = 0;
    else
        rss(i) = rss_err_nom(out.g.atm.ind_curr(i));
        rss2(i) = rss_err_mc(out.g.atm.ind_curr(i));
    end
    
%     ind_rss(i) = out.g.atm.ind_rss(i) / out.g.atm.rss_hist(i);
end

colors = get_plot_colors();

figure(); hold on
grid on
set(gca,'FontSize',14)
xline(out.tjett(1),'k--','LineWidth',1.5)
yyaxis left
plot(out.traj.t,out.g.atm.ind_curr,'LineWidth',1.5)
yline(mc.mcIndex,'k--','linewidth',1.5)
ylabel('GRAM Index');
yyaxis right
plot(out.traj.t, rss,'LineWidth',1.5)
yline(rss(mc.mcIndex),'--','color',colors(:,1),'LineWidth',1.5);
ylabel('Gram Index Total RSS Error');
% plot(out.traj.t,out.g.ha_err,'LineWidth',1.5)
% ylabel('Apoapsis Error (km)');
% ylim([-500 500])
xlim([0 300])
xlabel('Time (s)')


% figure(); hold on; grid on
% set(gca,'FontSize',14)
% yyaxis left
% plot(out.traj.t, rss,'LineWidth',1.5)
% ylabel('Gram Index Total RSS Error');
% yyaxis right
% plot(out.traj.t, out.g.atm.ind_rss, 'linewidth',2)
% xlabel('Time (s)')
% ylabel('ECF Minimum RSS Error (kg/m^3)')

%}

%% ensemble correlation filter - 'mc'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
% mc.mcIndex = randi(200);
mc.N = 1;
mc.mcIndex = 516;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = true;
gnc = exp_decay_ecrv(gnc);
gnc.pred_mode = 1;  % verbose

% gnc.atm_mode = uint8(1);    %dsf
% out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    %ecf
gnc.ecf.p_tol = 0.05;
gnc.ecf.mode = 0x81;
mc.debug = true;
out = run_dej_n(x0,gnc,aero,sim,mc);

%}

%% DI - MC
%% DI - 'nominal' sim
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 500;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = false;

% no nav errors
% gnc.pred_mode = 1;  % verbose
gnc.atm_mode = uint8(2);    %interpolator
% mc.debug = true;
out_nom = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

% nav errors
gnc = exp_decay_ecrv(gnc);
gnc.atm_mode = uint8(1);    %scale factor
out_K = run_dej_n(x0,gnc,aero,sim,mc);

% nav errors - interpolator
% gnc = exp_decay_ecrv(gnc);
gnc.atm_mode = uint8(2);    %interpolator
% mc.debug = true;
% gnc.pred_mode = 1;  % verbose
out_ecrv = run_dej_n(x0,gnc,aero,sim,mc);

%}

%% DI - MEJ
%{
[x0,aero,gnc,sim, mc] = base_venus_ac_mej(9, true);
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = false;
gnc.npc_mode = uint8(2);   %mej
mc.N = 500;
% mc.Index = 10;
gnc.atm_mode = uint8(2);    %interpolator
gnc = exp_decay_ecrv(gnc);

% sim.parMode = false;

% mc.debug = true;
% gnc.pred_mode = 1;
out_mej = run_dej_n(x0,gnc,aero,sim,mc);

%}


% histogram_plot(50, [out_K.haf_err, out_nom.haf_err, out_ecrv.haf_err, out_mej.haf_err], ... 
%     {'Density Factor, Nav:ECRV','Interpolator, Nav:Perfect', ... 
%     'Interpolator, Nav:ECRV','Interpolator, MEJ, Nav:ECRV'}, 2000)


%% interpolator - 'nominal' sim
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 10;
mc.N = length(mc.mcIndex);
sim.ignore_nominal = true;
mc.sigs = sigs_zero();
gnc.rss_flag = true;

% no nav errors
% gnc.pred_mode = 1;  % verbose
gnc.atm_mode = uint8(2);    %interpolator
% mc.debug = true;
out_nom = run_dej_n(x0,gnc,aero,sim,mc);

gnc = exp_decay_ecrv(gnc);
% mc.debug = true;
% gnc.pred_mode = 1;  % verbose
out_ecrv = run_dej_n(x0,gnc,aero,sim,mc);


figure(); hold on
set(gca,'FontSize',14)
grid on
yline(1, '--', 'Color', [0 0 0] + 0.5,'LineWidth',1.5)
p1=plot(out_ecrv.traj.t,out_ecrv.g.atm.rss_K,'b','LineWidth',2);
p0=plot(out_nom.traj.t,out_nom.g.atm.rss_K,'k--','LineWidth',2);
xlabel('Time (s)')
ylabel('RSS_\rho Error')
title('Density Interpolator RSS Errors vs. Time');
legend([p0 p1],'Perfect Nav','ECRV','location','best');

%}

%% nav errors - simulation
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 10;
mc.N = length(mc.mcIndex);
sim.ignore_nominal = true;
mc.sigs = sigs_zero();
gnc.rss_flag = true;
gnc_nom = gnc;

% no nav errors
% mc.debug = true;
out_nom = run_dej_n(x0,gnc_nom,aero,sim,mc);

% nav errors - exp decay ecrv
gnc = exp_decay_ecrv(gnc);
% mc.debug = true;
% sim.traj_rate = 1;
out_ecrv = run_dej_n(x0,gnc,aero,sim,mc);

keyboard;

figure(); hold on
plot(out_ecrv.traj.t,out_ecrv.g.atm.rss_nom,'k')
plot(out_ecrv.traj.t,out_ecrv.g.atm.rss_K,'b')
plot(out_ecrv.traj.t,out_ecrv.g.atm.rss_ecf,'r')



% figure(); hold on
% set(gca,'FontSize',16)
% plot(out_ecrv.traj.t, out_ecrv.g.K_true,'k','LineWidth',1.5)
% plot(out_ecrv.traj.t, out_ecrv.g.K_dens,'b','LineWidth',1.5)
% xlabel('Time (s)')
% ylabel('Density Scale Factor');
% legend('Truth','ECRV');



nom_err = nan(length(out_nom.traj.t),1);
ecrv_err = nan(length(out_ecrv.traj.t),1);
% exp_err = nan(length(out_exp.traj.t),1);
for i = 1:length(out_nom.traj.t)
    nom_err(i) = 100 * (out_nom.g.K_true(i) - out_nom.g.K_dens(i))/out_nom.g.K_true(i);
end
for i = 1:length(out_ecrv.traj.t)
    ecrv_err(i) = 100 * (out_ecrv.g.K_true(i) - out_ecrv.g.K_dens(i))/out_ecrv.g.K_true(i);
end
% for i = 1:length(out_exp.traj.t)
%     exp_err(i) = 100 * (out_exp.g.K_true(i) - out_exp.g.K_dens(i))/out_exp.g.K_true(i);
% end


figure(); hold on
grid on
title(['Atm index=' num2str(mc.mcIndex) ',ECRV: \tau=' num2str(gnc.nav.tau_ecrv) ...
    ', EXP: \tau =\{' num2str(gnc.nav.tau_r) ',' ...
    num2str(gnc.nav.tau_v) ',' ...
    num2str(gnc.nav.tau_a) '\}']);
set(gca,'FontSize',14)
plot(out_nom.traj.t, nom_err,'k','LineWidth',1.5)
plot(out_ecrv.traj.t, ecrv_err, 'b-.','LineWidth',1.5)
% plot(out_exp.traj.t, exp_err,'r-.','LineWidth',1.5)
xlim([0 120])
xlabel('Time (s)')
ylabel('Density Factor Error (%)');
legend(['Perfect Nav, \Delta r_a = ' num2str(round(out_nom.haf_err,3)) ' km'], ... 
    ['ECRV, \Delta r_a = ' num2str(round(out_ecrv.haf_err,3)) ' km'], ... 
    'location','se')

%}

%% nav error - mc
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 100;
mc.N = length(mc.mcIndex);
sim.ignore_nominal = true;
mc.sigs = sigs_zero();
gnc.rss_flag = true;
gnc_nom = gnc;

% no nav errors
% mc.debug = true;
% out_nom = run_dej_n(x0,gnc_nom,aero,sim,mc);

% nav errors - exp decay ecrv
gnc = exp_decay_ecrv(gnc);
% mc.debug = true;
% sim.traj_rate = 1;
out = run_dej_n(x0,gnc,aero,sim,mc);

for j = 1:mc.N
    for i = 1:length(out_nom.traj.t)
        nom_err(i,j) = 100 * (out_nom.g.K_true(i,j) - out_nom.g.K_dens(i,j))/out_nom.g.K_true(i,j);
    end
    for i = 1:length(out_ecrv.traj.t)
        ecrv_err(i,j) = 100 * (out_ecrv.g.K_true(i,j) - out_ecrv.g.K_dens(i,j))/out_ecrv.g.K_true(i,j);
    end
end

figure(); hold on
set(gca,'FontSize',16)
p0=plot(ecrv_err,'Color',0.5 + [0 0 0]);
p1=plot(nom_err,'Color','b');
[~,hobj,~,~] = legend([p0(1),p1(1)],'ECRV','Perfect Nav','LineWidth',1.5);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
xlabel('Time (s)')
ylabel('Density Scale Factor Error (%)')
%}

%% navigation errors - ecrv
%{
p = nav_msl();
p.rate = 20;
N = 100;
vSig = 3.7e-3 / 3;
rSig = 100 / 3;

% ecrv0 = zeros(9,1);

sx0 = rSig;   % position sigma each axis
sy0 = sx0;
sz0 = sx0;
svx0 = vSig;   % velocity sigma each acis
svy0 = svx0;
svz0 = svx0;
P0 = zeros(9,9);
P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs

ecrv0(1:3) = rSig;
ecrv0(4:6) = 0.1 / 3;%vSig;
ecrv0(7:9) = P0(7,7);
ecrv0 = ecrv0';

% p.P_SS = p.P_SS .^2;
p.P_SS = P0;
tau = [100 200 200 1000];  %ecrv, tau_r, tau_v, tau_a
% figure(1); hold on
fprintf('Simulating Hybrid ECRV Errors\n');
for i = 1:N
nav = simulate_nav(2, 0, 400, 1/p.rate, ...
    tau, p.P_SS, normrnd(0,ecrv0));

% keyboard;

% figure(); hold on
% set(gca,'FontSize',14)
% plot(nav.t, nav.x_ecrv(4,:))

% get vector norms each time step
for j = 1:length(nav.t)
    
%     rmagX(i,j) = nav.rva_err(1,j);
%     rmagY(i,j) = nav.rva_err(2,j);
%     rmagZ(i,j) = nav.rva_err(3,j);
%     
%     vmagX(i,j) = nav.rva_err(4,j);
%     vmagY(i,j) = nav.rva_err(5,j);
%     vmagZ(i,j) = nav.rva_err(6,j);
    
    
    rmagX(i,j) = nav.x_ecrv(1,j);
    rmagY(i,j) = nav.x_ecrv(2,j);
    rmagZ(i,j) = nav.x_ecrv(3,j);
    
    vmagX(i,j) = nav.x_ecrv(4,j);
    vmagY(i,j) = nav.x_ecrv(5,j);
    vmagZ(i,j) = nav.x_ecrv(6,j);
end

end %Ni

% figure(); hold on
% plot(nav.t,nav.rva_err(4:6,:))
% yline(3.7e-3)
% yline(-3.7e-3)
% return

% get mean, std of all cases at each time step
for j = 1:length(nav.t)    
    % r
    muRx(j) = mean(rmagX(:,j));
    stdRx(j) = 3*std(rmagX(:,j));
    muRy(j) = mean(rmagY(:,j));
    stdRy(j) = 3*std(rmagY(:,j));
    muRz(j) = mean(rmagZ(:,j));
    stdRz(j) = 3*std(rmagZ(:,j));
    
    % v
    muVx(j) = mean(vmagX(:,j));
    stdVx(j) = 3*std(vmagX(:,j));
    muVy(j) = mean(vmagY(:,j));
    stdVy(j) = 3*std(vmagY(:,j));
    muVz(j) = mean(vmagZ(:,j));
    stdVz(j) = 3*std(vmagZ(:,j));
    
%         muA(i,j) = norm(nav.x_ecrv(7:9,j));
end


% figure(); hold on
% plot(nav.t,nav.vmag)

% plot r
ind = randi([1,N],1);
figure(1); hold on
subplot(3,1,1); hold on
    plot(nav.t,rmagX,'Color',0.8+[0 0 0]);
    p0=plot(nav.t,rmagX(ind,:),'m','LineWidth',1.2);
    p1=plot(nav.t,stdRx,'b--','LineWidth',1.5);
    plot(nav.t,-stdRx,'b--','LineWidth',1.5);
    p2=yline(mean(stdRx),'r--','LineWidth',1.5);
    yline(-mean(stdRx),'r--','LineWidth',1.5);
    ylabel('v_x noise (m)')
subplot(3,1,2); hold on
    plot(nav.t,rmagY,'Color',0.8+[0 0 0]);
    plot(nav.t,rmagY(ind,:),'m','LineWidth',1.2);
    plot(nav.t,stdRy,'b--','LineWidth',1.5);
    plot(nav.t,-stdRy,'b--','LineWidth',1.5);
    yline(mean(stdRy),'r--','LineWidth',1.5);
    yline(-mean(stdRy),'r--','LineWidth',1.5);
    ylabel('r_y errors (m)')
subplot(3,1,3); hold on
    plot(nav.t,rmagZ,'Color',0.8+[0 0 0]);
    plot(nav.t,rmagZ(ind,:),'m','LineWidth',1.2);
    plot(nav.t,stdRz,'b--','LineWidth',1.5);
    plot(nav.t,-stdRz,'b--','LineWidth',1.5);
    yline(mean(stdRz),'r--','LineWidth',1.5);
    yline(-mean(stdRz),'r--','LineWidth',1.5);
    ylabel('r_z errors (m)')
xlabel('Time (s)','FontSize',14)
legend([p1,p2],'\mu \pm 3\sigma', ... 
    '3\sigma Mean', ... 
    'orientation','horizontal');



% plot v
ind = randi([1,N],1);
figure(2); hold on
subplot(3,1,1); hold on
    plot(nav.t,vmagX,'Color',0.8+[0 0 0]);
%     p0=plot(nav.t,muVx,'b','LineWidth',1.5);
    p0=plot(nav.t,vmagX(ind,:),'m','LineWidth',1.2);
    p1=plot(nav.t,stdVx,'b--','LineWidth',1.5);
    plot(nav.t,-stdVx,'b--','LineWidth',1.5);
    p2=yline(mean(stdVx),'r--','LineWidth',1.5);
    yline(-mean(stdVx),'r--','LineWidth',1.5);
%     p3=yline(3.7e-3,'k','LineWidth',1.5);
%     yline(-3.3e-3,'k','LineWidth',1.5);
    ylabel('v_x noise (m/s)')
subplot(3,1,2); hold on
    plot(nav.t,vmagY,'Color',0.8+[0 0 0]);
%     plot(nav.t,muVy,'b','LineWidth',1.5);
    plot(nav.t,vmagY(ind,:),'m','LineWidth',1.2);
    plot(nav.t,stdVy,'b--','LineWidth',1.5);
    plot(nav.t,-stdVy,'b--','LineWidth',1.5);
    yline(mean(stdVy),'r--','LineWidth',1.5);
    yline(-mean(stdVy),'r--','LineWidth',1.5);
%     yline(3.7e-3,'k','LineWidth',1.5);
%     yline(-3.3e-3,'k','LineWidth',1.5);
%     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
%         '3\sigma Mean','3\sigma Requirement','location','best');
    ylabel('v_y noise (m/s)')
subplot(3,1,3); hold on
    plot(nav.t,vmagZ,'Color',0.8+[0 0 0]);
%     plot(nav.t,muVz,'b','LineWidth',1.5);
    plot(nav.t,vmagZ(ind,:),'m','LineWidth',1.2);
    plot(nav.t,stdVz,'b--','LineWidth',1.5);
    plot(nav.t,-stdVz,'b--','LineWidth',1.5);
    yline(mean(stdVz),'r--','LineWidth',1.5);
    yline(-mean(stdVz),'r--','LineWidth',1.5);
%     yline(3.7e-3,'k','LineWidth',1.5);
%     yline(-3.3e-3,'k','LineWidth',1.5);
%     legend([p0,p1,p2,p3],'\mu','\mu \pm 3\sigma', ... 
%         '3\sigma Mean','3\sigma Requirement','location','best');
    ylabel('v_z noise (m/s)')
xlabel('Time (s)','FontSize',14)
legend([p1,p2],'\mu \pm 3\sigma', ... 
    '3\sigma Mean', ... 
    'orientation','horizontal');


%}


%% navigation errors - exponential decay
%{
% sx0 = 100 / 3;   % position sigma each axis
% sy0 = sx0;
% sz0 = sx0;
% svx0 = 3.7e-3 / 3;   % velocity sigma each acis
% svy0 = svx0;
% svz0 = svx0;
% P0 = zeros(9,9);
% P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
% P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs
% 
% % p = nav_msl();
% % P0 = p.P_SS;
% 
% nav = simulate_nav(3, 0, 500, 1/20, 0, P0, 0);
% 
% subplot(3,1,1)
% plot(nav.t,nav.rErr)
% ylabel('r_{noise} (m)','FontSize',14)
% subplot(3,1,2)
% plot(nav.t,nav.vErr)
% ylabel('v_{noise} (m/s)','FontSize',14);
% subplot(3,1,3)
% plot(nav.t,nav.aErr)
% ylabel('a_{noise} (m/s^2)','FontSize',14)
% xlabel('Time (s)');

%}

