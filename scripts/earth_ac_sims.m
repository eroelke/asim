 %% Evan Roelke
% Earth Aerocapture Simulations
% 

clear;clc
p = '..\earth_ac\dej_n\';


%% Monte Carlo Post Process
%{

% d = load([p '\atm_est\data\mc\' 'mc2k_50.4deg.mat']);
d = load([p '\atm_est\data\mc\mc1000_v10_2k_52deg.mat']);


latex_table(d);


errs = get_default_atm_errs(d);

% capture rate
haTgt = d.gnc.ha_tgt;
tols = 0:1:50;
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


% % ecf accuracy
get_ecf_accuracies(d);

%}


%% Monte Carlo: pert ON/OFF, nav errors
%{
clear;clc
[x0,aero,gnc,sim, mc] = base_earth_ac(true);
mc.N = 1000;
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

x0.fpa0 = -5.2;

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

save(['../earth_ac/dej_n/atm_est/data/mc/mc' ... 
    num2str(gnc.ha_tgt/1000) 'k_' num2str(abs(x0.fpa0*10)) 'deg.mat'])

%}


%% Post Process 'Nominal' Monte Carlo - dispersed GRAM index vs. err
%{
% d = load([p 'atm_est\data/entry_trades/v_10_2k_50.4deg_gramIdxTrade.mat']);
d = load([p 'atm_est\data/entry_trades/v_10_2k_52deg_gramIdxTrade.mat']);

inds = 1:1:1000;

err_k = d.out_k.haf_err;
err_di = d.out_di.haf_err;
err_ecf = d.out_ecf.haf_err;
err_ecfh = d.out_ecfh.haf_err;
err_ecf_p = d.out_ecf_p.haf_err;
err_ecfh_p = d.out_ecfh_p.haf_err;

mu_k = mean(err_k);             std_k = std(err_k);
mu_di = mean(err_di);           std_di = std(err_di);
mu_ecf = mean(err_ecf);         std_ecf = std(err_ecf);
mu_ecfh = mean(err_ecfh);       std_ecfh = std(err_ecfh);
mu_ecf_p = mean(err_ecf_p);     std_ecf_p = std(err_ecf_p);
mu_ecfh_p = mean(err_ecfh_p);   std_ecfh_p = std(err_ecfh_p);

labels = default_atm_labels();
haf_errs = [err_k err_di err_ecf err_ecfh err_ecf_p err_ecfh_p];

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
[x0,aero,gnc,sim, mc] = base_earth_ac(true);
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
gnc.guid_rate = 1;

ha_tgts = 2000;

efpas = [-5.04];

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
        out_ecf_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        err_ecf_p = out_ecf_p.haf_err;
        tjr_ecf_p = out_ecf_p.tjr(:,1);
        dvCirc_ecf_p = out_ecf_p.dv;
        % hybrid
        gnc.atm_mode = uint8(4);  %hyrbid
        gnc.ecf.mode = 1;	% no scaling
        out_ecfh_p = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator        
        err_ecfh_p = out_ecfh_p.haf_err;
        tjr_ecfh_p = out_ecfh_p.tjr(:,1);
        dvCirc_ecfh_p = out_ecfh_p.dv;
        
%         keyboard;
        save(['../earth_ac/dej_n/atm_est/data/entry_trades/perturbed/' ... 
            'v' num2str(round(x0.v0,0)) '_' ...
            num2str(gnc.ha_tgt/1000) 'k_' ... 
            num2str(abs(x0.fpa0*10)) 'deg_gramIdxTrade.mat']);

    end %j - ha tgt
end %i - efpa

keyboard
%}



%%  performance vs. efpa - 'nominal', perturbed
%{
clear;clc
[x0,aero,gnc,sim, mc] = base_earth_ac(true);
mc.N = 1;

low_ind = 20;
high_ind = 600;

inds = [low_ind high_ind];

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

ha_tgts = [2000];
efpas = linspace(-6,-4.8, 40);
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
end %j ,tgt

save(['..\earth_ac\dej_n\atm_est\data\entry_trades\perturbed\' ...
    num2str(gnc.ha_tgt / 1000) 'k_' ...
    'v' num2str(round(x0.v0,0)) ...
    'efpaTrade_GramIdx_' num2str(inds(k)) '.mat']);

end %k

keyboard
%}



%% add simulated variations to density profile
%{
if ~exist('atm','var')
    atm = load('./data/atm_data/atm_earth_gram2016.mat');
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
ylim([nom(end,1)/2000 nom(end,1)/1000]);
xlabel('K_{pert}')
ylabel('Altitude (km)')


figure(); hold on; grid on
title(['Amplitude=' num2str(round(A(ind),2)) ',period=2\pi/' num2str(period(ind)) ... 
    ' km, shift=' num2str(round(shift(ind)/pi,2)) '\pi'])
set(gca,'FontSize',16)
plot(K(:,ind),alt,'k','linewidth',1.5)
plot(Knew(:,ind),alt,'b','linewidth',1.5);
ylim([nom(end,1)/2000 nom(end,1)/1000]);
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
save('./data/atm_data/atm_earth_mc_pert.mat','nom_table','mc');

% ECF data array
atm = load('./data/atm_data/atm_earth_mc_pert.mat');
alt = atm.nom_table(:,1);
nom = atm.nom_table;
mc = atm.mc;
mcs = nan(1000,1001);
mcs(:,1) = alt;
for i = 1:1000 %prof
    mcs(:,i+1) = mc(i).table(:,2);
end

rss_error = zeros(1000,1);
for i = 1:1000 %profile
    disp = atm.mc(i).table(:,2);
    % compute rss error to nom
    rss = 0;
    for j = 1:1000 %altitude
        K(j,i) = atm.mc(i).table(j,2) / nom(j);
        err = (disp(j) - nom(j))^2;
        rss = rss + err;
    end
    rss_err(i) = sqrt(rss);
end

% 
save('./data/atm_data/atm_earth_all_pert.mat','mcs');

% keyboard
%}



%% performance vs. ha tgt - 'nominal'
%{
clear;clc
[x0,aero,gnc,sim, mc] = base_earth_ac(true);
mc.N = 1;

low_ind = 130;
high_ind = 20;

inds = [low_ind high_ind];

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

efpas = [-5.9, -6.2, -6.7];
% efpas = -5.5;

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

save(['../earth_ac/dej_n/atm_est/data/entry_trades/' ...
    num2str(abs(x0.fpa0*10)) 'deg_haTrade_gramIdx_' num2str(inds(k)) '.mat']);

end %k
%}


%% nominal performance vs. haTGt post process
%{
%  20 high, 130 low
% d = load([p 'atm_est\data\entry_trades\67deg_haTrade_gramIdx_20.mat']);
d = load([p 'atm_est\data\entry_trades\67deg_haTrade_gramIdx_130.mat']);

ind = 1;    %-5.9deg
% ind = 2;    %-6.2deg
% ind = 3;    %-6.7deg

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
clear;clc
[x0,aero,gnc,sim, mc] = base_earth_ac(true);
mc.N = 1;

low_ind = 130;
high_ind = 20;

inds = [low_ind high_ind];

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

ha_tgts = [2000 10000];
efpas = linspace(-7,-5.5, 40);
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

save(['../earth_ac/dej_n/atm_est/data/entry_trades/' ...
    'gramIdx_' num2str(inds(k)) '.mat']);

end %k

keyboard
%}


%% nominal performance vs. efpa post process
%{
% d = load([p 'atm_est\data\entry_trades\gramIdx_20.mat']);
% d = load([p 'atm_est\data\entry_trades\gramIdx_130.mat']);

d = load(['..\earth_ac\dej_n\atm_est\data\entry_trades\' ...
    '2k_v10_efpaTrade_GramIdx_600.mat']);

ind = 1;    %2000 km tgt
% ind = 2;  %10000 km tgt

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
% xlim([-5.7, -5.35])
ylim([-200 500])


%}


%% investigate GRAM indices
%{
clear; clc
atm = load('./data/atm_data/atm_earth_gram2016.mat');

alt = atm.nom_table(:,1);
nom = atm.nom_table(:,2);
disp = nan(1000,1000);
rss_err = nan(1000,1);
indices = 1:1:1000;
K = disp;
for i = 1:1000
    disp(:,i) = atm.mc(i).table(:,2);
    % compute rss error to nom
    rss = 0;
    for j = 1:1000
        K(j,i) = atm.mc(i).table(j,2) / nom(j);
        err = (disp(j,i) - nom(j))^2;
        rss = rss + err;
    end
    rss_err(i) = sqrt(rss);
end

min_ind = find(rss_err == min(rss_err));
max_ind = find(rss_err == max(rss_err));


figure(); hold on; grid on
title('EarthGRAM RSS Errors')
set(gca,'FontSize',14)
%histogram(rss_err,'NumBins',50)
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


%% Visualize lift/drag trajectory modifications
%{
clear; clc;
cd = 1;

LDs = [0 0.05 0.1 0.15 0.2];
betas = [50;100;125;250;500;];


[x0,aero,gnc,sim,mc] = base_earth_ac(true);
x0.v0 = 13;
x0.fpa0 = -5.6;
mc.flag = false;
x0.az0 = 0;
aero.rcs = 1.5;
aero.cds = cd;
aero.cls = 0;
sim.efpa_flag = false;
mc.debug = false;
gnc.atm_mode = uint8(1);    %exponential atm
gnc.guid_rate = 0;
gnc.mode = 0;   %dej
gnc.n = 0;
gnc.ha_tgt = 42164;
sim.traj_rate = 100;
sim.data_rate = 10;
sim.atm_mode = 1;
sim.t_max = 500;
gnc.tj0 = 125;
gnc.dtj_lim = 10;
sim.h_max = 125;
sim.h_min = 0;
sim.debug = false;

% drag
% labels = cell(length(betas),1);
% for i = 1:length(betas)
%     aero.m = betas(i) * aero.cds(1) * pi * aero.rcs(1)^2;
%     out(i) = run_dej_n(x0,gnc,aero,sim,mc);
%     
%     labels{i} = ['\beta = ' num2str(betas(i)) ' kg/m^2'];
% end
% 
% 
% figure(); hold on; grid on;
% title(['v0=' num2str(x0.v0) ', fpa=' num2str(x0.fpa0)]);
% set(gca,'FontSize',16);
% for i = 1:length(betas)
%     plot(out(i).traj.vel_pp_mag./1000,out(i).traj.alt./1000,'LineWidth',2)
% end
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')
% legend(labels,'location','nw');


% lift
Ni = length(LDs);
labels = cell(Ni,1);
for i = 1:Ni
    aero.cls = aero.cds * LDs(i);
    out(i) = run_dej_n(x0,gnc,aero,sim,mc);
    labels{i} = ['L/D = ' num2str(LDs(i))];
end
figure(); hold on; grid on;
title(['v0=' num2str(x0.v0) ', fpa=' num2str(x0.fpa0)]);
set(gca,'FontSize',16);
for i = 1:Ni
    plot(out(i).traj.vel_pp_mag./1000,out(i).traj.alt./1000,'LineWidth',2)
end
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')
legend(labels,'location','nw');

% keyboard
%}


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
[nom, mc] = load_atm_data(0, 'earth');
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



