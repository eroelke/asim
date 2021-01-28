%% venus_ac_atm.m
%   atmospheric estimation code/sims for venus aerocapture
% 
clear;clc;

%% Uniqueness of GRAM
atm = load('./data/atm_data/atm_venus_mc.mat');


%% DI/ECF hybrid - mc
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.N = 500;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = true;
% gnc.pred_mode = 1;  % verbose
gnc = exp_decay_ecrv(gnc);

gnc.atm_mode = uint8(1);    %dsf
out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(2);    %di
out_di = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(4);    %hybrid
gnc.p_tol = 0.1;
out_ecf = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator


save('../venus_ac/dej_n/atmospheric_estimation/data/hybrid/mc_test.mat');

for i = 1:mc.N
    if ~isnan(out_ecf.idj(i,1))
        ind_j(i) = out_ecf.g.atm.ind_curr(out_ecf.idj(i,1),i);
    else
        ind_j(i) = out_ecf.g.atm.ind_curr(out_ecf.idxend(i),i);
    end
end

figure(); hold on; grid on
histogram(ind_j, 'NumBins',200)
xline(100);

figure(); hold on; grid on
set(gca,'FontSize',14)
plot(ind_j, out_ecf.haf_err,'+','LineWidth',1.5)
xlabel('GRAM Index at Jettison')
ylabel('Apoapsis Error (km)')

figure(); hold on; grid on
set(gca,'FontSize',14)
plot(out_ecf.traj.t, out_ecf.g.atm.ind_curr,'LineWidth',1.2)
p1=plot(out_ecf.tjett(:,1), ind_j,'k*')
yline(100);
xlabel('Time (s)')
ylabel('GRAM Index')
legend([p1],'Jettison')
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
gnc.p_tol = 0.02;

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
gnc.p_tol = 0.02;
% mc.debug = true;

for i = 1:1000
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator
err(i) = out.haf_err;
t(:,i) = out.traj.t;
inds(:,i) = out.g.atm.ind_curr;
end

N = i;

% save('../venus_ac/dej_n/atmospheric_estimation/data/hybrid/gram_idx100_mc.mat')
for i = 1:N
%     idx0 = find(~isnan(inds(:,i)),1);
    ind_f(i) = inds(find(isnan(t(:,i)),1)-1,i);
end

% figure(); hold on; grid on
% set(gca,'FontSize',14)
% histogram(ind_f,'NumBins',100);
% xline(mc.mcIndex,'k--','Color',[0 0 0] + 0.5)

% figure(); hold on; grid on
% set(gca,'FontSize',14)
% plot(ind_f, err,'*')
% yline(0)
% xlabel('Final Ecf Index')
% ylabel('Apoapsis Error (km)')
% xline(mc.mcIndex,'k--','Color',[0 0 0] + 0.5)

%}

%% DI/ECF hybrid - nominal
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = 100;
mc.N = 1;
sim.ignore_nominal = true;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% gnc.pred_mode = 1;  % verbose
gnc = exp_decay_ecrv(gnc);


% out_k = run_dej_n(x0,gnc,aero,sim,mc);
gnc.atm_mode = uint8(4);    %hybrid
gnc.p_tol = 0.02;
mc.debug = true;
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator
%}

%% ensemble correlation filter - 'nominal'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
mc.mcIndex = randi(200) + 100;
mc.N = 1;
sim.ignore_nominal = true;
mc.sigs = zero_sigs();
gnc.rss_flag = true;
% gnc.pred_mode = 1;  % verbose

out_k = run_dej_n(x0,gnc,aero,sim,mc);
gnc.atm_mode = uint8(3);    %ecf
gnc.p_tol = 0.05;
% mc.debug = true;
out = run_dej_n(x0,gnc,aero,sim,mc);    %perfect nav, interpolator

figure(); hold on
grid on
set(gca,'FontSize',14)
yline(mc.mcIndex,'k--','LineWidth',1.5);
xline(out.tjett(1),'k--','LineWidth',1.5)
yyaxis left
plot(out.traj.t,out.g.atm.ind_curr - 1,'LineWidth',1.5)
ylabel('GRAM Index');
yyaxis right
plot(out.traj.t,out.g.ha_err,'LineWidth',1.5)
ylabel('Apoapsis Error (km)');
ylim([-500 500])
xlim([0 300])
xlabel('Time (s)')
%}

%% ensemble correlation filter - 'mc'
%{
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
% mc.mcIndex = randi(200);
mc.N = 500;
sim.ignore_nominal = true;
mc.sigs = default_sigs();
gnc.rss_flag = true;
gnc = exp_decay_ecrv(gnc);
% gnc.pred_mode = 1;  % verbose

gnc.atm_mode = uint8(1);    %dsf
out_k = run_dej_n(x0,gnc,aero,sim,mc);

gnc.atm_mode = uint8(3);    %ecf
gnc.p_tol = 0.05;
% mc.debug = true;
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
mc.mcIndex = 1:1000;
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
p.rate = 100;
N = 100;
sigReq = 3.7e-3;

% ecrv0 = zeros(9,1);

sx0 = 100 / 3;   % position sigma each axis
sy0 = sx0;
sz0 = sx0;
svx0 = sigReq / 3;   % velocity sigma each acis
svy0 = svx0;
svz0 = svx0;
P0 = zeros(9,9);
P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs

% p.P_SS = p.P_SS .^2;
p.P_SS = P0;
tau = [100 200 200 1000];  %ecrv, tau_r, tau_v, tau_a
% figure(1); hold on
for i = 1:N
ecrv0 = [normrnd(0,100,[3,1]); ...      %r_atm uncertainty
    normrnd(0,.1/3,[3,1]); ...   %v_atm uncertainty
    normrnd(0,P0(7,7),[3,1])];  

nav = simulate_nav(2, 0, 400, 1/p.rate, ...
    tau, p.P_SS, ecrv0);

nav.stdV = nav.stdV ./ (sigReq / 3);


% keyboard;

% figure(); hold on
% set(gca,'FontSize',14)
% plot(nav.t, nav.x_ecrv(4,:))

% get vector norms each time step
for j = 1:length(nav.t)
    
    vmagX(i,j) = nav.x_ecrv(4,j);
    vmagY(i,j) = nav.x_ecrv(5,j);
    vmagZ(i,j) = nav.x_ecrv(6,j);    
%     amag(i,j) = mean(nav.x_ecrv(7:9,j));
end

end %Ni

% figure(); hold on
% plot(nav.t,nav.rva_err(4:6,:))
% yline(3.7e-3)
% yline(-3.7e-3)
% return

% get mean, std of all cases at each time step
for j = 1:length(nav.t)    
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

ind = randi([1,N],1);
figure(1); hold on
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


function gnc = exp_decay_ecrv(gnc)
gnc.nav = nav_msl();
gnc.nav.mode = 5;
gnc.nav.rate = 100;  %hz
gnc.nav.tau_ecrv = 100;
gnc.nav.tau_r = 400;
gnc.nav.tau_v = 400;
gnc.nav.tau_a = 1000;
sx0 = 100 / 3;   % position sigma each axis
sy0 = sx0;
sz0 = sx0;
svx0 = 3.7e-3 / 3;   % velocity sigma each acis
svy0 = svx0;
svz0 = svx0;
% P0 = zeros(9,9);
P0 = diag([sx0^2, sy0^2, sz0^2, svx0^2, svy0^2, svz0^2]);   %covariance matrix
P0(7:9,7:9) = eye(3)*(9.81*.05*1e-6 / 3)^2;  %.05 micro-gs
gnc.nav.mode = 5;
gnc.nav.P_SS = P0;
gnc.nav.ercv0 = [normrnd(0,sx0,[3,1]); ...      %r_atm uncertainty
    normrnd(0,svx0,[3,1]); ...   %v_atm uncertainty
    normrnd(0,P0(7,7),[3,1])]; 
end
