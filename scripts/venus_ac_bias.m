% % venus_ac_biasing.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 2-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nominal Visualization
%{
clear;clc
% h vs v
[x0,aero,gnc,sim, mc] = base_venus_ac(false);

% reminder: set init_iters to 1!

tj0 = 110;
biasing = -500e3;
b21 = 3;
b31 = 10;
gnc.ha_tgt = 2000;
x0.v0 = 11;
sim.traj_rate = 100;

% gnc.n = 2;
m = [150 60 20];
rcs(1) = 1;
rcs(2) = sqrt( (m(2) * rcs(1)^2) / ...
    (m(1) * b21) );
rcs(3) = sqrt( (m(3) * rcs(1)^2) / ...
    (m(1) * b31) );
cd = 1.05;
b = nan(3,1);   %betas
br = nan(2,1);  %beta ratios
for i = 1:3
    b(i) = m(i) / (cd * pi * rcs(i)^2);
end
for i = 1:2
    br(i) = b(i+1)/b(i);
end
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br31 = b(3)/b(1);
aero.rn = 0.1;
aero.cls = 0;

sim.traj_rate = 100;
sim.efpa_flag = true;
mc.flag = false;
gnc.guid_rate = 0;

% beta 1 traj
% aero.m = m(1);
% aero.cds = cd;
% aero.rcs = [rcs(1)];
% gnc.n = 0;
% out1 = run_dej_n(x0,gnc,aero,sim,mc);
% fpa1 = out1.traj.gamma_pp(1) * 180/pi;

% beta 2 traj
% aero.m = m(2);
% aero.cds = cd;
% aero.rcs = [rcs(2)];
% gnc.n = 0;
% out2 = run_dej_n(x0,gnc,aero,sim,mc);
% fpa2 = out2.traj.gamma_pp(1) * 180/pi;

% beta 3 traj
% aero.m = m(3);
% aero.cds = cd;
% aero.rcs = [rcs(3)];
% gnc.n = 0;
% x0.fpa0 = fpa2;
% out3 = run_dej_n(x0,gnc,aero,sim,mc);
% fpa3 = out3.traj.gamma_pp(1) * 180/pi;

% n3 biased control


% n2 biased control
aero.m = [m(1) m(2) m(3)];
aero.cds = [cd cd cd];
aero.rcs = [rcs(1) rcs(2) rcs(3)];
aero.cls = [0 0 0];
gnc.n = 2;
gnc.guid_rate = 0.1;
gnc.npc_mode = uint8(2);
gnc.bias(1).tgt_ap = biasing;
gnc.tj0 = [tj0 160];
x0.fpa0 = -5.4;
sim.efpa_flag = false;

% sim.debug = true;
out2 = run_dej_n(x0,gnc,aero,sim,mc);

%single jettison bo bias
aero.m = [m(1) m(2)];
aero.cds = [cd cd];
aero.rcs = [rcs(1) rcs(2)];
aero.cls = [0 0];
gnc.n = 1;
gnc.guid_rate = 0.1;
gnc.npc_mode = uint8(1);
gnc.bias(1).tgt_ap = 0;
gnc.tj0 = tj0;
x0.fpa0 = -5.4;
sim.efpa_flag = false;
outSEJ = run_dej_n(x0,gnc,aero,sim,mc);

%n3 control
[aero, bijs] = get_n3_betas(3, 5, 10);
gnc.n = 3;
gnc.guid_rate = 0.1;
gnc.npc_mode = uint8(2);    %multi
gnc.bias(1).tgt_ap = biasing;
gnc.bias(2).tgt_ap = -100e3;
gnc.tj0 = [tj0 160 200];
x0.fpa0 = -5.4;
sim.efpa_flag = false;

out3 = run_dej_n(x0,gnc,aero,sim,mc);

colors = get_plot_colors();
% time vs. guidance apoapsis estimate
figure(); hold on
grid on;
set(gca,'FontSize',16)
title('-500km, -100 km bias, 2k target, v=11')
set(gca,'FontSize',16)
p3=plot(out3.traj.time, out3.g.dej_n.dr_ap_true./1000,'Color',colors(1,:),'LineWidth',2);
p2=plot(out2.traj.time, out2.g.dej_n.dr_ap_true./1000,'Color',colors(2,:),'LineWidth',2);
p1=plot(outSEJ.traj.time, outSEJ.g.dej_n.dr_ap./1000,'Color',colors(3,:),'LineWidth',2);
plot(outSEJ.t_jett(1), outSEJ.g.dej_n.dr_ap(outSEJ.idj(1))./1000, ... 
    'x','Color',colors(3,:),'LineWidth',1.5);
plot(out2.t_jett(1), out2.g.dej_n.dr_ap_true(out2.idj(1))./1000, ...
    'x','Color',colors(2,:),'LineWidth',1.5);
plot(out2.t_jett(2), out2.g.dej_n.dr_ap_true(out2.idj(2))./1000, ...
    'o','Color',colors(2,:),'LineWidth',1.5);
p3_1=plot(out3.t_jett(1), out3.g.dej_n.dr_ap_true(out3.idj(1))./1000, ...
    'x','Color',colors(1,:),'LineWidth',1.5);
p3_2=plot(out3.t_jett(2), out3.g.dej_n.dr_ap_true(out3.idj(2))./1000, ...
    'o','Color',colors(1,:),'LineWidth',1.5);
p3_3=plot(out3.t_jett(3), out3.g.dej_n.dr_ap_true(out3.idj(3))./1000, ...
    '^','Color',colors(1,:),'LineWidth',1.5);
% p3_1 = xline(out3.t_jett(1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
% p3_2 = xline(out3.t_jett(2),'Color',colors(1,:),'LineStyle','-.','LineWidth',1.5);
% p3_3 = xline(out3.t_jett(3),'Color',colors(1,:),'LineStyle',':','LineWidth',1.5);
% p2_1 = xline(out2.t_jett(1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
% p2_2 = xline(out2.t_jett(2),'Color',colors(2,:),'LineStyle','-.','LineWidth',1.5);
% p1_1 = xline(
yline(-100,'Color',[0 0 0] + 0.5,'LineWidth',1.5);
yline(-500,'Color',[0 0 0] + 0.5,'LineWidth',1.5);
legend([p3, p2, p1 p3_1, p3_2, p3_3], ... 
    '$\beta_{i+1}/\beta_i = \{3, 5, 10\}$', ... 
    '$\beta_{i+1}/\beta_i = \{3, 10\}$', ...
    '$\beta_2/\beta_1 = 3 (0\% Bias)$', ... 
    'Jettison 1','Jettison 2','Jettison 3', ...
    'Interpreter','latex', ...
    'location','ne')
ylabel('Apoapsis Error Estimate (km)')
xlabel('Time (s)')

keyboard;

% time vs. jettison time difference
figure(); hold on
grid on;
set(gca,'FontSize',16)
title('-500km, -100 km bias, 2k target, v=11')
set(gca,'FontSize',16)
p3=plot(out3.traj.time, out3.g.dej_n.tj ./out3.traj.time(out3.idxend),'Color',colors(1,:),'LineWidth',2);
p2=plot(out2.traj.time, out2.g.dej_n.tj./out2.traj.time(out2.idxend),'Color',colors(2,:),'LineWidth',2);
p1=plot(outSEJ.traj.time, outSEJ.g.dej_n.tj./outSEJ.traj.time(outSEJ.idxend),'Color',colors(3,:),'LineWidth',2);
% p3_1 = xline(out3.t_jett(1),'Color',colors(1,:),'LineStyle','--','LineWidth',1.5);
% p3_2 = xline(out3.t_jett(2),'Color',colors(1,:),'LineStyle','-.','LineWidth',1.5);
% p3_3 = xline(out3.t_jett(3),'Color',colors(1,:),'LineStyle',':','LineWidth',1.5);
% p2_1 = xline(out2.t_jett(1),'Color',colors(2,:),'LineStyle','--','LineWidth',1.5);
% p2_2 = xline(out2.t_jett(2),'Color',colors(2,:),'LineStyle','-.','LineWidth',1.5);
% p1_1 = xline(outSEJ.t_jett(1),'Color',colors(3,:),'LineStyle','--','LineWidth',1.5);
plot(outSEJ.t_jett(1), outSEJ.t_jett(1)/outSEJ.traj.time(outSEJ.idxend), ... 
    'x','Color',colors(3,:),'LineWidth',1.5);
plot(out2.t_jett(1), out2.t_jett(1)/out2.traj.time(out2.idxend), ...
    'x','Color',colors(2,:),'LineWidth',1.5);
plot(out2.t_jett(2), out2.t_jett(2)/out2.traj.time(out2.idxend), ...
    'o','Color',colors(2,:),'LineWidth',1.5);
p3_1=plot(out3.t_jett(1), out3.t_jett(1)/out3.traj.time(out3.idxend), ...
    'x','Color',colors(1,:),'LineWidth',1.5);
p3_2=plot(out3.t_jett(2), out3.t_jett(2)/out3.traj.time(out3.idxend), ...
    'o','Color',colors(1,:),'LineWidth',1.5);
p3_3=plot(out3.t_jett(3), out3.t_jett(3)/out3.traj.time(out3.idxend), ...
    '^','Color',colors(1,:),'LineWidth',1.5);
legend([p3, p2, p1, p3_1, p3_2, p3_3], ... 
    '$\beta_{i+1}/\beta_i = \{3, 5, 10\}$', ... 
    '$\beta_{i+1}/\beta_i = \{3, 10\}$', ...
    '$\beta_2/\beta_1 = 3$', ... 
    'Jettison 1', ...
    'Jettison 2', ...
    'Jettison 3', ...
    'Interpreter','latex', ...
    'location','ne')
ylabel('t_j^{est}/t_f')
xlabel('Time (s)')

keyboard;

% % alt vs vel
% figure(); hold on
% title('-500km bias, 2k target, v=11')
% set(gca,'FontSize',16)
% plot(out1.traj.vel_pp_mag/1000, out1.traj.alt/1000,'k','LineWidth',2)
% plot(out2.traj.vel_pp_mag/1000, out2.traj.alt/1000,'k--','LineWidth',2)
% plot(out3.traj.vel_pp_mag/1000, out3.traj.alt/1000,'k-.','LineWidth',2)
% plot(outMEJ.traj.vel_pp_mag/1000, outMEJ.traj.alt/1000,'b', 'LineWidth',2)
% plot(outSEJ.traj.vel_pp_mag/1000, outSEJ.traj.alt/1000, 'r-.','LineWidth',2)
% plot(outMEJ.traj.vel_pp_mag(outMEJ.idj(1))/1000, ...
%     outMEJ.traj.alt(outMEJ.idj(1))/1000, 'bo','LineWidth',2)
% plot(outMEJ.traj.vel_pp_mag(outMEJ.idj(2))/1000, ...
%     outMEJ.traj.alt(outMEJ.idj(2))/1000, 'bo','LineWidth',2)
% plot(outSEJ.traj.vel_pp_mag(outSEJ.idj(1))/1000, ...
%     outSEJ.traj.alt(outSEJ.idj(1))/1000, 'ro','LineWidth',2)
% legend(['\beta_1, EFPA = ', num2str(fpa1)], ...
%     ['\beta_2 = ' num2str(outMEJ.beta_ratios(1)) '\beta_1, EFPA = ', num2str(fpa2)], ...
%     ['\beta_3 = ' num2str(outMEJ.beta_ratios(2) * outMEJ.beta_ratios(1)) '\beta_1, EFPA = ', num2str(fpa3)], ...
%     ['\beta_{i+1}/\beta_i = \{' num2str(outMEJ.beta_ratios(1)) ... 
%      ', ' num2str(outMEJ.beta_ratios(2)) '\}, EFPA = ', num2str(x0.fpa0)], ...
%     ['SEJ, \beta_2/\beta_1 = ' num2str(outMEJ.beta_ratios(1)) ', EFPA = ', num2str(x0.fpa0)], ...
%     'location','best');
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')

keyboard;

%}

%% Bias Value Trade
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [400 2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 3;
v0 = 11;

for z = 1:Nz
biasing = ha_tgts(z) .* linspace(0.001, 0.5, 200) .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
% biasing = ha_tgts(z) * 0.25 * -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac(false);

x0.v0 = v0;
x0.fpa0 = -5.6;

gnc.ha_tgt = ha_tgts(z);
gnc.ha_tol = 0.1;
gnc.n = 2;
tj0s = [100 170];
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
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;
    
mc.debug = false;
mc.flag = false;

sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage bias, b2/b1 = %2.0f, ha* = %4.0f. ',br21, gnc.ha_tgt); print_current_time();
gnc.force_jett = true;
gnc.comp_curr = false;

% sim.nWorkers = 5;
gnc.tj0 = tj0s;
tjr = nan(length(biasing),gnc.n);
tj = tjr;
err = nan(length(biasing),1);

for i = 1:length(biasing)

gnc.bias(1).tgt_ap = biasing(i);
out = run_dej_n(x0,gnc,aero,sim,mc);

tj(i,:) = out.t_jett';
tjr(i,:) = out.idj/out.idxend;
err(i) = out.haf_err;

end %biasing
% betas = ['b_ijs' num2str(b21) '_' num2str(b31)];

if (ha_tgts(z) > 1000)
    tgt = [num2str(floor(ha_tgts(z)/1000)) 'k_'];
else
    tgt = [num2str(ha_tgts(z)) '_'];
end
% keyboard;
save(['..\venus_ac\dej_n\HYDRA\2stage_bias\bias_val\' tgt ... 
    'v' num2str(v0) '_' num2str(abs(x0.fpa0)) 'deg_b21_' num2str(b21s(j)) '.mat'])
fprintf('Finished 2stage bias sim. '); print_current_time();

end %betas

end %ha_tgts
keyboard
%}

%% Bias Value Post Processing
%{
% plot_tj_vs_bias(2, 400, 5.4, 3, 5);
% plot_tj_vs_bias(2, 400, 5.4, 9, 5);

% plot_tj_vs_bias(2, 400, 5.6, 3, 5);
% plot_tj_vs_bias(2, 400, 5.6, 3, 5);

% plot_tj_vs_bias(2, 2000, 5.4, 3, 5);
% plot_tj_vs_bias(2, 2000, 5.4, 9, 5);

% plot_tj_vs_bias(2, 2000, 5.6, 3, 5);
% plot_tj_vs_bias(2, 2000, 5.6, 3, 5);

%}

%% Biasing EFPA Trade Studies
%{
clearvars -except save_path ;clc

Nfpa = 50;

% ha_tgts = [400 2000 10000];
ha_tgts = [400 2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 2;
v0 = 11;

for z = 1:Nz
biasing = ha_tgts(z) .* [0.1 0.2 0.3 0.4 0.5] .* -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac(false);
x0.v0 = v0;
% efpas = [-5.35:-0.0025:-5.4, ...
%     -5.4:-0.01:-5.52, ...
%     -5.52:-0.0025:-5.6, ... 
%     -5.6:-0.05:-6];

% efpas = [-5.3:-0.0025:-5.46, ...
%     -5.46:-0.02:-5.8];
% efpas = -5.4;

efpas =linspace(-5.35 ,-5.85, Nfpa);
% efpas = -5.4;
Nk = length(efpas);

gnc.ha_tgt = ha_tgts(z);
gnc.ha_tol = 0.5;
gnc.n = 2;
tj0s = [100 170];
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

sim.traj_rate = 100;
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;
    
mc.debug = false;
mc.flag = false;

sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage bias, b2/b1 = %2.0f, ha* = %4.0f. ',br21, gnc.ha_tgt); print_current_time();
gnc.force_jett = true;
gnc.comp_curr = false;

% sim.nWorkers = 5;
gnc.tj0 = tj0s;
tjr = nan(length(biasing),Nk,gnc.n);
tj = tjr;
err = nan(length(biasing),Nk);

for k = 1:Nk
x0.fpa0 = efpas(k);

for i = 1:length(biasing)

gnc.bias(1).tgt_ap = biasing(i);
out = run_dej_n(x0,gnc,aero,sim,mc);

tj(i,k,:) = out.t_jett';
tjr(i,k,:) = out.idj/(out.idxend / sim.data_rate);
err(i,k) = out.haf_err;

end %biasing
% betas = ['b_ijs' num2str(b21) '_' num2str(b31)];
end %fpas

if (ha_tgts(z) > 1000)
    tgt = [num2str(floor(ha_tgts(z)/1000)) 'k_'];
else
    tgt = [num2str(ha_tgts(z)) '_'];
end
% keyboard;
save(['..\venus_ac\dej_n\HYDRA\2stage_bias\entry_trades\' tgt ... 
    'v' num2str(v0) '_b21_' num2str(b21s(j)) '.mat'])
fprintf('Finished 2stage bias sim. '); print_current_time();

end %betas

end %ha_tgts
keyboard

% est err vs time
% figure(); hold on
% grid on
% title('2k, v11, b21 = 3')
% set(gca,'FontSize',16);
% plot(out.traj.time,out.g.dej_n.dr_ap_true./1000,'b','LineWidth',2);
% p1=xline(out.traj.time(out.idj(1)),'k--','LineWidth',1.5);
% p2=xline(out.traj.time(out.idj(2)),'k-.','LineWidth',1.5);
% % yline(-100, '-','Color',[0,0,0]+0.5,'LineWidth',1.5)
% yline(-500, '-','Color',[0,0,0]+0.5,'LineWidth',1.5)
% ylim([-500 200]);
% xlim([0 out.traj.time(out.idxend)])
% xlabel('Time (s)')
% ylabel('Estimated Apoapsis Error (km)')
% legend([p1 p2],'Jettison 1','Jettison 2','location','ne')


%}

%% 2stage EFPA trade study post-processing
%{
clear;clc

tol = 5;
plot_tj_opt(2, 400, 3, tol, true);
plot_tj_opt(2, 400, 9, tol, true);
plot_tj_opt(2, 2000, 3, tol, true);
plot_tj_opt(2, 2000, 9, tol, true);



%}

%% Regular 2stage Monte Carlo
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [400];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 2;
v0 = 11;

for z = 1:Nz
% biasing = ha_tgts(z) .* [0.1 0.2 0.3 0.4 0.5] .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
biasing = ha_tgts(z) * [0.1 0.3] * -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac(true);
x0.v0 = v0;
efpas = [-5.4 -5.6];
Nk = length(efpas);

gnc.ha_tgt = ha_tgts(z);
gnc.ha_tol = 10;
gnc.n = 2;
tj0s = [100 170];
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

sim.traj_rate = 100;
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;
    
mc.debug = false;
mc.flag = true;
mc.N = 1000;

sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage bias MC, b2/b1 = %2.0f, ha* = %4.0f. ',br21, gnc.ha_tgt); print_current_time();
gnc.force_jett = true;
gnc.comp_curr = false;

% sim.nWorkers = 5;
gnc.tj0 = tj0s;
tjr = nan(length(biasing),Nk,gnc.n);
tj = tjr;
err = nan(length(biasing),Nk);

for k = 1:Nk
x0.fpa0 = efpas(k);

for i = 1:length(biasing)

gnc.bias(1).tgt_ap = biasing(i);
out = run_dej_n(x0,gnc,aero,sim,mc);

if (ha_tgts(z) > 1000)
    tgt = [num2str(floor(ha_tgts(z)/1000)) 'k_'];
else
    tgt = [num2str(ha_tgts(z)) '_'];
end
% keyboard;
save(['..\venus_ac\dej_n\HYDRA\2stage_bias\data\' tgt ... 
    'v' num2str(v0) '_b21=' num2str(b21s(j)) '_bias=' ... 
    num2str(abs(biasing(z)/1000)) '_set2.mat'])
fprintf('Finished 2stage bias MC. '); print_current_time();

end %biasing
% betas = ['b_ijs' num2str(b21) '_' num2str(b31)];
end %fpas

end %betas

end %ha_tgts
keyboard

%}

%% 2stage Monte Carlo as func of biasing percent
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [400];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 3;
v0 = 11;

for z = 1:Nz
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
x0.v0 = v0;
efpas = [-5.4];
Nk = length(efpas);

gnc.ha_tgt = ha_tgts(z);
gnc.ha_tol = 10;
gnc.n = 2;
tj0s = [100 170];
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);    %multi event jettison

for j = 1:length(b21s)
% biasing = ha_tgts(z) * ([0 1 5 10 15 20] ./ 100) * -1000;
% biasing = ha_tgts(z) * ([25 30 35 40] ./ 100) * -1000;
% biasing = ha_tgts(z) * ([0 1 2 3 4 5] ./ 100) * - 1000;

if (b21s(j) == 3)
    biasPercents = [0 1 5 10 15 20 25 30 35 40]';
elseif (b21s(j) == 9)
    biasPercents = [0 1 2 3 4 5 10 15 20]';
end
biasing = ha_tgts(z) * (biasPercents ./ 100) * - 1000;

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
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;
sim.ignore_nominal = true;
    
mc.debug = false;
mc.flag = true;
mc.N = 2500;

sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage bias val MC, b2/b1 = %2.0f, ha* = %4.0f. ',br21, gnc.ha_tgt); print_current_time();
gnc.force_jett = true;
gnc.comp_curr = false;

% sim.nWorkers = 5;
gnc.tj0 = tj0s;
err = nan(mc.N, length(biasing));
tj1 = err; tj2 = err;
tjr1 = err; tjr2 = err;

for k = 1:Nk
x0.fpa0 = efpas(k);

for i = 1:length(biasing)
gnc.bias(1).tgt_ap = biasing(i);
out = run_dej_n(x0,gnc,aero,sim,mc);
% fprintf('%g% done\n',100 * i / length(biasing));
err(:,i) = out.haf_err;
tj1(:,i) = out.tjett(:,1);
tj2(:,i) = out.tjett(:,2);
tjr1(:,i) = out.tjr(:,1);
tjr2(:,i) = out.tjr(:,2);

end %biasing
if (ha_tgts(z) > 1000)
    tgt = [num2str(floor(ha_tgts(z)/1000)) 'k_'];
else
    tgt = [num2str(ha_tgts(z)) '_'];
end

save(['..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\' tgt ... 
    'v' num2str(v0) '_' num2str(abs(x0.fpa0)) 'deg_b21=' num2str(b21s(j)) '_full.mat'])
fprintf('Finished 2stage bias val MC. '); print_current_time();

% betas = ['b_ijs' num2str(b21) '_' num2str(b31)];
end %fpas

end %betas

end %ha_tgts
keyboard
%}

%% Biasing Monte Carlo mean/std trends
%{
clear;clc
% kk -> biasing1
% zz -> biasing2

% b21 = 3;
% tgt = 400;
% if ~exist('bias','var')
%     [bias,err] = load_bias_val_data(tgt, b21, 2);
% end

b21 = 9;
tgt = 400;
p = '..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\';
% d = load([p 'mc2500_400_v11_5.4deg_b21=3.mat']);
% d2 = load([p '400_v11_5.4deg_b21=3_set2.mat']);
d = load([p '400_v11_5.4deg_b21=9_full.mat']);

fj = load(['..\venus_ac\dej_n\HYDRA\2stage\' ... 
    'data\400\v11_5.4deg_b_ijs' num2str(b21) '_10.mat']);
fj_ind = 'out';

bias1 = 100 * (([d.biasing] ./ 1000) ./ d.gnc.ha_tgt);
if exist('d2','var')
    bias2 = 100 * (([d2.biasing] ./ 1000) ./ d2.gnc.ha_tgt);
% bias = [bias2 bias1(4:end)]';
    bias = [bias1 bias2]';
else
    bias = bias1';
end
Ni = length(bias);

tols = 0:0.1:20;
Nj = length(tols);
pc = nan(Ni, Nj);
j1s = pc; j2s = pc;
pcFJ = nan(Nj,1);

meanHa = nan(Ni,1);
stdHa = meanHa;
entries = cell(2, 1);

for i = 1:Ni
%     for j = 1:Nj
%         pc(i,j) = get_percent_captured( ... 
%             d.err(:,i),d.gnc.ha_tgt, tols(j), d.mc.N);
% 
%         [j1s(i,j), j2s(i,j), ~] = ... 
%             count_jettisons(d.tj1(:,i), d.tj2(:,i));            
%     end
%     if (i <= length(bias1))
%         meanHa(i) = mean(d.err(:,i));
%         stdHa(i) = std(d.err(:,i));
%     else
%         ind = i - length(bias1);
%         for j = 1:size(d2.err,1)
%             if (abs(d2.err(j,ind) > 3200))
%                 d2.err(j,ind) = nan;
%             end  
%         end
%         
%         meanHa(i) = nanmean(d2.err(:,ind));
%         stdHa(i) = nanstd(d2.err(:,ind));
%     end

    for j = 1:size(d.err,1)
        if (abs(d.err(j,i) > prctile(d.err(:,i),99.9)))
            d.err(j,i) = nan;
        end
    end
    meanHa(i) = nanmean(d.err(:,i));
    stdHa(i) = nanstd(d.err(:,i));
    
    nCrashes(i) = length(find(d.err(:,i) <= -(tgt + 150)));
end

% meanHa(1) = mean(fj.(fj_ind).haf_err);
% stdHa(1) = std(fj.(fj_ind).haf_err);
% for j = 1:Nj
%     pcFJ(j) = get_percent_captured(fj.(fj_ind).haf_err, 2000, tols(j),d.mc.N);
% %     pcFJ(i) = get_percent_captured(fj.out4.haf_err, 2000,tols(i),d.mc.N);
% end

entries{1} = 'Free Jettison (0% Bias)';
entries{2} = 'Biased Jettison';
% for i = 1:Ni
%     entries{i + 1} = ['Stage 1 Bias: ' num2str(b1(i)) '% of Tgt'];
% end

% colors = get_plot_colors();

% mean & std apoapsis vs. biasing
figure(); hold on
grid on
title(['Mean ha vs. bias, b21 = ' num2str(b21) ', ha_{tgt} = ' num2str(tgt)]);
set(gca,'FontSize',16);
yyaxis left
% plot(b1, ones(Ni,1) * meanHa(1),'--','LineWidth',2);
plot(bias, ones(Ni,1) * meanHa(1),'--','LineWidth',2);
plot(bias, meanHa,'-o','LineWidth',2);
xlabel('Stage 1 Bias (% of Tgt)')
ylabel('Mean Apoapsis Error (km)')
yyaxis right
% plot(b1, ones(Ni,1) * stdHa(1),'--','LineWidth',2);
plot(bias, ones(Ni,1) * stdHa(1),'--','LineWidth',2);
plot(bias, stdHa,'-o','LineWidth',2);
ylabel('Apoapsis 1\sigma (km)');
legend(entries,'location','best');

%}

%% Biasing MC mean/std 400km trend investigations
%{
clear;clc
% kk -> biasing1
% zz -> biasing2

d = load('..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\400_v11_5.4deg_b21=3.mat');
fj = load('..\venus_ac\dej_n\HYDRA\2stage\data\400\v11_5.4deg_b_ijs3_10.mat');
fj_ind = 'out';
Ni = length(d.biasing);

pcBias = ((d.biasing/1000) ./ d.gnc.ha_tgt) .* 100; %percent of target

errNew = nan(1000,Ni);
for i = 1:1000
    for j = 1:Ni
        if (abs(d.err(i,j)) < 1000)
            errNew(i,j) = d.err(i,j);
        end
    end
end

for i = 1:Ni
for j = 1:d.mc.N
    sigma(i,j) = nanstd(d.err(1:j,i));
end
end

% std convergence
figure(); hold on
for i=1:Ni
    plot(sigma(i,:));
end

for i = 1:Ni
    dej_scatter_plot_bias(d.tjr1(:,i), d.tjr2(:,i), ... 
        d.err(:,i), 2, d.gnc.ha_tgt, d.b21, pcBias(i));
end
%}

%% Biasing Percent Monte Carlo post-process
%{

% 400, -5.4 deg, b21 = 3
dfree = load('..\venus_ac\dej_n\HYDRA\2stage\data\400\v11_5.4deg_b_ijs3_10.mat');
d = load('..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\400_v11_5.4deg_b21=3.mat');
% dfree = load('..\venus_ac\dej_n\HYDRA\2stage\data\400\v11_5.4deg_b_ijs9_10.mat');
% d = load('..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\400_v11_5.4deg_b21=9.mat');

pcBias = [0.01 0.05 0.1 0.15 0.2] .* 100;
tols = 0:1:50;
Ni = length(tols);
pc = nan(Ni,5);
for j = 1:length(pcBias)
    pc(:,j) = get_percent_captured(d.err(:,j), d.gnc.ha_tgt, tols, d.mc.N);
end
pcFree = get_percent_captured(dfree.out.haf_err, dfree.gnc.ha_tgt, tols, dfree.mc.N);


colors = get_plot_colors();
figure(); hold on
grid on;
set(gca,'FontSize',16)
plot(tols, pcFree,'k--','LineWidth',2)
entry = cell(length(pcBias) + 1,1);
entry{1} = 'Free Jettison';
for j = 1:length(pcBias)
    plot(tols, pc(:,j),'Color',colors(j,:),'LineWidth',2)
    entry{j+1} = ['Biasing = -' num2str(pcBias(j)) '%'];
end
yline(100,'k--','LineWidth',1.5)
legend(entry,'location','se')
xlabel('Tolerance(% of Target)')
ylabel('Percent Captured (%)')


% % 2k, -5.4 deg, b21 = 3,9
% dfree = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\v11_5.4deg_b_ijs3_10.mat');
% d = load('..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\2k_v11_5.4deg_b21=3.mat');
dfree = load('..\venus_ac\dej_n\HYDRA\2stage\data\2k\v11_5.4deg_b_ijs9_10.mat');
d = load('..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_bias_val\2k_v11_5.4deg_b21=9.mat');

pcBias = [0.01 0.05 0.1 0.15 0.2] .* 100;
tols = 0:1:20;
Ni = length(tols);
pc = nan(Ni,5);
for j = 1:length(pcBias)
    pc(:,j) = get_percent_captured(d.err(:,j), d.gnc.ha_tgt, tols, d.mc.N);
end
pcFree = get_percent_captured(dfree.out.haf_err, dfree.gnc.ha_tgt, tols, dfree.mc.N);


colors = get_plot_colors();
figure(); hold on
grid on;
set(gca,'FontSize',16)
y0=plot([0 20],[100 100],'k--','LineWidth',1.5);
p0=plot(tols, pcFree,'k--','LineWidth',2);
entry = cell(length(pcBias) + 1,1);
entry{1} = 'Free Jettison';
for j = 1:length(pcBias)
    p1(j)=plot(tols, pc(:,j),'Color',colors(j,:),'LineWidth',2);
    entry{j+1} = ['Biasing = -' num2str(pcBias(j)) '%'];
end
legend([p0 p1],entry,'location','se')
xlabel('Tolerance(% of Target)')
ylabel('Percent Captured (%)')



%}

%% 2stage Monte Carlo trade: capture rate vs. b2/b1
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000];
ha_tgts = [2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [1.1, 2, 3, 4, 5, 6, 7, 8, 9, 9.9];
% b21s = 2;
v0 = 11;

for z = 1:Nz
% biasing = ha_tgts(z) .* [0.1 0.2 0.3 0.4 0.5] .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
% biasing = ha_tgts(z) * linspace(0.01, .2, 20)' * -1000;
% biasing = ha_tgts(z) * [.005 0.01 0.05 0.1 0.15 0.2] * -1000;
biasing = [1 10 20 30] ./ 100;  %percent of target
% biasing = 0;

[x0,aero,gnc,sim, mc] = base_venus_ac(true);
x0.v0 = v0;
efpas = [-5.4];
Nk = length(efpas);

gnc.ha_tgt = ha_tgts(z);
gnc.ha_tol = 5;
gnc.n = 2;
tj0s = [100 170];
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);

mc.debug = false;
mc.flag = true;
mc.N = 2000;

Nj = length(b21s);

for i = 1:length(biasing)
    
err = nan(mc.N, Nj);
tj1 = err; tj2 = err;
tjr1 = err; tjr2 = err;
    
for j = 1:Nj
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
for w = 1:3
    b(w) = aero.m(w) / (aero.cds(1) * pi * aero.rcs(w)^2);
end
for w = 1:2
    br(w) = b(w+1)/b(w);
end
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br31 = b(3)/b(1);

aero.rn = 0.1;
aero.cds = [1.05 1.05 1.05];
aero.cls = [0 0 0];

sim.traj_rate = 100;
sim.data_rate = 1;
sim.planet = 'venus';
sim.efpa_flag = false;

sim.parMode = true;
gnc.dtj_lim = 10;
% gnc.rss_flag = true;

% x0.fpa0 = fpas;
fprintf('2stage bias val MC, b2/b1 = %2.0f, ha* = %4.0f. ',br21, gnc.ha_tgt); print_current_time();
gnc.force_jett = true;
gnc.comp_curr = false;

for k = 1:Nk
x0.fpa0 = efpas(k);
    
% sim.nWorkers = 5;
gnc.tj0 = tj0s;

gnc.bias(1).tgt_ap = ha_tgts(z) * biasing(i) * -1000;
out = run_dej_n(x0,gnc,aero,sim,mc);
err(:,j) = out.haf_err;
tj1(:,j) = out.tjett(:,1);
tj2(:,j) = out.tjett(:,2);
tjr1(:,j) = out.tjr(:,1);
tjr2(:,j) = out.tjr(:,2);

end %k, efpas

end %j, b21s

if (ha_tgts(z) > 1000)
    tgt = [num2str(floor(ha_tgts(z)/1000)) 'k_'];
else
    tgt = [num2str(ha_tgts(z)) '_'];
end

save(['..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc2k_b21\' tgt ... 
    'v' num2str(v0) '_' num2str(abs(x0.fpa0)) 'deg_bias=' num2str(biasing(i)*100) '%.mat'])
fprintf('Finished 2stage bias val MC. '); print_current_time();

end %biasing

end %ha_tgts
keyboard
%}

%% capture rate vs b21 post process
%{
clear;clc
% p = '..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_b21\';
p = '..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc2k_b21\';

% % -5.4deg
d0 = load(['..\venus_ac\dej_n\HYDRA\2stage_bias\data\mc_b21\2k_v11_5.4deg_bias=0%.mat']);
d1 = load([p '2k_v11_5.4deg_bias=1%.mat']);
d10 = load([p '2k_v11_5.4deg_bias=10%.mat']);
d20 = load([p '2k_v11_5.4deg_bias=20%.mat']);
d30 = load([p '2k_v11_5.4deg_bias=30%.mat']);
% 
% d0 = load([p '400_v11_5.4deg_bias=0%.mat']);
% d1 = load([p '400_v11_5.4deg_bias=1%.mat']);
% d10 = load([p '400_v11_5.4deg_bias=10%.mat']);
% d20 = load([p '400_v11_5.4deg_bias=20%.mat']);
% d30 = load([p '400_v11_5.4deg_bias=30%.mat']);

% % -5.6deg
% d0 = load([p '2k_v11_5.6deg_bias=0%.mat']);
% d1 = load([p '2k_v11_5.6deg_bias=1%.mat']);
% d10 = load([p '2k_v11_5.6deg_bias=10%.mat']);
% d20 = load([p '2k_v11_5.6deg_bias=20%.mat']);
% d30 = load([p '2k_v11_5.6deg_bias=30%.mat']);
% 
% d0 = load([p '400_v11_5.6deg_bias=0%.mat']);
% d1 = load([p '400_v11_5.6deg_bias=1%.mat']);
% d10 = load([p '400_v11_5.6deg_bias=10%.mat']);
% d20 = load([p '400_v11_5.6deg_bias=20%.mat']);
% d30 = load([p '400_v11_5.6deg_bias=30%.mat']);

% % % % dej scatter plots
% plot_bias_scatter(d1, 2);




% % % % all bias percents
tols = [1 10];
Ni = length(tols);
Nj = length(d1.b21s);
pc1 = nan(Ni,Nj); pc0 = pc1;
pc10 = pc1; pc20 = pc1; pc30 = pc1;
entries = cell(Nj,1);
for i = 1:Ni
    for j = 1:Nj
        pc0(i,j) = get_percent_captured(d0.err(:,j), d0.gnc.ha_tgt, tols(i), d0.mc.N);
        pc1(i,j) = get_percent_captured(d1.err(:,j), d1.gnc.ha_tgt, tols(i), d1.mc.N);
        pc10(i,j) = get_percent_captured(d10.err(:,j), d10.gnc.ha_tgt, tols(i), d10.mc.N);
        pc20(i,j) = get_percent_captured(d20.err(:,j), d20.gnc.ha_tgt, tols(i), d20.mc.N);
        pc30(i,j) = get_percent_captured(d30.err(:,j), d30.gnc.ha_tgt, tols(i), d30.mc.N);
        
        [pc0_1(j), pc0_2(j)] = ...
            count_jettisons(d0.tj1(:,j), d0.tj2(:,j));
        
        [pc1_1(j), pc1_2(j)] = ...
            count_jettisons(d1.tj1(:,j), d1.tj2(:,j));
        
        [pc10_1(j), pc10_2(j)] = ...
            count_jettisons(d10.tj1(:,j), d10.tj2(:,j));
        
        [pc20_1(j), pc20_2(j)] = ...
            count_jettisons(d20.tj1(:,j), d20.tj2(:,j));
        
        [pc30_1(j), pc30_2(j)] = ...
            count_jettisons(d30.tj1(:,j), d30.tj2(:,j));
    end
    entries{i} = [num2str(tols(i)) '% Tolerance'];
end


% specific tols
colors = get_plot_colors();

figure(); hold on
set(gca,'FontSize',16)
title(['h_a^* Tol, ' num2str(d0.gnc.ha_tgt) ' km tgt, ' num2str(d0.x0.fpa0) 'deg']);

p0=plot(d0.b21s, pc0(1,:),'-o','Color',colors(1,:),'LineWidth',2);
p1=plot(d1.b21s, pc1(1,:),'-o','Color',colors(2,:),'LineWidth',2);
p10=plot(d10.b21s, pc10(1,:),'-o','Color',colors(3,:),'LineWidth',2);
p20=plot(d20.b21s, pc20(1,:),'-o','Color',colors(4,:),'LineWidth',2);
p30=plot(d30.b21s, pc30(1,:),'-o','Color',colors(5,:),'LineWidth',2);

p0_2=plot(d0.b21s, pc0(2,:),'--o','Color',colors(1,:),'LineWidth',2);
plot(d1.b21s, pc1(2,:),'--o','Color',colors(2,:),'LineWidth',2);
plot(d10.b21s, pc10(2,:),'--o','Color',colors(3,:),'LineWidth',2);
plot(d20.b21s, pc20(2,:),'--o','Color',colors(4,:),'LineWidth',2);
plot(d30.b21s, pc30(2,:),'--o','Color',colors(5,:),'LineWidth',2);

legend([p0 p1 p10 p20 p30 p0_2], ...
    '0% Bias 1% Tol','1% Bias','10% Bias','20% Bias', ... 
    '30% Bias', ... 
    '0% Bias, 10% Tol', ...
    'location','sw')
% legend([p0 p0_2], ...
%     '1% Tol','10% Tol', ...
%     'location','sw')
xlabel('\beta_2/\beta_1');
ylabel('Percent Captured (%)');
xlim([1.1 9.9])
xticks(d1.b21s)
xticklabels({'1.1','2','3','4','5','6','7','8','9','9.9'})


% percent 2nd jett
colors = get_plot_colors();
figure(); hold on
set(gca,'FontSize',16)
title(['h_a^* Tol, ' num2str(d0.gnc.ha_tgt) ' km tgt, ' num2str(d0.x0.fpa0) 'deg']);
p0=plot(d0.b21s, pc0_2,'-o','Color',colors(1,:),'LineWidth',2);
p1=plot(d1.b21s, pc1_2,'-o','Color',colors(2,:),'LineWidth',2);
p10=plot(d10.b21s, pc10_2,'-o','Color',colors(3,:),'LineWidth',2);
p20=plot(d20.b21s, pc20_2,'-o','Color',colors(4,:),'LineWidth',2);
p30=plot(d30.b21s, pc30_2,'-o','Color',colors(5,:),'LineWidth',2);

% p0=plot(d0.b21s, pc0(1,:),'-o','Color',colors(1,:),'LineWidth',2);
% p1=plot(d1.b21s, pc1(1,:),'-o','Color',colors(2,:),'LineWidth',2);
% p10=plot(d10.b21s, pc10(1,:),'-o','Color',colors(3,:),'LineWidth',2);
% p20=plot(d20.b21s, pc20(1,:),'-o','Color',colors(4,:),'LineWidth',2);
% p30=plot(d30.b21s, pc30(1,:),'-o','Color',colors(5,:),'LineWidth',2);

legend([p0 p1 p10 p20 p30], ...
    '0% Bias','1% Bias','10% Bias','20% Bias', ... 
    '30% Bias', ... 
    'location','sw')
xlabel('\beta_2/\beta_1');
ylabel('Percent Stage 2 Jettison (%)');
xlim([1.1 9.9])
xticks(d1.b21s)
xticklabels({'1.1','2','3','4','5','6','7','8','9','9.9'})



% % % % single bias percentage
% tols = 0:1:100;
% Ni = length(tols);
% Nj = length(d1.b21s);
% pc1 = nan(Ni,Nj);
% 
% entries = cell(Nj,1);
% for i = 1:Ni
%     for j = 1:Nj
%         pc1(i,j) = get_percent_captured(d1.err(:,j), d1.gnc.ha_tgt, tols(i), 1000);
%         entries{j} = ['\beta_2/\beta_1 = ' num2str(d1.b21s(j))];
%     end
%     %entries{i} = [num2str(tols(i)) '% Tolerance'];
% end
% 
% figure(); hold on
% set(gca,'FontSize',14)
% plot(tols, pc1,'LineWidth',2)
% xlabel('\beta_2/\beta_1');
% ylabel('Percent Captured (%)');
% legend(entries,'location','se')
% xlim([0 20])
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 3-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bias Value Trade Studies
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(false);
gnc.n = 3;
br41 = 10;
% br21s = [3 3 9];
% br31s = [5 7 9.5];
br21s = [3 9];
br31s = [5 9.5];

x0.fpa0 = -5.6;
ha_tgts = [400 2000];

stage2biases = linspace(0.001, 0.6, 30);    % percentage of 1st stage bias

for jj = 1:length(ha_tgts)
% first stage bias
biasing1 = linspace(0.001, 0.2, 30)';
Nkk = length(biasing1);
% second stage

Nzz = length(stage2biases);
for zz = 1:Nzz
    biasing2(:,zz) = biasing1 * stage2biases(zz);
end

tjr = nan(Nkk, Nzz, gnc.n);
tj = tjr;
err = nan(Nkk, Nzz, 1);

for k = 1:length(br31s)
    
for kk = 1:Nkk

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

gnc.ha_tgt = ha_tgts(jj);
gnc.tj0 = [90 170 260];

sim.t_max = 1000;
% 
for x = 1:4
    b(x) = aero.m(x)/( aero.cds(1) * pi * aero.rcs(x)^2 );
end
br41 = b(4)/b(1);
br31 = b(3)/b(1);
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br43 = b(4)/b(3);
br42 = b(4)/b(2);

bijs = [num2str(br21) '_' num2str(br31) '_' num2str(br41)];

mc.flag = false;
mc.N = 1;
mc.sigs = zero_sigs();

sim.traj_rate = 250;
sim.planet = 'venus';
sim.efpa_flag = false;
% sim.debug = true;
gnc.dtj_lim = 10;

gnc.iters = uint8(1);
gnc.guid_rate = 0.5;
gnc.ha_tol = 0.1;
gnc.npc_mode = uint8(2); %multi stage
gnc.comp_curr = false;
gnc.force_jett = true;


t = datetime();
fprintf(['Starting 3stage bias sim, betas = ' ... 
    num2str(br21) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

for zz = 1:Nzz %stage 2 bias
    
gnc.bias(1).tgt_ap = gnc.ha_tgt * biasing1(kk) * -1000;
gnc.bias(2).tgt_ap = gnc.ha_tgt * biasing2(kk,zz) * -1000;

out = run_dej_n(x0,gnc,aero,sim,mc);

tj(kk,zz,:) = out.t_jett';
tjr(kk,zz,:) = out.t_jett./out.traj.time(out.idxend);
err(kk,zz) = out.haf_err;

end %zz biasing stage 2s
fprintf(['Finished 3stage bias sim, b_ijs = ' ... 
    num2str(bijs) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

end %kk biasing 1

if (gnc.ha_tgt < 1000)
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\bias_val\' num2str(floor(gnc.ha_tgt)) ... 
        '_v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_bijs=' num2str(bijs) '.mat'])
else
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\bias_val\' num2str(floor(gnc.ha_tgt/1000)) ... 
        'k_v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_bijs=' num2str(bijs) '.mat'])
end

end %k br31s

end %jj ha tgts
keyboard;

%}

%% Bias Val post-process
%{

tol = 10;
% -5.4deg
inds = [2 4 6 10 15];
plot_tj_vs_bias(3, 400, -5.4, [3 5], tol, inds, true, 0);
plot_tj_vs_bias(3, 2000, -5.4, [3 5], tol, inds, true, 0);

plot_tj_vs_bias(3, 400, -5.4,  [9 9.5], tol, inds, true, 1);
plot_tj_vs_bias(3, 2000, -5.4, [9 9.5], tol, inds, true, 1);

% -5.6deg
inds = [2 4 6 8 10 12];
plot_tj_vs_bias(3, 400, -5.6, [9 9.5], 1, inds, true, 1);

plot_tj_vs_bias(3, 2000, -5.6, [3 5], 1, inds, true, 1);
plot_tj_vs_bias(3, 2000, -5.6, [9 9.5], 1, inds, true, 0);




d = load('..\venus_ac\dej_n\HYDRA\3stage_bias\bias_val\2k_v11_5.4deg_bijs=3_5_10.mat');
bias1 = d.biasing1;
bias2 = d.biasing2; [Ni, Nj] = size(bias2);
err = d.err;
tj1 = nan(Ni,Nj,1); tj2 = tj1; tj3 = tj1;
tjr1 = tj1; tjr2 = tj1; tjr3 = tj1;
b1 = tj1;
b2 = tj1;

tol = 1;    %km
for i = 1:Ni
    for j = 1:Nj
        if (abs(err(i,j)) <= tol)            
            tjr1(i,j) = d.tjr(i,j,1);
            tjr2(i,j) = d.tjr(i,j,2);
            tjr3(i,j) = d.tjr(i,j,3);
            
            tj1(i,j) = d.tj(i,j,1);
            tj2(i,j) = d.tj(i,j,2);
            tj3(i,j) = d.tj(i,j,3);         
        end
        
        b1(:,j) = bias1;
        b2(i,j) = bias2(i,j)/bias1(i);
    end
end

% figure(); hold on;
% plot(bias2(:,2)./1000, tjr3(:,2))
% plot(bias2(:,4)./1000, tjr3(:,4))


inds = [2, 6, 10, 15];
figure(); hold on
colors = get_plot_colors();
set(gca,'FontSize',16)
entries = cell(length(inds),1);
for i = 1:length(inds)
    entries{i} = ['Stage 1 Bias=' ... 
        num2str(round(100 * bias1(inds(i))/1000 / d.gnc.ha_tgt,2)) '% Tgt'];
    
    plot(bias2(inds(i),:)./bias1(inds(i)), ...
     tjr2(inds(i),:),'Color',colors(i,:),'LineStyle','--','LineWidth',2);
    p(i) = plot(bias2(inds(i),:)./bias1(inds(i)), ...
        tjr3(inds(i),:),'Color',colors(i,:),'LineWidth',2);
end
xlabel('Stage 2 Bias (% of Stage 1)')
ylabel('t_{jett}/t_f')
legend(p,entries,'location','best');
title(['EFPA = ' num2str(d.x0.fpa0) ', b21 = ' num2str(d.br21) ...
    ', ' num2str(d.gnc.ha_tgt)  'km tgt'])

% figure(); hold on
% contour(bias1./1000, bias2(end,:)./1000, tjr3')



figure(); hold on
set(gca,'FontSize',16)
surf(b1./1000, bias2./1000, tjr3)
xlabel('Bias Stage 1 (km)')
ylabel('Bias Stage 2 (km)')
zlabel('tjr3')


% bias2 = f(bias1).. plot of tjrX vs bias2 (% of bias1)


%}

%% Biasing EFPA Trade Studies
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(false);
gnc.n = 3;
br41 = 10;
% br21s = [3 3 9];
% br31s = [5 7 9.5];
br21s = 3;
br31s = 5;

efpas = linspace(-5.35,-5.8,100);
Nkk = length(efpas);

ha_tgts = [2000];

for jj = 1:length(ha_tgts)
% first stage bias
biasing(:,1) = ha_tgts(jj) .* [0.05 0.2]' .* -1000;
% second stage bias
%biasing(:,2) = ha_tgts(jj) .* [0.1 0.2 0.1]' .* -1000;

stage2biases = linspace(0.01, 0.9, 50);
Nzz = length(stage2biases);

for k = 1:length(br31s)
    
for zz = 1:Nzz
biasing(:,2) = biasing(:,1) .* stage2biases(zz);

tjr = nan(length(biasing),Nkk,gnc.n);
tj = tjr;
err = nan(length(biasing),Nkk);

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

gnc.ha_tgt = ha_tgts(jj);
gnc.tj0 = [90 170 260];

sim.t_max = 1000;
% 
for x = 1:4
    b(x) = aero.m(x)/( aero.cds(1) * pi * aero.rcs(x)^2 );
end
br41 = b(4)/b(1);
br31 = b(3)/b(1);
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br43 = b(4)/b(3);
br42 = b(4)/b(2);

bijs = [num2str(br21) '_' num2str(br31) '_' num2str(br41)];

mc.flag = false;
mc.N = 1;
mc.sigs = zero_sigs();

sim.traj_rate = 250;
sim.planet = 'venus';
sim.efpa_flag = false;
% sim.debug = true;
gnc.dtj_lim = 10;

gnc.iters = uint8(1);
gnc.guid_rate = 0.5;
gnc.ha_tol = 0.1;
gnc.npc_mode = uint8(2); %multi stage
gnc.comp_curr = false;
gnc.force_jett = true;


t = datetime();
fprintf(['Starting 3stage bias sim, betas = ' ... 
    num2str(br21) ',' num2str(br31) ',' num2str(br41)]); print_current_time();
for kk = 1:length(efpas)
x0.fpa0 = efpas(kk);

for i = 1:length(biasing)

gnc.bias(1).tgt_ap = biasing(i,1);
gnc.bias(2).tgt_ap = biasing(i,2);
out = run_dej_n(x0,gnc,aero,sim,mc);

tj(i,kk,:) = out.t_jett';
tjr(i,kk,:) = out.idj/out.idxend;
err(i,kk) = out.haf_err;

end %biasing, i

end %kk fpas

end %zz biasing stage 2s
fprintf(['Finished 3stage bias sim, b_ijs = ' ... 
    num2str(bijs) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

if (gnc.ha_tgt < 1000)
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\entry_trades\' num2str(floor(gnc.ha_tgt)) ... 
        '_v' num2str(x0.v0) '_bijs=' num2str(bijs) '_biasSet' ... 
        '_1.mat'])
else
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\entry_trades\' num2str(floor(gnc.ha_tgt/1000)) ... 
        'k_v' num2str(x0.v0) '_bijs=' num2str(bijs) '_biasSet' ... 
        '_1.mat'])
end

end %k br31s

end %jj ha tgts
keyboard;


% est err vs time
figure(); hold on
title('2k, v11, b21 = 3')
set(gca,'FontSize',16);
plot(out.traj.time,out.g.dej_n.dr_ap_true./1000,'r','LineWidth',2);
p1=xline(out.traj.time(out.idj(1)),'k--','LineWidth',1.5);
p2=xline(out.traj.time(out.idj(2)),'k-.','LineWidth',1.5);
p3=xline(out.traj.time(out.idj(3)),'k:','LineWidth',1.5);
yline(-100, '-','Color',[0,0,0]+0.5,'LineWidth',1.5)
yline(-500, '-','Color',[0,0,0]+0.5,'LineWidth',1.5)
ylim([-500 200]);
xlim([0 out.traj.time(out.idxend)])
xlabel('Time (s)')
ylabel('Estimated Apoapsis Error (km)')
legend([p1 p2 p3],'Jettison 1','Jettison 2','Jettison 3','location','ne')
%}

%% 3stage EFPA trade post-processing
%{
clear;clc
tol = 75;

% plot_tj_opt(3, 400, [3 5], tol, true, false);
% % plot_tj_opt(3, 400, [3 7], tol, true, false);
% plot_tj_opt(3, 400, [9 9.5], tol, true, false);

plot_tj_opt(3, 2000, [3 5], tol, true, false);
% plot_tj_opt(3, 2000, [3 7], tol, true, false);
plot_tj_opt(3, 2000, [9 9.5], tol, true, false);

%}

%% 3stage - investigate expected_haf_err vs time
%{
% It was an issue with the apoapis error tolerance being too high

clear;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(false);
gnc.n = 3;
br21 = 3;
br31 = 5;
br41 = 10;

efpas = [-5.45 -5.5 -5.55];
Ni = length(efpas);
ha_tgts = 400;
    

biasing(2) = ha_tgts * [0.01]' .* -1000;
biasing(1) = ha_tgts * [0.1]' .* -1000;

tjr = nan(length(biasing),Ni,gnc.n);
tj = tjr;
err = nan(length(biasing),Ni);

aero.m = [150 60 40 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * br21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * br31) );
aero.rcs(4) = sqrt( (aero.m(4) * aero.rcs(1)^2) / ...
    (aero.m(1) * br41) );

aero.rn = 0.75/2;
aero.cds = [1.05 1.05 1.05 1.05];
aero.cls = [0 0 0 0];

gnc.ha_tgt = ha_tgts;
gnc.tj0 = [90 170 260];
sim.t_max = 1750;
% 
for x = 1:4
    b(x) = aero.m(x)/( aero.cds(1) * pi * aero.rcs(x)^2 );
end
br41 = b(4)/b(1);
br31 = b(3)/b(1);
br21 = b(2)/b(1);
br32 = b(3)/b(2);
br43 = b(4)/b(3);
br42 = b(4)/b(2);

bijs = [num2str(br21) '_' num2str(br31) '_' num2str(br41)];

mc.flag = false;
mc.N = 1;
mc.sigs = zero_sigs();

sim.traj_rate = 1000;
sim.planet = 'venus';
sim.efpa_flag = false;

% %% pdg algorithm
% gnc.mode = 3;   %pdg
% gnc.guid_rate = 50;
% 
% % sim.debug = true;
% for i = 1:Ni
%     fprintf('Case %i: ',i); print_current_time();
%     gnc.bias(1).tgt_ap = biasing(1);
%     gnc.bias(2).tgt_ap = biasing(2);
%     x0.fpa0 = efpas(i);
%     out(i) = run_dej_n(x0,gnc,aero,sim,mc);
% end
% 
% save('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/investigate_3stage.mat');
% keyboard;
% d = load('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/investigate_3stage.mat');
% Ni = length(d.efpas);
% out1 = d.out(1);
% out2 = d.out(2);
% out3 = d.out(3);
% 
% for i = 1:Ni
%     entries{i} = ['EFPA = ' num2str(d.efpas(i))];
% end
% 
% figure(); hold on
% title(['Biasing = [' num2str(d.biasing./1000) '] km'])
% set(gca,'FontSize',14)
% plot(out1.traj.time./out1.traj.time(out1.idxend), out1.g.pd.dr_ap_true./1000,'LineWidth',2)
% plot(out2.traj.time./out2.traj.time(out2.idxend), out2.g.pd.dr_ap_true./1000,'LineWidth',2)
% plot(out3.traj.time./out3.traj.time(out3.idxend), out3.g.pd.dr_ap_true./1000,'LineWidth',2)
% legend(entries)
% xlabel('t/t_f')
% ylabel('Estimated Apoapsis Error (km)')
% 
% figure(); hold on
% title(['Biasing = [' num2str(d.biasing./1000) '] km'])
% set(gca,'FontSize',14)
% plot(out1.traj.vel_pp_mag./1000, out1.traj.alt./1000,'LineWidth',2)
% plot(out2.traj.vel_pp_mag./1000, out2.traj.alt./1000,'LineWidth',2)
% plot(out3.traj.vel_pp_mag./1000, out3.traj.alt./1000,'LineWidth',2)
% legend(entries)
% xlabel('Velocity (km/s)')
% ylabel('Altitude (km)')

%% npc
gnc.mode = 4;   %npc
gnc.dtj_lim = 10;
gnc.ha_tol = 0.5;
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);
sim.traj_rate = 100;

% sim.debug = true;
for i = 1:Ni
    fprintf('Case %i: ',i); print_current_time();
    gnc.bias(1).tgt_ap = biasing(1);
    gnc.bias(2).tgt_ap = biasing(2);
    x0.fpa0 = efpas(i);
    out(i) = run_dej_n(x0,gnc,aero,sim,mc);
end
% % save('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/investigate_3stage_npc.mat');
% keyboard;
% % post process
d = load('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/investigate_3stage_npc.mat');
Ni = length(d.efpas);
out1 = d.out(1);
out2 = d.out(2);
out3 = d.out(3);
% 
% for i = 1:Ni
%     entries{i} = ['EFPA = ' num2str(d.efpas(i))];
% end
% entries{Ni + 1} = 'Jettison 1';
% entries{Ni + 2} = 'Jettison 2';
% entries{Ni + 3} = 'Jettison 3';
% 
% 
% 
% tf1 = out1.traj.time(out1.idxend);
% tf2 = out2.traj.time(out2.idxend);
% tf3 = out3.traj.time(out3.idxend);
% 
% figure(); hold on
% title(['Biasing = [' num2str(d.biasing./1000) '] km'])
% set(gca,'FontSize',14)
% plot(out1.traj.time./tf1, out1.g.dej_n.dr_ap_true./1000,'LineWidth',2)
% plot(out2.traj.time./tf2, out2.g.dej_n.dr_ap_true./1000,'LineWidth',2)
% plot(out3.traj.time./tf3, out3.g.dej_n.dr_ap_true./1000,'LineWidth',2)
% xline(out1.t_jett(1)./tf1,'k-','LineWidth',2)
% xline(out1.t_jett(2)./tf1,'k--','LineWidth',2)
% xline(out1.t_jett(3)./tf1,'k:','LineWidth',2)
% xline(out2.t_jett(1)./tf2,'k-','LineWidth',2)
% xline(out3.t_jett(1)./tf3,'k-','LineWidth',2)
% xline(out2.t_jett(2)./tf2,'k--','LineWidth',2)
% xline(out3.t_jett(2)./tf3,'k--','LineWidth',2)
% xline(out2.t_jett(3)./tf2,'k:','LineWidth',2)
% xline(out3.t_jett(3)./tf3,'k:','LineWidth',2)
% yline(d.biasing(1)/1000,'Color',[0 0 0] + 0.5,'LineWidth',2)
% yline(d.biasing(2)/1000,'Color',[0 0 0] + 0.5,'LineWidth',2)
% xlabel('t/t_f')
% ylabel('Estimated Apoapsis Error (km)')
% legend(entries,'location','se')
% ylim([-5 5])

%% manual
gnc.mode = uint8(6);   %manual
gnc.dtj_lim = 10;
gnc.guid_rate = 1;
sim.traj_rate = 1000;

x0.fpa0 = round(out3.traj.gamma_pp(1) * 180/pi, 3);
gnc.tj0 = out3.t_jett(1:d.gnc.n);

tj3s = linspace(d.gnc.tj0(2), 300,50);
Nj = length(tj3s);

% sim.debug = true;
% for i = 1:Ni
% fprintf('Case %i: ',i); print_current_time();
% gnc.bias(1).tgt_ap = biasing(1);
% gnc.bias(2).tgt_ap = biasing(2);
err = nan(Nj,1);
for j = 1:Nj
% x0.fpa0 = efpas(i);
    gnc.tj0(3) = tj3s(j);
    out = run_dej_n(x0,gnc,aero,sim,mc);
    err(j) = out.haf_err;
end
save('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/3stage_manual_3.mat');
keyboard;
% post process
d1 = load('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/3stage_manual_1.mat');
d2 = load('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/3stage_manual_2.mat');
d3 = load('../venus_ac/dej_n/HYDRA/3stage_bias/entry_trades/3stage_manual_3.mat');

entries{1} = ['EFPA = ' num2str(d1.x0.fpa0)];
entries{2} = ['EFPA = ' num2str(d2.x0.fpa0)];
entries{3} = ['EFPA = ' num2str(d3.x0.fpa0)];


figure(); hold on
plot(d1.tj3s, d1.err,'LineWidth',2)
plot(d2.tj3s, d2.err,'LineWidth',2)
plot(d3.tj3s, d3.err,'LineWidth',2)
yline(0);
% xlim([100 300]);
ylim([-10 200]);
legend(entries)

%}

%% 3stage - monte carlo
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac(true);
gnc.n = 3;
br41 = 10;
br21s = [3];
br31s = [5];

x0.fpa0 = -5.4;
ha_tgts = [2000];
% ha_tgts = 2000;

stage2biases = linspace(0.1, 90, 10) ./ 100;    % percentage of 1st stage

for jj = 1:length(ha_tgts)
% first stage bias
% biasing1 = linspace(0.001, 0.2, 30)';
biasing1 = [0.1, 1, 10, 20] ./100;  %percent of apoapsis tgt
Nkk = length(biasing1);
% second stage

Nzz = length(stage2biases);
for zz = 1:Nzz
    biasing2(:,zz) = biasing1 * stage2biases(zz);
end

mc.flag = true;
mc.N = 1000;
sim.parMode = true;
sim.nWorkers = 10;

tjr1 = nan(Nkk, Nzz, mc.N); tjr2 = tjr1; tjr3 = tjr1;
tj1 = tjr1; tj2 = tj1; tj3 = tj1;
err = nan(Nkk, Nzz, mc.N);

gnc.ha_tgt = ha_tgts(jj);
gnc.tj0 = [90 170 260];
sim.t_max = 1500;

for k = 1:length(br31s)
    
for kk = 1:Nkk

br21 = br21s(k);
br31 = br31s(k);

[aero, bijs] = get_n3_betas(aero, br21, br31, br41);

sim.traj_rate = 100;
sim.planet = 'venus';
sim.efpa_flag = false;
% sim.debug = true;
gnc.dtj_lim = 10;

gnc.iters = uint8(1);
gnc.guid_rate = 0.5;
gnc.ha_tol = 0.1;
gnc.npc_mode = uint8(2); %multi stage
gnc.comp_curr = false;
gnc.force_jett = true;


t = datetime();
fprintf(['Starting 3stage MC sim, betas = ' ... 
    num2str(br21) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

for zz = 1:Nzz %stage 2 bias
    
gnc.bias(1).tgt_ap = gnc.ha_tgt * biasing1(kk) * -1000; %m
gnc.bias(2).tgt_ap = gnc.ha_tgt * biasing2(kk,zz) * -1000; %m

out = run_dej_n(x0,gnc,aero,sim,mc);
% beep
err(kk,zz,:) = out.haf_err;
tj1(kk,zz,:) = out.tjett(:,1);
tj2(kk,zz,:) = out.tjett(:,2);
tj3(kk,zz,:) = out.tjett(:,3);
tjr1(kk,zz,:) = out.tjr(:,1);
tjr2(kk,zz,:) = out.tjr(:,2);
tjr3(kk,zz,:) = out.tjr(:,3);

end %zz biasing stage 2s
fprintf(['Finished 3stage MC sim, b_ijs = ' ... 
    num2str(bijs) ',' num2str(br31) ',' num2str(br41)]); print_current_time();

end %kk biasing 1

if (gnc.ha_tgt < 1000)
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\mc\' num2str(floor(gnc.ha_tgt)) ... 
        '_v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_bijs=' num2str(bijs) '_2.mat'])
else
    save(['..\venus_ac\dej_n\HYDRA\3stage_bias\mc\' num2str(floor(gnc.ha_tgt/1000)) ... 
        'k_v' num2str(x0.v0) '_' num2str(abs(x0.fpa0)) 'deg_bijs=' num2str(bijs) '_2.mat'])
end

end %k br31s

end %jj ha tgts
keyboard;

%}


%% 3stage MC post processing
%{
clear;clc
% kk -> biasing1
% zz -> biasing2
% d = load('..\venus_ac\dej_n\HYDRA\3stage_bias\mc\2k_v11_5.4deg_bijs=3_5_10.mat');
d = load('..\venus_ac\dej_n\HYDRA\3stage_bias\mc\2k_v11_5.4deg_bijs=9_9.5_10.mat');


p = '../venus_ac/dej_n/HYDRA/';
p3 = [p '3stage/data/'];
temp = load([p3 '2k/v11_5.4deg_b_ijs3_5_10.mat']);
fj.out1 = temp.out;
temp = load([p3 '2k/v11_5.4deg_b_ijs5_7_10.mat']);
fj.out2 = temp.out;
temp = load([p3 '2k/v11_5.4deg_b_ijs7_9_10.mat']);
fj.out3 = temp.out;
temp = load([p3 '2k/v11_5.4deg_b_ijs9_9.5_10.mat']);
fj.out4 = temp.out;


% 400 FREE JETTISON MISSING
% d = load('..\venus_ac\dej_n\HYDRA\3stage_bias\mc\400_v11_5.4deg_bijs=3_5_10.mat');
% d = load('..\venus_ac\dej_n\HYDRA\3stage_bias\mc\2k_v11_5.4deg_bijs=9_9.5_10.mat');

% a = 'out1';   %b21=3
a = 'out4'; %b21 = 9

tols = 0:0.1:20;
Ni = length(tols);
pc = nan(d.Nkk, d.Nzz, Ni);
j1s = pc; j2s = pc; j3s = pc;
pcFJ = nan(Ni,1);

meanHa = nan(d.Nkk,d.Nzz);
stdHa = meanHa;
entries = cell(d.Nkk + 1,1);

b1 = d.biasing1 .* -100; % percent of tgt
b2 = d.biasing2 .* -100; % percent of tgt
for kk = 1:d.Nkk
    for zz = 1:d.Nzz
        for i = 1:Ni
            pc(kk,zz, i) = get_percent_captured( ... 
                d.err(kk,zz,:),d.gnc.ha_tgt, tols(i), 1000);
            
            [j1s(kk,zz,i), j2s(kk,zz,i), j3s(kk,zz,i)] = ... 
                count_jettisons(d.tj1(kk,zz,:), ... 
                d.tj2(kk,zz,:), d.tj3(kk,zz,:));            
        end
        
        meanHa(kk,zz) = mean(d.err(kk,zz,:));
        stdHa(kk,zz) = std(d.err(kk,zz,:));
    end
end

for i = 1:Ni
    pcFJ(i) = get_percent_captured(fj.(a).haf_err, 2000,tols(i),1000);
end

j1FJ = nan(6,1); j2FJ = j1FJ; j3FJ = j1FJ;
for i = 1:6
    [j1FJ(i), j2FJ(i), j3FJ(i)] = count_jettisons(... 
            fj.(a).tjett(:,1), fj.(a).tjett(:,2), fj.(a).tjett(:,3));
end

entries{1} = 'Free Jettison (0% Bias)';
for i = 1:6
    entries{i + 1} = ['Stage 1 Bias: ' num2str(b1(i)) '% of Tgt'];
end

colors = get_plot_colors();

% mean apoapsis vs. biasing
figure(); hold on
grid on
title(['Mean apoapsis vs. biasing, b21 = ' num2str(d.br21)]);
set(gca,'FontSize',16);
plot(b2(:,end), ones(6,1) * mean(fj.(a).haf_err),'k--','LineWidth',2);
plot(b2, meanHa,'-o','LineWidth',2);
xlabel('Stage 2 Bias (% of Tgt)')
ylabel('Mean Apoapsis Error (km)')
legend(entries,'location','best');
% ylim([-50, mean(fj.(a).haf_err) + 10]);

% apoapsis std vs. biasing
figure(); hold on
grid on
title(['apoapsis std vs. biasing, b21 = ' num2str(d.br21)]);
set(gca,'FontSize',16);
plot(b2(:,end), ones(6,1) * std(fj.(a).haf_err),'k--','LineWidth',2);
plot(b2, stdHa,'-o','LineWidth',2)
xlabel('Stage 2 Bias (% of Tgt)')
ylabel('Apoapsis 1\sigma (km)')
legend(entries,'location','best');
% ylim([std(fj.(a).haf_err) - 5, 70]);



% ind = d.Nkk;
ind = 6;
entries = cell(d.Nzz + 1,1);
entries{1} = 'Free Jettison (0% Bias)';
for z = 1:d.Nzz
    entries{z + 1} = ['Stage 2 Bias = ' num2str(b2(ind,z)) '% Tgt'];
end

% percent captured vs. tol at specific stage 1 bias
figure(); hold on
title(['cdf, stage1 bias = ' num2str(b1(ind)) '%']);
set(gca,'FontSize',16);
plot(tols, pcFJ, 'k--','LineWidth',2);
plot(tols, reshape(pc(ind,:,:),[d.Nkk Ni]),'LineWidth',2);
xlim([0 20]);
xlabel('Apoapsis Altitude Tolerance (% of Tgt)')
ylabel('Percent Captured (%)')
legend(entries,'location','se');


% 2nd stage jettison rate vs. tol at specific stage 1 bias
figure(); hold on
title(['stage 3 jettison rate, b21 = ' num2str(d.br21)]);
set(gca,'FontSize',16);
% plot(b2(:,1), j3FJ, 'k--','LineWidth',2);
yline(j2FJ(1), 'k--','LineWidth',2)
plot(b2, j2s(:,:,1),'-o','LineWidth',2);
xlabel('Stage 2 Bias (% of Tgt)')
ylabel('Second Stage Jettison Rate (%)')
legend(entries,'location','se');

% 3rd stage jettison rate vs. tol at specific stage 1 bias
figure(); hold on
title(['stage 3 jettison rate, b21 = ' num2str(d.br21)]);
set(gca,'FontSize',16);
% plot(b2(:,1), j3FJ, 'k--','LineWidth',2);
yline(j3FJ(1), 'k--','LineWidth',2)
plot(b2, j3s(:,:,1),'-o','LineWidth',2);
xlabel('Stage 2 Bias (% of Tgt)')
ylabel('Third Stage Jettison Rate (%)')
legend(entries,'location','se');


% % percent captured vs. stage 2 bias, specific tol
% tol = [1 5];
% % tol = 1;
% ind = d.Nkk;
% 
% for j = 1:length(tol)
%     tInd(j) = find(tols == tol(j),1);
% end
% 
% 
% 
% entries = cell(d.Nzz + 1,1);
% entries{1} = 'Free Jettison (0% Bias)';
% for z = 1:d.Nzz
%     entries{z + 1} = ['Stage 2 Bias = ' num2str(b2(ind,z)) '% Tgt'];
% end
% 
% figure(); hold on
% title(['cdf, stage1 bias = ' num2str(b1(ind)) '%']);
% set(gca,'FontSize',16);
% % plot(b2, pcFJ, 'k--','LineWidth',2);
% plot(b1, pc(:,:,tInd(1)),'LineWidth',2);
% % plot(b1, pc(:,:,tInd(2)),'--','LineWidth',2);
% % xlim([0 20]);
% xlabel('Stage 1 Bias (% of Tgt)')
% ylabel('Percent Captured (%)')
% legend(entries,'location','se');

%}


