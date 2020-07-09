% % venus_ac_biasing.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 2-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nominal Visualization
%{
clear;clc
% h vs v
[x0,aero,gnc,sim, mc] = base_venus_ac();

biasing = -500e3;
b21 = 4;
b31 = 10;
gnc.ha_tgt = 2000;
x0.v0 = 11;
sim.traj_rate = 1000;

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

sim.traj_rate = 200;
sim.efpa_flag = true;
mc.flag = false;
gnc.guid_rate = 0;

% beta 1 traj
aero.m = m(1);
aero.cds = cd;
aero.rcs = [rcs(1)];
gnc.n = 0;
out1 = run_dej_n(x0,gnc,aero,sim,mc);
fpa1 = out1.traj.gamma_pp(1) * 180/pi;

% beta 2 traj
aero.m = m(2);
aero.cds = cd;
aero.rcs = [rcs(2)];
gnc.n = 0;
out2 = run_dej_n(x0,gnc,aero,sim,mc);
fpa2 = out2.traj.gamma_pp(1) * 180/pi;

% beta 3 traj
aero.m = m(3);
aero.cds = cd;
aero.rcs = [rcs(3)];
gnc.n = 0;
x0.fpa0 = fpa2;
out3 = run_dej_n(x0,gnc,aero,sim,mc);
fpa3 = out3.traj.gamma_pp(1) * 180/pi;

% biased control
aero.m = [m(1) m(2) m(3)];
aero.cds = [cd cd cd];
aero.rcs = [rcs(1) rcs(2) rcs(3)];
aero.cl = [0 0 0];
gnc.n = 2;
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(2);
gnc.bias(1).tgt_ap = biasing;
gnc.tj0 = [90 160];
x0.fpa0 = -5.5;
sim.efpa_flag = false;
outMEJ = run_dej_n(x0,gnc,aero,sim,mc);

%single jettison
aero.m = [m(1) m(2)];
aero.cds = [cd cd];
aero.rcs = [rcs(1) rcs(2)];
aero.cl = [0 0];
gnc.n = 1;
gnc.guid_rate = 0.5;
gnc.npc_mode = uint8(1);
gnc.bias(1).tgt_ap = 0;
gnc.tj0 = [90];
x0.fpa0 = -5.5;
sim.efpa_flag = false;
outSEJ = run_dej_n(x0,gnc,aero,sim,mc);

keyboard;

% alt vs vel
figure(); hold on
title('-500km bias, 2k target, v=11')
set(gca,'FontSize',16)
plot(out1.traj.vel_pp_mag/1000, out1.traj.alt/1000,'k','LineWidth',2)
plot(out2.traj.vel_pp_mag/1000, out2.traj.alt/1000,'k--','LineWidth',2)
plot(out3.traj.vel_pp_mag/1000, out3.traj.alt/1000,'k-.','LineWidth',2)
plot(outMEJ.traj.vel_pp_mag/1000, outMEJ.traj.alt/1000,'b', 'LineWidth',2)
plot(outSEJ.traj.vel_pp_mag/1000, outSEJ.traj.alt/1000, 'r-.','LineWidth',2)
plot(outMEJ.traj.vel_pp_mag(outMEJ.idj(1))/1000, ...
    outMEJ.traj.alt(outMEJ.idj(1))/1000, 'bo','LineWidth',2)
plot(outMEJ.traj.vel_pp_mag(outMEJ.idj(2))/1000, ...
    outMEJ.traj.alt(outMEJ.idj(2))/1000, 'bo','LineWidth',2)
plot(outSEJ.traj.vel_pp_mag(outSEJ.idj(1))/1000, ...
    outSEJ.traj.alt(outSEJ.idj(1))/1000, 'ro','LineWidth',2)
legend(['\beta_1, EFPA = ', num2str(fpa1)], ...
    ['\beta_2 = ' num2str(outMEJ.beta_ratios(1)) '\beta_1, EFPA = ', num2str(fpa2)], ...
    ['\beta_3 = ' num2str(outMEJ.beta_ratios(2) * outMEJ.beta_ratios(1)) '\beta_1, EFPA = ', num2str(fpa3)], ...
    ['\beta_{i+1}/\beta_i = \{' num2str(outMEJ.beta_ratios(1)) ... 
     ', ' num2str(outMEJ.beta_ratios(2)) '\}, EFPA = ', num2str(x0.fpa0)], ...
    ['SEJ, \beta_2/\beta_1 = ' num2str(outMEJ.beta_ratios(1)) ', EFPA = ', num2str(x0.fpa0)], ...
    'location','best');
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')

% time vs. apoapsis
figure(); hold on
title('-500km bias, 2k target, v=11')
set(gca,'FontSize',16)
plot(outMEJ.traj.time, outMEJ.traj.h_ap./1000,'LineWidth',2)
plot(outSEJ.traj.time, outSEJ.traj.h_ap./1000, 'LineWidth',2)
ylim([1800 10000])
set(gca, 'YScale', 'log')
legend('MEJ','SEJ')

%}

%% Bias Value Trade
% %{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [400 2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 3;
v0 = 11;
x0.fpa0 = -5.6;

for z = 1:Nz
biasing = ha_tgts(z) .* linspace(0.001, 0.8, 50) .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
% biasing = ha_tgts(z) * 0.25 * -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac();
x0.v0 = v0;

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

sim.traj_rate = 500;
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
% %{
d1 = load('..\venus_ac\dej_n\HYDRA\2stage_bias\bias_val\2k_v11_5.4deg_b21_3.mat');
d2 = load('..\venus_ac\dej_n\HYDRA\2stage_bias\bias_val\2k_v11_5.4deg_b21_9.mat');
Ni = length(d1.biasing);
tjr1_1 = nan(Ni,1); tjr1_2 = nan(Ni,1);
tjr2_1 = nan(Ni,1); tjr2_2 = nan(Ni,1);
tol = 1;
for i = 1:Ni
    if (abs(d1.err(i)) <= tol)
        tjr1_1(i) = d1.tjr(i,1);
        tjr1_2(i) = d1.tjr(i,2);
    end
    
    if (abs(d2.err(i)) <= tol)
        tjr2_1(i) = d2.tjr(i,1);
        tjr2_2(i) = d2.tjr(i,2);
    end
end

figure(); hold on
set(gca,'FontSize',16)
plot(100*(d1.biasing./1000)/d1.gnc.ha_tgt, tjr1_1,'b--','LineWidth',2);
p1 = plot(100*(d1.biasing./1000)/d1.gnc.ha_tgt, tjr1_2,'b','LineWidth',2);
plot(100*(d2.biasing./1000)/d2.gnc.ha_tgt, tjr2_1,'r--','LineWidth',2);
p2 = plot(100*(d2.biasing./1000)/d2.gnc.ha_tgt, tjr2_2,'r','LineWidth',2);
xlabel('Stage 1 Biasing (% of Target)')
ylabel('t_{jett}^2/t_f')
legend([p1 p2],'\beta_2/\beta_1 = 3','\beta_2/\beta_1 = 9','location','nw')
title(['2000 km, EFPA=-5.4'])

%}

%% Biasing Trade Studies
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [400 2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 2;
v0 = 11;

for z = 1:Nz
biasing = ha_tgts(z) .* [0.1 0.2 0.3 0.4 0.5] .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
% biasing = ha_tgts(z) * 0.25 * -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac();
x0.v0 = v0;
% efpas = [-5.35:-0.0025:-5.4, ...
%     -5.4:-0.01:-5.52, ...
%     -5.52:-0.0025:-5.6, ... 
%     -5.6:-0.05:-6];

% efpas = [-5.3:-0.0025:-5.46, ...
%     -5.46:-0.02:-5.8];
% efpas = -5.4;

efpas =linspace(-5.35 ,-5.85, 30);
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

sim.traj_rate = 500;
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
tjr(i,k,:) = out.idj/out.idxend;
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

%% 2stage trade study post-processing
%{
clear;clc

tol = 5;
plot_tj_opt(2, 400, 3, tol, true);
plot_tj_opt(2, 400, 9, tol, true);
plot_tj_opt(2, 2000, 3, tol, true);
plot_tj_opt(2, 2000, 9, tol, true);



%}

%% 2stage Monte Carlo
%{
clearvars -except save_path ;clc

% ha_tgts = [400 2000 10000];
ha_tgts = [2000];
Nz = length(ha_tgts);
b31 = 10;
b21s = [3 9];
% b21s = 2;
v0 = 11;

for z = 1:Nz
% biasing = ha_tgts(z) .* [0.1 0.2 0.3 0.4 0.5] .* -1000;
% biasing = ha_tgts(z) .* [0.4 0.5] .* -1000;
biasing = ha_tgts(z) * [0.1 0.3] * -1000;

[x0,aero,gnc,sim, mc] = base_venus_ac();
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
    num2str(abs(biasing(z)/1000)) '.mat'])
fprintf('Finished 2stage bias MC. '); print_current_time();

end %biasing
% betas = ['b_ijs' num2str(b21) '_' num2str(b31)];
end %fpas

end %betas

end %ha_tgts
keyboard

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 3-STAGE JETTISON %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bias Value Trade Studies
%{
clearvars -except save_path ;clc
[x0,aero,gnc,sim, mc] = base_venus_ac();
gnc.n = 3;
br41 = 10;
% br21s = [3 3 9];
% br31s = [5 7 9.5];
br21s = 9;
br31s = 9.5;

x0.fpa0 = -5.4;
ha_tgts = [2000];

for jj = 1:length(ha_tgts)
% first stage bias
biasing1 = ha_tgts(jj) .* linspace(0.001, 0.6, 20)' .* -1000;
Nkk = length(biasing1);
% second stage
stage2biases = linspace(0.001, 1, 50);    % percentage of 1st
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

gnc.bias(1).tgt_ap = biasing1(kk);
gnc.bias(2).tgt_ap = biasing2(kk,zz);

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
[x0,aero,gnc,sim, mc] = base_venus_ac();
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

%% 3stage post-processing
%{
clear;clc

tol = 5;
plot_tj_opt(3, 400, [3 5], tol, true);
plot_tj_opt(3, 400, [3 7], tol, true);
plot_tj_opt(3, 400, [9 9.5], tol, true);

plot_tj_opt(3, 2000, [3 5], tol, true);
plot_tj_opt(3, 2000, [3 7], tol, true);
plot_tj_opt(3, 2000, [9 9.5], tol, true);

%}

%% 3stage - investigate expected_haf_err vs time
%{
% It was an issue with the apoapis error tolerance being too high

clear;clc
[x0,aero,gnc,sim, mc] = base_venus_ac();
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




