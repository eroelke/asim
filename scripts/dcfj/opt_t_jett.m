% find optimal jettison time with ac_venus_manual
% high traj rate and guidance rate
clear;clc;
% fpas = linspace(-5.2,-5.6,20);
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';

% temp = load('C:\Users\Evan Roelke\Documents\research\asim\venus_ac\ac_guid_algorithms\DCF\opt_t.mat');
% t0 = temp.t_min;
% for i = 1:length(t0)
%    times(:,i) = (t0(i)-5):0.1:(t0(i)+5); 
% end

% dt = 0.1;
% t1 = 25;
% t2 = 170;
% times = t1:dt:t2;   % jettison time events
% times = times';

% 1000Hz precision optimal jettison time
% temp = load('.\venus_ac\ac_guid_algorithms\NPC\nominal_errs.mat');
temp = load([save_path 'dej_n/DCFJ/data/opt_t_jett_beta10.mat']);
fpas = temp.fpas;

ha_tgts = [2000 10000];  % apoapse targets

Ni = length(fpas);
% Nj = length(times);
% Nj = 100;
% times1 = nan(Ni,Nj);
% times2 = times1;

haf2k = nan(Ni,1);
err2k = haf2k;
haf10k = haf2k;
err10k = haf2k;

t2k = nan(Ni,1);
t10k = t2k;

guid_rate = 0.5;
traj_rate = 1000;  %4 decimal precision (5 causes memory issues)
tol_ha = 10;
n = 1;
tj0 = 85;
% m = [80 40];
% rcs = [0.75 0.175];
m = [150 60];
rcs = [1 0.2];
iters = uint8(1);
root_mode = uint8(1);   %newton
mc.flag = uint8(0); %nominal only

[x0,aero,gnc,sim,mc] = base_jsr();
sim.efpa_flag = false;
mc.flag = false;
sim.traj_rate = 500;
sim.data_rate = 10;
gnc.npc_mode = uint8(1);    %newton
gnc.guid_rate = 0.2;

for i = 1:Ni
    fprintf('Run %i/%i\n',i,Ni)
    
    x0.fpa0 = fpas(i);
    
    gnc.ha_tgt = 2000;
    out1 = run_dej_n(x0,gnc,aero,sim,mc);
    
    gnc.ha_tgt = 10000;
    out2 = run_dej_n(x0,gnc,aero,sim,mc);
%     tic
    % run NPC at very large rate
%     out1 = dej_n_auto_venus(fpas(i), ha_tgts(1), ... 
%         n, tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
%     out2 = dej_n_auto_venus(fpas(i), ha_tgts(2), ... 
%         n, tj0, rcs, m, root_mode,guid_rate,mc,iters,traj_rate);
%     toc
%     out1 = sej_auto_venus(fpas(i),ha_tgts(1),10,100,0.25);
%     out2 = sej_auto_venus(fpas(i),ha_tgts(2),10,100,0.25);
    
%     times1(i,:) = round(linspace(out1.t_jett - 5,out1.t_jett + 5,Nj),2);
%     times2(i,:) = round(linspace(out2.t_jett - 5,out2.t_jett + 5,Nj),2);
    
% %     for j = 1:Nj
%         % run sim, jettison at time j, traj rate 1000 Hz
%         out1 = ac_venus_manual( fpas(i), times1(i,j) );
%         out2 = ac_venus_manual( fpas(i), times2(i,j) );
%         
%         haf1(i,j) = out1.haf;    %km
%         err1(i,j) = abs(haf1(i,j) - ha_tgts(1));  %km, err 2000km tgt
%         
%         haf2(i,j) = out2.haf;
%         err2(i,j) = abs(haf2(i,j) - ha_tgts(2));  %km, err 1000km tgt
% %     end
        
    % 2000 km optimal
    haf2k(i) = out1.haf;
    t2k(i) = out1.t_jett;
    err2k(i) = out1.haf_err;

    % 10000 km optimal
    haf10k(i) = out2.haf;
    t10k(i) = out2.t_jett;
    err10k(i) = out2.haf_err;
    
    % find smallest apoapse errors for this fpa, 2000 km
%     min_ind1 = find(err1(i,:) == min(err1(i,:)),1);
%     t2k(i) = times1(i,min_ind1);  % jettison time of min error    
    
    % find smallest apoapse errors for this fpa, 10000 km
%     min_ind2 = find(err2(i,:) == min(err2(i,:)),1);
%     t_min2(i) = times2(i,min_ind2);  % jettison time of min error
    
end
%

save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';
save([save_path '\dej_n\DCFJ\data\tj_optErrs_beta10.mat'], ... 
    't2k','t10k','err2k','err10k','fpas','ha_tgts');

temp = load([save_path '\dej_n\DCFJ\data\tj_optErrs_beta10.mat']);
keyboard

% 1000Hz precision optimal jettison time (errs < 1km)
% temp = load('.\venus_ac\ac_guid_algorithms\NPC\nominal_errs.mat');
% fpas = temp.fpas;
% t2k = temp.t_2k;
% t10k = temp.t_10k;
% err2k = temp.err_2k;

figure(); hold on
set(gca,'FontSize',14)
plot(fpas,t2k,'LineWidth',2);
plot(fpas,t10k,'LineWidth',2);
xlabel('FPA_{atm} (deg)')
ylabel('Minimal Error Jettison Time (sec)')
legend('h_{tgt} = 2000 km','h_{tgt} = 10000 km','location','nw')

figure(); hold on; grid on; box on
set(gca,'FontSize',14)
plot(fpas,t2k,'LineWidth',2)
legend('h_{a,tgt} = 2000km')
xlabel('FPA_{atm} (deg)')
ylabel('Optimal Jettison Time (sec)')





%% different target apoapses
Ni = length(fpas);  % num loops
t_jett = [t2k t10k];
Nk = size(t_jett,2);
g_earth = 9.81;
g1 = 0.3;               % first measurement g-load
delta_t = 10;   % 10 seconds between measurements
traj_rate = 100;
ha_tgt = [2000 10000];
% preallocate
haf = nan(Ni,Nk);
err = haf;

g1_idx = zeros(Ni,Nk);   % index at which g1 measurement occurs
time = nan(Ni,1000,Nk);
decel = time;
g2 = nan(Ni,Nk);       % g2 values (g's)
ttogo = nan(Ni,Nk);     % Time-to-go from g2 time to jettison (sec)
yfit = nan(Ni,Nk);
rsq = nan(Nk,1);
p = nan(4,Nk);
nan_ind = rsq;
tic
veh.m = aero.m;
veh.r = aero.rcs;
veh.cd = aero.cds(1);
veh.rn = aero.rn;

for k = 1:Nk
for i = 1:Ni
    out = ac_venus_manual( fpas(i),t_jett(i,k), traj_rate, ha_tgt(k), veh );
    
    % find last index
%     if isnan(out.traj.alt(end))
%         idxend = find(isnan(out.traj.alt),1) - 1;
%     else
%         idxend = length(out.traj.alt);
%     end
    idxend = out.idxend;
    
    haf(i,k) = out.haf;         % resulting orbital apoapsis, km
    err(i,k) = haf(i,k) - ha_tgt(k);   % delta r_a
    
    % get decel and time measurements
    for j = 1:idxend
        time(i, j, k) = out.traj.time(j);      % time vector
        decel(i, j, k) = out.traj.g_loading(j);    % deceleration vector
        % Determines index at which deceleration surpasses g1 value, and
        % sets g2 value
        if (decel(i, j,k)*g_earth >= g1 && g1_idx(i,k) == 0)
            g1_idx(i,k) = j;
        end
    end 
end
toc
%% Data Processing
% Parses through deceleration curves to obtain g2 and t-go values

for i = 1:Ni
    t1 = time(i, g1_idx(i,k));    % Time at g1 (s)
    t2 = t1 + delta_t;          % Time at g2 (s)
    g2_idx = find(round(time(i, :,k), 2) == round(t2, 2),1); % Index of g2 value
    g2(i,k) = decel(i, g2_idx)*9.81;
    ttogo(i,k) = t_jett(i,k) - t2;
    
    if ttogo(i,k) < 0
        g2(i,k) = nan;
        ttogo(i,k) = nan;
    end
    
%     if i > 1 && ttogo(i,j,k) < ttogo(i-1,j,k)
%         g2(i,j,k) = nan;
%         ttogo(i,j,k) = nan;
%     end
    
end

ind = find(isnan(ttogo(:,k)) == 1,1)-1;
if isempty(ind)
    nan_ind(k) = Ni;
else
    nan_ind(k) = ind;
end
% Least squares fit:
p(:,k) = polyfit(g2(:,k), ttogo(:,k), 3)';  % Least-squares polynomial variables (ttogo as func of g2)
yfit(:,k) = polyval(p(:,k), g2(:,k));   % Linear least-squares y vals
yresid = ttogo(:,k) - yfit(:,k);
SSresid = sum(yresid.^2);
SStotal = (length(ttogo(:,k))-1)*var(ttogo(:,k));
rsq(k) = 1 - SSresid/SStotal;  % R^2 value of fit
end %k

figure(1); hold on;
grid on;
set(gca, 'FontSize', 14);
p1 = plot(tgo_2k, g2_2k/9.81, 'o', 'LineWidth', 2);
p2 = plot(yfit_2k, g2_2k/9.81, '-k', 'LineWidth', 2);
p3 = plot(tgo_10k, g2_10k/9.81, 'o', 'LineWidth', 2);
p4 = plot(yfit_10k, g2_10k/9.81, '-k', 'LineWidth', 2);
xlabel('t_{go} (s)');
ylabel('g_{2} Measurement (Earth g)');
legend([p1 p3],{[num2str(ha_tgt(1)) ' km Target'],[num2str(ha_tgt(2)) ' km Target']}, ... 
    'location','ne')
% text(40,1.3,['R^{2} = ' num2str(rsq)]);
% text(40,1.35,['\Delta t = ' num2str(delta_t)]);
dim = [0.6 0.65 1 0.1];
str = {['R^{2}_{2k} = ' num2str(rsq(1))], ... 
    ['R^{2}_{10k} = ' num2str(rsq(2))], ... 
    ['\Delta t = ' num2str(delta_t) ' s']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% saveas(gca,'C:\Users\Evan Roelke\Documents\research\venus_aerocapture\figures\DCF\curve_ha_tgt.fig')

save([save_path 'dej_n\DCFJ\data\curve_ha_tgt.mat'], ... 
    'ttogo', 'g2', 'rsq', 'delta_t', 'ha_tgt','p','yfit');

% fix up data
clear;clc;
save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';
load([save_path 'dej_n\DCFJ\data\curve_ha_tgt.mat']);
sI = 1;
ttogo(1:sI,:) = [];
g2(1:sI,:) = [];
g2_10k = g2(:,2);
tgo_10k = ttogo(:,2);

p_10k = polyfit(g2_10k, tgo_10k, 3)';  % Least-squares polynomial variables (ttogo as func of g2)
yfit_10k = polyval(p_10k, g2_10k);   % Linear least-squares y vals
yresid = tgo_10k - yfit_10k;
SSresid = sum(yresid.^2);
SStotal = (length(tgo_10k-1)*var(tgo_10k));
rsq_10k = 1 - SSresid/SStotal;  % R^2 value of fit


temp = load([save_path 'dej_n\DCFJ\data\curve_dt.mat']);
g2_2k = temp.g2(2:end,1);
tgo_2k = temp.ttogo(2:end,1);
p_2k = polyfit(g2_2k, tgo_2k, 3)';
yfit_2k = polyval(p_2k, g2_2k);
yresid = tgo_2k - yfit_2k;
SSresid = sum(yresid.^2);
SStotal = (length(tgo_2k-1)*var(tgo_2k));
rsq_2k = 1 - SSresid/SStotal;  % R^2 value of fit

save('data/dcfj_data/curve_ha_tgt.mat','tgo_2k','g2_2k','p_2k','tgo_10k','g2_10k','tgo_10k','p_10k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELTA T TRADE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% as func of delta t
% clear;clc
% temp = load('.\venus_ac\ac_guid_algorithms\NPC\nominal_errs.mat');
% fpas = temp.fpas;
% t2k = temp.t_2k;
% t10k = temp.t_10k;
% err2k = temp.err_2k;

Ni = length(fpas);  % num loops
t_jett = t2k;
g_earth = 9.81;
g1 = 0.3;               % first measurement g-load
delta_t = [10 15 20];   % 10 seconds between measurements
traj_rate = 100;
ha_tgt = 2000;
Nj = length(delta_t);
% preallocate
haf = nan(Ni);
err = haf;

g1_idx = zeros(Ni,1);   % index at which g1 measurement occurs
time = nan(Ni,1000);
decel = time;
g2 = nan(Ni,Nj);       % g2 values (g's)
ttogo = nan(Ni,Nj);     % Time-to-go from g2 time to jettison (sec)
yfit = nan(Ni,Nj);
rsq = nan(Nj,1);
p = nan(4,Nj);
for i = 1:Ni
    out = ac_venus_manual( fpas(i),t_jett(i), traj_rate, ha_tgt, veh );
    
    % find last index
%     if isnan(out.traj.alt(end))
%         idxend = find(isnan(out.traj.alt),1) - 1;
%     else
%         idxend = length(out.traj.alt);
%     end
    idxend = out.idxend;
    
    haf(i) = out.haf;         % resulting orbital apoapsis, km
    err(i) = haf(i) - ha_tgt;   % delta r_a
    
    % get decel and time measurements
    for k = 1:idxend
        time(i, k) = out.traj.time(k);      % time vector
        decel(i, k) = out.traj.g_loading(k);    % deceleration vector
        % Determines index at which deceleration surpasses g1 value, and
        % sets g2 value
        if (decel(i, k)*g_earth >= g1 && g1_idx(i) == 0)
            g1_idx(i) = k;
        end
    end 
end

%% Data Processing
% Parses through deceleration curves to obtain g2 and t-go values
for j = 1:Nj

for i = 1:Ni
    t1 = time(i, g1_idx(i));    % Time at g1 (s)
    t2 = t1 + delta_t(j);          % Time at g2 (s)
    g2_idx = find(round(time(i,:), 2) == round(t2, 2),1); % Index of g2 value
    g2(i,j) = decel(i, g2_idx)*9.81;
    ttogo(i,j) = t_jett(i) - t2;
    
    if ttogo(i,j) < 0
        g2(i,j) = nan;
        ttogo(i,j) = nan;
    end
    
%     if i > 1 && ttogo(i,j,k) < ttogo(i-1,j,k)
%         g2(i,j,k) = nan;
%         ttogo(i,j,k) = nan;
%     end
    
end

% ind = find(isnan(ttogo(:,j)) == 1,1)-1;
% if isempty(ind)
%     nan_ind(j) = Ni;
% else
%     nan_ind(j) = ind;
% end
% Least squares fit:
p(:,j) = polyfit(g2(:,j), ttogo(:,j), 3)';  % Least-squares polynomial variables (ttogo as func of g2)
yfit(:,j) = polyval(p(:,j), g2(:,j));   % Linear least-squares y vals
yresid = ttogo(:,j) - yfit(:,j);
SSresid = sum(yresid.^2);
SStotal = (length(ttogo(:,j))-1)*var(ttogo(:,j));
rsq(j) = 1 - SSresid/SStotal;  % R^2 value of fit
end %j

figure(); hold on;
grid on;
set(gca, 'FontSize', 14);
p1 = plot(ttogo(:,1), g2(:,1)/9.81, 'o', 'LineWidth', 2);
p2 = plot(yfit(:,1), g2(:,1)/9.81, '-k', 'LineWidth', 2);
p3 = plot(ttogo(:,2), g2(:,2)/9.81, 'o', 'LineWidth', 2);
p4 = plot(yfit(:,2), g2(:,2)/9.81, '-k', 'LineWidth', 2);
p5 = plot(ttogo(:,3), g2(:,3)/9.81, 'o', 'LineWidth', 2);
p6 = plot(yfit(:,3), g2(:,3)/9.81, '-k', 'LineWidth', 2);
xlabel('t_{go} (s)');
ylabel('g_{2} Measurement (Earth g)');
legend([p1 p3 p5],{['\Delta t = ' num2str(delta_t(1)) 's'], ... 
    ['\Delta t = ' num2str(delta_t(2)) 's'], ...
    ['\Delta t = ' num2str(delta_t(3)) 's']}, ...
    'location','ne')
dim = [0.6 0.6 1 0.1];
str = {['R^{2}_{10s} = ' num2str(rsq(1))], ... 
    ['R^{2}_{15s} = ' num2str(rsq(2))], ... 
    ['R^{2}_{20s} = ' num2str(rsq(3))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');


save([save_path 'dej_n\DCFJ\data\curve_dt.mat'], ... 
    'ttogo', 'g2', 'rsq', 'delta_t', 'ha_tgt','p','yfit');


% % fix up data
% save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';
% load([save_path 'dej_n\DCFJ\data\curve_dt.mat']);
% sI = 3;
% ttogo(1:sI,:) = [];
% g2(1:sI,:) = [];



