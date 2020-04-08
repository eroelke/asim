% DCF set up
% determine 20 evenly spaced FPA values within a corridor to get as close
% as possible to the target apoapsis
clear; clc;

fpa_corr = [-5.3 -5.7];    %deg

g_earth = 9.81;

Ni = 20;
fpas = linspace(fpa_corr(1),fpa_corr(2),Ni)';  % fpas to investigate
ha_tgt = 2000;

N_mc = 1;
% Ni = length(ha_tgts);

ha_tol = 10e3;  % orbit tolerance, 10km

% v_atm = 11e3; %m/s
% h_atm = 150e3;  %m
% m1 = 80;
% m2 = 40;
[x0,aero,gnc,sim,mc] = base_jsr();
sim.efpa_flag = false;
mc.flag = false;
sim.traj_rate = 500;
sim.data_rate = 10;
gnc.npc_mode = uint8(1);    %newton
gnc.guid_rate = 0.2;

% iter = uint8(10);

t_steps = [10 20];
t_incs = [0.01 2];

t_jett = nan(Ni,1);
haf_err = t_jett;
haf = t_jett;

g_venus = 8.870;    % graviational acceleration at Venus surface, m/s2
g_earth = 9.81;



delta_t = 10;           % time difference between g1 and g2 measurements

g1 = 0.3;               % first measurement g-load
g1_idx = zeros(Ni,length(delta_t));   % index at which g1 measurement occurs
time = nan(Ni,1000);
decel = time;
g2 = nan(Ni, length(delta_t));       % g2 values (g's)
ttogo = nan(Ni, length(delta_t));     % Time-to-go from g2 time to jettison (sec)

% temp = load('C:\Users\Evan Roelke\Documents\research\venus_ac\dej_n\DCFJ\data\opt_t_jett_betaRatios.mat');
% tj_opts = temp.tj
tol = 10;
veh.m = [150 60];
veh.r = [1 0.2];
% veh.m = [80 40];
% veh.r = [0.75 .175];
veh.rn = 0.175/2;
veh.cd = 1.05;
b1 = veh.m(1)/(veh.cd * pi * veh.r(1)^2);
b2 = veh.m(2)/(veh.cd * pi * veh.r(2)^2);
% 
% % find jettison time with NPC
for i = 1:Ni
    
%     tj0 = 160-5*(i-1);
    gnc.tj0 = 160 - 5*(i-1);
%     [t_jett(i),~,haf_err(i),out] = find_opt_tj(ha_tgt,fpas(i),tj0,tol,veh);
%     x0.fpa0 = fpas(i);
    %run numerical predictor-corrector
    x0.fpa0 = fpas(i);
    out = run_dej_n(x0,gnc,aero,sim,mc);
%     
%     % find jettison time
%     ind = find(out.g.dej_n.tjett(1) == 1,1); % find index of jettison event
%     if isempty(ind)
%         fprintf('No Jettison\n')
%     else
%         t_jett(i) = out.traj.time(ind);   % jettison timing
%     end
%     
%     % find last index
%     if isnan(out.traj.alt(end))
%         idxend = find(isnan(out.traj.alt),1) - 1;
%     else
%         idxend = length(out.traj.alt);
%     end
    
    haf(i) = out.haf;         % resulting orbital apoapsis
    t_jett(i) = out.t_jett;
%     err(i) = haf(i) - ha_tgt;   % delta r_a
    
    t_jett(17) = mean([t_jett(16) t_jett(18)]);
    % get decel and time measurements
    for j = 1:out.idxend
        time(i, j) = out.traj.time(j);      % time vector
        decel(i, j) = out.traj.g_loading(j);    % deceleration vector
        % Determines index at which deceleration surpasses g1 value, and
        % sets g2 value
        if (decel(i, j)*g_earth >= g1 && g1_idx(i) == 0)
            g1_idx(i) = j;
        end
    end
end

% figure(); hold on
% plot(t_jett,err./1000)
% xlabel('Jettison Time (s)')
% ylabel('Apoapse Err (km)')

% keyboard

% show jettison times
% [t_jett err]

%% investigate jettison times using manual jettison
clearvars -except t_jett Ni fpas ha_tgt delta_t

Nj = 500;
t_min = nan(Ni,Nj);  % minimum jettison times
t_diff = t_min;
errs = t_min;

g_earth = 9.81;
delta_t = 10;           % time difference between g1 and g2 measurements

g1 = 0.3;               % first measurement g-load
g1_idx = zeros(Ni,length(delta_t));   % index at which g1 measurement occurs
time = nan(Ni,1000);
decel = time;
g2 = nan(Ni, length(delta_t));       % g2 values (g's)
ttogo = nan(Ni, length(delta_t));     % Time-to-go from g2 time to jettison (sec)

haf_err = t_jett; tjr = t_jett;

traj_rate = 200;

ha_tgt = 2000;

veh.m = [150 60];
% veh.m = [80 40];
veh.r = [1 0.2];
veh.rn = 0.175/2;
veh.cd = 1.05;
tj0 = 170;
tol = 10;



for i = 1:Ni
    
%     times = linspace(t_jett(i)-(i*5), t_jett(i)+(i*5),Nj);
%     for j = 1:Nj
%         if times(j) <= 0
%             times(j) = [];
%             break;
%         end
%         out = ac_venus_manual( fpas(i),times(j),traj_rate, ha_tgt );
%         % get apoapse error
%         errs(i,j) = abs(out.haf - ha_tgt);    % apoapsis error, m
%     end
    tj_init = tj0 - 5*(i-1);
    [t_jett(i),tjr(i),haf_err(i),out] = find_opt_tj(ha_tgt,fpas(i),tj_init,tol,veh);
    
    for j = 1:out.idxend
        time(i, j) = out.traj.time(j);      % time vector
        decel(i, j) = out.traj.g_loading(j);    % deceleration vector
        % Determines index at which deceleration surpasses g1 value, and
        % sets g2 value
        if (decel(i, j)*g_earth >= g1 && g1_idx(i) == 0)
            g1_idx(i) = j;
        end
    end
    
    % find minimum error index
%     min_ind = find(errs(i,j) == min(errs(i,j)));
%     t_min(:,j) = times(min_ind);  % exact jettison time, s
%     t_diff(:,j) = t_min(:,j) - t_jett(i);   % difference in jettison time
    
end

save_path = 'C:\Users\Evan Roelke\Documents\research\venus_ac\';
save([save_path 'dej_n/DCFJ/data/opt_t_jett_beta10.mat'], ... 
    't_jett','fpas','haf_err','g1_idx','time','decel','g1');
% keyboard;

% figure(); hold on
% set(gca,'FontSize',14)
% yyaxis left
% plot(fpas, t_jett,'LineWidth',2)
% ylabel('Jettison Time (s)')
% yyaxis right
% plot(fpas, haf_err,'LineWidth',2)
% ylabel('Apoapsis Error (km)')
% ylim([-5 5])
% legend('\beta_2/\beta_1 = 10','location','nw')
% xlabel('Entry Flight-Path Angle (deg)')

% keyboard;

%% Data Processing
% Parses through deceleration curves to obtain g2 and t-go values

% load('../venus_ac/dcfj/data/opt_t_jett.mat');
Ni = length(fpas);
delta_t = 10;

for i = 1:Ni
    t1 = time(i, g1_idx(i));    % Time at g1 (s)
    t2 = t1 + delta_t;          % Time at g2 (s)
    g2_idx = find(round(time(i, :), 2) == round(t2, 2),1); % Index of g2 value
    g2(i) = decel(i, g2_idx)*9.81;
    ttogo(i) = t_jett(i) - t2;
    
    if ttogo(i) < 0
        g2(i) = [];
        ttogo(i) = [];
    end
end
% Least squares fit:
p = polyfit(g2(4:end), ttogo(4:end), 3);  % Least-squares polynomial variables (ttogo as func of g2)
yfit = polyval(p, g2(4:end));   % Linear least-squares y vals
yresid = ttogo(4:end) - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(ttogo(4:end))-1)*var(ttogo(4:end));
rsq = 1 - SSresid/SStotal;  % R^2 value of fit

figure(1); hold on; 
grid on;
set(gca, 'FontSize', 14);
plot(ttogo(4:end), g2(4:end)/9.81, 'o', 'LineWidth', 2);
plot(yfit, g2(4:end)/9.81, '-k', 'LineWidth', 2);
xlabel('t_{go} (s)');
ylabel('g_{2} Measurement (Earth g)');
% text(40,1.3,['R^{2} = ' num2str(rsq)]);
% text(40,1.35,['\Delta t = ' num2str(delta_t)]);
dim = [0.7 0.5 1 0.1];
str = {['R^{2} = ' num2str(rsq)],['\Delta t = ' num2str(delta_t)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

save([save_path 'dej_n/DCFJ/data/VenusDCF_dt10_bratio10'],'p')






